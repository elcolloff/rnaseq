# rustar-aligner integration: observed differences vs STAR

Running notes captured while wiring up the experimental `--use_rustar_star`
path in PR #1855. The intent here is to track every divergence we observe so
that nothing surprises us at review and so we can file targeted upstream
issues at https://github.com/scverse/rustar-aligner (or fix things on our
side) rather than discovering them in production. This will be cleaned up
before merge.

The verification setup: standard `-profile test,docker` on the
`nf-dev-rnaseq` VM (36 CPU / 69 GB), back-to-back STAR and rustar runs,
identical inputs.

## Verified

### Wall-time and RAM (test profile, one tile of yeast + GFP)

From `pipeline_info/execution_trace_*.txt`, comparing the per-task medians:

| Process              | n | Wall (s) STAR → rustar | Peak RSS (GB) STAR → rustar |
|----------------------|---|------------------------|------------------------------|
| `STAR_GENOMEGENERATE` / `RUSTAR_GENOMEGENERATE` | 1 | 0.3 → 0.3 | 0.01 → 0.02 |
| `STAR_ALIGN` / `RUSTAR_ALIGN`                   | 5 | 68.0 → 33.8 | 0.92 → 0.12 |

Caveat: this is on the tiny test genome (a yeast subset plus GFP transgene)
with ≤10 k reads per sample, run inside Docker. The absolute numbers say
nothing about human-scale performance. Re-running on the `test_full`
samplesheet on AWS is a follow-up.

### Mapping rate (per `Log.final.out`)

| Sample              | STAR  | rustar | Δ (pp) |
|---------------------|-------|--------|--------|
| RAP1_IAA_30M_REP1   | 90.44 | 90.23  | -0.21  |
| RAP1_UNINDUCED_REP1 | 95.96 | 95.88  | -0.08  |
| RAP1_UNINDUCED_REP2 | 95.85 | 95.80  | -0.05  |
| WT_REP1             | 88.99 | 88.81  | -0.18  |
| WT_REP2             | 89.54 | 89.39  | -0.15  |

All within ±0.25 pp of STAR. Consistent with what rustar reports upstream
on its yeast 10 k-read benchmark.

### Quantification concordance (per-sample Pearson on merged Salmon matrices)

| Sample              | gene_tpm | gene_counts |
|---------------------|----------|-------------|
| RAP1_IAA_30M_REP1   | 0.996808 | 0.999848    |
| RAP1_UNINDUCED_REP1 | 0.999673 | 0.999904    |
| RAP1_UNINDUCED_REP2 | 0.999746 | 0.999906    |
| WT_REP1             | 0.995496 | 0.999890    |
| WT_REP2             | **0.985040** | 0.999842 |

`gene_counts` (raw `NumReads`) is essentially identical across both
runs. `gene_tpm` is also very close on three samples but diverges
materially on `WT_REP2`, with `RAP1_IAA_30M_REP1` and `WT_REP1` showing
the same effect at smaller magnitude. The two single-end samples
(`RAP1_UNINDUCED_REP1/2`) are clean.

This is **not** sample-specific. See
[`rustar_investigation_wt_rep2.md`](rustar_investigation_wt_rep2.md)
for the deep dive - short version: rustar v0.1.0's
`Aligned.toTranscriptome.out.bam` doesn't populate mate-pair fields
(`RNEXT` / `PNEXT` / `TLEN`) or set the proper-pair flag on paired-end
records, so Salmon can't infer a fragment-length distribution and
falls back to its default prior (mean 250, SD 25). That distorts
`EffectiveLength` for short transcripts, which is what we see in TPM.
The hit is bigger on `WT_REP2` because of how its mapped reads
distribute across short vs long transcripts.

This is the headline bug to file upstream once the report is in
shape; everything below is secondary.

## Module-level workarounds we had to add

These are deltas baked into `modules/local/rustar_align/` so the rustar
modules slot into the existing `ALIGN_STAR` subworkflow without
collateral damage. They are not user-visible. The goal is to keep them
small and clearly marked so they can be retired as rustar tightens its
STAR compatibility.

### `--limitGenomeGenerateRAM` is not accepted

STAR exposes `--limitGenomeGenerateRAM`; the upstream `STAR_GENOMEGENERATE`
module derives a value from `task.memory` and passes it. rustar v0.1.0
rejects this flag at startup (`error: unexpected argument
'--limitGenomeGenerateRAM' found`).

`modules/local/rustar_align/genomegenerate/main.nf` therefore omits the
flag and relies on rustar's built-in memory management. We should
check whether this matters on full-size genomes.

### `--outFileNamePrefix` ending in `.` is treated as a directory

STAR treats `--outFileNamePrefix SAMPLE.` as a literal string prefix and
writes `SAMPLE.Aligned.out.bam`, `SAMPLE.Log.final.out`, etc. side by
side in the work directory.

rustar v0.1.0 instead interprets the same value as a directory name and
writes bare-named files inside it:

```
SAMPLE./
  Aligned.out.bam
  Aligned.toTranscriptome.out.bam
  Log.final.out
  SJ.out.tab
  SJ.pass1.out.tab
```

`modules/local/rustar_align/align/main.nf` post-processes by flattening
that directory back into STAR-style prefixed filenames so the downstream
emit globs (`*Log.final.out`, etc.) still match. Worth filing upstream.

### `Log.out` and `Log.progress.out` are not written

STAR emits three log files: `Log.final.out` (summary stats, MultiQC
input), `Log.out` (verbose run log) and `Log.progress.out` (per-chunk
progress). rustar v0.1.0 only writes `Log.final.out`.

Marked `Log.out` / `Log.progress.out` as `optional: true` outputs in
`RUSTAR_ALIGN`. Nothing in the pipeline currently consumes them, but if
that changes we'll need to re-evaluate.

### Extra `SJ.pass1.out.tab` is emitted

rustar writes both `SJ.out.tab` and `SJ.pass1.out.tab` (the two-pass
intermediate). STAR keeps the intermediate inside `<prefix>_STARpass1/`
rather than at the top level. Currently the rustar one is caught by the
existing `*.tab` glob and silently emitted - harmless but unusual.

### Version reporting

The rustar container (`ghcr.io/scverse/rustar-aligner:dev` on debian-slim)
does not bundle `samtools` or `gawk`, which are present in the STAR Wave
container. STAR_GENOMEGENERATE uses `samtools faidx` + `gawk` to
auto-compute `--genomeSAindexNbases`.

To avoid adding a `samtools`/`gawk` dependency to the rustar image,
`RUSTAR_GENOMEGENERATE` does the same heuristic in Groovy from the
on-disk FASTA size. The approximation is well inside the floor() of
`log2(len)/2 - 1` so the chosen index size matches.

`RUSTAR_ALIGN` emits only the `rustar-aligner` version through the
topic-based versions channel - no `samtools` / `gawk` emissions.

## Nextflow-side, not rustar's fault, but bites us anyway

### Boolean CLI flags get coerced to the string `"true"`

`--use_rustar_star`, `--use_rustar_star=true`, and
`--use_rustar_star true` all fail nf-schema validation with `Value is
[string] but should be [boolean]` on Nextflow 26.04 + nf-schema 2.6.1.
This is not rustar-specific; the same error occurs for
`--use_parabricks_star`. A YAML params file works:

```yaml
use_rustar_star: true
outdir: results-rustar
```

then `nextflow run ... -params-file rustar.params.yml`. Worth raising
upstream (Nextflow / nf-schema), separately from rustar.

## Still to verify

- Full-size run on the `test_full` samplesheet (GRCh37, larger reads) to
  produce performance and concordance numbers that map to user
  expectations. The test-profile numbers above are not load-bearing.
- Whether the `--limitGenomeGenerateRAM` omission matters at human-genome
  scale.
- Whether rustar's `--quantTranscriptomeSAMoutput BanSingleEnd` matches
  STAR's interpretation byte-for-byte. Almost certainly fine, but worth
  a glance once the paired-end mate-field bug is fixed.

## Tracked upstream

A categorical BAM-level comparison is in
[`rustar_bam_comparison.md`](rustar_bam_comparison.md). Punch list (file
ordering loosely by severity):

- **Paired-end transcriptome BAM missing mate fields** ->
  [scverse/rustar-aligner#22](https://github.com/scverse/rustar-aligner/issues/22)
  (filed). Full evidence in
  [`rustar_investigation_wt_rep2.md`](rustar_investigation_wt_rep2.md).
  Salmon falls back to its default fragment-length prior, distorting
  paired-end TPMs.
- **`NM` tag replaced with `nM`, semantics changed** (high). rustar
  emits `nM` everywhere STAR emits `NM`. STAR's `NM` is the SAM-spec
  edit distance (mismatches + indels); rustar's `nM` counts only
  substitutions. ~2% of identical-CIGAR records disagree. Any tool
  reading `NM:i:` (samtools stats, Picard, MultiQC parsers) silently
  gets nothing.
- **`XS` tag never emitted** despite `--outSAMstrandField intronMotif`
  (high). Breaks StringTie, Cufflinks, rseqc's `infer_experiment.py`.
- **GTF splice junctions not seeded into pass 1** (medium). ~50% of
  splices dropped; `Number of splices: Annotated (sjdb) = 0` in every
  `Log.final.out`. Manifests as per-read CIGAR collapse from spliced
  to plain `M` on ~70% of CIGAR diffs in WT_REP2.
- **~17% more secondary alignments than STAR** (medium, possibly
  intentional). `NH` tail extends to 20 vs STAR's 7 on the same data;
  no `outFilterMultimapScoreRange`-equivalent appears to be applied.
  Worth documenting upstream if intentional.
- **`@PG` header is content-free** (low). Just `ID:rustar-aligner`,
  no `PN`/`VN`/`CL`. Provenance loss.
- **Transcriptome BAM missing per-record `RG:Z:`** (low). Sample
  association lost on transcriptome records. Genome BAM is fine.
- **`AS` tag drifts by ±2-5 score units on identical CIGAR** (low),
  864 records. rustar's README promises byte-equivalent `AS` on
  unique mappers, so this is a small contract break.
