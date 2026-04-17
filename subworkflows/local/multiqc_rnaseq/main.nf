//
// MultiQC report assembly for nf-core/rnaseq.
//

include { MULTIQC                 } from '../../../modules/nf-core/multiqc'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc    } from '../../nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../utils_nfcore_rnaseq_pipeline'
include { multiqcNameReplacements } from '../utils_nfcore_rnaseq_pipeline'
include { multiqcSampleMergeYaml  } from '../utils_nfcore_rnaseq_pipeline'

workflow MULTIQC_RNASEQ {

    take:
    ch_multiqc_files           // channel: [ val(meta), path(file_or_file_list) ]
    ch_fastq                   // channel: [ val(meta), [ reads ] ]
    ch_collated_versions       // channel: path(versions yaml)
    samplesheet_path           // path: pipeline input samplesheet
    samplesheet_schema         // path: samplesheet JSON schema
    mqc_default_config         // path: pipeline-bundled MultiQC config
    mqc_custom_config          // path (or []): optional user MultiQC config
    mqc_logo                   // path (or []): optional custom logo
    methods_description_yml    // path: methods-description YAML template
    skip_quantification_merge  // boolean
    ch_expected_count          // channel: [ id, groupKey(id, n) ] per sample

    main:

    // Per-run table_sample_merge config: only PE samples from the samplesheet
    // get their _1/_2 rows grouped in the General Stats table.
    ch_mqc_dynamic_config = channel.of(multiqcSampleMergeYaml(samplesheet_path, samplesheet_schema))
        .collectFile(name: 'multiqc_sample_merge.yml')

    // Workflow summary and methods description rendered as MultiQC sections.
    ch_workflow_summary = channel.value(
        paramsSummaryMultiqc(
            paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        )
    ).collectFile(name: 'workflow_summary_mqc.yaml')

    ch_methods_description = channel.value(
        methodsDescriptionText(methods_description_yml)
    ).collectFile(name: 'methods_description_mqc.yaml')

    // Per-sample MultiQC swaps the full collated versions yaml (which only closes
    // after every task emits into the `versions` topic — the blocker that
    // undermines progressive closure) for a manifest-only stub. The stub keeps
    // per-sample multiqc_software_versions.txt present in each report so the
    // file structure matches dev; its contents are .nftignored.
    ch_static_versions = channel.value(
        "Workflow:\n    ${workflow.manifest.name}: ${workflow.manifest.version}\n    Nextflow: ${workflow.nextflow.version.toString()}\n"
    ).collectFile(name: 'nf_core_rnaseq_software_mqc_versions.yml')

    ch_static_globals = ch_workflow_summary.mix(ch_methods_description)

    // --replace-names TSV so MultiQC uses sample IDs rather than FASTQ basenames.
    ch_name_replacements = multiqcNameReplacements(ch_fastq)

    if (skip_quantification_merge) {
        // One MultiQC report per sample. Each per-sample item carries a
        // caller-supplied groupKey so groupTuple closes that sample's group
        // as soon as its expected files arrive (see perSampleMultiqcExpectedCount).
        // Combined with the ch_static_* globals below — which don't source from
        // ch_multiqc_files and so don't wait for the whole run to close — each
        // sample's MULTIQC fires ASAP instead of waiting for the slowest sample.

        // Value-channel lookup id -> groupKey so combine broadcasts the keys to
        // every per-sample emission without re-materialising the upstream queue.
        ch_sample_keys = ch_expected_count
            .reduce([:]) { acc, row -> acc + [(row[0]): row[1]] }
            .first()

        // Globals available in finite time (value channels + manifest-only versions).
        // Dynamic globals (DESEQ2, fail_mapped_samples_mqc, fail_strand_check_mqc)
        // and the full collated versions yaml are only used in the merged branch.
        ch_per_sample_globals = ch_static_globals.mix(ch_static_versions).collect()

        ch_multiqc_input = ch_multiqc_files
            .filter { meta, _file -> meta.id != null }
            .combine(ch_sample_keys)
            .map { meta, f, keys -> [keys[meta.id] ?: groupKey(meta.id, 0), f] }
            .groupTuple(remainder: true)
            .map { key, files ->
                // Compare key.size against the tuple count (not flat file count):
                // perSampleMultiqcExpectedCount predicts contributor tuples, some of
                // which (e.g. DUPRADAR's pair of _mqc.txt files) are list-valued.
                def id = key.toString()
                if (key.size > 0 && files.size() != key.size) {
                    log.warn "[nf-core/rnaseq] MultiQC per-sample contributor count drift for '${id}': expected ${key.size}, got ${files.size()}. Update perSampleMultiqcExpectedCount() to match the current ch_multiqc_files contributors."
                }
                [id, files.flatten()]
            }
            .combine(ch_per_sample_globals.toList())
            .combine(ch_mqc_dynamic_config)
            .map { id, sample_files, global_files, dyn ->
                [
                    [id: id],
                    sample_files + (global_files ?: []),
                    [mqc_default_config, dyn, mqc_custom_config].findAll { it },
                    mqc_logo,
                    [],  // no replace_names — each report contains one sample's files
                    [],
                ]
            }
    } else {
        // One merged MultiQC report. 'multiqc_report' is a sentinel meta.id
        // used by conf/modules/multiqc.config to pick the merged output
        // path/prefix. Wrap the collected file list in a 1-tuple so
        // .combine() doesn't spread it across the downstream closure args.
        ch_all_files = ch_multiqc_files
            .map { _meta, f -> f }
            .mix(ch_static_globals)
            .mix(ch_collated_versions)
            .collect()
            .map { files -> [files] }

        ch_multiqc_input = ch_all_files
            .combine(ch_name_replacements.ifEmpty([]).toList())
            .combine(ch_mqc_dynamic_config)
            .map { files, replace_names, dyn ->
                [
                    [id: 'multiqc_report'],
                    files,
                    [mqc_default_config, dyn, mqc_custom_config].findAll { it },
                    mqc_logo,
                    replace_names ?: [],
                    [],
                ]
            }
    }

    MULTIQC(ch_multiqc_input)

    emit:
    report = MULTIQC.out.report.map { _meta, report -> report }
}
