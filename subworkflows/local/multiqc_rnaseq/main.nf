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
    ch_expected_count          // channel: [ id, groupKey(id, n), n ] per sample

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

    // Everything MultiQC will ingest: the collected per-sample QC outputs
    // plus the global context files, tagged with an empty meta.
    ch_multiqc_all = ch_multiqc_files.mix(
        ch_workflow_summary
            .mix(ch_collated_versions)
            .mix(ch_methods_description)
            .map { f -> [[:], f] }
    )

    // --replace-names TSV so MultiQC uses sample IDs rather than FASTQ basenames.
    ch_name_replacements = multiqcNameReplacements(ch_fastq)

    if (skip_quantification_merge) {
        // One MultiQC report per sample. Split incoming files into per-sample
        // and global buckets, then attach the global bucket to every sample.
        //
        // Each per-sample file is tagged with a per-sample groupKey supplied
        // by the caller so groupTuple can close each sample's group as soon
        // as its expected items arrive, instead of waiting for the slowest
        // sample in the run to release the upstream ch_multiqc_files channel.
        ch_branched = ch_multiqc_all
            .branch { meta, _file ->
                per_sample: meta.id != null
                global: true
            }

        ch_global_files = ch_branched.global
            .map { _meta, f -> f }
            .collect()

        // Build a single value-channel lookup map from id -> [groupKey, expected_count].
        // Using a value channel avoids the Cartesian-broadcast behaviour we saw with
        // `.combine(by:)` on queue-to-queue, which duplicated left-side items per
        // subscriber.
        ch_expected_count_map = ch_expected_count
            .map { id, key, n -> [(id): [key, n]] }
            .reduce([:]) { acc, entry -> acc + entry }

        ch_multiqc_input = ch_branched.per_sample
            .combine(ch_expected_count_map)
            .map { meta, f, count_map ->
                def entry = count_map[meta.id]
                def key   = entry ? entry[0] : groupKey(meta.id, 0)
                def n     = entry ? entry[1] : 0
                [key, f, n]
            }
            .groupTuple()
            .map { key, files, ns ->
                // Flatten because some upstream contributors emit bundled
                // [meta, list_of_files] tuples (e.g. RUSTQC) and others emit
                // [meta, single_file]; both forms land in the grouped list.
                def id = key.toString()
                def flat = files.collectMany { it instanceof List ? it : [it] }
                def expected = ns ? ns[0] : 0
                if (expected > 0 && flat.size() != expected) {
                    log.warn "[nf-core/rnaseq] MultiQC per-sample file count drift for '${id}': expected ${expected}, got ${flat.size()}. Update perSampleMultiqcExpectedCount() to match the current ch_multiqc_files contributors."
                }
                [id, flat]
            }
            .combine(ch_global_files.toList())
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
        ch_all_files = ch_multiqc_all
            .map { _meta, f -> f }
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
