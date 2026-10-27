#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Joint PureTarget pipeline — run-level steps:
 *   pt_samplesSummary   runInfo/<yyMMdd>_<sheet>/<sheet>.pureTarget.samples.summary.txt
 *   pt_cellInventory    QC/cellInventory/: every uBAM of each SMRT Cell
 *   ptRunQcReport       reports/<sheet>.pureTarget.runQC.{html,json}: both panels, every cell
 * Scripts: lrsCustomTools summaryScripts/ptcp/ (params.ptcpReportScripts).
 */


// Every sample sheet row: panel, testlist, analysed or no uBAM found, uBAM count and size
process pt_samplesSummary {
    label "low"

    publishDir (
        path: {"${params.outputDirBase}/runInfo/${params.dateStamp}_${params.ssBase}/"},
        mode: 'copy')

    input:
    val(text)

    output:
    path("${params.ssBase}.pureTarget.samples.summary.txt")

    script:
    def safe = text.toString().replace("'", "'\\''")
    """
    printf '%s\\n' '${safe}' > ${params.ssBase}.pureTarget.samples.summary.txt
    """
}

// movies: comma-separated movie names; runDirs: run folders to search
process pt_cellInventory {
    label "low"
    conda "${params.condaEnvs.ptcp}"     // a known python3 (the script is stdlib only)

    publishDir (
        path: {"${params.runRoot}/QC/cellInventory/"},
        mode: 'copy')

    input:
    tuple val(movies), val(runDirs)

    output:
    path("${params.ssBase}.cellInventory.json"), emit: json
    path("${params.ssBase}.cellInventory.tsv"),  emit: tsv

    script:
    """
    python3 ${params.ptcpReportScripts}/pt_cell_inventory.py \
    --movies ${movies} \
    --runs ${runDirs.collect { "'${it}'" }.join(' ')} \
    --out ${params.ssBase}.cellInventory.json \
    --tsv ${params.ssBase}.cellInventory.tsv
    """
}

// QC only, no results. Either panel's QC may be empty (no samples of that panel).
process ptRunQcReport {
    label "low"
    conda "${params.condaEnvs.ptcp}"

    publishDir (
        path: {"${params.runRoot}/reports/"},
        mode: 'copy',
        saveAs: { fn -> file(fn).name })

    input:
    path(carrierQc, stageAs: 'carrier_qc/*')     // qc.<npn>.json + qc.aggregate.json, or none
    path(repeatQc,  stageAs: 'repeat_qc/*')      // idem
    path(inventory, stageAs: 'inventory/*')      // <sheet>.cellInventory.json
    val(sampleRows)                              // [[npn, panel, sex M|K, testlist, material, cells, runs, barcodes], ...]
    val(runMeta)

    output:
    path("reports/*.html"), emit: html
    path("reports/*.json"), emit: json

    script:
    def tsv      = (["id\tpanel\tsex\ttestlist\tmaterial\tcells\truns\tbarcodes"] + sampleRows.collect { it.collect { v -> v ?: '' }.join('\t') }).join('\n')
    def metaJson = groovy.json.JsonOutput.toJson(runMeta).replace("'", "'\\''")
    def thr      = params.ptcpReportThresholds ? "--thresholds ${params.ptcpReportThresholds}" : ""
    def ref      = params.ptcpReportReference  ? "--reference ${params.ptcpReportReference}" : ""
    def scripts  = params.ptcpReportScripts
    """
    printf '%s\\n' '${tsv}' > samples.tsv
    printf '%s\\n' '${metaJson}' > run_meta.json

    python3 ${scripts}/pt_run_qc_collector.py \
    --carrier-qc-dir carrier_qc \
    --repeat-qc-dir repeat_qc \
    --targets-bed ${params.ptcpQcBed} \
    --samples samples.tsv \
    --cell-inventory inventory/${params.ssBase}.cellInventory.json \
    --run-meta run_meta.json \
    ${thr} ${ref} \
    --template ${scripts}/pt_run_qc_report_template.html \
    --outdir reports
    """
}
