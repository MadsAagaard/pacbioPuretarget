#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * =============================================================================
 * pureTarget v2 — repeat expansion panel branch (modules)
 * =============================================================================
 * Copied verbatim from MadsAagaard/pureTarget dev (commit 971d6e7,
 * modules/dnaModules.nf): the nine processes the v1 workflow runs. Commands,
 * envs, file names and publish subfolders are unchanged (PLAN.md D4); every
 * publishDir goes through params.outBase(meta), which now returns the run's
 * repeatExpansion/ subtree. Changes, all marked "KG v2":
 *   - trgt5_pureTarget_plots: the '|| echo WARN' lines started on a new line
 *     after '-o ...pdf', a bash syntax error that failed every task; moved
 *     onto the '-o' line.
 *   - trgt4_pureTarget_plots: '|| echo WARN' added per plot, so one failing
 *     locus does not fail the task now that processes run under bash -eo
 *     pipefail (v1 on slurm ran bash -l without -e: same tolerance).
 * Not copied (not run by v1): pbmm2_align, methylationSegm, qualimap,
 * mosdepthROI, whatsHap_stats, cramino, nanoStat, multiQC, multiQC_ALL,
 * vntyper2, advntr.
 * =============================================================================
 */


process create_fofn {
    publishDir path: {"${params.outBase(meta)}/documents/"}, mode: 'copy',pattern: '*.fofn'  // KG v2: named path, Nextflow >= 25.10 does not resolve meta in a positional closure

    input:
    tuple val(meta), path(data) //ubam

    output:
    tuple val(meta), path("${meta.id}.fofn")
    script:
    """
    `realpath ${data} > ${meta.id}.fofn`
    """
}

process inputFiles_symlinks_ubam{
    label "low"
    publishDir path: {"${params.outBase(meta)}/inputSymlinks/"}, mode: 'symlink', pattern: '*.{bam,pbi}'  // KG v2: named path, Nextflow >= 25.10 does not resolve meta in a positional closure

    input:
    tuple val(meta), path(data)   
    //data (default): 0:ubam, 1:ubam pbi

    output:
    tuple val(meta), path(data)

    script:
    """
    """
}

process pbmm2_align_mergedData {
    label "medium"
    tag "$meta.id"
    conda "${params.condaEnvs.pbmm2}"

    publishDir path: {"${params.outBase(meta)}/alignments/"}, mode: 'copy', pattern: '*.{bam,bai}'  // KG v2: named path, Nextflow >= 25.10 does not resolve meta in a positional closure

    input:
    tuple val(meta), path(fofn)
    
    output:
    tuple val(meta), path("${meta.id}.${params.genomeVersion}.${params.readSet}.bam"), path("${meta.id}.${params.genomeVersion}.${params.readSet}*bai"),  emit: bam
    
    script:
    """
    pbmm2 align \
    --preset HIFI \
    --sort \
    --num-threads ${task.cpus} \
    --bam-index BAI \
    --sample ${meta.id} \
    ${params.genomeMmi} \
    ${fofn} \
    ${meta.id}.${params.genomeVersion}.${params.readSet}.bam
    """
}

process trgt4_pureTarget{
    tag "$meta.id"
    label "medium"
    conda "${params.condaEnvs.trgt4}"

    publishDir path: {"${params.outBase(meta)}/TRGT/bam/"}, mode: 'copy', pattern: "*.trgt4.ba*"  // KG v2: named path, Nextflow >= 25.10 does not resolve meta in a positional closure
    publishDir path: {"${params.outBase(meta)}/TRGT/vcf/"}, mode: 'copy', pattern: "*.trgt4.vcf.*"  // KG v2: named path, Nextflow >= 25.10 does not resolve meta in a positional closure

/*
    publishDir {params.groupedOutput ? "${outputDir}/TRGT/pureTarget/bam" : "${outputDir}/${meta.id}/TRGT/pureTarget/bam/"}, mode: 'copy', pattern: "*.sorted.ba*"

    publishDir {params.groupedOutput ? "${outputDir}/TRGT/pureTarget/perSampleVcf/" : "${outputDir}/${meta.id}/TRGT/pureTarget/perSampleVcf/"}, mode: 'copy', pattern: "*.sorted.vcf.*"
*/


    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path("${meta.id}.${params.genomeVersion}.${params.readSet}.trgt4.bam"), path("${meta.id}.${params.genomeVersion}.${params.readSet}.trgt4.bam.bai"), path("${meta.id}.${params.genomeVersion}.${params.readSet}.trgt4.vcf.gz"), path("${meta.id}.${params.genomeVersion}.${params.readSet}.trgt4.vcf.gz.tbi"),emit: trgt_full
    tuple val(meta),path ("*.trgt4.*")
    path("${meta.id}.${params.genomeVersion}.${params.readSet}.trgt4.vcf.gz"),emit: trgt_vcf
    script:

    //def sex=meta.sex=="male" ? "--karyotype XY" : "--karyotype XX"
    def karyotype=(meta.sex=="male"||meta.gender=="M") ? "--karyotype XY" : "--karyotype XX"
    """
    trgt genotype \
    --preset targeted \
    --genome ${params.genomeFasta} \
    --repeats ${params.trDiseaseCatalog} \
    $karyotype \
    --reads ${data[0]} \
    --output-prefix ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt4.STRchive

    bcftools sort -Oz -o ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt4.vcf.gz ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt4.STRchive.vcf.gz 
    bcftools index -t ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt4.vcf.gz

    samtools sort -o ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt4.bam ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt4.STRchive.spanning.bam
    samtools index ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt4.bam
    """
}

process trgt5_pureTarget{
    tag "$meta.id"
    label "medium"
    conda "${params.condaEnvs.trgt5}"

    publishDir path: {"${params.outBase(meta)}/TRGT5/bam/"}, mode: 'copy', pattern: "*.trgt5.ba*"  // KG v2: named path, Nextflow >= 25.10 does not resolve meta in a positional closure
    publishDir path: {"${params.outBase(meta)}/TRGT5/vcf/"}, mode: 'copy', pattern: "*.trgt5.vcf.*"  // KG v2: named path, Nextflow >= 25.10 does not resolve meta in a positional closure


    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path("${meta.id}.${params.genomeVersion}.${params.readSet}.trgt5.bam"), path("${meta.id}.${params.genomeVersion}.${params.readSet}.trgt5.bam.bai"), path("${meta.id}.${params.genomeVersion}.${params.readSet}.trgt5.vcf.gz"), path("${meta.id}.${params.genomeVersion}.${params.readSet}.trgt5.vcf.gz.tbi"),emit: trgt_full
    tuple val(meta),path ("*.trgt5.*")
    path("${meta.id}.${params.genomeVersion}.${params.readSet}.trgt5.vcf.gz"),emit: trgt_vcf
    script:

    //def sex=meta.sex=="male" ? "--karyotype XY" : "--karyotype XX"
    def karyotype=(meta.sex=="male"||meta.gender=="M") ? "--karyotype XY" : "--karyotype XX"
    """
    trgt genotype \
    --preset targeted \
    --genome ${params.genomeFasta} \
    --repeats ${params.trDiseaseCatalog} \
    $karyotype \
    --reads ${data[0]} \
    --output-prefix ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt5.STRchive

    bcftools sort -Oz -o ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt5.vcf.gz ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt5.STRchive.vcf.gz 
    bcftools index -t ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt5.vcf.gz

    samtools sort -o ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt5.bam ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt5.STRchive.spanning.bam
    samtools index ${meta.id}.${params.genomeVersion}.${params.readSet}.trgt5.bam
    """
}

process trgt4_pureTarget_plots{
    tag "$meta.id"
    label "medium"
    conda "${params.condaEnvs.trgt4}"
    publishDir path: {"${params.outBase(meta)}/TRGT/plots/${meta.id}/"}, mode: 'copy',  pattern: "*.{pdf,png,svg}"  // KG v2: named path, Nextflow >= 25.10 does not resolve meta in a positional closure

    input:
    tuple val(meta), val(data)
    
    output:
    tuple val(meta), path("*.{pdf,png,svg}")
    script:

    def geneList = params.puretargetPlotGenes.join(' ')
    """
    for gene in ${geneList}; do
    trgt plot \
    --genome ${params.genomeFasta} \
    --repeats ${params.trDiseaseCatalog} \
    --vcf ${data.vcf} \
    --spanning-reads ${data.bam} \
    --repeat-id \$gene \
    --squished \
    -o ${meta.id}.${params.genomeVersion}.${params.readSet}.\$gene.alleleSquished.pdf || echo "WARN: squished plot failed for \$gene" >&2  # KG v2: tolerant under bash -e, as v1 was on slurm

    trgt plot \
    --genome ${params.genomeFasta} \
    --repeats ${params.trDiseaseCatalog} \
    --vcf ${data.vcf} \
    --spanning-reads ${data.bam} \
    --repeat-id \$gene \
    --plot-type waterfall \
    -o ${meta.id}.${params.genomeVersion}.${params.readSet}.\$gene.waterfall.pdf || echo "WARN: waterfall plot failed for \$gene" >&2  # KG v2: tolerant under bash -e, as v1 was on slurm
    done
    """
}

process trgt4_pureTarget_plots_meth{
    tag "$meta.id"
    label "medium"
    conda "${params.condaEnvs.trgt4}"

    publishDir path: {"${params.outBase(meta)}/TRGT/methylationPlots/${meta.id}/"}, mode: 'copy',  pattern: "*.{pdf,png,svg}"  // KG v2: named path, Nextflow >= 25.10 does not resolve meta in a positional closure

    input:
    tuple val(meta), val(data)
    
    output:
    tuple val(meta), path("*.{pdf,png,svg}")
    script:

    """
    trgt plot \
    --genome ${params.genomeFasta} \
    --repeats ${params.trDiseaseCatalog} \
    --vcf ${data.vcf} \
    --spanning-reads ${data.bam} \
    --repeat-id FXS_FMR1 \
    --show meth \
    --squished \
    --max-allele-reads 75 \
    -o FXS_FMR1.${meta.id}.${params.genomeVersion}.${params.readSet}.METH.alleleSquished.pdf

    trgt plot \
    --genome ${params.genomeFasta} \
    --repeats ${params.trDiseaseCatalog} \
    --vcf ${data.vcf} \
    --spanning-reads ${data.bam} \
    --repeat-id FXS_FMR1 \
    --plot-type waterfall \
    --show meth \
    --max-allele-reads 75 \
    -o FXS_FMR1.${meta.id}.${params.genomeVersion}.${params.readSet}.METH.waterfall.pdf

    """
}

process trgt5_pureTarget_plots{
    tag "$meta.id"
    label "medium"
    conda "${params.condaEnvs.trgt5}"

    publishDir path: {"${params.outBase(meta)}/TRGT5/plots/${meta.id}/"}, mode: 'copy',  pattern: "*.{pdf,png,svg}"  // KG v2: named path, Nextflow >= 25.10 does not resolve meta in a positional closure


    input:
    tuple val(meta), val(data)
    
    output:
    tuple val(meta), path("*.{pdf,png,svg}")
    
    script:
    def geneList = params.puretargetPlotGenes.join(' ')
    """
    for gene in ${geneList}; do
    trgt plot \
    --genome ${params.genomeFasta} \
    --repeats ${params.trDiseaseCatalog} \
    --vcf ${data.vcf} \
    --spanning-reads ${data.bam} \
    --repeat-id \$gene \
    --squished \
    -o ${meta.id}.${params.genomeVersion}.${params.readSet}.\$gene.alleleSquished.pdf || echo "WARN: squished plot failed for \$gene" >&2  # KG v2: was on its own line (bash syntax error)

    trgt plot \
    --genome ${params.genomeFasta} \
    --repeats ${params.trDiseaseCatalog} \
    --vcf ${data.vcf} \
    --spanning-reads ${data.bam} \
    --repeat-id \$gene \
    --plot-type waterfall \
    -o ${meta.id}.${params.genomeVersion}.${params.readSet}.\$gene.waterfall.pdf || echo "WARN: waterfall plot failed for \$gene" >&2  # KG v2: was on its own line (bash syntax error)
    done
    """
}

process methylationBW{
    tag "$meta.id"
    errorStrategy 'ignore'
    cpus 20
    conda "${params.condaEnvs.pbCPGtools}"

    publishDir path: {"${params.outBase(meta)}/methylation/BigWigBed/"}, mode: 'copy',  pattern: "*.methylation.{hap1,hap2,combined}.*"  // KG v2: named path, Nextflow >= 25.10 does not resolve meta in a positional closure


    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path("${meta.id}.${params.genomeVersion}.${params.readSet}.methylation*")
    
    script:
    """
    aligned_bam_to_cpg_scores \
    --bam ${data[0]} \
    --output-prefix ${meta.id}.${params.genomeVersion}.${params.readSet}.methylation
    """
}

///////////////////////////////////////////////////
/////// ------- JOINT PIPELINE ADDITION ------- ///
///////////////////////////////////////////////////
// Not in pureTarget v1. ptcp-qc 1.0.0 for the repeat panel, same targets BED as
// the carrier panel, on the repeat branch's HiFi + fail alignment. Mapped BAMs
// only, --no-trgt --no-paraphase (the repeat TRGT VCFs use the STRchive catalog,
// whose TRIDs are not the BED's). Not --only-groups: in testing it changed
// total_reads. Sample name = NPN (one repeat analysis per NPN).
// Output: <run>/repeatExpansion/QC/ptcpQC/

process ptcpQC_repeat {
    label "intermediate"
    conda "${params.condaEnvs.ptcp}"

    publishDir (
        path: {"${params.panelOut.repeat}/QC/ptcpQC/"},
        mode: 'copy',
        saveAs: { fn -> file(fn).name })

    input:
    tuple val(records), path(files)      // records: [[npn, bamName, baiName], ...]

    output:
    path("ptcp_qc/qc.*.json"),                            emit: json
    path("ptcp_qc/qc.*.csv"),                             emit: csv
    path("${params.ssBase}.repeat.ptcpQC.manifest.json"), emit: manifest

    script:
    def manifestName = "${params.ssBase}.repeat.ptcpQC.manifest.json"
    // ptcp-qc reads the sample name from the file name (<sample>.mapped.bam)
    def links    = []
    def manifest = [ sample_names: [], mapped_bams: [] ]
    records.each { r ->
        def (npn, bam, bai) = r
        links << "ln -s ../${bam} qc_input/${npn}.mapped.bam"
        links << "ln -s ../${bai} qc_input/${npn}.mapped.bam.bai"
        manifest.sample_names << npn
        manifest.mapped_bams  << "qc_input/${npn}.mapped.bam".toString()
    }
    def manifestJson = groovy.json.JsonOutput.toJson(manifest)
    """
    mkdir -p ptcp_qc qc_input
    ${links.join('\n    ')}

    printf '%s\\n' '${manifestJson}' > ${manifestName}

    ptcp-qc -v analyze \
    --manifest ${manifestName} \
    --targets-bed ${params.ptcpQcBed} \
    --output-prefix ptcp_qc/qc \
    --threads ${task.cpus} \
    --no-trgt \
    --no-paraphase \
    --json

    python3 ${params.ptcpScripts}/ptcpqc_json_to_csv.py --inputs ptcp_qc/qc*.json \
    -o ptcp_qc/qc.csv
    """
}


// Joint-pipeline addition: one report per repeat sample (repeatExpansion/reports/),
// built by lrsCustomTools summaryScripts/ptcp/pt_repeat_report_collector.py.
// TRGT 5 genotype per locus (TRGT 4 alongside), the sample's own ptcp-qc QC, and
// links to its plots / VCFs / BAMs relative to reports/. Plot links come from the
// plot names the plot steps actually produced (a failed plot is not linked).
process ptRepeatReports {
    label "low"
    conda "${params.condaEnvs.ptcp}"

    publishDir (
        path: {"${params.panelOut.repeat}/reports/"},
        mode: 'copy',
        saveAs: { fn -> file(fn).name })

    input:
    val(sampleRows)                                  // [[npn, file_id, sex M|K, testlist, cells, runs], ...]
    path(trgt5Vcfs, stageAs: 'trgt5/*')
    path(trgt4Vcfs, stageAs: 'trgt4/*')
    path(qcJson,    stageAs: 'qc/*')                 // qc.<npn>.json + qc.aggregate.json
    val(plotRows)                                    // [[file_id, kind, file name], ...]
    val(runMeta)

    output:
    path("reports/*.html"), emit: html
    path("reports/*.json"), emit: json

    script:
    def tsv      = (["id\tfile_id\tsex\ttestlist\tcells\truns"] + sampleRows.collect { it.collect { v -> v ?: '' }.join('\t') }).join('\n')
    def plots    = (["file_id\tkind\tfile"] + plotRows.collect { it.join('\t') }).join('\n')
    def metaJson = groovy.json.JsonOutput.toJson(runMeta)
    def thr      = params.ptcpReportThresholds ? "--thresholds ${params.ptcpReportThresholds}" : ""
    def scripts  = params.ptcpReportScripts
    // quoted heredocs: no shell expansion, and no 128 kB limit on a single argument
    // (the plot list has ~150 names per sample)
    """
    cat > samples.tsv <<'PT_EOF'
${tsv}
PT_EOF
    cat > plots.tsv <<'PT_EOF'
${plots}
PT_EOF
    cat > run_meta.json <<'PT_EOF'
${metaJson}
PT_EOF

    python3 ${scripts}/pt_repeat_report_collector.py \
    --samples samples.tsv \
    --qc-dir qc \
    --trgt5-dir trgt5 \
    --trgt4-dir trgt4 \
    --plots plots.tsv \
    --targets-bed ${params.ptcpQcBed} \
    --run-meta run_meta.json \
    ${thr} \
    --template ${scripts}/pt_repeat_sample_report_template.html \
    --outdir reports
    """
}
