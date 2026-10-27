#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * REPEAT — PureTarget repeat expansion panel, the pureTarget v1 analysis
 * (MadsAagaard/pureTarget dev 971d6e7, main.nf workflow block), unchanged:
 *   uBAM symlinks + fofn -> pbmm2 on HiFi + fail together -> TRGT 4 and TRGT 5
 *   (STRchive KGVejle catalog) -> plots (TRGT4 allele/waterfall + FMR1
 *   methylation, TRGT5 allele/waterfall) -> methylation bigwigs,
 * then ptcp-qc across all repeat samples of the run (joint pipeline addition).
 * create_fofn and inputFiles_symlinks_ubam are aliased: the carrier branch
 * has processes with the same names.
 */

include {
        create_fofn              as repeat_create_fofn;
        inputFiles_symlinks_ubam as repeat_inputFiles_symlinks_ubam;
        pbmm2_align_mergedData;
        trgt4_pureTarget;
        trgt5_pureTarget;
        trgt4_pureTarget_plots;
        trgt4_pureTarget_plots_meth;
        trgt5_pureTarget_plots;
        methylationBW;
        ptcpQC_repeat;
        } from "../modules/repeatModules.nf"


workflow REPEAT {

    take:
    finalUbamInput          // tuple(meta, [uBAMs]), meta.panel == 'repeat'

    main:

    repeat_inputFiles_symlinks_ubam(finalUbamInput)
    repeat_create_fofn(finalUbamInput)
    pbmm2_align_mergedData(repeat_create_fofn.out)

    pbmm2_align_mergedData.out.bam
        .map { meta, bam, bai -> tuple(meta, [bam, bai]) }
        .set { alignedFinal }

    trgt4_pureTarget(alignedFinal)
    trgt5_pureTarget(alignedFinal)

    trgt4_pureTarget.out.trgt_full
        .map { meta, bam, bai, vcf, tbi -> tuple(meta, [bam: bam, bai: bai, vcf: vcf, tbi: tbi]) }
        .set { trgt4_plot_ch }

    trgt4_pureTarget_plots(trgt4_plot_ch)
    trgt4_pureTarget_plots_meth(trgt4_plot_ch)

    trgt5_pureTarget.out.trgt_full
        .map { meta, bam, bai, vcf, tbi -> tuple(meta, [bam: bam, bai: bai, vcf: vcf, tbi: tbi]) }
        .set { trgt5_plot_ch }

    trgt5_pureTarget_plots(trgt5_plot_ch)

    methylationBW(alignedFinal)

    // ptcp-qc, all repeat samples together (sample name = NPN)
    pbmm2_align_mergedData.out.bam
        .map { meta, bam, bai -> tuple([meta.npn, bam.name, bai.name], [bam, bai]) }
        .collect(flat: false)                   // emits nothing without repeat samples
        .map { rows ->
            def sorted = rows.sort { it[0][0] }
            tuple(sorted.collect { it[0] }, sorted.collect { it[1] }.flatten())
        }
        .set { qc_input_ch }

    ptcpQC_repeat(qc_input_ch)

    // names of the plots each plot step produced, for the sample reports' links
    // (joint pipeline addition): [file_id, kind, file name]
    trgt4_pureTarget_plots.out.map { meta, f -> [meta.id, 'trgt4', f] }
        .mix(trgt5_pureTarget_plots.out.map { meta, f -> [meta.id, 'trgt5', f] },
             trgt4_pureTarget_plots_meth.out.map { meta, f -> [meta.id, 'trgt4_meth', f] })
        .flatMap { id, kind, f -> (f instanceof List ? f : [f]).collect { p -> [id, kind, p.name] } }
        .set { plot_names_ch }

    emit:
    trgt5     = trgt5_pureTarget.out.trgt_full    // tuple(meta, spanBam, spanBai, vcf.gz, tbi)
    trgt4     = trgt4_pureTarget.out.trgt_full    // tuple(meta, spanBam, spanBai, vcf.gz, tbi)
    qcJson    = ptcpQC_repeat.out.json            // qc.<npn>.json + qc.aggregate.json (one list)
    plotNames = plot_names_ch                     // [file_id, kind, file name]
}
