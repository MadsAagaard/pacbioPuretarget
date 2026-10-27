#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * PTCP — carrier panel analyses (PTCP main.wdl scatter), germline read-set setup:
 *   TRGT        HiFi BAM + --fail-reads (genotype, plots)
 *   Paraphase   HiFi BAM (+havanno) -> F8 inversion
 *   Sawfish     HiFi BAM, HBA region
 *   ptcp-qc and the repeat-reads extract: HiFi + fail merged (allBam)
 * then ptcp-qc across all carrier samples of the run (PTCP common.wdl ptcp_qc;
 * was PTCP_QC.nf in PTCP v1, moved here unchanged).
 */

include {
        samtools_merge_allReads;
        trgt5_ptcp_repeatReads;
        trgt5_ptcp;
        trgt5_ptcp_plots;
        paraphase_ptcp;
        f8inversion_ptcp;
        sawFish2_ptcp_HBA;
        ptcpQC;
        } from "../modules/ptcpModules.nf"


workflow PTCP {

    take:
    aligned     // tuple(meta, [hifiBam, hifiBai, failBam, failBai])

    main:

    // -------------------------------------------------------------------------
    // HiFi + fail merged, for ptcp-qc and the repeat-reads extract.
    // A sample without a fail BAM (--hifiReads, or no fail_reads uBAMs) uses
    // its HiFi BAM as is.
    // -------------------------------------------------------------------------
    aligned
        .branch { meta, aln ->
            withFail : aln.failBam
            hifiOnly : true
        }
        .set { alignedBySet }

    samtools_merge_allReads(alignedBySet.withFail)

    alignedBySet.withFail
        .join(samtools_merge_allReads.out.bam)
        .map { meta, aln, allBam, allBai -> tuple(meta, aln + [allBam: allBam, allBai: allBai]) }
        .mix(
            alignedBySet.hifiOnly
                .map { meta, aln -> tuple(meta, aln + [allBam: aln.hifiBam, allBai: aln.hifiBai]) }
        )
        .set { alignedAll }

    // -------------------------------------------------------------------------
    // Repeat expansions  (sex-dependent: --karyotype)
    // -------------------------------------------------------------------------
    trgt5_ptcp_repeatReads(alignedAll)
    trgt5_ptcp(aligned)

    trgt5_ptcp.out.trgt_full
        .map { meta, vcf, spanBam, spanBai ->
            tuple(meta, [vcf: vcf, spanBam: spanBam, spanBai: spanBai])
        }
        .set { trgt_ch }

    // PTCP trgt.wdl runs four plot sets per sample.
    def trgtPlotModes = [
        ['motifs', 'allele'],
        ['meth',   'allele'],
        ['motifs', 'waterfall'],
        ['meth',   'waterfall'],
    ]
    trgt5_ptcp_plots(trgt_ch, trgtPlotModes)

    // -------------------------------------------------------------------------
    // Paraphase -> F8 inversion  (F8 is sex-dependent: --sex)
    // -------------------------------------------------------------------------
    paraphase_ptcp(aligned)

    paraphase_ptcp.out.paraphase
        .map { meta, bam, bai, json ->
            tuple(meta, [paraphaseBam: bam, paraphaseBai: bai, paraphaseJson: json])
        }
        .set { paraphase_ch }

    f8inversion_ptcp(paraphase_ch)

    // -------------------------------------------------------------------------
    // HBA structural variants
    // -------------------------------------------------------------------------
    sawFish2_ptcp_HBA(aligned)

    // -------------------------------------------------------------------------
    // Everything ptcp-qc needs, one record per sample
    // -------------------------------------------------------------------------
    alignedAll
        .join(trgt_ch)
        .join(paraphase_ch)
        .join(f8inversion_ptcp.out.f8)
        .join(sawFish2_ptcp_HBA.out.sv_vcf)
        .map { meta, aln, trgt, para, f8Vcf, f8Json, svVcf, svTbi ->
            tuple(meta, aln + trgt + para + [
                f8Vcf      : f8Vcf,
                f8Json     : f8Json,
                sawfishVcf : svVcf,
                sawfishTbi : svTbi
            ])
        }
        .set { qcInput_ch }

    // -------------------------------------------------------------------------
    // ptcp-qc, all carrier samples together. The manifest is built from the
    // joined per-sample channel, so sample order is identical in every list.
    // -------------------------------------------------------------------------
    qcInput_ch
        .map { meta, d ->
            // staged file names; ptcpQC links them under the names ptcp-qc expects
            def record = [
                id            : meta.id,
                mappedBam     : d.allBam.name,     // HiFi + fail (see samtools_merge_allReads)
                mappedBai     : d.allBai.name,
                trgtVcf       : d.vcf.name,
                trgtSpanBam   : d.spanBam.name,
                trgtSpanBai   : d.spanBai.name,
                paraphaseBam  : d.paraphaseBam.name,
                paraphaseBai  : d.paraphaseBai.name,
                paraphaseJson : d.paraphaseJson.name,
                f8Json        : d.f8Json.name,
                sawfishVcf    : d.sawfishVcf.name,
                sawfishTbi    : d.sawfishTbi.name,
            ]
            // indexes are staged next to their files so ptcp-qc can open them
            def files = [ d.allBam, d.allBai, d.vcf, d.spanBam, d.spanBai,
                          d.paraphaseBam, d.paraphaseBai, d.paraphaseJson,
                          d.f8Json, d.sawfishVcf, d.sawfishTbi ]
            tuple(record, files)
        }
        .collect(flat: false)                   // emits nothing if no sample survived
        .map { rows ->
            def sorted = rows.sort { it[0].id }
            tuple(sorted.collect { it[0] }, sorted.collect { it[1] }.flatten())
        }
        .set { ptcpQC_input_ch }

    ptcpQC(ptcpQC_input_ch)

    emit:
    qcInput      = qcInput_ch      // tuple(meta, [hifiBam, hifiBai, failBam, failBai, allBam, allBai, vcf, spanBam, spanBai, paraphaseBam, paraphaseBai, paraphaseJson, f8Vcf, f8Json, sawfishVcf, sawfishTbi])
    trgtVcf      = trgt_ch
    paraphaseDir = paraphase_ptcp.out.paraphase_dir
    havanno      = paraphase_ptcp.out.havanno
    qcJson       = ptcpQC.out.json     // qc.<npn>.json + qc.aggregate.json (one list)
}
