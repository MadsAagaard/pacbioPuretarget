#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

import java.util.Locale

/*
 * =============================================================================
 * KG Vejle — joint PureTarget pipeline: repeat expansion panel + carrier panel
 * =============================================================================
 * One sample sheet, one run. Every sample sheet row gets its panel from its
 * testlist (conf/pureTarget.config) and runs that panel's analysis:
 *   repeat   pureTarget v1 analysis (subworkflows/REPEAT.nf)
 *   carrier  PTCP v1 analysis       (subworkflows/PREPROCESS.nf + PTCP.nf)
 * Both get ptcp-qc; the run gets the cell inventory and one run QC report.
 * Single-sample analysis only (intRef etc. in LabWare sheets are not used).
 * =============================================================================
 */

log.info """\
======================================================
Clinical Genetics Vejle: PacBio PureTarget (repeat + carrier panels)
======================================================
Samplesheet   : ${params.samplesheet}
Read mode     : ${params.hifiReads ? 'HiFi only' : 'HiFi + fail'}
uBAM search   : ${params.input ?: params.dataArchive}
Run output    : ${params.runRoot}
lrsCustomTools: ${params.lrsCustomToolsBase}
Reports       : ${params.reports}
"""

if (!params.samplesheet) {
    exit 1, "USER INPUT ERROR: --samplesheet is required (SMRT Link sample sheet or LabWare metadata)."
}


///////////////////////////////////////////////////
/////// ------- PANEL FROM TESTLIST ------- ///////
///////////////////////////////////////////////////
// Leading 'SL-' dropped, '_' read as '-', case ignored (LabWare SL-PT-KENDT =
// SMRT Link / uBAM PT-KENDT). null = not a PureTarget testlist.
def normTestlist = { String tl -> (tl ?: '').trim().replaceFirst(/^(?i)SL-/, '').replace('_', '-').toUpperCase() }
def carrierExtra = params.ptCarrierTestlists.collect { normTestlist(it) }
def panelOf = { String tl ->
    def t = normTestlist(tl)
    if (t.startsWith(params.ptCarrierPrefix.toUpperCase()) || t in carrierExtra) return 'carrier'
    if (t.startsWith(params.ptRepeatPrefix.toUpperCase()))                        return 'repeat'
    return null
}


///////////////////////////////////////////////////
/////// ------- SAMPLESHEET PARSING ------- ///////
///////////////////////////////////////////////////
//   SMRT Link CSV: every line after 'Bio Sample Name,Plate Well,Adapter,Adapter2';
//   column 1 is <npn>_<material>_<testlist>_<gender>.
//   LabWare metadata: first tab column rekv_npn_material_testlist_gender_proband_intRef[_expectedCount],
//   anchored from the right on the proband marker (T/F), as pacbio.ingestMetadata.sh reads it.
def SS_HEADER = 'Bio Sample Name,Plate Well,Adapter,Adapter2'

def sexFromGender = { gender, sampleId ->
    def g = (gender ?: '').toString().trim().toLowerCase()
    if (g in ['m', 'male'])        return 'male'
    if (g in ['k', 'f', 'female']) return 'female'
    throw new IllegalArgumentException("Sample ${sampleId}: gender field is '${gender}', expected M or K.")
}

// <npn>_<material>_<testlist>_<gender>; testlist may contain '_'. null if not that shape.
def parseSampleName = { String name ->
    def t = name.split('_', -1) as List
    if (t.size() < 4) return null
    def r = [ npn: t[0], material: t[1], testlist: t[2..-2].join('_'), gender: t[-1] ]
    (r.npn && r.testlist) ? r : null
}

def parseLabwareLine = { String line ->
    def f = line.split('\t')[0].trim().tokenize('_')
    def n = f.size()
    def g = (n >= 7 && f[n - 2] in ['T', 'F']) ? n - 3          // ... gender proband intRef
          : (n >= 8 && f[n - 3] in ['T', 'F']) ? n - 4          // ... gender proband intRef expectedCount
          : null
    (g == null || g < 4) ? null : [ npn: f[1], material: f[2], testlist: f[3..g - 1].join('_'), gender: f[g] ]
}

def ssLines   = file(params.samplesheet, checkIfExists: true).readLines()*.trim().findAll { it }
def headerIdx = ssLines.findIndexOf { it == SS_HEADER }
def ssFormat  = headerIdx >= 0 ? 'SMRT Link CSV' : 'LabWare metadata'
def rawRows   = headerIdx >= 0
    ? ssLines.drop(headerIdx + 1).collect { [it, parseSampleName(it.split(/\s*,\s*/, -1)[0])] }
    : ssLines.findAll { !it.startsWith('#') }.collect { [it, parseLabwareLine(it)] }

def badLines = rawRows.findAll { it[1] == null }*.get(0)
if (badLines) {
    exit 1, "USER INPUT ERROR: sample line(s) in ${params.samplesheet} could not be read (${ssFormat}):\n    ${badLines.take(5).join('\n    ')}"
}

// One analysis per NPN + panel. File names as in v1 / PTCP inside the panel's
// subtree: repeat <npn>.<testlist>, carrier <npn>.
def ssAll = rawRows*.get(1).collect { r ->
    def panel = panelOf(r.testlist)
    def tlUbam = r.testlist.replaceFirst(/^(?i)SL-/, '')
    r + [ panel: panel,
          key  : "${r.npn}.${panel}".toString(),
          sex  : panel ? sexFromGender(r.gender, r.npn) : null,
          id   : panel == 'repeat' ? "${r.npn}.${tlUbam}".toString() : r.npn ]
}
def notPT = ssAll.findAll { !it.panel }
if (notPT) log.warn "Skipped, not a PureTarget testlist: ${notPT.collect { "${it.npn} (${it.testlist})" }.join(', ')}"

def ssRows = ssAll.findAll { it.panel }.unique { it.key }
def dupes  = ssAll.findAll { it.panel }.size() - ssRows.size()
if (dupes) log.warn "${dupes} sample sheet row(s) repeat an NPN already in the same panel: analysed once, with all its uBAMs."
log.info "Samplesheet (${ssFormat}): ${ssRows.count { it.panel == 'carrier' }} carrier, ${ssRows.count { it.panel == 'repeat' }} repeat"


///////////////////////////////////////////////////
/////// ------- INPUT DISCOVERY ------- ///////////
///////////////////////////////////////////////////
// uBAM name: <npn>_<material>_<testlist>_<gender>.<movie>.<readset>.<barcode>.bam
// in <run>/<hifi_reads|failed_reads>/ (archive) or anywhere under --input.
// A sample takes every uBAM with its NPN and a testlist of its panel, from any
// run (a resequenced sample gets all its data).
def root     = params.input ?: params.dataArchive
def inputBam = params.hifiReads
    ? (params.input ? "${root}/**/*.hifi_reads.*.bam" : "${root}/**/hifi_reads/*.hifi_reads.*.bam")
    : (params.input ? "${root}/**/*.bam"              : "${root}/**/*_reads/*.bam")

channel.fromPath(inputBam, followLinks: true)
    .map { bam -> tuple(parseSampleName(bam.name.tokenize('.')[0]), bam) }
    .filter { sn, bam -> sn && panelOf(sn.testlist) }
    .map { sn, bam -> tuple("${sn.npn}.${panelOf(sn.testlist)}".toString(), bam) }
    .set { ubam_by_key }

channel.fromList(ssRows)
    .map { row -> tuple(row.key, row) }
    .combine(ubam_by_key, by: 0)
    .map { key, meta, bam -> tuple(meta, bam) }
    .groupTuple(sort: true)
    .set { finalUbamInput }

// uBAM name tokens: [sample, movie, readset, barcode]
def cellsOf    = { bams -> bams.collect { it.name.tokenize('.')[1] }.unique().sort() }
def runsOf     = { bams -> bams.collect { it.parent.parent.name }.unique().sort() }
def barcodesOf = { bams -> bams.collect { it.name.tokenize('.').size() > 3 ? it.name.tokenize('.')[3] : '' }.findAll { it }.unique().sort() }

// One list of every analysed sample: [npn, panel, sex M|K, testlist, material, cells, runs, barcodes, nBams, GB]
finalUbamInput
    .map { meta, bams -> [ meta.npn, meta.panel, meta.sex == 'male' ? 'M' : 'K', meta.testlist, meta.material,
                           cellsOf(bams).join(','), runsOf(bams).join(','), barcodesOf(bams).join(','),
                           bams.size(), String.format(Locale.US, "%.2f", (bams.sum { it.size() } as long) / 1024.0 ** 3) ] }
    .collect(flat: false)
    .ifEmpty([])
    .set { samples_ch }

// sheet rows without uBAMs
def missingOf = { rows -> def found = rows.collect { "${it[0]}.${it[1]}".toString() }; ssRows.findAll { !(it.key in found) } }

// SMRT Cells and run folders of the analysed uBAMs, for the cell inventory
finalUbamInput
    .flatMap { meta, bams -> bams.collect { [ it.name.tokenize('.')[1], it.parent.parent.toString() ] } }
    .unique()
    .collect(flat: false)
    .map { pairs -> tuple(pairs*.get(0).unique().sort().join(','), pairs*.get(1).unique().sort()) }
    .set { cell_inventory_in_ch }

finalUbamInput
    .branch { meta, bams ->
        carrier : meta.panel == 'carrier'
        repeat  : meta.panel == 'repeat'
    }
    .set { byPanel }


///////////////////////////////////////////////////
/////// ------- MODULES / SUBWORKFLOWS ------- ////
///////////////////////////////////////////////////

include { PREPROCESS }  from './subworkflows/PREPROCESS.nf'     // carrier: PTCP v1
include { PTCP }        from './subworkflows/PTCP.nf'
include { ptcpReports } from './modules/ptcpModules.nf'
include { REPEAT }      from './subworkflows/REPEAT.nf'         // repeat: pureTarget v1
include { pt_samplesSummary; pt_cellInventory; ptRunQcReport } from './modules/runModules.nf'


///////////////////////////////////////////////////
/////// ------- MAIN WORKFLOW ------- /////////////
///////////////////////////////////////////////////

workflow {

    // ---- runInfo: every sheet row, analysed or not -----------------------------------
    samples_ch
        .map { rows ->
            def miss = missingOf(rows)
            if (miss) log.warn "Sample sheet rows with no matching uBAM: ${miss.collect { "${it.npn} (${it.panel}, ${it.testlist})" }.join(', ')}"
            ( ["sample\tpanel\ttestlist\tstatus\tubams\ttotal_gb\tcells"] +
              rows.collect { "${it[0]}\t${it[1]}\t${it[3]}\tanalysed\t${it[8]}\t${it[9]}\t${it[5]}" } +
              miss.collect { "${it.npn}\t${it.panel}\t${it.testlist}\tno uBAM found\t0\t0.00\t" } +
              notPT.collect { "${it.npn}\t-\t${it.testlist}\tskipped: not a PureTarget testlist\t0\t0.00\t" } ).join('\n')
        }
        | pt_samplesSummary

    // ---- carrier panel: PTCP v1 --------------------------------------------------------
    PREPROCESS(byPanel.carrier)
    PTCP(PREPROCESS.out.alignedFinal)

    // ---- repeat expansion panel: pureTarget v1 ----------------------------------------
    REPEAT(byPanel.repeat)

    // ---- reports ----------------------------------------------------------------------
    if (params.reports) {

        // carrier sample reports + carrier run QC (carrierPanel/reports/), as PTCP v1
        samples_ch
            .map { rows -> rows.findAll { it[1] == 'carrier' }.collect { [it[0], it[2], it[3], it[5], it[6]] } }   // id, sex, testlist, cells, runs
            .set { carrier_rows_ch }

        samples_ch
            .map { rows -> rows.findAll { it[1] == 'carrier' } }
            .map { rows ->
                [ samplesheet           : params.ssBase,
                  samplesheet_format    : ssFormat,
                  samples_in_sheet      : ssRows.count { it.panel == 'carrier' },
                  runfolder             : rows.collect { it[6].tokenize(',') }.flatten().unique().sort().join(', '),
                  run_output            : params.panelOut.carrier,
                  read_mode             : params.hifiReads ? 'HiFi only' : 'HiFi + fail',
                  ptcp_version          : params.ptcpVersion,
                  genome_version_tag    : params.genomeVersion,
                  tag_hifi              : params.tagHifi,
                  tag_str               : params.strTag,
                  link_prefix           : '../',
                  lrscustomtools_release: file(params.lrsCustomToolsBase).name,
                  min_gb                : 0,
                  tools                 : params.ptcpToolVersions,
                  smn_model             : params.ptcpSmnRegression ? file(params.ptcpSmnRegression).name : null,
                  missing               : missingOf(rows).findAll { it.panel == 'carrier' }.collect { [id: it.npn, testlist: it.testlist] },
                  dropped               : [] ]
            }
            .set { carrier_meta_ch }

        ptcpReports(
            PTCP.out.qcJson,
            PTCP.out.havanno.map { meta, f -> f }.collect().ifEmpty([]),
            carrier_rows_ch,
            carrier_meta_ch)

        // everything on each SMRT Cell, and the run QC report for both panels (reports/)
        pt_cellInventory(cell_inventory_in_ch)

        samples_ch
            .map { rows ->
                [ samplesheet              : params.ssBase,
                  samplesheet_format       : ssFormat,
                  samples_in_sheet_by_panel: [carrier: ssRows.count { it.panel == 'carrier' }, repeat: ssRows.count { it.panel == 'repeat' }],
                  runfolder                : rows.collect { it[6].tokenize(',') }.flatten().unique().sort().join(', '),
                  run_output               : params.runRoot,
                  read_mode                : params.hifiReads ? 'HiFi only' : 'HiFi + fail',
                  ptcp_version             : params.ptcpVersion,
                  lrscustomtools_release   : file(params.lrsCustomToolsBase).name,
                  panel_rules              : [carrier_prefix: params.ptCarrierPrefix, carrier_testlists: params.ptCarrierTestlists,
                                              repeat_prefix: params.ptRepeatPrefix],
                  panel_dirs               : params.panelDirs,
                  tools                    : params.ptcpToolVersions,
                  tools_repeat             : params.repeatToolVersions,
                  missing                  : missingOf(rows).collect { [id: it.npn, panel: it.panel, testlist: it.testlist] } ]
            }
            .set { run_meta_ch }

        // a panel without samples gives [] (its ptcp-qc never runs)
        ptRunQcReport(
            PTCP.out.qcJson.ifEmpty([]),
            REPEAT.out.qcJson.ifEmpty([]),
            pt_cellInventory.out.json,
            samples_ch.map { rows -> rows.collect { it.take(8) } },
            run_meta_ch)
    }
}


workflow.onComplete {
    def failed = (workflow.stats.failedCount ?: 0) + (workflow.stats.ignoredCount ?: 0)
    if (failed > 0) {
        def f = file("${params.outputDirBase}/runInfo/${params.dateStamp}_${params.ssBase}/FAILED_TASKS.pureTarget.txt")
        f.parent.mkdirs()
        f.text = """\
            Run          : ${workflow.runName}
            Failed       : ${workflow.stats.failedCount}
            Ignored      : ${workflow.stats.ignoredCount}
            Error mode   : ${params.errorMode}
            Run output   : ${params.runRoot}
            Per-task detail: trace file under NextflowReports/. Output of the
            affected samples is INCOMPLETE and must not be reported without review.
            """.stripIndent()
        log.warn "*** ${failed} task(s) failed or were ignored — see ${f} ***"
    }
}
