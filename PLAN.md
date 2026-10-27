# pureTarget v2: repeat panel + carrier panel in one pipeline (plan, 26 Oct 2026; updated 27 Oct)

## Goal
One Nextflow pipeline for a PureTarget SMRT Cell that holds both repeat expansion panel and carrier panel (PTCP) samples. One sample sheet, one run, one run QC report covering every sample on every cell. The analysis of each panel stays as validated:
- repeat: the pureTarget v1 analysis (MadsAagaard/pureTarget, dev, lrsShared v1);
- carrier: the PTCP v1 analysis (project docs `claude/pacbioPTCP_v1/`, validated on `261016_10ptcpSamples`).

Draft code: project docs `claude/pureTarget_v2/` until it is pushed (D1).

## Why
- PacBio (reply 26 Oct): pool at least 8 samples per cell. Our run QC report showed 5 or fewer per cell, because a report only knows its own sample sheet. Repeat and carrier samples share cells, so per-cell QC needs both panels and everything else on the cell.
- pureTarget v1 has no testlist matching (it joins uBAMs on NPN only and searches the whole PureTarget archive). On a mixed sheet it analyses the carrier samples too, and a patient on both panels gets both sets of uBAMs merged into one repeat analysis.
- PacBio's QC suggestions (cell metadata, barcodes, DNA source, panel, bases on target, share of cell, target percentages) all need the cell-level view.

## Decisions (27 Oct; simplified 27 Oct after review)
- D1 Home: a NEW repo for the joint pipeline; the working PTCP (pacbioPTCP_v1) and pureTarget (MadsAagaard/pureTarget) pipelines stay as they are until the joint one is validated. PTCP's code is the carrier branch, pureTarget v1's the repeat branch. Layout and steps: IMPLEMENTATION.md.
- D2 Panel from the testlist, by prefix (leading `SL-` dropped, `_` = `-`, case ignored): carrier = `PTCP*` (any suffix) or one of `params.ptCarrierTestlists` (now `PT-KENDT`, the implementation testlist); repeat = `PT-*` (any suffix). Any other testlist is not PureTarget: its sheet rows are skipped with a warning and its uBAMs are never picked up.
- D3 One analysis per NPN + panel, single-sample only (intRef etc. not used). A sample takes every uBAM with its NPN and a testlist of its panel, from any run: a resequenced sample gets all its data, also under another suffix. A second sheet row for the same NPN + panel is analysed once (warning). File names as before inside the panel's subtree (repeat `<npn>.<testlist>` with the sheet row's testlist, carrier `<npn>`).
- D4 Repeat analysis unchanged: pbmm2 on HiFi + fail together (`--preset HIFI --sort`), TRGT 4 and TRGT 5 with the STRchive KGVejle catalog, plots, methylation. Only change: processes run under `bash -eo pipefail` (failure detection, not results). Harmonising with the carrier read setup is a later, separately validated change.
- D4b Proposed, not done: `pbmm2 --unmapped` in the repeat alignment, so ptcp-qc can count unmapped reads (now "n/a").
- D5 Carrier analysis unchanged: PTCP v1 modules (ptcpModules.nf, lrsFunctions.nf) and PREPROCESS.nf byte-identical; PTCP.nf = PTCP v1 PTCP.nf + the PTCP_QC.nf wiring moved in.
- D6 Output: `<root>/pureTargetJoint/<sheet>_<readSet>/` with `repeatExpansion/` (v1 layout + `QC/ptcpQC/`), `carrierPanel/` (PTCP layout incl. `QC/ptcpQC/`, `reports/`), `QC/cellInventory/`, `reports/`. Root `/lnx01_data2/shared/patients/pacbioLRS/<genome>/<year>/pureTargetJoint/`, next to the production roots `pureTarget/` and `pureTargetCarrier/`. runInfo: one summary file (every sheet row: panel, testlist, analysed / no uBAM found, uBAM count, GB, cells) + Nextflow reports + FAILED_TASKS.
- D7 QC: ptcp-qc 1.0.0 for both panels, one run per panel, same BED (carrier targets + all 37 repeat panel 2.0 loci). Repeat: mapped BAMs only, `--no-trgt --no-paraphase`, not `--only-groups` (see the test below).
- D8 Cell inventory: every uBAM of every SMRT Cell the analysed samples came from, reads/bases from the `.pbi`. Script in lrsCustomTools `summaryScripts/ptcp/pt_cell_inventory.py`.
- D9 Supersedes lrsShared PLAN.md D14 (each pipeline takes only its own samples from a mixed run).
- D10 Change control: a new pipeline version for both panels; validated (phase 5) before clinical use.
- D11 errorMode: the lrsShared default `strict`; a failed task stops the run, rerun with `-resume`.
- D12 Minimal gates (27 Oct review): no launch checks of scripts or resources (a missing file fails its task), no size gate, no `--panels`. Remaining checks: `--samplesheet` given, sheet lines readable, gender M/K. `--reports false` skips all reports (e.g. before an lrsCustomTools release has the report scripts).

## ptcp-qc on repeat samples (tested 26 Oct)
ptcp-qc 1.0.0 (release binary, MD5 checked) on two synthetic repeat-panel BAMs (reads at the 37 TR loci, random off-target reads, some rq 0.95 reads, unmapped reads), manifest with `sample_names` and `mapped_bams` only:
- runs without `--no-paraphase` / `--no-trgt`; logs "Skipping TRGT / Paraphase / Sawfish analysis ... files not available" and writes `qc.<sample>.json` + `qc.aggregate.json`;
- gives stats (total, on/off-target, unmapped, each split HiFi/fail by rq), per-locus reads (39 groups with reads), sex ratios, aggregate locus_stats;
- the HiFi/fail split comes from the rq tag, so it works on pureTarget's HiFi + fail merged alignment;
- unmapped reads are only counted if the BAM has them (D4b).
Not yet tested: TRGT VCFs from the STRchive catalog in the manifest (TRIDs may not match the BED names). Not needed for QC; callability for repeat loci can come from the TRGT VCF directly.
27 Oct, further tests (same binary, synthetic data):
- `--only-groups <the 39 TR groups>` restricts the coverage percentages to the repeat loci ("Percent of Targets >=10-fold" 100% instead of 85%: the carrier targets counted as empty targets), and moves reads on carrier targets to off-target. But it also changed total_reads: of 100 reads placed on HBB and GBA1, 10 vanished from every count. Not explained, so not used; the repeat-only target metrics are computed in the run QC report instead (phase 4).
- PacBio's `ptcpqc_json_to_csv.py` works on repeat-only JSONs (paraphase/trgt CSVs header-only, coverage and aggregate filled).
- A read without the `np` tag makes ptcp-qc stop for that sample (error, no JSON). Real HiFi/fail reads always carry `np`; noted for the failure checks.
- The repeat QC task's shell (links, manifest, ptcp-qc, CSV) ran under `bash -eo pipefail` on two synthetic repeat samples: exit 0, all outputs.

## Phases
1. Front end + config (draft 26 Oct; simplified 27 Oct): both sheet formats (from PTCP), panel by testlist prefix, NPN + panel matching, one runInfo summary, cell inventory.
2. Branches (draft 27 Oct, this folder):
   - carrier: PTCP v1 files copied unchanged: `modules/ptcpModules.nf`, `modules/lrsFunctions.nf`, `subworkflows/PREPROCESS.nf`; `PTCP.nf` with PTCP_QC.nf's wiring moved in; plus the carrier reports (`ptcpReports`, carrier rows only, to `carrierPanel/reports/`). Transcribed from the project docs: check byte identity on the server (IMPLEMENTATION.md).
   - repeat: `modules/repeatModules.nf` = the nine processes v1 runs, extracted verbatim from dnaModules.nf (dev 971d6e7; diffed per process, identical apart from the marked lines); `subworkflows/REPEAT.nf` = the v1 workflow block. create_fofn / inputFiles_symlinks_ubam aliased as `repeat_*` (same names exist in the PTCP modules) and run on the head node. Two marked fixes: the trgt5 plot syntax error, and `|| echo` on the TRGT4 plot loop (tolerant under `-e`, as v1 was on slurm without `-e`).
   - repeat QC: `ptcpQC_repeat` appended to `modules/repeatModules.nf` (marked as a joint-pipeline addition), wired in `REPEAT.nf`, to `repeatExpansion/QC/ptcpQC/`.
   - uBAM glob: archive `**/*_reads/*.bam` (as PTCP) instead of v1's `**/*.bam`; same files for the archive layout (`hifi_reads/`, `failed_reads/`).
3. QC: ptcp-qc per panel; off-target metrics per sample (count, median length, AT/GC skew), excluding fixed fragments outside the BED (chr16:112,769-117,220; chr19:44,517,839-44,526,570 once confirmed); bases per target (mosdepth `--by` on the ptcp-qc BED).
4. Reports:
   - Done 27 Oct (draft, tested on synthetic data): run QC report for both panels, `<run>/reports/<sheet>.pureTarget.runQC.html` (process `ptRunQcReport`, modules/runModules.nf). Code in lrsCustomTools `summaryScripts/ptcp/`: `pt_run_qc_collector.py` (imports ptcp_report_collector.py, so a carrier sample is judged exactly as in its PTCP report), `src/pt_run_body.html` -> `pt_run_qc_report_template.html` (built by `src/build_templates.py`, now three templates), `tests/make_mock_pt_run.py` + `tests/test_pt_run_qc.py` (8 tests). Project copy: `claude/pacbioPTCP_v1/report/` (+ `examples/DEMO_mixed.pureTarget.runQC.html`). Not yet in lrsCustomTools: push / PR pending.
     Sections: 01 run integrity per panel (sheet vs analysed vs ptcp-qc, missing, sex, ptcp-qc versions, cells with < 8 barcodes, unassigned, repeat loci not captured); 02 SMRT Cells from the cell inventory (every barcode: this sheet by NPN, others by testlist only; HiFi reads/bases, unassigned, share of cell vs even); 03 sample QC both panels (verdict first, sex, reads, share of cell, on/off-target, unmapped, lowest target of its own panel, read Q); 04 carrier targets (reads + share % + callability); 05 repeat loci (39 TR groups of the BED, reads + share %, loci with run median < fail flagged and not judged).
     Thresholds: `qc` (carrier, as PTCP), `qc_repeat` (defaults = carrier limits), `cell` (min_barcodes 8, unassigned 10%, low share 0.5x even): placeholders until approved.
   - The carrier sample reports and the carrier-only run QC (`carrierPanel/reports/`) stay as PTCP v1; the carrier-only run QC can be dropped once the new one is accepted.
   - Still open: bases on target, off-target metrics (phase 3), PacBio's DNA source field (needs the material code table); repeat sample reports (need approved repeat classes per locus); a repeat X-dose reference (now shown as a raw ratio, no verdict).
5. Validation: 261016_10ptcpSamples (carrier output identical to the PTCP v1 run), one production repeat run (repeat VCFs identical to production pureTarget), one mixed cell.

## Server checks and test runs
See IMPLEMENTATION.md (identity checks against PTCP v1 and pureTarget dev, prerequisites, test runs).

## PacBio's report suggestions → source (status 27 Oct)
- Samples in this report out of all on the cell; total HiFi reads per cell; each sample's share: cell inventory (D8).
- Barcode: uBAM name (4th dot token); SMRT Link sheet also has well and adapters.
- DNA source: material code from the sheet (parsed); needs the lab's code table. Extraction method/date: not in either sheet format (LabWare export change).
- PT version / panel: panel from the testlist map; panel/kit version from config.
- Bases on target: mosdepth on the ptcp-qc BED (phase 3).
- Carrier target percentages: shown in the run QC report (section 04), and the same for the repeat loci (section 05). Done.
- Samples on the cell / HiFi reads per cell / share of cell / barcode / panel: run QC report sections 02-03. Done.
- Bases on target: phase 3. DNA source: material code shown; text needs the lab's code table.

## Open questions
1. Lab: final testlist names per panel before the first clinical run (placeholders in use).
2. Lab: material code → DNA source table; can extraction method/date be added to the LabWare export?
3. PacBio: chr16:112,769-117,220 fixed fragment (see below); full list of guide cut sites; Paraphase 3.5.0 vs 3.4.0; smn_homology training JSONs.

## Found in v1 (pureTarget dev)
- `trgt5_pureTarget_plots`: the `|| echo "WARN..."` lines start on a new line after `-o ...pdf`, a bash syntax error (`bash -n`: unexpected token `||`). Every TRGT5 plot task fails; with errorMode cohort it is only recorded in FAILED_TASKS. Fix: end the `-o` line with ` \` or move `||` onto it.
- uBAMs matched on NPN only, across the whole PureTarget archive (see "Why").
- Under the slurm profile v1 processes ran `bash -l` without `-e`: a failing command mid-script went unnoticed if a later one succeeded. v2 runs every process with `-eo pipefail`.
- SMRT Link CSV only (PTCP reads both formats).

## Off-target findings (22-26 Oct, for the QC design)
- Off-target reads are random genomic fragments; their amount tracks input DNA damage (113291162253: 497k off-target, median 2.5 kb, AT-biased; normal cell-mate 113623114462: 32k, median 3.7 kb, flat). PacBio agrees: damaged ends + incomplete dephosphorylation; suggests the protocol's DNA repair step.
- Fixed two-cut fragment chr16:112,769-117,220 (reads 4,442-4,453 bp) in 10/10 samples, 591-2,269 reads, about 3x the hba_4.2_3.7 window; outside the ptcp-qc BED, the Paraphase HBA extract region (168,656-177,541) and PTCP's Sawfish HBA window (126,123-209,491). Asked PacBio whether it is designed and what for.
- Weaker fixed fragment chr19:44,517,839-44,526,570 (about 8.73 kb, about 160 reads); chr11:61,798,785 (about 3.52 kb, 28 reads) looks like a guide off-target pair.

## Sandbox test (27 Oct, after the repo was created)
Nextflow 25.04.8 and 25.10.0 in a sandbox mirroring the server layout (lrsShared v1, an lrsCustomTools checkout, a mock archive with `.pbi` files). Analysis processes stubbed (tools not installed); bookkeeping, the cell inventory, both report collectors and the ptcp-qc process scripts ran for real (ptcp-qc itself faked). Mixed, carrier-only and repeat-only sheets: exit 0, the IMPLEMENTATION.md section 5 checks hold. Fixes from it:
- repeatModules.nf: the v1 `publishDir {closure}, ...` form cannot see `meta` from Nextflow 25.10.0 on (works up to 25.04.8); now `publishDir path: {closure}, ...`, marked "KG v2". pureTarget dev has the same form.
- main.nf: the samples summary also lists rows skipped as not PureTarget.
