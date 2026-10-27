# KG Vejle PacBio PureTarget — repeat expansion panel + carrier panel (joint pipeline)

One Nextflow (DSL2) pipeline for a PureTarget run that holds both repeat
expansion panel and carrier panel (PTCP) samples, usually on the same SMRT Cells.

- One sample sheet, one run, single-sample analysis. Each sample gets its panel from its testlist:
  carrier `PTCP*` (and the implementation testlist `PT-KENDT`), repeat `PT-*`
  (`conf/pureTarget.config`). Other testlists are skipped. A sample takes every uBAM with its
  NPN and a testlist of its panel, from any run (resequenced samples get all their data).
- **Repeat panel**: the pureTarget v1 analysis, unchanged (pbmm2 on HiFi + fail,
  TRGT 4 + 5 with the STRchive KGVejle catalog, plots, methylation), plus a sample
  report per repeat sample (TRGT 5 genotype per locus with TRGT 4 alongside, collapsible
  locus cards, links to the plots).
- **Carrier panel**: the PTCP v1 analysis, unchanged (PacBio PTCP 3.2.1 commands,
  KG Vejle read setup; Paraphase, F8 inversion, Sawfish HBA, TRGT 5.1, ptcp-qc,
  carrier sample reports).
- **QC for both panels**: ptcp-qc 1.0.0 per panel, the cell inventory (everything on
  each SMRT Cell) and one run QC report for the whole run.

Design decisions and status: `PLAN.md`.

## Run

```
nextflow run <repo> -r <tag> -profile slurm --samplesheet <sheet> [--input <run folder>] [--outdir <test base>]
```

Sample sheet: the SMRT Link CSV used to sequence the run, or LabWare metadata
(`rekv_npn_material_testlist_gender_proband_intRef[_expectedCount]`). Gender M or K.

Options:

```
--samplesheet  [path]  SMRT Link CSV or LabWare metadata. Required
--input        [path]  folder to search for uBAMs. Default: the PureTarget archive
--hifiReads            HiFi reads only (default: HiFi + fail)
--outdir       [path]  output base (<outdir>_<readSet>)
--reports      [bool]  reports, cell inventory, run QC report (default true)
--ptcpSmnRegression [path]   SMN homology model (optional)
--ptcpReportThresholds [path] approved QC thresholds (optional; placeholders without)
--ptcpReportReference  [path] carrier target-share reference (optional)
--lrsCustomToolsBase   [path] lrsCustomTools checkout (default: releases/current)
--lrsSharedDir         [path] shared config version (default /data/configs/lrsShared/v1)
--errorMode    [str]   strict (default) | cohort
-profile slurm         KG Vejle SLURM cluster
```

## Output

```
<root>/pureTargetJoint/<sheet>_<readSet>/
    repeatExpansion/     repeat panel, pureTarget v1 layout (alignments/, TRGT/, TRGT5/, methylation/, ...)
                         + QC/ptcpQC/ (ptcp-qc for the repeat samples, on ptcp-qc.repeatCatalog.bed: a target per
                           TRGT catalog locus, built each run from params.trDiseaseCatalog)
                         + reports/<npn>.repeat.sampleReport.html (links relative, so the run folder can move)
    carrierPanel/        carrier panel, PTCP v1 layout (alignments/, repeatExpansions/, specialAnalysis/,
                         structuralVariants/, QC/ptcpQC/, reports/ = sample reports + carrier run QC)
    QC/cellInventory/    every uBAM of each SMRT Cell: barcode, testlist, reads, bases
    reports/             <sheet>.pureTarget.runQC.html: run QC for both panels and every cell
<root>/pureTargetJoint/runInfo/<yyMMdd>_<sheet>/   <sheet>.pureTarget.samples.summary.txt (every sheet row:
                         analysed / no uBAM found), NextflowReports, FAILED_TASKS.pureTarget.txt
```

`<root>` = `/lnx01_data2/shared/patients/pacbioLRS/<genome>/<year>`.

## Dependencies

- Shared config layer `lrsShared/v1` (`lrsShared.config`, `finalize.config`).
- lrsCustomTools release: `ptcp/pacbio/` (PacBio PTCP scripts), `resources/ptcp/v3.2.1/`,
  `summaryScripts/ptcp/` (report collectors and templates, `pt_cell_inventory.py`).
  Nothing is checked at launch: a missing script or resource fails the task that uses it
  (`--reports false` runs the analysis without the report scripts).
- Conda envs: lrsShared `condaEnvs` (pbmm2, trgt4, trgt5, trgt51, paraphase35, pbCPGtools)
  plus `ptcp` (`env/build_ptcp_env.sh`) and `sawfish_221`.
- samtools and bcftools on PATH.
