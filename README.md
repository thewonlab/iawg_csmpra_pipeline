# iawg_csmpra_pipeline

A small R workflow for the CSMPRA **custom** reanalysis, adapted from `src/IAWG/A5.NANA_mpralm.Rmd`. This is separate from the smaller published summary-statistics dataset. It generates analysis tables, not IGVF submission metadata, and makes no network requests or uploads.

## Requirements

R with `mpra` and `limma`; `jsonlite` is needed only for the optional sequence-design mapping.

```r
install.packages("BiocManager")
BiocManager::install(c("mpra", "limma"))
install.packages("jsonlite")
```

The default inputs are IGVF Reporter Experiment files named
`Reporter_Experiment_Long_veh.tsv.gz` and `Reporter_Experiment_Long_stim.tsv.gz`.
Each row is an oligo/replicate pair. The script pivots `dna_counts` and
`rna_counts` into the model's wide matrices; it does not use normalized counts
or the supplied log2FoldChange. Required replicates are `rep1`–`rep11` and
`rep14`–`rep16`. Duplicate pairs and invalid counts are rejected. Missing rows
and blank count cells remain missing, never zero-filled. Include the Scrambled
negative controls as well as tested `_A1_`/`_A2_` oligos.

For legacy wide matrices, add `--input-format matrix`; filenames are
`ENIGMA.Veh.activity.csv` and `ENIGMA.Stim.activity.csv`, with unique oligo names
in the first column and `dna1`–`dna11`, `dna14`–`dna16` plus corresponding RNA
columns. Both input formats represent post-QC oligo counts, not raw barcode
counts. Data files are not included; the underlying public study is GEO
**GSE289858**.

## Run

From this directory, with inputs placed in `data/`:

```bash
Rscript run_analysis.R --counts-dir data --output-dir results --input-format reporter
```

Reporter Experiment input and both conditions are the defaults. The same count files can support Published and Custom analyses; this script produces only the Custom results.

Both conditions run by default. Use `--condition veh` or `--condition stim` for one condition. Paths can point outside this directory. Existing condition output directories are never overwritten; choose a new output directory for a new run.

To add genomic ALT/REF effects, supply a **verified** IGVF sequence-design TSV (plain or gzip):

```bash
Rscript run_analysis.R --counts-dir data --output-dir results \
  --design data/MPRA_sequence_design.tsv.gz
```

The mapping needs `name`, `SPDI`, and `allele` columns. SPDI and allele use single-item JSON arrays; allele is `ref` or `alt`. Both original A1/A2 names must map to the same SPDI and opposite genomic alleles. Without that evidence, `SPDI` and `logAlt_Ref` in the allelic output remain missing. The script trusts the supplied design's genomic annotation; reference/sequence validation is a separate preparation step.

## Models and corrections

- **Allelic:** the original `mpralm` paired model is retained: `normalize=TRUE`, `aggregate="none"`, `model_type="corr_groups"`, with matched replicate blocks. The raw coefficient is **A2 minus A1**. `logAlt_Ref` is derived from the verified A1/A2 allele mapping, never from a presumed dbSNP allele order.
- **Activity:** the original intercept model is retained: `normalize=FALSE`, `aggregate="sum"`, `model_type="indep_groups"`. Effects are model estimates on the source log2 RNA/DNA scale, not effects centered on the negative-control median and not published MPRAnalyze alpha values.
- The activity threshold is the median, across negative controls, of each control's median log2(RNA/DNA). Nonfinite ratios are omitted. A negative or undefined threshold stops the run for review.
- `treat(lfc=threshold)` and BH FDR are retained. An activity call additionally requires a **positive effect above the same threshold**. The original extra `log2(threshold)` in the label expression is removed. `treat` tests effect magnitude, so its P-value should not be described as a one-sided activity test.
- VEH and STIM have separate inputs, fits, thresholds, FDR calculations, and output directories. The historical duplicated cleaned STIM files are not used.
- Both raw P-values and BH FDR are saved, as are fitted models and R/package versions. All fitted rows are exported; tables are not restricted to significant results.
- `emVar` requires an allelic call and at least one active allele for the same original variant identifier. This rule uses the explicit positive activity call; it is not necessarily identical to historical FDR-only emVar labels.

## Outputs

Each of `results/veh/` and `results/stim/` contains:

| File | Contents |
|---|---|
| `custom_activity.tsv` | Per-oligo activity effect, raw P, FDR, activity call and optional SPDI |
| `custom_allelic.tsv` | Per-pair A2−A1 effect, raw P, FDR, allelic call and optional genomic ALT/REF effect |
| `custom_emvar.tsv` | Variant-level activity, allelic and emVar calls |
| `allele_pairs.tsv` | Exact A1/A2 names used in the paired model |
| `unpaired_oligos.tsv` | Oligos retained for activity but unavailable for paired allelic analysis |
| `excluded_oligos.tsv` | All-missing DNA/RNA oligos excluded before fitting |
| `run_summary.tsv` | Input/design checksums, threshold and row counts |
| `model_fits.rds` | Fitted model objects for later statistics export |
| `sessionInfo.txt` | R and package versions |

Outputs use original oligo/variant identifiers as primary keys. Mapping disagreements with historical cleaned tables must be reconciled before combining with them. This workflow does not update those tables or the published submission candidate.

## GitHub contents

Upload this code directory (`run_analysis.R`, `README.md`, `.gitignore`). Input matrices, sequence designs and result/model files are excluded by default. No machine-specific paths or credentials are embedded in the code.

## Validation

The Reporter Experiment reader was checked against the complete original VEH/STIM
wide matrices, including missing values. This checks input equivalence, not a
new full-dataset model fit.

Checked with R 4.5.0 on 2026-10-12 using a temporary subset of 400 paired variants plus the source negative controls per condition. Both VEH and STIM completed with 800 activity rows and 400 allelic rows. Raw P/FDR ranges and all 399 mapped ALT/REF orientations per condition passed checks; the unmapped pair retained missing genomic orientation. This was an execution check, not a full-dataset refit or a reproduction claim for historical statistics. No original result files were overwritten.
