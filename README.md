# iawg_csmpra_pipeline

CSMPRA activity and allelic-effect analysis using `mpralm` and `limma` (GEO: GSE289858). VEH and STIM are analyzed separately.

## Requirements

R packages: `mpra`, `limma`, and `jsonlite` for optional REF/ALT mapping.

## Inputs

Place these IGVF Reporter Experiment files in `data/`:

- `Reporter_Experiment_Long_veh.tsv.gz`
- `Reporter_Experiment_Long_stim.tsv.gz`

Required columns: `oligo_name`, `replicate`, `dna_counts`, `rna_counts`. Include Scrambled controls and replicates `rep1`–`rep11`, `rep14`–`rep16`. Raw counts are used; missing values are preserved.

## Run

```bash
Rscript run_analysis.R --counts-dir data --output-dir results
```

Add `--design data/MPRA_sequence_design.tsv.gz` to map A2−A1 effects to genomic ALT−REF using verified sequence-design annotations. Use `--condition veh` or `--condition stim` to run one condition.

## Outputs

`results/veh/` and `results/stim/` contain `custom_activity.tsv`, `custom_allelic.tsv`, and `custom_emvar.tsv`, plus model fits and run records. Tables include effect estimates, raw P-values, and BH FDR. An emVar requires a significant allelic effect and at least one active allele.
