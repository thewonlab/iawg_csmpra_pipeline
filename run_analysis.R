#!/usr/bin/env Rscript
# CSMPRA custom summary statistics. Original model: src/IAWG/A5.NANA_mpralm.Rmd.
# Raw allelic effects are A2 minus A1; genomic orientation requires a design map.

usage <- function() {
  cat(paste0(
    "Usage: Rscript run_analysis.R --counts-dir DIR --output-dir DIR\n",
    "       [--input-format matrix|reporter] [--condition both|veh|stim] [--design MPRA_sequence_design.tsv.gz]\n",
    "\nInputs: matrix ENIGMA.Veh/Stim.activity.csv; reporter Reporter_Experiment_Long_veh/stim.tsv.gz\n",
    "Outputs: separate condition directories; existing directories are never overwritten.\n"))
}

parse_args <- function(args) {
  if (!length(args) || "--help" %in% args) { usage(); quit(status = 0) }
  opts <- list(condition = "both", design = NULL, `input-format` = "reporter")
  allowed <- c("counts-dir", "output-dir", "condition", "design", "input-format")
  if (length(args) %% 2L) stop("Each option needs a value; see --help")
  for (i in seq(1, length(args), by = 2)) {
    key <- sub("^--", "", args[i])
    if (!startsWith(args[i], "--") || !key %in% allowed) stop("Unknown option: ", args[i])
    opts[[key]] <- args[i + 1]
  }
  if (is.null(opts[["counts-dir"]]) || is.null(opts[["output-dir"]])) stop("--counts-dir and --output-dir are required")
  if (!opts$condition %in% c("both", "veh", "stim")) stop("Invalid --condition")
  if (!opts[["input-format"]] %in% c("matrix", "reporter")) stop("Invalid --input-format")
  opts
}

input_name <- function(condition, format) {
  if (format == "reporter") paste0("Reporter_Experiment_Long_", condition, ".tsv.gz") else
    paste0("ENIGMA.", if (condition == "veh") "Veh" else "Stim", ".activity.csv")
}

read_counts <- function(input, format) {
  if (format == "matrix") {
    counts <- read.csv(input, check.names = FALSE, stringsAsFactors = FALSE)
    names(counts)[1] <- "oligo_name"
    return(counts)
  }
  con <- gzfile(input, open = "rt")
  on.exit(close(con))
  long <- read.delim(con, check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN"))
  required <- c("oligo_name", "replicate", "dna_counts", "rna_counts")
  if (!all(required %in% names(long))) stop("Reporter experiment needs: ", paste(required, collapse = ", "))
  if (!nrow(long) || anyNA(long$oligo_name) || any(!nzchar(long$oligo_name))) stop("Empty/missing oligo names")
  reps <- paste0("rep", c(1:11, 14:16))
  if (anyNA(long$replicate) || !setequal(unique(long$replicate), reps)) stop("Expected replicates rep1-rep11 and rep14-rep16")
  if (anyDuplicated(long[, c("oligo_name", "replicate")])) stop("Duplicate oligo/replicate rows")
  for (col in c("dna_counts", "rna_counts")) {
    x <- long[[col]]
    if (!is.numeric(x) || any(!is.na(x) & (!is.finite(x) | x < 0))) stop("Invalid raw counts: ", col)
  }
  oligos <- unique(long$oligo_name)
  counts <- data.frame(oligo_name = oligos, stringsAsFactors = FALSE)
  # Missing rows and missing count cells stay NA. Never use normalized counts.
  for (kind in c("dna", "rna")) {
    mat <- matrix(NA_real_, nrow = length(oligos), ncol = length(reps))
    mat[cbind(match(long$oligo_name, oligos), match(long$replicate, reps))] <- long[[paste0(kind, "_counts")]]
    for (j in seq_along(reps)) counts[[paste0(kind, sub("^rep", "", reps[j]))]] <- mat[, j]
  }
  counts
}

read_design <- function(path) {
  if (is.null(path)) return(NULL)
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite for --design")
  con <- if (grepl("\\.gz$", path)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con))
  d <- read.delim(con, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("name", "SPDI", "allele") %in% names(d))) stop("Design needs name, SPDI, allele")
  if (anyDuplicated(d$name)) stop("Duplicate oligo names in design")
  scalar <- function(x) {
    if (is.na(x) || !nzchar(x)) return(NA_character_)
    value <- jsonlite::fromJSON(x)
    if (length(value) != 1L) return(NA_character_)
    as.character(value)
  }
  d$SPDI <- vapply(d$SPDI, scalar, character(1))
  d$allele <- vapply(d$allele, scalar, character(1))
  d[, c("name", "SPDI", "allele")]
}

attach_orientation <- function(result, pairs, design) {
  result$SPDI <- NA_character_
  result$logAlt_Ref <- NA_real_
  result$orientation_status <- "no_verified_design_mapping"
  if (is.null(design)) return(result)
  a1 <- design[match(pairs$A1_name, design$name), ]
  a2 <- design[match(pairs$A2_name, design$name), ]
  valid <- !is.na(a1$SPDI) & !is.na(a2$SPDI) & a1$SPDI == a2$SPDI &
    ((a1$allele == "ref" & a2$allele == "alt") | (a1$allele == "alt" & a2$allele == "ref"))
  valid[is.na(valid)] <- FALSE
  pos <- match(result$variant_id, pairs$variant_id)
  good <- which(valid[pos])
  result$SPDI[good] <- a1$SPDI[pos[good]]
  result$logAlt_Ref[good] <- result$logFC_allelic[good] * ifelse(a2$allele[pos[good]] == "alt", 1, -1)
  result$orientation_status[good] <- "verified_design_alleles"
  result
}

run_condition <- function(condition, opts, design) {
  input <- file.path(opts[["counts-dir"]], input_name(condition, opts[["input-format"]]))
  out <- file.path(opts[["output-dir"]], condition)
  counts <- read_counts(input, opts[["input-format"]])
  if (anyNA(counts$oligo_name) || anyDuplicated(counts$oligo_name)) stop("Missing/duplicate oligo names")
  reps <- c(1:11, 14:16)
  dna_cols <- paste0("dna", reps); rna_cols <- paste0("rna", reps)
  if (!all(c(dna_cols, rna_cols) %in% names(counts))) stop("Missing expected count columns: ", input)
  raw <- as.matrix(counts[, c(dna_cols, rna_cols)])
  if (!is.numeric(raw) || any(!is.na(raw) & (!is.finite(raw) | raw < 0))) stop("Invalid count values")
  # Source matrices are post-QC oligo counts. Keep missing values; do not impute zeros.
  neg <- grepl("Scrambled", counts$oligo_name, fixed = TRUE)
  if (!any(neg)) stop("No Scrambled controls found")
  logrd <- log2(as.matrix(counts[neg, rna_cols]) / as.matrix(counts[neg, dna_cols]))
  logrd[!is.finite(logrd)] <- NA_real_
  control_medians <- apply(logrd, 1, median, na.rm = TRUE)
  threshold <- median(control_medians, na.rm = TRUE)
  if (!is.finite(threshold) || threshold < 0) stop("Negative-control threshold must be finite and nonnegative for treat()")
  v <- counts[!neg, ]
  if (any(!grepl("_A[12]_", v$oligo_name))) stop("Unexpected variant name; expected _A1_ or _A2_")
  v$variant_id <- sub("_A[12]_", "_", v$oligo_name)
  v$allele <- ifelse(grepl("_A1_", v$oligo_name), "A1", "A2")
  if (anyDuplicated(paste(v$variant_id, v$allele))) stop("Duplicate variant/allele pairs")
  # All-missing oligos cannot be fitted. Record them instead of silently imputing.
  usable <- rowSums(!is.na(v[, dna_cols])) > 0 & rowSums(!is.na(v[, rna_cols])) > 0
  excluded <- data.frame(oligo_name = v$oligo_name[!usable], reason = rep("all_missing_DNA_or_RNA", sum(!usable)))
  v <- v[usable, ]
  a1 <- v[v$allele == "A1", ]; a2 <- v[v$allele == "A2", ]
  keys <- a1$variant_id[a1$variant_id %in% a2$variant_id]
  a1 <- a1[match(keys, a1$variant_id), ]; a2 <- a2[match(keys, a2$variant_id), ]
  if (length(keys) < 10) stop("Too few complete allele pairs to fit mpralm")
  pairs <- data.frame(variant_id = keys, A1_name = a1$oligo_name, A2_name = a2$oligo_name)
  # Column order is A1 replicates, then A2 replicates. Coefficient 2 is A2-A1.
  dna <- cbind(as.matrix(a1[, dna_cols]), as.matrix(a2[, dna_cols]))
  rna <- cbind(as.matrix(a1[, rna_cols]), as.matrix(a2[, rna_cols]))
  dimnames(dna) <- dimnames(rna) <- list(keys, NULL)
  group <- rep(c(0, 1), each = length(reps))
  model <- model.matrix(~ group)
  fit_allelic <- mpra::mpralm(mpra::MPRASet(DNA = dna, RNA = rna, eid = keys),
    design = model, aggregate = "none", normalize = TRUE,
    block = rep(seq_along(reps), 2), model_type = "corr_groups", plot = FALSE)
  tab <- limma::topTable(fit_allelic, coef = 2, number = Inf, sort.by = "none")
  allelic <- data.frame(variant_id = rownames(tab), logFC_allelic = tab$logFC,
    pval_allelic = tab$P.Value, fdr_allelic = p.adjust(tab$P.Value, "BH"))
  allelic$allelic <- allelic$fdr_allelic < 0.05
  allelic <- attach_orientation(allelic, pairs, design)
  dna <- as.matrix(v[, dna_cols]); rna <- as.matrix(v[, rna_cols])
  dimnames(dna) <- dimnames(rna) <- list(v$oligo_name, NULL)
  fit_activity <- mpra::mpralm(mpra::MPRASet(DNA = dna, RNA = rna, eid = v$oligo_name),
    design = matrix(1, nrow = length(reps), ncol = 1, dimnames = list(NULL, "intercept")),
    aggregate = "sum", normalize = FALSE, model_type = "indep_groups", plot = FALSE)
  tested_activity <- limma::treat(fit_activity, lfc = threshold)
  tab <- limma::topTreat(tested_activity, coef = 1, number = Inf, sort.by = "none")
  activity <- data.frame(oligo_name = rownames(tab), logFC_active = tab$logFC,
    pval_active = tab$P.Value, fdr_active = p.adjust(tab$P.Value, "BH"))
  idx <- match(activity$oligo_name, v$oligo_name)
  activity$variant_id <- v$variant_id[idx]; activity$allele <- v$allele[idx]
  # threshold already is log2(RNA/DNA). Do not take log2 a second time.
  # treat() tests magnitude; positive direction is required for an activity call.
  activity$active <- activity$fdr_active < 0.05 & activity$logFC_active > threshold
  activity$SPDI <- if (is.null(design)) NA_character_ else design$SPDI[match(activity$oligo_name, design$name)]
  any_active <- tapply(activity$active, activity$variant_id, function(x) {
    if (any(x %in% TRUE)) TRUE else if (anyNA(x)) NA else FALSE
  })
  calls <- allelic[, c("variant_id", "SPDI", "allelic")]
  calls$active <- unname(any_active[calls$variant_id])
  calls$emVar <- calls$allelic & calls$active
  unpaired <- v[!v$variant_id %in% keys, c("oligo_name", "variant_id", "allele")]
  dir.create(out, recursive = TRUE)
  write_tsv <- function(x, name) write.table(x, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  write_tsv(activity, "custom_activity.tsv")
  write_tsv(allelic, "custom_allelic.tsv")
  write_tsv(calls, "custom_emvar.tsv")
  write_tsv(pairs, "allele_pairs.tsv")
  write_tsv(unpaired, "unpaired_oligos.tsv")
  write_tsv(excluded, "excluded_oligos.tsv")
  provenance <- data.frame(condition = condition, input_format = opts[["input-format"]], input = normalizePath(input),
    input_md5 = unname(tools::md5sum(input)),
    design = if (is.null(opts$design)) NA_character_ else normalizePath(opts$design),
    design_md5 = if (is.null(opts$design)) NA_character_ else unname(tools::md5sum(opts$design)),
    negative_control_log2_threshold = threshold, controls = sum(neg),
    activity_rows = nrow(activity), allelic_rows = nrow(allelic),
    oriented_rows = sum(!is.na(allelic$logAlt_Ref)))
  write_tsv(provenance, "run_summary.tsv")
  saveRDS(list(allelic = fit_allelic, activity = fit_activity, activity_treat = tested_activity), file.path(out, "model_fits.rds"))
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  message(condition, ": ", nrow(activity), " activity rows; ", nrow(allelic), " allelic rows; written to ", out)
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  for (p in c("mpra", "limma")) if (!requireNamespace(p, quietly = TRUE)) stop("Missing R package: ", p)
  suppressPackageStartupMessages(library(mpra))
  conditions <- if (opts$condition == "both") c("veh", "stim") else opts$condition
  # Preflight all destinations before starting either condition.
  for (c in conditions) {
    if (dir.exists(file.path(opts[["output-dir"]], c))) stop("Output directory already exists: ", c)
    name <- input_name(c, opts[["input-format"]])
    if (!file.exists(file.path(opts[["counts-dir"]], name))) stop("Input file not found: ", name)
  }
  design <- read_design(opts$design)
  for (c in conditions) run_condition(c, opts, design)
}

if (sys.nframe() == 0) main()
