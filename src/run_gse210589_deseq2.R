#!/usr/bin/env Rscript
# Frozen GSE210589 rin-1 versus WT perturbation analysis with transcript-level tximport.

suppressPackageStartupMessages({
  library(DESeq2)
  library(tximport)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("usage: run_gse210589_deseq2.R SAMPLES QUANT_ROOT TX2GENE RIN_BINDING OUTPUT_DIR SEED")
}
samples_path <- args[[1]]
quant_root <- args[[2]]
tx2gene_path <- args[[3]]
binding_path <- args[[4]]
output_dir <- args[[5]]
seed <- as.integer(args[[6]])
set.seed(seed)

samples <- read_csv(samples_path, show_col_types = FALSE)
samples <- samples[samples$inclusion_status == "admitted", ]
if (nrow(samples) != 6 || any(table(samples$condition) != 3)) {
  stop("Frozen design requires exactly three WT and three rin-1 biological replicates")
}
files <- file.path(quant_root, samples$sra_run, "quant.sf")
complete <- file.path(quant_root, samples$sra_run, "COMPLETE")
if (!all(file.exists(files)) || !all(file.exists(complete))) {
  stop("All six version-locked Salmon quantifications must be complete")
}
names(files) <- samples$sra_run
tx2gene <- read_tsv(tx2gene_path, show_col_types = FALSE)
txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "no",
  dropInfReps = TRUE
)

coldata <- data.frame(
  row.names = samples$sra_run,
  condition = relevel(factor(samples$condition), ref = "WT_34DPA")
)
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~condition)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("condition", "rin1_34DPA", "WT_34DPA"), alpha = 0.05)
res_table <- as.data.frame(res)
res_table$gene_id <- rownames(res_table)
res_table <- res_table[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
res_table$differential <- !is.na(res_table$padj) & res_table$padj < 0.05 & abs(res_table$log2FoldChange) >= 1
res_table <- res_table[order(res_table$padj, na.last = TRUE), ]

binding <- read_csv(binding_path, show_col_types = FALSE)
binding <- binding[binding$admission_status == "admitted", ]
evidence <- merge(binding, res_table, by = "gene_id", all.x = TRUE, sort = FALSE)
evidence$perturbation_direction <- ifelse(
  is.na(evidence$log2FoldChange), "not_measured",
  ifelse(evidence$log2FoldChange > 0, "higher_in_rin1", "lower_in_rin1")
)
evidence$evidence_grade <- ifelse(
  evidence$differential %in% TRUE,
  "binding_plus_independent_perturbation_expression",
  "binding_only_or_no_significant_perturbation"
)
evidence$claim_boundary <- paste(
  "RIN occupancy and an independently generated rin-1 expression change support a context-limited candidate;",
  "regulatory sign, direct causality, and mechanism are not established."
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(res_table, file.path(output_dir, "GSE210589_rin1_vs_WT_DESeq2.csv"))
write_csv(evidence, file.path(output_dir, "GSE210589_RIN_binding_perturbation_evidence.csv"))
vsd <- vst(dds, blind = TRUE)
sample_cor <- cor(assay(vsd), method = "spearman")
write_csv(data.frame(sample_id = rownames(sample_cor), sample_cor, check.names = FALSE),
          file.path(output_dir, "GSE210589_sample_spearman.csv"))

summary_json <- sprintf(
  paste0(
    "{\n",
    "  \"source_id\": \"GSE210589\",\n",
    "  \"contrast\": \"rin-1_34DPA_vs_WT_34DPA\",\n",
    "  \"biological_replicates_per_group\": 3,\n",
    "  \"salmon_version\": \"1.10.2\",\n",
    "  \"DESeq2_version\": \"%s\",\n",
    "  \"tximport_version\": \"%s\",\n",
    "  \"genes_tested\": %d,\n",
    "  \"genes_differential_padj_lt_0_05_abs_lfc_ge_1\": %d,\n",
    "  \"admitted_signature_RIN_binding_genes\": %d,\n",
    "  \"binding_plus_perturbation_genes\": %d,\n",
    "  \"random_seed\": %d,\n",
    "  \"claim_boundary\": \"Perturbation-expression support is context-limited and is not by itself a causal regulatory sign.\"\n",
    "}"
  ),
  as.character(packageVersion("DESeq2")),
  as.character(packageVersion("tximport")),
  nrow(res_table),
  sum(res_table$differential, na.rm = TRUE),
  nrow(evidence),
  sum(evidence$evidence_grade == "binding_plus_independent_perturbation_expression", na.rm = TRUE),
  seed
)
writeLines(summary_json, file.path(output_dir, "GSE210589_deseq2_summary.json"))
cat(summary_json, "\n")
