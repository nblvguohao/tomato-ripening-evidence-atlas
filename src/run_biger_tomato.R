#!/usr/bin/env Rscript
# BiGER baseline adapter for tomato signed effect ranks.
# BiGER itself aggregates unsigned importance. This adapter gives it absolute
# within-study ranks, then adds a clearly separate precision-weighted direction
# layer; it must not be called a native signed BiGER model.

args <- commandArgs(trailingOnly = TRUE)
value <- function(flag, required = TRUE) {
  index <- match(flag, args)
  if (is.na(index)) {
    if (required) stop(sprintf("Missing %s", flag))
    return(NULL)
  }
  args[[index + 1]]
}

suppressPackageStartupMessages(library(BiGER))

contrasts_path <- value("--contrasts")
genes_path <- value("--gene-list")
output_dir <- value("--output-dir")
external_contrasts_path <- value("--external-contrasts", required = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

raw <- read.csv(contrasts_path, check.names = FALSE)
keep <- read.csv(genes_path, check.names = FALSE)$gene_id
raw <- raw[raw$gene_id %in% keep & is.finite(raw$signed_effect_rank), ]

# Collapse shared-publication strata before all fit/holdout operations.
units <- aggregate(signed_effect_rank ~ independence_group + gene_id, raw, mean)
names(units)[names(units) == "signed_effect_rank"] <- "effect"
context_source <- unique(raw[, c("independence_group", "early_ordinal", "late_ordinal")])
context_source$transition <- paste(context_source$early_ordinal, context_source$late_ordinal, sep = "->")
context_map <- setNames(context_source$transition, context_source$independence_group)

fit_one <- function(train, target_context = NULL) {
  genes <- sort(unique(train$gene_id))
  studies <- sort(unique(train$independence_group))
  matrix <- matrix(NA_real_, nrow = length(genes), ncol = length(studies), dimnames = list(genes, studies))
  for (study in studies) {
    subset <- train[train$independence_group == study, ]
    matrix[subset$gene_id, study] <- subset$effect
  }
  rank_matrix <- apply(abs(matrix), 2, function(x) rank(-x, na.last = "keep", ties.method = "min"))
  if (is.null(dim(rank_matrix))) rank_matrix <- matrix(rank_matrix, ncol = 1, dimnames = dimnames(matrix))
  n_ranked <- colSums(is.finite(rank_matrix))
  fit <- BiGER(r = rank_matrix, n_r = n_ranked, n_u = rep(0, length(studies)), return_mu_s2 = TRUE)
  precision <- fit$sigma2_inv
  signs <- sign(matrix)
  signed_weight <- sweep(signs, 2, precision, "*")
  denominator <- rowSums(sweep(1 * is.finite(signs), 2, precision, "*"), na.rm = TRUE)
  global_direction <- rowSums(signed_weight, na.rm = TRUE) / denominator
  global_direction[!is.finite(global_direction)] <- 0
  direction <- global_direction
  # Context extension: only the held study's known stage transition is used.
  # A context direction is shrunk toward global direction when it has few
  # precision-weighted supporting studies, preventing a single study from
  # dominating the signed prediction.
  if (!is.null(target_context)) {
    same_context <- unname(context_map[studies]) == target_context
    if (any(same_context)) {
      context_weight <- precision[same_context]
      context_signs <- signs[, same_context, drop = FALSE]
      context_denominator <- rowSums(sweep(1 * is.finite(context_signs), 2, context_weight, "*"), na.rm = TRUE)
      context_direction <- rowSums(sweep(context_signs, 2, context_weight, "*"), na.rm = TRUE) / context_denominator
      context_direction[!is.finite(context_direction)] <- global_direction[!is.finite(context_direction)]
      shrinkage <- context_denominator / (context_denominator + 2.0)
      direction <- global_direction + shrinkage * (context_direction - global_direction)
    }
  }
  # BiGER's latent mu is ordered importance, not a non-negative magnitude.
  # Convert only by a monotone rank transform before attaching direction.
  importance <- rank(fit$mu, ties.method = "average") / length(fit$mu)
  data.frame(gene_id = genes, biger_mu = fit$mu, biger_importance_percentile = importance, direction_vote = direction,
             biger_signed_score = importance * sign(direction), direction_probability = abs(direction),
             stringsAsFactors = FALSE)
}

metric <- function(prediction, held, top_k) {
  joined <- merge(prediction, held, by = "gene_id")
  rho <- cor(joined$biger_signed_score, joined$effect, method = "spearman")
  concordance <- mean(sign(joined$biger_signed_score) == sign(joined$effect))
  k <- min(top_k, nrow(joined))
  p_top <- joined$gene_id[order(abs(joined$biger_signed_score), decreasing = TRUE)[seq_len(k)]]
  h_top <- joined$gene_id[order(abs(joined$effect), decreasing = TRUE)[seq_len(k)]]
  jaccard <- length(intersect(p_top, h_top)) / length(union(p_top, h_top))
  data.frame(gene_count = nrow(joined), spearman_effect_rank = rho, direction_concordance = concordance,
             top_k = k, top_k_jaccard = jaccard)
}

full <- fit_one(units)
write.csv(full, file.path(output_dir, "biger_abs_rank_signed_wrapper_scores.csv"), row.names = FALSE)

if (!is.null(external_contrasts_path)) {
  external_raw <- read.csv(external_contrasts_path, check.names = FALSE)
  external_raw <- external_raw[external_raw$gene_id %in% keep & is.finite(external_raw$signed_effect_rank), ]
  external <- aggregate(signed_effect_rank ~ independence_group + gene_id, external_raw, mean)
  names(external)[names(external) == "signed_effect_rank"] <- "effect"
  external_rows <- list()
  for (external_group in sort(unique(external$independence_group))) {
    held <- external[external$independence_group == external_group, c("gene_id", "effect")]
    result <- metric(full, held, 200)
    result$held_out_independence_group <- external_group
    result$method <- "biger_abs_rank_signed_wrapper_frozen_external"
    external_rows[[length(external_rows) + 1]] <- result
  }
  write.csv(do.call(rbind, external_rows), file.path(output_dir, "biger_frozen_external_validation.csv"), row.names = FALSE)
}

rows <- list()
for (held_group in sort(unique(units$independence_group))) {
  fit <- fit_one(units[units$independence_group != held_group, ], target_context = context_map[[held_group]])
  held <- units[units$independence_group == held_group, c("gene_id", "effect")]
  result <- metric(fit, held, 200)
  result$held_out_independence_group <- held_group
  result$method <- "context_biger_signed_wrapper"
  rows[[length(rows) + 1]] <- result
}
benchmark <- do.call(rbind, rows)
benchmark <- benchmark[, c("held_out_independence_group", "method", "gene_count", "spearman_effect_rank", "direction_concordance", "top_k", "top_k_jaccard")]
write.csv(benchmark, file.path(output_dir, "context_biger_signed_wrapper_loco.csv"), row.names = FALSE)
writeLines(sprintf("median_loco_spearman=%.8f", median(benchmark$spearman_effect_rank)), file.path(output_dir, "summary.txt"))
