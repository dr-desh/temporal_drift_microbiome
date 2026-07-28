# ============================================================
# 02_analysis_lme.R
# Analysis 2 — Linear Mixed Effects Model
# treatment × timepoint interaction per OTU
# ============================================================

library(data.table)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(pbapply)

OUT <- file.path(getwd(), "outputs")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ── SECTION 1: Load metadata ─────────────────────────────────
cat("\n[1/5] Loading metadata...\n")

meta_raw <- fread(file.path(OUT, "merged_metadata_clean.csv"), header = FALSE)
col_names    <- as.character(meta_raw[1, ])
col_names[1] <- "sample_uid"
setnames(meta_raw, col_names)
meta <- meta_raw[-1]

if ("sample_id_2" %in% names(meta) & !"sample_uid" %in% names(meta)) {
  setnames(meta, "sample_id_2", "sample_uid")
}

meta[, treatment  := factor(treatment, levels = c("control", "fiber"))]
meta[, timepoint  := factor(timepoint,  levels = c("before", "after"))]
meta[, study      := as.factor(study)]
meta[, subject_id := as.factor(subject_id)]

cat("  Rows:", nrow(meta), "| Treatment:", paste(names(table(meta$treatment)), collapse="/"),
    "| Timepoint:", paste(names(table(meta$timepoint)), collapse="/"), "\n")


# ── SECTION 2: Select top 200 OTUs by SHAP ───────────────────
cat("\n[2/5] Loading SHAP rankings, selecting top 200 OTUs...\n")

shap     <- fread(file.path(OUT, "delta_feature_importances.csv"))
top200   <- shap[order(-mean_abs_shap)][1:200, OTU_ID]
cat("  Top 200 OTUs selected. Example:", head(top200, 2), "\n")


# ── SECTION 3: Load CLR matrix, subset, pivot to long ────────
cat("\n[3/5] Loading CLR matrix (2-3 min)...\n")

clr_wide <- fread(file.path(OUT, "preprocessed_otu_clr.csv"))
cat("  Loaded:", nrow(clr_wide), "samples x", ncol(clr_wide)-1, "OTUs\n")

# Keep sample_uid + top 200 OTU columns only
keep_cols <- c("sample_uid", intersect(top200, names(clr_wide)))
cat("  Top 200 OTUs found in CLR matrix:", length(keep_cols)-1, "\n")
clr_sub <- clr_wide[, ..keep_cols]

cat("  Pivoting to long format...\n")
clr_long <- melt(clr_sub,
                 id.vars       = "sample_uid",
                 variable.name = "OTU_ID",
                 value.name    = "clr_abundance")
clr_long[, sample_uid := as.character(sample_uid)]
clr_long[, OTU_ID     := as.character(OTU_ID)]
cat("  Long format rows:", nrow(clr_long), "\n")


# ── SECTION 4: Join metadata ──────────────────────────────────
cat("\n[4/5] Joining metadata...\n")

meta[, sample_uid := as.character(sample_uid)]
lme_data <- merge(clr_long, meta, by = "sample_uid", all.x = FALSE)
lme_data <- lme_data[!is.na(clr_abundance) & !is.na(treatment) &
                       !is.na(timepoint) & !is.na(subject_id) & !is.na(study)]

cat("  Merged rows   :", nrow(lme_data), "\n")
cat("  Unique OTUs   :", uniqueN(lme_data$OTU_ID), "\n")
cat("  Unique subjects:", uniqueN(lme_data$subject_id), "\n")

# ── SECTION 5: Run LME per OTU ───────────────────────────────
cat("\n[5/5] Running mixed effects models for 200 OTUs...\n")
cat("  Model: clr_abundance ~ treatment * timepoint + (1|subject_id) + (1|study)\n")
cat("  Estimated runtime: 30-60 min\n\n")

otu_list <- unique(lme_data$OTU_ID)

run_lme <- function(otu) {
  dat <- lme_data[OTU_ID == otu]
  tryCatch({
    fit <- lmer(clr_abundance ~ treatment * timepoint + (1 | subject_id) + (1 | study),
                data    = dat,
                REML    = FALSE,
                control = lmerControl(optimizer = "bobyqa",
                                      optCtrl   = list(maxfun = 2e5)))
    coefs   <- as.data.frame(coef(summary(fit)))
    int_row <- coefs[grep("treatmentfiber:timepointafter", rownames(coefs)), ]
    if (nrow(int_row) == 0) {
      return(data.table(OTU_ID=otu, estimate=NA, se=NA, t_value=NA, p_value=NA, converged=FALSE))
    }
    data.table(OTU_ID   = otu,
               estimate = int_row[1, "Estimate"],
               se       = int_row[1, "Std. Error"],
               t_value  = int_row[1, "t value"],
               p_value  = int_row[1, "Pr(>|t|)"],
               converged = TRUE)
  }, error = function(e) {
    data.table(OTU_ID=otu, estimate=NA, se=NA, t_value=NA, p_value=NA, converged=FALSE)
  })
}

results_list <- pblapply(otu_list, run_lme)
results      <- rbindlist(results_list)

# ── FDR correction ────────────────────────────────────────────
results_clean <- results[converged == TRUE & !is.na(p_value)]
results_clean[, q_value    := p.adjust(p_value, method = "BH")]
results_clean[, significant := q_value < 0.05]

# ── Summary ───────────────────────────────────────────────────
n_sig   <- sum(results_clean$significant)
n_total <- nrow(results_clean)

cat("\n========== ANALYSIS 2 COMPLETE ==========\n")
cat(sprintf("  Models converged      : %d / %d\n", n_total, length(otu_list)))
cat(sprintf("  Significant (FDR<0.05): %d / %d\n", n_sig, n_total))
cat(sprintf("  Proportion significant: %.1f%%\n", 100 * n_sig / n_total))
cat("\n  Top 20 by interaction effect size:\n")
print(results_clean[order(-abs(estimate))][1:20,
                                           .(OTU_ID, estimate, se, t_value, p_value, q_value, significant)])

# ── Save ──────────────────────────────────────────────────────
fwrite(results_clean[order(q_value)],
       file.path(OUT, "lme_interaction_results.csv"))
cat(sprintf("\nSaved: lme_interaction_results.csv (%d OTUs)\n", n_total))



# ── SECTION 6: Taxonomy annotation of significant OTUs ───────
cat("\n[6/6] Annotating significant OTUs with taxonomy...\n")

tax <- fread(file.path(OUT, "taxonomy_table.csv"))
cat("  Taxonomy table:", nrow(tax), "OTUs\n")
cat("  Columns:", paste(names(tax), collapse=", "), "\n")

# Merge significant results with taxonomy
sig_otus <- results_clean[significant == TRUE][order(q_value)]

tax_col <- names(tax)[1]   # OTU_ID column
setnames(tax, tax_col, "OTU_ID")

sig_annotated <- merge(sig_otus, tax, by = "OTU_ID", all.x = TRUE)
sig_annotated <- sig_annotated[order(q_value)]

cat("\n  Significant OTUs with taxonomy (top 20):\n")
print(sig_annotated[1:20, .(OTU_ID, estimate, q_value, taxonomy)])

# Save annotated results
fwrite(sig_annotated, file.path(OUT, "lme_significant_otus_annotated.csv"))
cat(sprintf("\nSaved: lme_significant_otus_annotated.csv (%d significant OTUs)\n", 
            nrow(sig_annotated)))
