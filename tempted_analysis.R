# ============================================================
# 03_analysis_tempted.R
# Analysis 3 — TEMPTED temporal tensor decomposition
# Input : preprocessed_otu_clr.csv
#         merged_metadata_clean.csv
# Output: tempted_subject_loadings.csv
#         tempted_feature_loadings.csv
# ============================================================

library(tempted)
library(data.table)
library(ggplot2)
library(pROC)

OUT <- file.path(getwd(), "outputs")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ── SECTION 1: Load metadata ──────────────────────────────────
cat("\n[1/4] Loading metadata...\n")

meta_raw     <- fread(file.path(OUT, "merged_metadata_clean.csv"), header = FALSE)
col_names    <- as.character(meta_raw[1, ])
col_names[1] <- "sample_uid"
setnames(meta_raw, col_names)
meta <- meta_raw[-1]

# Keep one row per subject — use after timepoint for treatment label
meta_subject <- meta[timepoint == "after", 
                     .(subject_id, treatment, study)] |> unique()

cat("  Subjects:", nrow(meta_subject), "\n")
cat("  Treatment:\n")
print(table(meta_subject$treatment))


# ── SECTION 2: Load CLR matrix ────────────────────────────────
cat("\n[2/4] Loading CLR matrix...\n")

clr_wide <- fread(file.path(OUT, "preprocessed_otu_clr.csv"))
cat("  Loaded:", nrow(clr_wide), "samples x", ncol(clr_wide)-1, "OTUs\n")

# Use top 200 SHAP OTUs to keep computation manageable
shap   <- fread(file.path(OUT, "delta_feature_importances.csv"))
top200 <- shap[order(-mean_abs_shap)][1:200, OTU_ID]
top200_present <- intersect(top200, names(clr_wide))
cat("  Top 200 OTUs present in CLR:", length(top200_present), "\n")

keep_cols <- c("sample_uid", top200_present)
clr_sub   <- clr_wide[, ..keep_cols]

# ── SECTION 3: Prepare data for TEMPTED ──────────────────────
cat("\n[3/4] Preparing data for TEMPTED...\n")

# Align meta to clr_sub sample order
clr_samples  <- clr_sub$sample_uid
meta_aligned <- meta[match(clr_samples, sample_uid)]

# Identify subjects with both timepoints
subject_timepoint_counts <- meta_aligned[, .(n_timepoints = uniqueN(timepoint)), 
                                         by = subject_id]
valid_subjects <- subject_timepoint_counts[n_timepoints == 2, subject_id]
invalid_subjects <- subject_timepoint_counts[n_timepoints != 2, subject_id]

cat("  Removing", length(invalid_subjects), "subjects with single timepoint:",
    paste(invalid_subjects, collapse=", "), "\n")
cat("  Keeping", length(valid_subjects), "subjects with both timepoints\n")

# Filter meta and clr to valid subjects only
meta_valid   <- meta_aligned[subject_id %in% valid_subjects]
valid_samples <- meta_valid$sample_uid
clr_valid    <- clr_sub[sample_uid %in% valid_samples]

# Realign after filter
clr_valid_ordered  <- clr_valid[match(meta_valid$sample_uid, sample_uid)]

# Build feature matrix — rows = samples, cols = OTUs
feat_mat <- as.matrix(clr_valid_ordered[, -1])
rownames(feat_mat) <- clr_valid_ordered$sample_uid

# Build timepoint and subjectID vectors
timepoint_vec <- ifelse(meta_valid$timepoint == "before", 0, 1)
subject_vec   <- as.character(meta_valid$subject_id)

cat("  Feature matrix   :", nrow(feat_mat), "samples x", ncol(feat_mat), "OTUs\n")
cat("  Timepoint (0/1)  :", paste(names(table(timepoint_vec)), 
                                  table(timepoint_vec), sep="=", collapse=", "), "\n")
cat("  Unique subjects  :", length(unique(subject_vec)), "\n")


# ── SECTION 4: Run TEMPTED ────────────────────────────────────
cat("\n[4/4] Running tempted_all()...\n")
cat("  Data is already CLR-transformed — setting transform='none'\n")
cat("  This may take 10-30 min...\n")

set.seed(42)
tempted_result <- tempted_all(
  featuretable = feat_mat,
  timepoint    = timepoint_vec,
  subjectID    = subject_vec,
  r            = 2,
  smooth       = 1e-5,
  transform    = "none",
  pseudo       = NULL
)

cat("\n  tempted_all complete.\n")
cat("  Result components:", paste(names(tempted_result), collapse=", "), "\n")

# ── Extract subject loadings ──────────────────────────────────
cat("\nExtracting subject loadings...\n")



subj_load <- as.data.table(tempted_result$A_hat, keep.rownames = "subject_id")

# Join treatment label
meta_subject[, subject_id := as.character(subject_id)]
subj_load[,   subject_id := as.character(subject_id)]
subj_load <- merge(subj_load, meta_subject, by = "subject_id", all.x = TRUE)

cat("  Subject loadings shape:", nrow(subj_load), "x", ncol(subj_load), "\n")
cat("  Treatment in loadings:\n")
print(table(subj_load$treatment))

# ── Classify fiber vs control ─────────────────────────────────
cat("\nClassifying fiber vs control using subject loadings...\n")

subj_load[, treatment_bin := ifelse(treatment == "fiber", 1, 0)]

roc1 <- roc(subj_load$treatment_bin, subj_load$PC1, quiet = TRUE)
roc2 <- roc(subj_load$treatment_bin, subj_load$PC2, quiet = TRUE)
cat(sprintf("  Component 1 (PC1) AUC: %.4f\n", auc(roc1)))
cat(sprintf("  Component 2 (PC2) AUC: %.4f\n", auc(roc2)))

cat("\n  PC1 mean by treatment:\n")
print(subj_load[, .(mean_PC1 = mean(PC1), sd_PC1 = sd(PC1),
                    mean_PC2 = mean(PC2), sd_PC2 = sd(PC2)), by = treatment])

# ── Extract feature loadings ──────────────────────────────────
feat_load <- as.data.table(tempted_result$B_hat, keep.rownames = "OTU_ID")
shap[, rank := 1:.N]
feat_load <- merge(feat_load, shap[, .(OTU_ID, mean_abs_shap, rank)],
                   by = "OTU_ID", all.x = TRUE)
feat_load <- feat_load[order(-abs(PC1))]

cat("\n  Top 20 OTUs by PC1 loading:\n")
print(feat_load[1:20, .(OTU_ID, PC1, PC2, mean_abs_shap, rank)])

# ── Save outputs ──────────────────────────────────────────────
fwrite(subj_load, file.path(OUT, "tempted_subject_loadings.csv"))
fwrite(feat_load, file.path(OUT, "tempted_feature_loadings.csv"))
cat("\nSaved: tempted_subject_loadings.csv\n")
cat("Saved: tempted_feature_loadings.csv\n")

# ── Plot ──────────────────────────────────────────────────────
p <- ggplot(subj_load, aes(x = PC1, y = PC2, color = treatment)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("control" = "#4575b4", "fiber" = "#d73027")) +
  labs(title = "TEMPTED Subject Loadings",
       x = "Component 1 (PC1)", y = "Component 2 (PC2)",
       color = "Treatment") +
  theme_bw()

ggsave(file.path(OUT, "tempted_subject_loadings_plot.png"),
       p, width = 7, height = 5, dpi = 150)
cat("Saved: tempted_subject_loadings_plot.png\n")

# ── Final summary ─────────────────────────────────────────────
cat("\n========== ANALYSIS 3 COMPLETE ==========\n")
cat(sprintf("  Subjects analysed : %d\n", nrow(subj_load)))
cat(sprintf("  PC1 AUC           : %.4f\n", auc(roc1)))
cat(sprintf("  PC2 AUC           : %.4f\n", auc(roc2)))
cat("  Variance explained:\n")
cat(sprintf("  Component 1 R²    : %.4f\n", tempted_result$r_square[1]))
cat(sprintf("  Component 2 R²    : %.4f\n", tempted_result$r_square[2]))
cat(sprintf("  Cumulative R²     : %.4f\n", tempted_result$accum_r_square[2]))
cat("=========================================\n")