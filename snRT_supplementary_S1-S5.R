###############################################################################
# SUPPLEMENTARY ANALYSES S1–S5  (append to the END of snRT.R)
#
# Run AFTER the pseudobulk block of snRT.R, in the SAME R session, so that
# these objects are in memory:  ifnb, pb, lut, de_all, cts, fit_dds,
#                               out_dir, PADJ, LFC, MIN_REPS_PER_SIDE
# Quick check before running:   stopifnot(exists("pb"), exists("lut"),
#                                         exists("de_all"), exists("ifnb"))
#
# Produces, in <out_dir>/supplementary/ :
#   S1  within-celltype (composition-adjusted) baseline distance   (R2 #2)
#   S2  community-coloured sample graph + membership + stability   (R2 #2)
#   S3  permutation null for the metformin reversal fraction       (R2 #3)
#   S4  QC / UMAP / marker / provenance from the annotated object  (R1 #2, R2 #1)
#   S5  segmentation-check & inflammation-score correlation stubs  (R2 #5)
###############################################################################
suppressPackageStartupMessages({
  library(Matrix); library(igraph); library(dplyr); library(ggplot2)
})
sup_dir <- file.path(out_dir, "supplementary")
dir.create(sup_dir, showWarnings = FALSE, recursive = TRUE)

## paths / keys used below (define defensively in case only pseudobulk was run)
if (!exists("cond_key"))
  cond_key <- c(Day0_Control="Day 0", Day14_Control="Day 14",
                Day14_Metformin="Day 14 metformin",
                Day21_Control="Day 21", Day21_Metformin="Day 21 metformin")
# <-- EDIT this to your imaging output if different:
IMAGING_DIST_CSV <- "/Users/mac/Desktop/itegrinmetprev/image analysis/results_graph/distance_to_baseline.csv"

.log2cpm <- function(M){                       # M = genes x samples counts
  cs <- Matrix::colSums(M); cs[cs == 0] <- 1
  log2(sweep(as.matrix(M), 2, cs, "/") * 1e6 + 1)
}
.top_hvg <- function(L, n = 2000){             # top-variance genes of a log matrix
  v <- apply(L, 1, var); names(sort(v, decreasing = TRUE))[seq_len(min(n, sum(v > 0)))]
}

###############################################################################
## S1. Within-cell-type baseline distance (composition-adjusted) — R2 #2
##  The main-text transcriptomic distance is whole-sample pseudobulk, so the
##  large composition shift can inflate it. Here the distance is computed
##  WITHIN each major cell type (per sample), then averaged first within
##  (condition x cell type) and then across cell types, so abundant cell types
##  and composition changes do not dominate. If Day21_Metformin is still the
##  disease condition closest to baseline, the convergence claim is robust to
##  the composition confound.
###############################################################################
tryCatch({
  base_sg <- "Day0_Control"
  # cell types with >=2 Day0 samples (needed to define a baseline correlation)
  d0_counts <- table(lut$celltype[lut$Subgroup == base_sg])
  ct_ok <- names(d0_counts)[d0_counts >= MIN_REPS_PER_SIDE]
  ct_ok <- intersect(cts, ct_ok)

  rows <- list()
  for (ct in ct_ok) {
    cols <- which(lut$celltype == ct)
    if (length(cols) < 4) next                 # need a few samples to be meaningful
    M  <- pb[, cols, drop = FALSE]
    colnames(M) <- lut$sample[cols]
    sg <- as.character(lut$Subgroup[cols])
    L  <- .log2cpm(M)
    hv <- .top_hvg(L, 2000); L <- L[hv, , drop = FALSE]
    C  <- suppressWarnings(cor(L))             # sample x sample within this cell type
    d0 <- which(sg == base_sg)
    if (length(d0) < 2) next
    for (j in seq_len(ncol(L))) {
      ref <- setdiff(d0, j)                    # leave-one-out for Day0 samples
      rows[[length(rows) + 1]] <- data.frame(
        celltype = ct, sample = colnames(L)[j], Subgroup = sg[j],
        dist = 1 - mean(C[j, ref]), stringsAsFactors = FALSE)
    }
  }
  per_samp <- dplyr::bind_rows(rows)
  write.csv(per_samp, file.path(sup_dir, "S1_within_celltype_distance_per_sample.csv"),
            row.names = FALSE)

  # average within (condition x cell type), then across cell types (unweighted)
  by_cs <- per_samp %>% dplyr::group_by(Subgroup, celltype) %>%
    dplyr::summarise(d = mean(dist), .groups = "drop")
  within_ct <- by_cs %>% dplyr::group_by(Subgroup) %>%
    dplyr::summarise(within_celltype_dist = mean(d),
                     sd = sd(d), n_celltypes = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(condition = cond_key[Subgroup])

  # compare with the whole-sample distance (from the integration output, if present)
  ws_path <- file.path(out_dir, "cross_modal_integration.csv")
  if (file.exists(ws_path)) {
    ws <- read.csv(ws_path)[, c("Subgroup","transcriptomic_dist")]
    within_ct <- dplyr::left_join(within_ct, ws, by = "Subgroup")
  }
  within_ct <- within_ct %>% dplyr::arrange(within_celltype_dist)
  write.csv(within_ct, file.path(sup_dir, "S1_within_celltype_distance_by_condition.csv"),
            row.names = FALSE)
  cat("\n== S1 within-cell-type distance to Day 0 (composition-adjusted) ==\n")
  print(as.data.frame(within_ct))
  dgroups <- within_ct %>% dplyr::filter(Subgroup != base_sg)
  closest <- dgroups$condition[which.min(dgroups$within_celltype_dist)]
  cat(sprintf("Closest disease condition to baseline (within-cell-type): %s\n", closest))

  ggplot(within_ct, aes(reorder(condition, within_celltype_dist), within_celltype_dist)) +
    geom_col(fill = "#4C72B0", width = .65) + coord_flip() +
    labs(x = NULL, y = "Within-cell-type distance to Day 0 (mean over cell types)",
         title = "S1 Composition-adjusted transcriptomic distance to baseline") +
    theme_bw(12)
  ggsave(file.path(sup_dir, "S1_within_celltype_distance.png"), width = 7, height = 4.5, dpi = 200)
}, error = function(e) message("S1 failed: ", conditionMessage(e)))

###############################################################################
## S2. Community-coloured sample-similarity graph + membership + stability — R2 #2
##  Recolour the sample graph by DETECTED community (not condition) so the
##  Day21_Metformin membership claim is verifiable, and test stability of the
##  Day21_Metformin<->Day0 co-clustering across kNN values and seeds.
###############################################################################
tryCatch({
  # whole-sample pseudobulk (collapse cell types), log2 CPM, HVGs
  samp <- lut$sample
  Msamp <- t(rowsum(t(as.matrix(pb)), samp))           # genes x sample
  s2c <- lut %>% dplyr::distinct(sample, Subgroup)
  s2c <- s2c[match(colnames(Msamp), s2c$sample), ]
  L <- .log2cpm(Msamp); L <- L[.top_hvg(L, 2000), , drop = FALSE]
  Cs <- cor(L)                                          # sample x sample

  build_graph <- function(C, k){
    n <- ncol(C); A <- matrix(0, n, n, dimnames = dimnames(C))
    for (i in seq_len(n)) { nb <- order(C[i, ], decreasing = TRUE)[2:(k + 1)]; A[i, nb] <- C[i, nb] }
    A <- pmax(A, t(A))
    igraph::graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
  }

  KNN_MAIN <- 4
  g <- build_graph(Cs, KNN_MAIN)
  V(g)$condition <- as.character(s2c$Subgroup)
  set.seed(1); comm <- igraph::cluster_louvain(g, weights = E(g)$weight)
  V(g)$community <- igraph::membership(comm)
  memb <- data.frame(sample = V(g)$name, condition = V(g)$condition,
                     community = V(g)$community)
  write.csv(memb, file.path(sup_dir, "S2_graph_communities.csv"), row.names = FALSE)

  pal <- setNames(grDevices::hcl.colors(max(V(g)$community), "Dark3"),
                  sort(unique(V(g)$community)))
  png(file.path(sup_dir, "S2_graph_by_community.png"), 1500, 1150, res = 200)
  set.seed(1); lay <- igraph::layout_with_fr(g)
  plot(g, layout = lay, vertex.color = pal[as.character(V(g)$community)],
       vertex.size = 10, vertex.label = V(g)$condition, vertex.label.cex = 0.5,
       edge.width = 1 + 3 * E(g)$weight,
       main = "Sample-similarity graph coloured by Louvain community")
  legend("bottomleft", legend = paste("community", names(pal)), col = pal,
         pch = 19, cex = 0.7, bty = "n")
  dev.off()

  # stability: does Day21_Metformin share Day0's community? across kNN and seeds
  d0_samp <- s2c$sample[s2c$Subgroup == "Day0_Control"]
  d21m     <- s2c$sample[s2c$Subgroup == "Day21_Metformin"]
  stab <- list()
  for (k in 3:6) for (seed in 1:10) {
    gg <- build_graph(Cs, k); set.seed(seed)
    m <- igraph::membership(igraph::cluster_louvain(gg, weights = E(gg)$weight))
    d0_comm <- as.integer(names(sort(table(m[d0_samp]), decreasing = TRUE))[1])
    frac <- mean(m[d21m] == d0_comm)
    stab[[length(stab)+1]] <- data.frame(knn = k, seed = seed, frac_D21met_with_D0 = frac)
  }
  stab <- dplyr::bind_rows(stab)
  write.csv(stab, file.path(sup_dir, "S2_community_stability.csv"), row.names = FALSE)
  cat(sprintf("\n== S2 stability: mean fraction of Day21_Metformin samples in Day0's community = %.2f (over kNN 3-6, 10 seeds) ==\n",
              mean(stab$frac_D21met_with_D0)))
}, error = function(e) message("S2 failed: ", conditionMessage(e)))

###############################################################################
## S3. Permutation null for the metformin reversal fraction — R2 #3
##  Observed: within Day21, fit ~Treatment (MET vs CTR) per cell type; reversal
##  fraction = fraction of disease genes moved opposite by metformin. Null:
##  shuffle Treatment labels among Day21 samples (within cell type), refit,
##  recompute. Runtime scales with N_PERM x (#cell types) DESeq2 fits.
###############################################################################
tryCatch({
  N_PERM <- 200                     # <- runtime knob (200 ~ a few minutes)
  # disease gene set + sign, from the existing Day21-vs-Day0 contrast
  dis <- de_all %>% dplyr::filter(contrast == "disease_Day21_vs_Day0",
                                  !is.na(padj), padj < PADJ, abs(log2FoldChange) > LFC) %>%
    dplyr::select(celltype, gene, disease_lfc = log2FoldChange)

  # per-cell-type Day21 units with both arms >= MIN_REPS
  fit_met_day21 <- function(ct, treatment_vec = NULL){
    sel <- which(lut$celltype == ct & lut$Group == "Day21")
    cd  <- droplevels(lut[sel, ])
    if (!is.null(treatment_vec)) cd$Treatment <- treatment_vec
    cd$Treatment <- factor(as.character(cd$Treatment), levels = c("CTR","MET"))
    if (any(table(cd$Treatment) < MIN_REPS_PER_SIDE) || nlevels(cd$Treatment) < 2) return(NULL)
    dds <- fit_dds(pb[, sel, drop = FALSE], cd, ~ Treatment)
    if (is.null(dds)) return(NULL)
    r <- as.data.frame(DESeq2::results(dds, contrast = c("Treatment","MET","CTR")))
    data.frame(gene = rownames(r), met_lfc = r$log2FoldChange)
  }
  reversal_fraction <- function(met_by_ct){
    tot <- 0; rev <- 0
    for (ct in names(met_by_ct)) {
      m <- met_by_ct[[ct]]; d <- dis[dis$celltype == ct, ]
      j <- merge(d, m, by = "gene")
      j <- j[is.finite(j$met_lfc) & j$met_lfc != 0, ]
      tot <- tot + nrow(j); rev <- rev + sum(sign(j$met_lfc) != sign(j$disease_lfc))
    }
    if (tot == 0) NA_real_ else rev / tot
  }
  ct_test <- intersect(cts, unique(dis$celltype))
  obs_by_ct <- setNames(lapply(ct_test, fit_met_day21), ct_test)
  obs_by_ct <- obs_by_ct[!vapply(obs_by_ct, is.null, logical(1))]
  obs <- reversal_fraction(obs_by_ct)

  set.seed(1)
  null <- numeric(N_PERM)
  for (p in seq_len(N_PERM)) {
    perm_by_ct <- list()
    for (ct in names(obs_by_ct)) {
      sel <- which(lut$celltype == ct & lut$Group == "Day21")
      tv  <- sample(as.character(lut$Treatment[sel]))       # shuffle within cell type
      perm_by_ct[[ct]] <- fit_met_day21(ct, treatment_vec = tv)
    }
    perm_by_ct <- perm_by_ct[!vapply(perm_by_ct, is.null, logical(1))]
    null[p] <- reversal_fraction(perm_by_ct)
    if (p %% 25 == 0) cat(sprintf("  S3 permutation %d/%d\n", p, N_PERM))
  }
  null <- null[is.finite(null)]
  p_perm <- (1 + sum(null >= obs)) / (1 + length(null))
  write.csv(data.frame(permutation = seq_along(null), reversal_fraction = null),
            file.path(sup_dir, "S3_permutation_null.csv"), row.names = FALSE)
  cat(sprintf("\n== S3 reversal fraction: observed = %.3f;  null mean = %.3f;  p_perm = %.4f (N=%d) ==\n",
              obs, mean(null), p_perm, length(null)))
  ggplot(data.frame(x = null), aes(x)) +
    geom_histogram(bins = 30, fill = "grey70", colour = "white") +
    geom_vline(xintercept = obs, colour = "#C44E52", linewidth = 1) +
    labs(x = "Reversal fraction under permuted treatment labels", y = "count",
         title = sprintf("S3 Permutation null (observed = %.2f, p = %.3f)", obs, p_perm)) +
    theme_bw(12)
  ggsave(file.path(sup_dir, "S3_permutation_null.png"), width = 7, height = 4.5, dpi = 200)
}, error = function(e) message("S3 failed: ", conditionMessage(e)))

###############################################################################
## S4. QC / composition / UMAP / marker evidence / versions — R1 #2, R2 #1
##  Only what the ANNOTATED object can supply. Upstream provenance (platform,
##  pre-QC counts, QC thresholds, doublet method, normalization/integration,
##  reference genome) must be reported from the original processing pipeline.
###############################################################################
tryCatch({
  md <- ifnb@meta.data
  # nuclei per cell type x sample (post-QC) and x condition
  write.csv(as.data.frame.matrix(table(md$celltype, md$Sample)),
            file.path(sup_dir, "S4_nuclei_per_celltype_sample.csv"))
  write.csv(as.data.frame.matrix(table(md$celltype, md$Subgroup)),
            file.path(sup_dir, "S4_nuclei_per_celltype_condition.csv"))
  # per-cell QC summaries if present
  qc_cols <- intersect(c("nCount_RNA","nFeature_RNA","percent.mt"), colnames(md))
  if (length(qc_cols))
    write.csv(aggregate(md[qc_cols], list(Sample = md$Sample),
                        function(x) round(c(median = median(x), min = min(x), max = max(x)), 1)),
              file.path(sup_dir, "S4_qc_summary_by_sample.csv"), row.names = FALSE)

  # UMAPs (only if a umap reduction exists)
  if ("umap" %in% names(ifnb@reductions) || "umap" %in% Seurat::Reductions(ifnb)) {
    ggsave(file.path(sup_dir, "S4_umap_celltype.png"),
           Seurat::DimPlot(ifnb, group.by = "celltype", label = TRUE, repel = TRUE) + ggplot2::ggtitle("Cell type"),
           width = 8, height = 6, dpi = 200)
    ggsave(file.path(sup_dir, "S4_umap_condition.png"),
           Seurat::DimPlot(ifnb, group.by = "Subgroup") + ggplot2::ggtitle("Condition"),
           width = 8, height = 6, dpi = 200)
  } else message("S4: no 'umap' reduction in object — run RunUMAP() upstream to add S4 UMAPs.")

  # marker evidence DotPlot — EDIT this list to your canonical markers
  markers <- c("Tnnt2","Myh6","Actc1",          # cardiomyocytes
               "Pecam1","Cdh5","Vwf",            # endothelium
               "Col1a1","Dcn","Pdgfra",          # fibroblasts
               "Ptprc","Cd68","Lyz2","Cd3e","Cd79a", # immune (mac/T/B)
               "Rgs5","Pdgfrb","Kcnj8",          # pericytes
               "Msln","Wt1",                     # mesothelial
               "Npr3","Pln")                     # endocardial
  markers <- intersect(markers, rownames(ifnb))
  if (length(markers) >= 3) {
    ggsave(file.path(sup_dir, "S4_marker_dotplot.png"),
           Seurat::DotPlot(ifnb, features = markers, group.by = "celltype") +
             ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)),
           width = 11, height = 6, dpi = 200)
  }
  # package versions
  writeLines(capture.output(sessionInfo()), file.path(sup_dir, "S4_sessionInfo.txt"))
  cat("\n== S4 QC/UMAP/markers written. NOTE: report platform, reference genome,\n",
      "pre-QC counts, QC thresholds, doublet method, normalization/integration\n",
      "from the ORIGINAL pipeline — these are not recoverable from the annotated object. ==\n", sep="")
}, error = function(e) message("S4 failed: ", conditionMessage(e)))

###############################################################################
## S5. Segmentation validation & inflammation-score correlation — R2 #5
##  These need inputs you create outside R:
##   (a) manual_counts.csv : columns  image, region(inflamed/quiet), manual, quPath
##   (b) inflammation_scores.csv : columns  condition, inflammation_score  (blinded)
##  Place them in <out_dir>/supplementary/ and re-run; both steps auto-skip if absent.
###############################################################################
tryCatch({
  # (a) segmentation vs manual counts
  mc_path <- file.path(sup_dir, "manual_counts.csv")
  if (file.exists(mc_path)) {
    mc <- read.csv(mc_path)
    r  <- cor(mc$manual, mc$quPath)
    bias <- mean(mc$quPath - mc$manual)
    write.csv(data.frame(pearson_r = r, mean_bias = bias, n = nrow(mc)),
              file.path(sup_dir, "S5_segmentation_validation.csv"), row.names = FALSE)
    ggplot(mc, aes(manual, quPath, colour = region)) + geom_point(size = 2) +
      geom_abline(slope = 1, intercept = 0, linetype = 2) +
      labs(title = sprintf("S5 QuPath vs manual nuclei (r = %.2f)", r)) + theme_bw(12)
    ggsave(file.path(sup_dir, "S5_segmentation_validation.png"), width = 6, height = 5, dpi = 200)
    cat(sprintf("\n== S5a segmentation: r = %.2f, mean bias = %.1f ==\n", r, bias))
  } else message("S5a: manual_counts.csv not found — skipping segmentation validation.")

  # (b) inflammation score vs the two baseline-distances
  infl_path <- file.path(sup_dir, "inflammation_scores.csv")
  if (file.exists(infl_path) && file.exists(IMAGING_DIST_CSV)) {
    infl <- read.csv(infl_path)                                  # condition, inflammation_score
    img  <- read.csv(IMAGING_DIST_CSV)[, c("condition","dist_to_baseline")]
    m <- dplyr::inner_join(infl, img, by = "condition")
    ws_path <- file.path(out_dir, "cross_modal_integration.csv")
    if (file.exists(ws_path)) {
      ws <- read.csv(ws_path)[, c("condition","transcriptomic_dist")]
      m <- dplyr::left_join(m, ws, by = "condition")
    }
    rr <- c(spatial  = suppressWarnings(cor(m$inflammation_score, m$dist_to_baseline, method = "spearman")),
            transcriptomic = if ("transcriptomic_dist" %in% names(m))
              suppressWarnings(cor(m$inflammation_score, m$transcriptomic_dist, method = "spearman")) else NA)
    write.csv(m, file.path(sup_dir, "S5_inflammation_vs_distance.csv"), row.names = FALSE)
    cat("\n== S5b inflammation-score Spearman vs distances ==\n"); print(round(rr, 3))
  } else message("S5b: inflammation_scores.csv or imaging CSV not found — skipping severity anchor.")
}, error = function(e) message("S5 failed: ", conditionMessage(e)))

cat("\n==== SUPPLEMENTARY S1–S5 COMPLETE ====\nOutputs in: ", sup_dir, "\n", sep = "")
###############################################################################
