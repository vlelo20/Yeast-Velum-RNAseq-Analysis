## ============================================================
# BINF 6110 - Assignment 2
# Author: Vian Lelo Last edited: March 4th, 2026
# Description: Differential expression analysis of S. cerevisiae
#              flor yeast velum development (Early / Thin / Mature)
#              BioProject: PRJNA592304
# ============================================================

# ── 0. Package installation (run once) ──────────────────────
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "tximport", "GenomicFeatures",
                       "txdbmaker", "AnnotationDbi",
                       "clusterProfiler", "org.Sc.sgd.db"),
                     ask = FALSE, update = FALSE)

install.packages(c("tidyverse", "pheatmap", "RColorBrewer",
                   "ggrepel", "patchwork"), repos = "https://cloud.r-project.org")

# ── 1. Load libraries ────────────────────────────────────────
suppressPackageStartupMessages({
  library(DESeq2)
  library(tximport)
  library(GenomicFeatures)
  library(txdbmaker)
  library(AnnotationDbi)
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(patchwork)
  library(clusterProfiler)
  library(org.Sc.sgd.db)
})

# ── 2. Sample metadata ───────────────────────────────────────
# samplesheet.csv must have columns: sample, stage
# stage values: Early, Thin, Mature  (9 samples total)
samples <- read.csv("samplesheet.csv", stringsAsFactors = FALSE)

# Ensure stage is an ordered factor for sensible contrasts
samples$stage <- factor(samples$stage,
                        levels = c("Early", "Thin", "Mature"))

# Sanity check
stopifnot(all(c("sample", "stage") %in% colnames(samples)))
message("Samples loaded: ", nrow(samples))
print(table(samples$stage))

# ── 3. Transcript-to-gene mapping ───────────────────────────
txdb <- makeTxDbFromGFF(
  "reference/Saccharomyces_cerevisiae.R64-1-1.110.gtf.gz"
)

k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME") %>%
  dplyr::rename(tx = TXNAME, gene = GENEID) %>%
  dplyr::distinct()

message("tx2gene rows: ", nrow(tx2gene))
head(tx2gene)   # e.g. tx ~ "YAL001C_mRNA", gene ~ "YAL001C"

# ── 4. Salmon file paths ─────────────────────────────────────
files <- file.path("results/salmon_quant", samples$sample, "quant.sf")
names(files) <- samples$sample

missing <- files[!file.exists(files)]
if (length(missing) > 0) stop("Missing quant.sf files:\n", paste(missing, collapse = "\n"))
message("All quant.sf files found.")

# ── 5. tximport (gene-level, raw counts for DESeq2) ──────────
txi <- tximport(files,
                type                = "salmon",
                tx2gene             = tx2gene,
                countsFromAbundance = "no")   # use inferredLibraryType counts

# ── 6. DESeq2 dataset ────────────────────────────────────────
dds <- DESeqDataSetFromTximport(txi,
                                colData = samples %>% column_to_rownames("sample"),
                                design  = ~ stage)

# Pre-filter: keep genes with >= 10 counts in >= 3 samples
keep <- rowSums(counts(dds) >= 10) >= 3
dds  <- dds[keep, ]
message("Genes after filtering: ", nrow(dds))

# ── 7. Run DESeq2 ────────────────────────────────────────────
dds <- DESeq(dds)
resultsNames(dds)

# ── 8. Extract pairwise contrasts ───────────────────────────
res_thin_vs_early   <- results(dds, contrast = c("stage", "Thin",   "Early"), alpha = 0.05)
res_mature_vs_early <- results(dds, contrast = c("stage", "Mature", "Early"), alpha = 0.05)
res_mature_vs_thin  <- results(dds, contrast = c("stage", "Mature", "Thin"),  alpha = 0.05)

# Apply lfcShrink for better LFC estimates (apeglm requires coef, use ashr for contrasts)
BiocManager::install("ashr", ask = FALSE, update = FALSE)
res_thin_vs_early   <- lfcShrink(dds, contrast = c("stage","Thin",  "Early"),
                                 res = res_thin_vs_early,   type = "ashr")
res_mature_vs_early <- lfcShrink(dds, contrast = c("stage","Mature","Early"),
                                 res = res_mature_vs_early, type = "ashr")
res_mature_vs_thin  <- lfcShrink(dds, contrast = c("stage","Mature","Thin"),
                                 res = res_mature_vs_thin,  type = "ashr")

# Summary
message("\n── Thin vs Early ──");   summary(res_thin_vs_early,   alpha = 0.05)
message("\n── Mature vs Early ──"); summary(res_mature_vs_early, alpha = 0.05)
message("\n── Mature vs Thin ──");  summary(res_mature_vs_thin,  alpha = 0.05)

# ── 9. Save DE tables ────────────────────────────────────────
dir.create("results/DE", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

write.csv(as.data.frame(res_thin_vs_early),
          "results/DE/DE_thin_vs_early.csv",   row.names = TRUE)
write.csv(as.data.frame(res_mature_vs_early),
          "results/DE/DE_mature_vs_early.csv", row.names = TRUE)
write.csv(as.data.frame(res_mature_vs_thin),
          "results/DE/DE_mature_vs_thin.csv",  row.names = TRUE)

# ── 10. Variance-stabilising transform ──────────────────────
vsd <- vst(dds, blind = FALSE)

# ── 11. PCA plot ─────────────────────────────────────────────
stage_palette <- c("Early" = "#4E79A7", "Thin" = "#F28E2B", "Mature" = "#59A14F")

pca_data   <- plotPCA(vsd, intgroup = "stage", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = stage, label = name)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = stage_palette, name = "Stage") +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("Principal Component Analysis",
          subtitle = "VST-normalised counts — coloured by velum stage") +
  coord_fixed() +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4)
  )

ggsave("results/figures/PCA_stage.png", p_pca,
       width = 6, height = 5, dpi = 300, bg = "white")

# ── 12. Sample distance heatmap ──────────────────────────────
sampleDists       <- dist(t(assay(vsd)))
sampleDistMatrix  <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- colnames(vsd)

annotation_col <- data.frame(Stage = vsd$stage)
rownames(annotation_col) <- colnames(vsd)

ann_colors <- list(Stage = stage_palette)

pheatmap(sampleDistMatrix,
         annotation_col           = annotation_col,
         annotation_colors        = ann_colors,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         color                    = colorRampPalette(rev(brewer.pal(9, "Blues")))(100),
         border_color             = NA,
         show_rownames            = TRUE,
         show_colnames            = FALSE,
         fontsize                 = 10,
         main                     = "Sample-to-Sample Euclidean Distances (VST)",
         filename                 = "results/figures/sample_distance_heatmap.png",
         width = 7, height = 6)

# ── 13. Volcano plots — all three contrasts ──────────────────

# Helper function: build a tidy DEG data frame from a results object
make_res_df <- function(res_obj) {
  as.data.frame(res_obj) %>%
    rownames_to_column("gene") %>%
    filter(!is.na(padj)) %>%
    mutate(sig = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE                               ~ "NS"
    ))
}

# Helper function: build one volcano ggplot object
make_volcano <- function(res_df, title_str) {
  top_labels <- res_df %>%
    filter(sig != "NS") %>%
    group_by(sig) %>%
    slice_min(padj, n = 10) %>%
    ungroup()
  
  ggplot(res_df, aes(log2FoldChange, -log10(padj), color = sig)) +
    geom_point(alpha = 0.45, size = 1.1) +
    geom_text_repel(data = top_labels, aes(label = gene),
                    size = 2.6, max.overlaps = 20, show.legend = FALSE,
                    box.padding = 0.4, segment.color = "grey55",
                    segment.size = 0.3) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed",
               color = "grey40", linewidth = 0.45) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey40", linewidth = 0.45) +
    scale_color_manual(
      values = c("Up" = "#D62728", "Down" = "#1F77B4", "NS" = "grey75"),
      labels = c("Up"   = paste0("Up (n=",   sum(res_df$sig == "Up"),   ")"),
                 "Down" = paste0("Down (n=", sum(res_df$sig == "Down"), ")"),
                 "NS"   = "NS"),
      name = NULL
    ) +
    labs(title    = title_str,
         subtitle = paste0("|LFC| > 1, FDR < 0.05"),
         x = expression(log[2]~fold~change),
         y = expression(-log[10]~adjusted~italic(p))) +
    theme_classic(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(size = 8, color = "grey45"),
      legend.position  = c(0.82, 0.12),
      legend.background = element_rect(fill = "white", color = "grey85", linewidth = 0.3),
      legend.text      = element_text(size = 8),
      panel.grid.major = element_line(color = "grey93", linewidth = 0.35)
    )
}

# ── Build data frames for each contrast ──────────────────────
df_thin_vs_early   <- make_res_df(res_thin_vs_early)
df_mature_vs_early <- make_res_df(res_mature_vs_early)
df_mature_vs_thin  <- make_res_df(res_mature_vs_thin)

# ── Individual volcano plots ──────────────────────────────────
p_v1 <- make_volcano(df_thin_vs_early,   "Thin vs Early")
p_v2 <- make_volcano(df_mature_vs_early, "Mature vs Early")
p_v3 <- make_volcano(df_mature_vs_thin,  "Mature vs Thin")

ggsave("results/figures/volcano_thin_vs_early.png",   p_v1,
       width = 6, height = 5, dpi = 300, bg = "white")
ggsave("results/figures/volcano_mature_vs_early.png", p_v2,
       width = 6, height = 5, dpi = 300, bg = "white")
ggsave("results/figures/volcano_mature_vs_thin.png",  p_v3,
       width = 6, height = 5, dpi = 300, bg = "white")

# ── Combined 3-panel figure ───────────────────────────────────
p_volcano_panel <- p_v1 + p_v2 + p_v3 +
  plot_layout(ncol = 3) +
  plot_annotation(
    title   = "Differential Expression Across Velum Development Stages",
    subtitle = "Flor yeast S. cerevisiae — DESeq2 with ashr LFC shrinkage",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5)
    )
  )

ggsave("results/figures/volcano_all_contrasts.png", p_volcano_panel,
       width = 16, height = 5.5, dpi = 300, bg = "white")

message("✓ Volcano plots saved (3 individual + 1 combined panel)")


# ── 14. Top-50 DEG heatmap (Mature vs Early) ─────────────────
top50_genes <- res_df %>%
  filter(sig != "NS") %>%
  arrange(padj) %>%
  slice(1:50) %>%
  pull(gene)

mat <- assay(vsd)[top50_genes, ]
mat <- mat - rowMeans(mat)   # Centre rows

annotation_hm <- data.frame(Stage = vsd$stage)
rownames(annotation_hm) <- colnames(mat)

pheatmap(mat,
         annotation_col  = annotation_hm,
         annotation_colors = ann_colors,
         color           = colorRampPalette(c("#2166AC","white","#D6604D"))(100),
         border_color    = NA,
         show_rownames   = TRUE,
         show_colnames   = TRUE,
         fontsize_row    = 7,
         fontsize_col    = 9,
         scale           = "none",
         clustering_method = "ward.D2",
         main            = "Top 50 DEGs: Mature vs Early (row-centred VST)",
         filename        = "results/figures/heatmap_top50_mature_vs_early.png",
         width = 8, height = 9)

# ── 15. GO Biological Process ORA ───────────────────────────
sig_genes <- res_df %>% filter(sig != "NS") %>% pull(gene)

ego <- enrichGO(gene          = sig_genes,
                OrgDb         = org.Sc.sgd.db,
                keyType       = "ORF",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = FALSE)

write.csv(as.data.frame(ego),
          "results/DE/GO_BP_ORA_mature_vs_early.csv", row.names = FALSE)

p_go <- dotplot(ego, showCategory = 20, font.size = 10) +
  scale_color_gradientn(colors = c("#D62728", "#FF7F0E", "#AEC7E8"),
                        name   = "Adjusted\np-value") +
  labs(title    = "GO Biological Process Enrichment",
       subtitle = "Over-representation analysis: Mature vs Early") +
  theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    axis.text.y   = element_text(size = 9)
  )

ggsave("results/figures/GO_ORA_mature_vs_early.png", p_go,
       width = 8, height = 7, dpi = 300, bg = "white")

message("\n✓ Analysis complete. Outputs in results/DE/ and results/figures/")




# ── 16. GO ORA — all three contrasts + combined panel ────────

run_ORA <- function(res_df, label) {
  genes <- res_df %>% filter(sig != "NS") %>% pull(gene)
  ego   <- enrichGO(gene          = genes,
                    OrgDb         = org.Sc.sgd.db,
                    keyType       = "ORF",
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.2,
                    readable      = FALSE)
  write.csv(as.data.frame(ego),
            paste0("results/DE/GO_BP_ORA_", label, ".csv"),
            row.names = FALSE)
  ego
}

ego_tve <- run_ORA(df_thin_vs_early,   "thin_vs_early")
ego_mve <- run_ORA(df_mature_vs_early, "mature_vs_early")
ego_mvt <- run_ORA(df_mature_vs_thin,  "mature_vs_thin")

make_go_plot <- function(ego_obj, subtitle_str, n = 15) {
  dotplot(ego_obj, showCategory = n, font.size = 9) +
    scale_color_gradientn(
      colors = c("#D62728", "#FF7F0E", "#AEC7E8"),
      name   = "adj. p"
    ) +
    labs(title = subtitle_str) +
    theme_classic(base_size = 10) +
    theme(
      plot.title      = element_text(face = "bold", size = 11),
      axis.text.y     = element_text(size = 8),
      legend.key.size = unit(0.4, "cm"),
      legend.position = "right",
      legend.text     = element_text(size = 7),
      legend.title    = element_text(size = 8)
    )
}

p_go1 <- make_go_plot(ego_tve, "Thin vs Early")
p_go2 <- make_go_plot(ego_mve, "Mature vs Early")
p_go3 <- make_go_plot(ego_mvt, "Mature vs Thin")

# Save individuals
ggsave("results/figures/GO_ORA_thin_vs_early.png",   p_go1,
       width = 7, height = 6, dpi = 300, bg = "white")
ggsave("results/figures/GO_ORA_mature_vs_thin.png",  p_go3,
       width = 7, height = 6, dpi = 300, bg = "white")

# Combined 3-panel
p_go_panel <- p_go1 + p_go2 + p_go3 +
  plot_layout(ncol = 3) +
  plot_annotation(
    title    = "GO Biological Process Enrichment Across Velum Development Contrasts",
    subtitle = "ORA (clusterProfiler): top 15 BP terms per contrast, BH FDR < 0.05",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 9, color = "grey40", hjust = 0.5)
    )
  )

ggsave("results/figures/GO_ORA_all_contrasts.png", p_go_panel,
       width = 22, height = 8, dpi = 300, bg = "white")

ggsave("results/figures/GO_ORA_all_contrasts.png", p_go_panel,
       width = 20, height = 8, dpi = 300, bg = "white")

message("✓ GO ORA complete — individual + combined panel saved.")

# ── 17. Expression trajectory plot — selected biologically justified genes ────
#
# Gene selection rationale:
#   FLO11  – master regulator of velum adhesion; abolishes biofilm when deleted
#            (Brückner & Mösch, 2012; Mardanov et al., 2020)
#   CYC1   – iso-1 cytochrome c; O2/HAP-activated electron carrier; mechanistic
#            link to "aerobic respiration" GO enrichment (Guarente et al., 1983)
#   QCR2   – Complex III subunit; qcr2Δ impairs growth on non-fermentable carbon
#            and reduces ethanol tolerance (SGD); links to ETC GO terms
#   ALD4   – mitochondrial ALDH; induced by ethanol in flor yeasts to channel
#            acetaldehyde into TCA; directly relevant to wine acetaldehyde
#            accumulation (Lasanta et al., 2003; Navarro-Aviño et al., 1999)
#   ALD6   – cytosolic ALDH; major acetaldehyde-to-acetate route with ALD4
#            (Navarro-Aviño et al., 1999)

genes_of_interest <- c("FLO11", "CYC1", "QCR2", "ALD4", "ALD6")

# Check which are present in the VST matrix (gene IDs must match your annotation)
# If your IDs are ORF names (e.g. YOR374W for ALD4), use those instead:
orf_map <- c(
  FLO11 = "YHR211W",
  CYC1  = "YJR048W",
  QCR2  = "YPR191W",
  ALD4  = "YOR374W",
  ALD6  = "YPL061W"
)

# Filter to genes actually present in the dataset
present <- orf_map[orf_map %in% rownames(assay(vsd))]
missing <- orf_map[!orf_map %in% rownames(assay(vsd))]
if (length(missing) > 0) {
  message("Note: these ORFs not found in VST matrix: ", paste(names(missing), collapse = ", "))
}

# Extract VST counts for selected genes
vst_mat <- assay(vsd)[present, ]

# Build tidy long-format data frame
traj_df <- as.data.frame(t(vst_mat)) %>%
  rownames_to_column("sample") %>%
  left_join(
    as.data.frame(colData(vsd)) %>% rownames_to_column("sample") %>% select(sample, stage),
    by = "sample"
  ) %>%
  pivot_longer(cols = -c(sample, stage), names_to = "orf", values_to = "vst") %>%
  mutate(
    gene  = names(present)[match(orf, present)],  # map ORF back to gene name
    stage = factor(stage, levels = c("Early", "Thin", "Mature"))
  )

# Compute mean + SE per gene x stage for the line
traj_summary <- traj_df %>%
  group_by(gene, stage) %>%
  summarise(mean_vst = mean(vst),
            se_vst   = sd(vst) / sqrt(n()),
            .groups  = "drop")

# Colour palette: group by biological role
gene_colors <- c(
  FLO11 = "#8C4B8C",   # purple  — adhesion
  CYC1  = "#E07B39",   # orange  — respiratory chain
  QCR2  = "#D4A537",   # amber   — respiratory chain
  ALD4  = "#3A7DC9",   # blue    — aldehyde metabolism
  ALD6  = "#4EAD81"    # green   — aldehyde metabolism
)

gene_labels <- c(
  FLO11 = "FLO11 — cell adhesion",
  CYC1  = "CYC1 — cytochrome c",
  QCR2  = "QCR2 — Complex III",
  ALD4  = "ALD4 — mito. ALDH",
  ALD6  = "ALD6 — cyto. ALDH"
)

p_traj <- ggplot(traj_summary, aes(x = stage, y = mean_vst,
                                   group = gene, color = gene)) +
  # Individual replicate points (jittered slightly)
  geom_jitter(data = traj_df,
              aes(x = stage, y = vst, color = gene),
              width = 0.06, size = 2, alpha = 0.45, inherit.aes = FALSE) +
  # Mean line
  geom_line(linewidth = 1.1) +
  # Mean point
  geom_point(aes(y = mean_vst), size = 3.5, shape = 21,
             fill = "white", stroke = 1.4) +
  # SE ribbon
  geom_ribbon(aes(ymin = mean_vst - se_vst,
                  ymax = mean_vst + se_vst,
                  fill = gene),
              alpha = 0.12, color = NA) +
  scale_color_manual(values = gene_colors, labels = gene_labels, name = "Gene") +
  scale_fill_manual(values  = gene_colors, labels = gene_labels, name = "Gene") +
  facet_wrap(~ gene, scales = "free_y", ncol = 5,
             labeller = labeller(gene = gene_labels)) +
  labs(
    title    = "Expression Trajectories of Biologically Selected Genes",
    subtitle = "VST-normalised counts — mean ± SE with individual replicates",
    x        = "Velum stage",
    y        = "VST-normalised expression"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 9, color = "grey40"),
    strip.text      = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "grey95", color = NA),
    legend.position = "none",   # facet labels are sufficient
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.4),
    axis.text.x     = element_text(size = 10)
  )

ggsave("results/figures/gene_trajectories.png", p_traj,
       width = 14, height = 4, dpi = 300, bg = "white")

message("✓ Gene trajectory plot saved: results/figures/gene_trajectories.png")


# ── 18. Shared DEGs across contrasts — UpSet plot ────────────────────────────
library(UpSetR)

# Build named lists of significant DEG IDs per contrast
sig_tve <- df_thin_vs_early   %>% filter(sig != "NS") %>% pull(gene)
sig_mve <- df_mature_vs_early %>% filter(sig != "NS") %>% pull(gene)
sig_mvt <- df_mature_vs_thin  %>% filter(sig != "NS") %>% pull(gene)

upset_list <- list(
  "Thin vs Early"   = sig_tve,
  "Mature vs Early" = sig_mve,
  "Mature vs Thin"  = sig_mvt
)

# UpSet plot
png("results/figures/DEG_upset.png",
    width = 2400, height = 1800, res = 300)

upset(
  fromList(upset_list),
  sets            = c("Thin vs Early", "Mature vs Early", "Mature vs Thin"),
  order.by        = "freq",
  keep.order      = TRUE,
  sets.bar.color  = c("#4E9A8C", "#D62728", "#3A7DC9"),
  main.bar.color  = "grey30",
  text.scale      = c(1.4, 1.2, 1.2, 1.0, 1.3, 1.1),
  mb.ratio        = c(0.6, 0.4),
  point.size      = 3.5,
  line.size       = 0.8
)

dev.off()

message("✓ UpSet plot saved: results/figures/DEG_upset.png")

# Quick count of intersections
length(intersect(intersect(sig_tve, sig_mve), sig_mvt))  # all three
length(intersect(sig_mve, sig_mvt))                       # Mature vs Early + Mature vs Thin only
length(intersect(sig_tve, sig_mve))                       # Thin vs Early + Mature vs Early only


# ── Session & package version info ───────────────────────────
session_info <- sessionInfo()

# Key packages to report
pkgs <- c("DESeq2", "tximport", "txdbmaker", "clusterProfiler", 
          "org.Sc.sgd.db", "ggplot2", "patchwork", "pheatmap", 
          "ggrepel", "ashr", "UpSetR", "dplyr", "tidyr", 
          "tibble", "BiocGenerics")

# Extract versions
pkg_versions <- sapply(pkgs, function(p) {
  if (p %in% names(session_info$otherPkgs)) {
    session_info$otherPkgs[[p]]$Version
  } else if (p %in% names(session_info$loadedOnly)) {
    session_info$loadedOnly[[p]]$Version
  } else {
    "not loaded"
  }
})

# Print cleanly
cat("R version:", R.version$major, ".", R.version$minor, "\n")
cat("\nPackage versions:\n")
for (pkg in names(pkg_versions)) {
  cat(sprintf("  %-20s %s\n", pkg, pkg_versions[pkg]))
}


# Rank all genes by Wald statistic from Mature vs Early
ranked_genes <- res_mature_vs_early %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  filter(!is.na(stat)) %>%
  arrange(desc(stat)) %>%
  { setNames(.$stat, .$gene) }

gsea_res <- gseGO(geneList     = ranked_genes,
                  OrgDb        = org.Sc.sgd.db,
                  keyType      = "ORF",
                  ont          = "BP",
                  minGSSize    = 10,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  nPermSimple  = 1000,
                  verbose      = FALSE)

write.csv(as.data.frame(gsea_res),
          "results/DE/GSEA_BP_mature_vs_early.csv", row.names = FALSE)

library(enrichplot)
p_gsea <- dotplot(gsea_res, showCategory = 15, split = ".sign") +
  facet_grid(. ~ .sign)
ggsave("results/figures/GSEA_BP_mature_vs_early.png", p_gsea,
       width = 12, height = 7, dpi = 300, bg = "white")