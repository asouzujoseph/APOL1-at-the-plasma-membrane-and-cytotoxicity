## =============================================================================
## 0. Libraries
## =============================================================================
library(DESeq2)
library(tidyverse)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)
library(VennDiagram)
library(enrichplot)
library(ReactomePA)
library(clusterProfiler)
library(GSVA)
library(AnnotationDbi)
library(dplyr)
library(limma)
library(msigdbr)
library(pheatmap)
library(RColorBrewer)
library(ggVennDiagram)
library(ggvenn)
library(GSEABase)

## =============================================================================
## 1. Paths and input
## =============================================================================
counts_file <- "C:/Users/User/Documents/Damola/counts/gene_counts.txt"
meta_file   <- "C:/Users/User/Documents/Damola/counts/samples.csv"
out_dir     <- "C:/Users/User/Documents/Damola/RNA_Seq2"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(out_dir)

## =============================================================================
## 2. Load counts and metadata
## =============================================================================
fc <- read.delim(counts_file, comment.char = "#", check.names = FALSE)
count_mat <- fc[, 7:ncol(fc)]
rownames(count_mat) <- fc$Geneid

colnames(count_mat) <- sub(".*/", "", colnames(count_mat))
colnames(count_mat) <- sub("\\.bam$", "", colnames(count_mat))

count_mat <- as.matrix(count_mat)
storage.mode(count_mat) <- "integer"
write.csv(count_mat, "gene_counts_formatted.csv")

meta <- read.csv(meta_file, stringsAsFactors = FALSE)
meta <- meta[match(colnames(count_mat), meta$sample), ]
rownames(meta) <- meta$sample
stopifnot(all(colnames(count_mat) == meta$sample))

## doxy only
doxy_samples <- meta %>% filter(treatment == "doxy") %>% select(sample, genotype, replicate)
count_doxy <- count_mat[, doxy_samples$sample]
stopifnot(all(colnames(count_doxy) == doxy_samples$sample))
write.csv(doxy_samples, file = "samples.csv", row.names = FALSE)

doxy_samples$genotype  <- factor(doxy_samples$genotype,  levels = c("EV", "G0", "G1", "G2"))
doxy_samples$replicate <- factor(doxy_samples$replicate)

dds <- DESeqDataSetFromMatrix(countData = count_doxy,colData  = doxy_samples, design   = ~ replicate + genotype)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
saveRDS(dds, "dds_full.rds")

## =============================================================================
## 3. Helper: annotate ENSEMBL → SYMBOL
## =============================================================================
annotate_results <- function(res) {
  res_df  <- as.data.frame(res) %>% rownames_to_column("gene_id")
  gene_ids <- gsub("\\..*", "", res_df$gene_id)
  
  annot <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = gene_ids,
    columns = "SYMBOL",
    keytype = "ENSEMBL"
  )
  
  res_annot <- res_df %>%
    mutate(ENSEMBL = gene_ids) %>%
    left_join(annot, by = "ENSEMBL")
  
  res_annot
}

## =============================================================================
## 4. Helper: volcano
## =============================================================================
plot_volcano <- function(res_annot, title, file) {
  volc <- res_annot %>%
    mutate(
      neglog10padj = -log10(padj),
      sig = case_when(
        padj < 0.05 & abs(log2FoldChange) >= 0.58 ~ "Significant",
        TRUE ~ "Not significant"
      )
    )
  
  df <- volc %>% filter(sig == "Significant")
  write.csv(df, paste0(title, "_filtered.csv"), row.names = FALSE)
  
  y_cap <- 50
  volc$neglog10padj_capped <- pmin(volc$neglog10padj, y_cap)
  
  label_genes <- volc %>%
    filter(sig == "Significant", !is.na(SYMBOL)) %>%
    arrange(padj) %>%
    slice(1:20)
  
  p <- ggplot(volc, aes(x = log2FoldChange, y = neglog10padj_capped)) +
    geom_point(aes(color = sig), alpha = 0.7, size = 1.8) +
    scale_color_manual(values = c("Significant" = "red", "Not significant" = "grey70")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
    geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "grey") +
    geom_text_repel(
      data = label_genes,
      aes(label = SYMBOL),
      size = 3,
      max.overlaps = Inf
    ) +
    labs(
      x = "log2 Fold Change",
      y = expression(-log[10]("adjusted p-value")),
      title = title,
      color = "Significance"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  
  ggsave(file, p, width = 6, height = 5, dpi = 300)
}

#ChatGPT revised for the plots
plot_volcano <- function(res_annot, title, file) {
  volc <- res_annot %>%
    mutate(
      neglog10padj = -log10(padj),
      sig = case_when(
        padj < 0.05 & log2FoldChange >= 0.58  ~ "Upregulated",
        padj < 0.05 & log2FoldChange <= -0.58 ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    )
  
  df <- volc %>% filter(sig != "Not significant")
  write.csv(df, paste0(title, "_filtered.csv"), row.names = FALSE)
  
  y_cap <- 50
  volc$neglog10padj_capped <- pmin(volc$neglog10padj, y_cap)
  
  label_genes <- volc %>%
    filter(sig != "Not significant", !is.na(SYMBOL)) %>%
    arrange(padj) %>%
    slice(1:20)
  
  p <- ggplot(volc, aes(x = log2FoldChange, y = neglog10padj_capped)) +
    geom_point(aes(color = sig), alpha = 0.7, size = 1.8) +
    scale_color_manual(values = c(
      "Upregulated" = "red",
      "Downregulated" = "blue",
      "Not significant" = "grey70"
    )) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
    geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "grey") +
    geom_text_repel(
      data = label_genes,
      aes(label = SYMBOL),
      size = 3,
      max.overlaps = Inf
    ) +
    labs(
      x = "log2 Fold Change",
      y = expression(-log[10]("adjusted p-value")),
      title = title,
      color = "Regulation"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  
  ggsave(file, p, width = 6, height = 5, dpi = 600, device = pdf)
}
## =============================================================================
## 5. Helper: PCA
## =============================================================================
plot_pca <- function(dds, title, file) {
  vsd <- vst(dds, blind = FALSE)
  pca <- prcomp(t(assay(vsd)))
  
  pc_df <- as.data.frame(pca$x) %>%
    rownames_to_column("sample")
  
  meta_df <- as.data.frame(colData(dds))
  meta_df$sample <- rownames(meta_df)
  pc_df <- left_join(pc_df, meta_df, by = "sample")
  
  percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
  
  p <- ggplot(pc_df, aes(
    x = PC1, y = PC2,
    color = replicate,
    shape = genotype
  )) +
    geom_point(size = 3, alpha = 0.9) +
    labs(
      x = paste0("PC1 (", percentVar[1], "%)"),
      y = paste0("PC2 (", percentVar[2], "%)"),
      color = "Replicate",
      shape = "Genotype",
      title = title
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  
  ggsave(file, p, width = 6, height = 5, dpi = 600, device = pdf)
}

## =============================================================================
## 6. Helper: heatmap of DEGs
## =============================================================================
plot_heatmap <- function(dds, res_annot, title, file, n_genes = 50) {
  top_sig <- res_annot %>%
    filter(!is.na(SYMBOL)) %>%
    filter(padj < 0.05, abs(log2FoldChange) >= 0.58) %>%
    head(n_genes)
  
  vsd <- vst(dds, blind = FALSE)
  vsd_mat <- assay(vsd)
  top_sig <- top_sig[top_sig$gene_id %in% rownames(vsd_mat), ]
  
  plot_mat <- vsd_mat[top_sig$gene_id, ]
  rownames(plot_mat) <- top_sig$SYMBOL
  
  anno_col <- as.data.frame(colData(dds)[, "genotype", drop = FALSE])
  anno_col$genotype <- droplevels(anno_col$genotype)
  
  full_colors <- c(
    "EV" = "green",
    "G0" = "black",
    "G1" = "red",
    "G2" = "blue"
  )
  present_genotypes <- levels(anno_col$genotype)
  anno_colors <- list(genotype = full_colors[present_genotypes])
  
  pdf(file, width = 800, height = 1000)  #, res = 120
  pheatmap(
    plot_mat,
    main = title,
    annotation_col = anno_col,
    annotation_colors = anno_colors,
    scale = "row",
    clustering_method = "ward.D2",
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    show_colnames = FALSE,
    fontsize_row = 8
  )
  dev.off()
}

## =============================================================================
## 7. Differential analysis
## =============================================================================
plot_pca(dds, "Doxycycline treated samples", "pca_doxy_only_0.58.pdf")

res_G0_vs_EV <- results(dds, contrast = c("genotype", "G0", "EV"))
res_G1_vs_EV <- results(dds, contrast = c("genotype", "G1", "EV"))
res_G2_vs_EV <- results(dds, contrast = c("genotype", "G2", "EV"))

res_G0_EV <- lfcShrink(dds, coef = "genotype_G0_vs_EV", res = res_G0_vs_EV, type = "apeglm")
res_G1_EV <- lfcShrink(dds, coef = "genotype_G1_vs_EV", res = res_G1_vs_EV, type = "apeglm")
res_G2_EV <- lfcShrink(dds, coef = "genotype_G2_vs_EV", res = res_G2_vs_EV, type = "apeglm")

res_G0_EV_annot <- annotate_results(res_G0_EV)
res_G1_EV_annot <- annotate_results(res_G1_EV)
res_G2_EV_annot <- annotate_results(res_G2_EV)

write.csv(res_G0_EV_annot, "doxy_G0_vs_EV_results_0.58.csv", row.names = FALSE)
write.csv(res_G1_EV_annot, "doxy_G1_vs_EV_results_0.58.csv", row.names = FALSE)
write.csv(res_G2_EV_annot, "doxy_G2_vs_EV_results_0.58.csv", row.names = FALSE)

plot_volcano(res_G0_EV_annot, "G0 vs EV", "volcano_doxy_G0_vs_EV_0.58.pdf")
plot_volcano(res_G1_EV_annot, "G1 vs EV", "volcano_doxy_G1_vs_EV_0.58.pdf")
plot_volcano(res_G2_EV_annot, "G2 vs EV", "volcano_doxy_G2_vs_EV_0.58.pdf")

## =============================================================================
## 8. Venn (SYMBOL-based)
## =============================================================================
sig_G0_sym <- res_G0_EV_annot %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 0.58) %>%
  pull(SYMBOL) %>%
  na.omit()

sig_G1_sym <- res_G1_EV_annot %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 0.58) %>%
  pull(SYMBOL) %>%
  na.omit()

sig_G2_sym <- res_G2_EV_annot %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 0.58) %>%
  pull(SYMBOL) %>%
  na.omit()

write.csv(sig_G0_sym, "list_of_genes_in_Venn_plot_for_G0_vs_EV.csv", row.names = FALSE)
write.csv(sig_G1_sym, "list_of_genes_in_Venn_plot_for_G1_vs_EV.csv", row.names = FALSE)
write.csv(sig_G2_sym, "list_of_genes_in_Venn_plot_for_G2_vs_EV.csv", row.names = FALSE)

gene_list <- list("G0 vs EV"  = unique(sig_G0_sym),
                  "G1 vs EV" = unique(sig_G1_sym),
                  "G2 vs EV" = unique(sig_G2_sym))

p_venn <- ggVennDiagram(gene_list, label_alpha = 0) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none") +
  labs(title = "Overlap of DEGs in the doxycycline treated samples",
       subtitle = "Significant genes")
ggsave("venn_doxy_treated.pdf", p_venn, width = 7, height = 6, dpi = 600, )

genotype_colors <- c("G0 vs EV" = "grey", "G1 vs EV" = "red", "G2 vs EV" = "blue")
p_ggvenn <- ggvenn(
  gene_list,
  fill_color   = genotype_colors,
  stroke_size  = 0.5,
  set_name_size = 5,
  text_size    = 5
) +
  ggtitle("Overlap of DEGs in the doxycycline treated samples")
ggsave("venn_doxy_treated_ggvenn.pdf", p_ggvenn, width = 7, height = 6, device = pdf)

## =============================================================================
## 9. Heatmaps per contrast
## =============================================================================
dds_G2_EV <- dds[, dds$genotype %in% c("G2", "EV")]
dds_G2_EV$genotype  <- droplevels(dds_G2_EV$genotype)
dds_G2_EV$replicate <- droplevels(dds_G2_EV$replicate)
plot_heatmap(dds_G2_EV, res_G2_EV_annot, "DEGs G2_vs_EV", "heatmap_G2_vs_EV.pdf")

dds_G1_EV <- dds[, dds$genotype %in% c("G1", "EV")]
dds_G1_EV$genotype  <- droplevels(dds_G1_EV$genotype)
dds_G1_EV$replicate <- droplevels(dds_G1_EV$replicate)
plot_heatmap(dds_G1_EV, res_G1_EV_annot, "DEGs G1_vs_EV", "heatmap_G1_vs_EV.pdf")

dds_G0_EV <- dds[, dds$genotype %in% c("G0", "EV")]
dds_G0_EV$genotype  <- droplevels(dds_G0_EV$genotype)
dds_G0_EV$replicate <- droplevels(dds_G0_EV$replicate)
plot_heatmap(dds_G0_EV, res_G0_EV_annot, "DEGs G0_vs_EV", "heatmap_G0_vs_EV.pdf")

## =============================================================================
## 10. Genotype progression heatmap (ENSEMBL-based)
## =============================================================================
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)
rownames(mat) <- sub("\\..*", "", rownames(mat))

res_G0_df <- as.data.frame(res_G0_EV_annot)
res_G1_df <- as.data.frame(res_G1_EV_annot)
res_G2_df <- as.data.frame(res_G2_EV_annot)

res_G0_df$gene_id <- res_G0_df$ENSEMBL
res_G1_df$gene_id <- res_G1_df$ENSEMBL
res_G2_df$gene_id <- res_G2_df$ENSEMBL

sig_G0_ens <- res_G0_df %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 0.58)
sig_G1_ens <- res_G1_df %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 0.58)
sig_G2_ens <- res_G2_df %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 0.58)

res_all <- bind_rows(sig_G0_ens, sig_G1_ens, sig_G2_ens)
top100 <- res_all %>% arrange(padj) %>% head(100)
top100_unique <- top100[!duplicated(top100$SYMBOL), ]

heatmap_mat <- mat[top100_unique$gene_id, ]
rownames(heatmap_mat) <- top100_unique$SYMBOL

desired_order <- c("G2", "G1", "G0", "EV")
sample_info <- as.data.frame(colData(dds))
ordered_samples <- rownames(sample_info)[order(factor(sample_info$genotype,
                                                      levels = desired_order))]
heatmap_mat <- heatmap_mat[, ordered_samples]

anno_col <- as.data.frame(colData(dds)[, "genotype", drop = FALSE])
anno_col <- anno_col[ordered_samples, , drop = FALSE]

full_colors <- c("EV" = "green", "G0" = "black", "G1" = "red", "G2" = "blue")
present_genotypes <- levels(anno_col$genotype)
anno_colors <- list(genotype = full_colors[present_genotypes])

pdf("genotype_progression_heatmap.pdf", width = 7.07, height = 7.87)
pheatmap(
  heatmap_mat,
  scale           = "row",
  cluster_cols    = FALSE,
  cluster_rows    = TRUE,
  annotation_col  = anno_col,
  annotation_colors = anno_colors,
  clustering_method = "ward.D2",
  color           = colorRampPalette(c("blue", "white", "red"))(100),
  show_colnames   = FALSE,
  fontsize_row    = 8,
#  main            = "Progression of DEGs in Doxycycline treated samples"
  main            = ""
)
dev.off()

## =============================================================================
## 11. Reactome GSVA (ENSEMBL-only)
## =============================================================================
reactome_df <- msigdbr(
  species       = "Homo sapiens",
  collection    = "C2",
  subcollection = "CP:REACTOME"
)
reactome_gs_ens <- split(reactome_df$ensembl_gene, reactome_df$gs_name)

vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)
rownames(mat) <- sub("\\..*", "", rownames(mat))

param    <- ssgseaParam(mat, reactome_gs_ens)
gsva_all <- GSVA::gsva(param)

plot_gsva_heatmap <- function(gsva_results, title, file) {
  pdf(file, width = 1000, height = 800) #, res = 120
  pathway_vars <- apply(gsva_results, 1, var)
  top_n <- min(30, nrow(gsva_results))
  top_pathways <- gsva_results[order(pathway_vars, decreasing = TRUE)[1:top_n], ]
  pheatmap(
    top_pathways,
    main  = title,
    scale = "row",
    color = colorRampPalette(c("blue", "white", "red"))(100),
    fontsize_row = 7
  )
  dev.off()
}
###GSVA = Gene set variation analysis
gsva_G0 <- gsva_all[, dds$genotype == "G0", drop = FALSE]
gsva_G1 <- gsva_all[, dds$genotype == "G1", drop = FALSE]
gsva_G2 <- gsva_all[, dds$genotype == "G2", drop = FALSE]
gsva_EV <- gsva_all[, dds$genotype == "EV", drop = FALSE]
## shows the activity level of pathways in each sample
plot_gsva_heatmap(gsva_EV, "Top Reactome Pathways: EV Activity", "gsva_heatmap_EV.pdf")
plot_gsva_heatmap(gsva_G2, "Top Reactome Pathways: G2 Activity", "gsva_heatmap_G2.pdf")
plot_gsva_heatmap(gsva_G1, "Top Reactome Pathways: G1 Activity", "gsva_heatmap_G1.pdf")
plot_gsva_heatmap(gsva_G0, "Top Reactome Pathways: G0 Activity", "gsva_heatmap_G0.pdf")

## =============================================================================
## 12. Enrichment analysis (GO BP/MF, Reactome, KEGG, GSEA)
## =============================================================================

## 12.1 ENSEMBL → ENTREZ converter (vector input)
convert_ids <- function(ensembl_ids) {
  ensembl_ids <- sub("\\..*", "", ensembl_ids)
  conv <- suppressWarnings(
    bitr(
      ensembl_ids,
      fromType = "ENSEMBL",
      toType   = "ENTREZID",
      OrgDb    = org.Hs.eg.db
    )
  )
  unique(na.omit(conv$ENTREZID))
}

## 12.2 Universe in ENTREZ
universe_ens <- rownames(mat)
valid_ens    <- universe_ens[universe_ens %in% keys(org.Hs.eg.db, keytype = "ENSEMBL")]
universe_entrez <- convert_ids(valid_ens)

## 12.3 DEG lists in ENSEMBL for ORA
sig_G0_ens_ids <- res_G0_EV_annot %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 0.58) %>%
  pull(ENSEMBL) %>%
  na.omit()

sig_G1_ens_ids <- res_G1_EV_annot %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 0.58) %>%
  pull(ENSEMBL) %>%
  na.omit()

sig_G2_ens_ids <- res_G2_EV_annot %>%
  filter(padj < 0.05, abs(log2FoldChange) >= 0.58) %>%
  pull(ENSEMBL) %>%
  na.omit()

entrez_G0 <- convert_ids(sig_G0_ens_ids)
entrez_G1 <- convert_ids(sig_G1_ens_ids)
entrez_G2 <- convert_ids(sig_G2_ens_ids)

## 12.4 GO BP
BF_G0 <- enrichGO(
  gene         = entrez_G0,
  universe     = universe_entrez,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",
  pvalueCutoff = 0.05,
  readable     = TRUE
)

BF_G1 <- enrichGO(
  gene         = entrez_G1,
  universe     = universe_entrez,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",
  pvalueCutoff = 0.05,
  readable     = TRUE
)

BF_G2 <- enrichGO(
  gene         = entrez_G2,
  universe     = universe_entrez,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",
  pvalueCutoff = 0.05,
  readable     = TRUE
)

dotplot(BF_G0, showCategory = 15) + ggtitle("GO BP – G0 vs EV")  ##few genes
p_BF1 <- dotplot(BF_G1, showCategory = 15) + ggtitle("GO BP – G1 vs EV") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.y  = element_text(size = 5, face = "bold"),
    axis.text.x  = element_text(size = 6),
    plot.title   = element_text(size = 10, hjust = 0.5),
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 7)
  )
ggsave("BF_G1_dotplot.pdf", p_BF1, width = 6, height = 5, dpi = 300)
p_BF2 <- dotplot(BF_G2, showCategory = 15) + ggtitle("GO BP – G2 vs EV") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.y  = element_text(size = 5, face = "bold"),
    axis.text.x  = element_text(size = 6),
    plot.title   = element_text(size = 10, hjust = 0.5),
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 7)
  )
ggsave("BF_G2_dotplot.pdf", p_BF2, width = 6, height = 5, dpi = 300)

## 12.5 GO MF
MF_G0 <- enrichGO(
  gene         = entrez_G0,
  universe     = universe_entrez,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "MF",
  pvalueCutoff = 0.05,
  readable     = TRUE
)

MF_G1 <- enrichGO(
  gene         = entrez_G1,
  universe     = universe_entrez,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "MF",
  pvalueCutoff = 0.05,
  readable     = TRUE
)

MF_G2 <- enrichGO(
  gene         = entrez_G2,
  universe     = universe_entrez,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  ont          = "MF",
  pvalueCutoff = 0.05,
  readable     = TRUE
)

dotplot(MF_G0, showCategory = 15) + ggtitle("GO MF – G0 vs EV")   ##few genes
p_MF1 <- dotplot(MF_G1, showCategory = 15) + ggtitle("GO MF – G1 vs EV") + 
  theme_bw(base_size = 8) +
  theme(
    axis.text.y  = element_text(size = 5, face = "bold"),
    axis.text.x  = element_text(size = 6),
    plot.title   = element_text(size = 10, hjust = 0.5),
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 7)
  )
ggsave("MF_G1_dotplot.pdf", p_MF1, width = 6, height = 5, dpi = 300)
p_MF2 <- dotplot(MF_G2, showCategory = 15) + ggtitle("GO MF – G2 vs EV") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.y  = element_text(size = 5, face = "bold"),
    axis.text.x  = element_text(size = 6),
    plot.title   = element_text(size = 10, hjust = 0.5),
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 7)
  )
ggsave("MF_G2_dotplot.pdf", p_MF2, width = 6, height = 5, dpi = 300)


## 12.8 GSEA (ReactomePA::gsePathway) using ENSEMBL IDs directly
## Convert ENSEMBL → ENTREZ for GSEA
convert_ids_df <- function(ensembl_ids) {
  ensembl_ids <- sub("\\..*", "", ensembl_ids)
  
  mapping <- suppressWarnings(
    bitr(
      ensembl_ids,
      fromType = "ENSEMBL",
      toType   = "ENTREZID",
      OrgDb    = org.Hs.eg.db
    )
  )
  
  mapping <- mapping[!is.na(mapping$ENTREZID), ]
  mapping <- mapping[!duplicated(mapping$ENSEMBL), ]
  
  mapping
}



## Prepare ranked vector for GSEA
prepare_ranks <- function(res_df) {
  
  res_df$ENSEMBL_clean <- sub("\\..*", "", res_df$gene_id)
  
  mapping <- convert_ids_df(res_df$ENSEMBL_clean)
  
  df2 <- res_df %>%
    inner_join(mapping, by = c("ENSEMBL_clean" = "ENSEMBL")) %>%
    distinct(ENTREZID, .keep_all = TRUE)
  
  ranks <- df2$stat
  names(ranks) <- df2$ENTREZID
  
  sort(ranks, decreasing = TRUE)
}

## Use unshrunken results for GSEA
res_G0_df <- as.data.frame(res_G0_vs_EV) %>% rownames_to_column("gene_id")
res_G1_df <- as.data.frame(res_G1_vs_EV) %>% rownames_to_column("gene_id")
res_G2_df <- as.data.frame(res_G2_vs_EV) %>% rownames_to_column("gene_id")

rank_G0 <- prepare_ranks(res_G0_df)
rank_G1 <- prepare_ranks(res_G1_df)
rank_G2 <- prepare_ranks(res_G2_df)

gsea_G0 <- gsePathway(rank_G0, organism = "human", pvalueCutoff = 0.25)
gsea_G1 <- gsePathway(rank_G1, organism = "human", pvalueCutoff = 0.25)
gsea_G2 <- gsePathway(rank_G2, organism = "human", pvalueCutoff = 0.25)

## Visualizaton
p_G0 <- ridgeplot(gsea_G0, showCategory = 20) +
  ggtitle("Reactome GSEA – G0 vs EV") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.y  = element_text(size = 3.5, face = "bold"),
    axis.text.x  = element_text(size = 6),
    plot.title   = element_text(size = 10, hjust = 0.5),
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 7)
  )
ggsave("GSEA_G0_ridgeplot.pdf", p_G0, width = 6, height = 5, dpi = 300)

p_G1 <- ridgeplot(gsea_G1, showCategory = 20) +
  ggtitle("Reactome GSEA – G1 vs EV") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.y  = element_text(size = 3.5, face = "bold"),
    axis.text.x  = element_text(size = 6),
    plot.title   = element_text(size = 10, hjust = 0.5),
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 7)
  )
ggsave("GSEA_G1_ridgeplot.pdf", p_G1, width = 6, height = 5, dpi = 300)

p_G2 <- ridgeplot(gsea_G2, showCategory = 20) +
  ggtitle("Reactome GSEA – G2 vs EV") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.y  = element_text(size = 3.5, face = "bold"),
    axis.text.x  = element_text(size = 6),
    plot.title   = element_text(size = 10, hjust = 0.5),
    legend.text  = element_text(size = 6),
    legend.title = element_text(size = 7)
  )
ggsave("GSEA_G2_ridgeplot.pdf", p_G2, width = 6, height = 5, dpi = 300)

