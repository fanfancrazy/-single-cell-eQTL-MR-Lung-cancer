# ============================================================
# Single-cell RNA-seq analysis pipeline with:
# Seurat preprocessing, visualization,
# scTenifoldKnk gene knockout simulation,
# differential regulation analysis,
# and functional enrichment (GO / KEGG / Reactome / GSEA)
# ============================================================

# ========================
# 1. Load required packages
# ========================

library(Seurat)
library(scTenifoldKnk)
library(Matrix)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(openxlsx)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(forcats)
library(enrichplot)
library(parallel)
library(ReactomePA)
library(msigdbr)

# ========================
# 2. System settings
# ========================

available_cores <- detectCores() - 1

# ========================
# 3. Working directory
# ========================

setwd("E:/Bio_analysis/SClungcancer")

dir.create("alltissue", showWarnings = FALSE)

# ========================
# 4. Plot saving function
# ========================

save_plot <- function(plot, filename, w = 6, h = 5){
  ggsave(paste0("alltissue/", filename, ".pdf"), plot,
         width = w, height = h, device = "pdf")
  
  ggsave(paste0("alltissue/", filename, ".jpg"), plot,
         width = w, height = h, dpi = 600)
}

# ========================
# 5. Load Seurat object
# ========================

seurat_obj <- readRDS(
  "GSE162498_NSCLC_CD3_4tumors_4Juxta_2Juxta.Rds"
)

seurat_obj <- UpdateSeuratObject(seurat_obj)
DefaultAssay(seurat_obj) <- "RNA"

# ========================
# 6. Cell annotation
# ========================

seurat_obj$celltype_major <- case_when(
  
  seurat_obj$predicted.id %in% c(
    "CD4 CD69", "CD4 GZMA", "CD4 memory", "CD4 Naive",
    "Stressed CD4 memory", "Tfh 1", "Tfh 2"
  ) ~ "CD4 T",
  
  seurat_obj$predicted.id %in% c(
    "CD8 effectors", "CD8 GZMH", "CD8 GZMK 1",
    "CD8 GZMK 2", "CD8 LAYN", "CD8 ZNF683"
  ) ~ "CD8 T",
  
  seurat_obj$predicted.id %in% c(
    "Tregs 1", "Tregs 2", "Tregs 3"
  ) ~ "Treg",
  
  seurat_obj$predicted.id %in% c("GDT cells", "MAIT") ~ "Innate T",
  
  seurat_obj$predicted.id == "IFN" ~ "IFN",
  
  seurat_obj$predicted.id %in% c("Cycling 1", "Cycling 2") ~ "Cycling",
  
  TRUE ~ "Other"
)

# ========================
# 7. Basic visualization (UMAP / Feature / Violin)
# ========================

p_umap <- DimPlot(
  seurat_obj,
  group.by = "celltype_major",
  label = TRUE,
  repel = TRUE
)

save_plot(p_umap, "UMAP_celltype_major")

p_feature <- FeaturePlot(
  seurat_obj,
  features = "RNASET2"
)

save_plot(p_feature, "Feature_RNASET2")

p_vln <- VlnPlot(
  seurat_obj,
  features = "RNASET2",
  group.by = "celltype_major"
)

save_plot(p_vln, "Vln_RNASET2")

# ========================
# 8. Blood subset analysis
# ========================

target_cells <- subset(
  seurat_obj,
  subset = tissue == "blood" & celltype_major == "CD4 T"
)

# Rebuild Seurat object to avoid assay errors
counts_data <- GetAssayData(target_cells, assay = "RNA", layer = "counts")
meta_data <- target_cells@meta.data

rm(target_cells)
gc()

target_cells <- CreateSeuratObject(
  counts = counts_data,
  meta.data = meta_data
)

# ========================
# 9. HVG selection
# ========================

target_cells <- NormalizeData(target_cells)
target_cells <- FindVariableFeatures(target_cells, nfeatures = 5000)

hvg_genes <- VariableFeatures(target_cells)

# ========================
# 10. scTenifoldKnk analysis
# ========================

target_gene <- "RNASET2"

counts_matrix <- GetAssayData(target_cells, assay = "RNA", layer = "counts")

genes_to_keep <- rowSums(counts_matrix > 0) >= 0
filtered_counts <- counts_matrix[genes_to_keep, ]

final_gene_list <- intersect(hvg_genes, rownames(filtered_counts))

if (!target_gene %in% rownames(filtered_counts)) {
  stop("Target gene not found in filtered matrix.")
}

if (!target_gene %in% final_gene_list) {
  final_gene_list <- c(final_gene_list, target_gene)
}

reduced_counts <- filtered_counts[final_gene_list, ]
reduced_counts <- as(reduced_counts, "dgCMatrix")

set.seed(123)

ko_result <- scTenifoldKnk(
  countMatrix = reduced_counts,
  gKO = target_gene,
  nc_nNet = 10,
  nc_nCells = 500,
  td_K = 5,
  nCores = max(1, floor(available_cores * 0.8))
)

saveRDS(
  ko_result,
  paste0("KO_Result_", target_gene, ".rds")
)

# ========================
# 11. Differential regulation analysis
# ========================

base_dir <- paste0(target_gene, "_KO_results")

dir.create(base_dir, showWarnings = FALSE)
dir.create(file.path(base_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(base_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(base_dir, "GSEA"), showWarnings = FALSE)

diff_reg <- ko_result$diffRegulation

sig_genes <- subset(diff_reg, p.adj < 0.05)

write.csv(sig_genes,
          file.path(base_dir, "tables", "sig_genes.csv"),
          row.names = FALSE)

write.csv(diff_reg,
          file.path(base_dir, "tables", "all_genes.csv"),
          row.names = FALSE)

# ========================
# 12. Volcano plot
# ========================

diff_reg$log2FC <- log2(diff_reg$FC + 1e-6)
diff_reg$change <- "NOT"
diff_reg$change[diff_reg$p.adj < 0.05 & diff_reg$log2FC > 0.25] <- "UP"
diff_reg$change[diff_reg$p.adj < 0.05 & diff_reg$log2FC < -0.25] <- "DOWN"

top_genes <- head(diff_reg[diff_reg$change != "NOT", ], 3)$gene

pdf(file.path(base_dir, "plots", "volcano.pdf"), 8, 7)

ggplot(diff_reg, aes(log2FC, -log10(p.adj))) +
  geom_point(aes(color = change), alpha = 0.7) +
  geom_text_repel(data = subset(diff_reg, gene %in% top_genes),
                  aes(label = gene)) +
  theme_classic()

dev.off()

# ========================
# 13. Functional enrichment analysis
# ========================

entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = toupper(sig_genes$gene),
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

entrez_ids <- na.omit(unique(entrez_ids))

go <- enrichGO(entrez_ids, org.Hs.eg.db, ont = "BP")
kegg <- enrichKEGG(entrez_ids, organism = "hsa")
reactome <- enrichPathway(entrez_ids, organism = "human")

msig_h <- msigdbr(
  species = "Homo sapiens",
  category = "H"
) %>%
  dplyr::select(gs_name, entrez_gene)

hallmark <- enricher(entrez_ids, TERM2GENE = msig_h)

# ========================
# 14. GSEA
# ========================

gene_rank <- diff_reg$Z
names(gene_rank) <- toupper(diff_reg$gene)
gene_rank <- sort(gene_rank, decreasing = TRUE)

gsea_res <- GSEA(
  gene_rank,
  TERM2GENE = msig_h,
  pvalueCutoff = 0.3
)

cat("Analysis completed.\n")