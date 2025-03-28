# FindMarkers, calculate All differential expression genes between Inulin_AD and Ctrl_AD in each cell type in interbrain
# run GSEA: GO-BP, REACTOME, KEGG

remove(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(DoubletFinder)
library(sctransform)
if(!require(AnnotationDbi))BiocManager::install('AnnotationDbi')
if(!require(org.Mm.eg.db))BiocManager::install('org.Mm.eg.db')
if(!require(clusterProfiler))BiocManager::install('clusterProfiler')
if(!require(tidyr))install.packages('tidyr', update = F, ask = F)
if(!require(tidyverse))install.packages('tidyverse', update = F, ask = F)
if(!require(ggrepel))install.packages('ggrepel', update = F, ask = F)
if(!require(hdf5r))install.packages('hdf5r', update = F, ask = F)
if(!require(pheatmap))install.packages('pheatmap', update = F, ask = F)
if(!require(msigdbr))install.packages('msigdbr')

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(qs))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(msigdbr))
suppressPackageStartupMessages(library(GSVA))

# suppressPackageStartupMessages(ReactomePA)


suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(stringr))


opt <- list(opd = "/home/FAD_IB_Integrate6/")

Obj_Integrate <- readRDS(paste0(opt$opd, "FAD_IB_Integrated6_res0.1_2023-10-20.rds"))

gg_title <- "interbrain"

Idents(Obj_Integrate) <- "seurat_clusters"
Obj_Integrate <- RenameIdents(Obj_Integrate, '0' = "Oligo", '1' = "Neuron",
                              '2' =  "Neuron", '3' = "Astrocyte", '4' = "Microglia", '5' = "Oligo", '6' = "OPC",
                              '7' = "Neuron", '8' = "Neuron", '9' = "Neuron", '10' = "Pericyte",
                              '11' = "Unknown", '12' = "Unknown")
Obj_Integrate$celltype <- Idents(Obj_Integrate)
Idents(Obj_Integrate) <- "celltype"

Obj_Integrate$celltype.group <- paste(Obj_Integrate$celltype, Obj_Integrate$group, sep = "_")
Idents(Obj_Integrate) <- "celltype.group"

cellfordeg <- unique(Obj_Integrate$celltype)
Obj_Integrate$celltype <- factor(x = Obj_Integrate$celltype,
                                 levels = c("Neuron", "Astrocyte", "Microglia", "Oligo", "OPC", "Pericyte", "Unknown"))

cellfordeg <- c("Neuron", "Astrocyte", "Microglia", "Oligo", "OPC", "Pericyte", "Unknown")

# -------------------------------------------------------------------------------------------------------------------------------------------------------
# All differential expression genes between groups in each cell type
# deg list for later use

DEGs <- list()
for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(Obj_Integrate,
                         logfc.threshold = 0,
                         ident.1 = paste0(cellfordeg[i],"_MBLI"),
                         ident.2 = paste0(cellfordeg[i],"_MB93"),
                         test.use = "wilcox", assay = "RNA",
                         verbose = T)
  DEGs[[cellfordeg[i]]] = CELLDEG
}
names(DEGs) <- cellfordeg
saveRDS(DEGs, paste0(opt$opd, file="FindMarkers.wilcox_IB6.celltype_AllDEGs_LIvs93.rds"))

gg_title <- "IB6_FindMarkers.wilcox_celltype_AllDEG_GSEA"

# -------------------------------------------------------------------------------------------------------------------------------------------------------
#  载入需要的数据库
msigdbr_species()
mouse <- msigdbr(species = "Mus musculus")
# table(mouse$gs_subcat)
mouse_GO_bp = msigdbr(species = "Mus musculus",
                      category = "C5", #GO在C5
                      subcategory = "GO:BP") %>%
  dplyr::select(gs_name,gene_symbol)#这里可以选择gene symbol，也可以选择ID，根据自己数据需求来，主要为了方便

mouse_CP_REACTOME = msigdbr(species = "Mus musculus",
                            category = "C2", #REACTOME在C2
                            subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name,gene_symbol)

mouse_CP_KEGG = msigdbr(species = "Mus musculus",
                        category = "C2", #KEGG在C2
                        subcategory = "CP:KEGG") %>%
  dplyr::select(gs_name,gene_symbol)

# -------------------------------------------------------------------------------------------------------------------------------------------------------
#  run GSEA
GSEA_for_each_celltype <- function(DEGs_celltype, name_celltype, TERM2GENE = TERM2GENE, data_source){
  gsea_input <- DEGs_celltype$avg_log2FC
  names(gsea_input) <- rownames(DEGs_celltype)
  gsea_input <- sort(gsea_input, decreasing = T)
  str(gsea_input)
  gse <- GSEA(gsea_input,
              TERM2GENE = TERM2GENE,
              pvalueCutoff = 0.5,
              pAdjustMethod = "BH",
              eps = 0)
 # return(gse)
  return(gse@result)
}
# saveRDS(GSEA_for_each_celltype, paste0(opt$opd, gg_title, ".rds"))  # save gse

# return gse@result
DEG_GO_BP <- list()
for (i in 1:length(DEGs)) {
  GSEA <- GSEA_for_each_celltype(DEGs_celltype=DEGs[[i]], name_celltype=names(DEGs[i]), TERM2GENE=mouse_GO_bp, data_source="mouse_GO_BP")
  DEG_GO_BP[[names(DEGs[i])]] <- GSEA
}
saveRDS(DEG_GO_BP, paste0(opt$opd, gg_title, "_mouse_GO_BP", ".rds"))

#
DEG_CP_REACTOME <- list()
for (i in 1:length(DEGs)) {
  GSEA <- GSEA_for_each_celltype(DEGs_celltype=DEGs[[i]], name_celltype=names(DEGs[i]), TERM2GENE=mouse_CP_REACTOME, data_source="mouse_CP_REACTOME")
  DEG_CP_REACTOME[[names(DEGs[i])]] <- GSEA
}
saveRDS(DEG_CP_REACTOME, paste0(opt$opd, gg_title, "_mouse_CP_REACTOME", ".rds"))

#
DEG_CP_KEGG <- list()
for (i in 1:length(DEGs)) {
  GSEA <- GSEA_for_each_celltype(DEGs_celltype=DEGs[[i]], name_celltype=names(DEGs[i]), TERM2GENE=mouse_CP_KEGG, data_source="mouse_CP_KEGG")
  DEG_CP_KEGG[[names(DEGs[i])]] <- GSEA
}
saveRDS(DEG_CP_KEGG, paste0(opt$opd, gg_title, "_mouse_CP_KEGG", ".rds"))

dat <- DEG_GO_BP
# 去掉Description列 REACTOME_ 前缀, 全部小写
for (i in 1:length(dat)) {
  dat[[i]]$Description <- gsub('REACTOME_', '', dat[[i]]$Description)  
  dat[[i]]$Description <- tolower(dat[[i]]$Description) 
}

# --------------------------------------------------------------------------------------------------------------------------------------------
## visualization of enriched pathways, dotplot - geneRatio
# @result, S3: data.frame
# GSEA基因集前缀分为: REACTOME, KEGG, GOBP ...

GSEA_list <- REACTOME; region = "BS"; opt_dir = opt$opd; data.base = "REACTOME"

GSEA_list <- Filter(function(dat) nrow(dat) > 0, GSEA_list) # 去除空的数据框

GSEA_list <- lapply(GSEA_list, function(dat){
  dat$gene_count <- sapply(dat$core_enrichment, function(genes) {
    genes_with_spaces <- gsub("/", " ", genes)               # 用空格分割基因名
    gene_vector <- unlist(strsplit(genes_with_spaces, " "))  # 计算基因数量
    length(gene_vector)
  })
  
  dat %>% 
    mutate(Description = str_replace_all(Description, "_", " ")) %>%   # 去除 Description列 字符连接符
    mutate(geneRatio = gene_count/setSize) %>%
    filter(p.adjust < 0.05)
    
  })

# 上调
up_GSEA_list <- lapply(GSEA_list, function(dat){
  dat[dat$enrichmentScore > 0, ]
})
up_GSEA_list <- Filter(function(dat) nrow(dat) > 0, up_GSEA_list)   # 去除空的数据框

# 下调
down_GSEA_list <- lapply(GSEA_list, function(dat){
  dat[dat$enrichmentScore < 0, ]
})
down_GSEA_list <- Filter(function(dat) nrow(dat) > 0, down_GSEA_list)   # 去除空的数据框


# plot_results <- list()

GSEA_Dotplot <- function(dat, region, data.base, celltype, regulation) {
    
  dat$Description <- gsub(rm_prefix, "", dat$Description)
  dat$Description <- gsub("_", " ", dat$Description)  # 去除 Description列 字符连接符

  dat$Description <- tolower(dat$Description)
  
  dat$Description <- factor(dat$Description, 
                                   levels = dat$Description[order(dat$geneRatio)])

  gsea_plot <- ggplot(data = dat, 
                      aes(x = geneRatio, y = Description, color = p.adjust, size = gene_count)) + 
    geom_point(alpha=0.8) +
    guides(size = guide_legend(order = 2),
           colour = guide_legend(order = 1)) +  
    theme(axis.title.y.left = element_blank()) 
  
#  plot_results[[names(GSEA_list[i])]] <- gsea_plot
#  return(plot_results)
  ggsave(paste0(opt_dir, region, "_", regulation, "_", data.base, "_", celltype, ".png"), 
         width = 7, height = 1+0.2*nrow(dat))

}

dat <- up_GSEA_list; regulation <- "up"
dat <- down_GSEA_list; regulation <- "down"

# 数据框为空时，会提示错误
for (i in 1:length(dat)) {
  GSEA_plot <- GSEA_Dotplot(dat = dat[[i]], region, regulation, data.base, celltype = names(dat[i]))  
}

## ------------------------------
# GSEA result dot plot, up&down
for (i in 1:length(REACTOME)) {  
  # 对每个数据框添加 sign 列  
  REACTOME[[i]]$sign <- ifelse(REACTOME[[i]]$enrichmentScore > 0, "up-regulated", "down-regulated")  
}  

dat <- REACTOME[[1]]

dat$gene_count <- sapply(dat$core_enrichment, function(genes) {
  genes_with_spaces <- gsub("/", " ", genes)               # 用空格分割基因名
  gene_vector <- unlist(strsplit(genes_with_spaces, " "))  # 计算基因数量
  length(gene_vector)
})
dat$geneRatio <- dat$gene_count/dat$setSize

dat <- dat[dat$p.adjust < 0.05, ]
dat$Description <- gsub("_", " ", dat$Description)

library(ggplot2)
library(forcats)
## from Tommy's code
p <- ggplot(dat, aes(x = geneRatio, y = fct_reorder(Description, geneRatio))) + 
  geom_point(aes(size = geneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.10), low="red") +
  ylab(NULL) +
  ggtitle("Reactome pathway enrichment")

p + facet_grid(.~sign)
ggsave("FB_Neuron_REACTOME_dotplot.png",dpi = 300)

# --------------------------------------------------------------------------------------------------------------------------------------------
# Visualization genes of enriched pathway, Gene-Concept Network
# S4, DOSE::gseaResult
          
Visul_GSEA <- function(Database, region, path, data.base, rm_prefix, celltype){
  
  # differential expression genes between groups in each cell type, FindMarkers
  # deg list
  deg.ct <- deg[[i]] # 1
  deg.ct$symbol <- rownames(deg.ct)
  geneList <- deg.ct$avg_log2FC
  names(geneList) <- rownames(deg.ct)
  
  Database@result$Description <- gsub(rm_prefix, "", Database@result$Description)

  Database@result$Description <- gsub('_', ' ', Database@result$Description)   # 去除 Description列 字符连接符
  Database@result$Description <- tolower(Database@result$Description) 

  cnetplot(Database, color.params = list(foldChange = geneList),
                 circular = TRUE)
  ggsave(paste0(path, region, "_", data.base, "_", celltype, "_circular_cnetplot.pdf"))
  
  if(F){
  ## enrichment map
  
  GSEA_enri <- pairwise_termsim(Database)
  emapplot(GSEA_enri, layout="kk")
  ggsave(paste0(path, region, "_", data.base, "_", celltype, "_emapplott.pdf"),
         height = 10, width = 10, dpi = 300)
  
  # ridgeline plot for GSEA result
  p <- ridgeplot(Database)
  p$data %>% ggplot(aes_string(x = value, y = "category")) +
    
  ggsave(paste0(path, region, "_", data.base, "_", celltype, "_ridgeplot.pdf"),
         height = 10, width = 10, dpi = 300)
  }
}

for (i in 1:length(deg)) {              # 确保deg列表与基因集列表一致，细胞类型
  IB <- Visul_GSEA(Database = REACTOME[[i]], data.base = "REACTOME",   # 1, 2      # Database = KEGG[[i]]
                   region = "FB", path = opt$opd_FB, 
                   rm_prefix = paste0(data.base, "_"), celltype = names(REACTOME[i]))  # 3
}

