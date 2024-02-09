#!/usr/bin/Rscript

#load custom functions & packages
source("./customFunctions.R")
library(singleseqgset)
library(circlize)

#############################################################
#######   subset and recluster on tumor/fibroblast   ########
#############################################################

#Load in processed data and complete analysis on all cells
seu.obj <- readRDS(file = "../output/s3/naive6_QCfilter_2000Feats_res0.8_dims45_dist0.35_neigh40_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "majorID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "freqID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "id")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "tumorO")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident_2", metaAdd = "name")
outName <- "tumor_naive6"

sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "name", metaAdd = "colz")

#subset on tumor cells
seu.obj <- subset(seu.obj,
                  subset = 
                  majorID ==  "tumor" | majorID ==  "cyclingTumor")

#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj, outDir = "../output/s2/", subName = "tumor_QCfiltered_3000", preSub = T, nfeatures = 3000,
                      vars.to.regress = "percent.mt"
                       )

# seu.obj <- readRDS(file = "../output/s2/tumor_QCfiltered_2000_S2.rds")
clusTree(seu.obj = seu.obj, dout = "../output/clustree/", outName = "tumor_QCfiltered_3000", test_dims = "40", algorithm = 3, prefix = "integrated_snn_res.")

#complete dim reduction
seu.obj <- readRDS("../output/s2/tumor_QCfiltered_3000_S2.rds")
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = "tumor_QCfiltered_3000", final.dims = 40, final.res = 0.5, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.5, n.neighbors = 50, assay = "integrated", saveRDS = F,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

#check QC parameters
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("../output/", outName, "/", outName, "_QC_feats.png", sep = ""), width = 9, height = 3)


#remove low qualtiy and suspected doublet clusters
seu.obj.sub <- subset(seu.obj,invert = T,
                  subset = 
                  clusterID_sub ==  "10" | clusterID_sub ==  "4")

seu.obj.sub$clusterID_sub <- NULL

# seu.obj <- readRDS(file = "../output/s2/tumor_QCfiltered_2000_S2.rds")
clusTree(seu.obj = seu.obj.sub, dout = "../output/clustree/", outName = "tumor_QCfiltered_2_3000", test_dims = "40", algorithm = 3, prefix = "integrated_snn_res.")

#repeat dim reduction
seu.obj <- dataVisUMAP(seu.obj = seu.obj.sub, outDir = "../output/s3/", outName = "tumor_QCfiltered_2_3000", final.dims = 40, final.res = 0.5, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.5, n.neighbors = 50, assay = "integrated", saveRDS = F,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

#check QC parameters
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("../output/", outName, "/", outName, "_QC_feats.png", sep = ""), width = 9, height = 3)


#remvoe the last low QC cells
seu.obj.sub <- subset(seu.obj,invert = T,
                  subset = 
                  clusterID_sub ==  "10")

seu.obj.sub$clusterID_sub <- NULL

# seu.obj <- readRDS(file = "../output/s2/tumor_QCfiltered_2000_S2.rds")
clusTree(seu.obj = seu.obj.sub, dout = "../output/clustree/", outName = "tumor_QCfiltered_2_3000", test_dims = "40", algorithm = 3, prefix = "integrated_snn_res.")

#repeat dim reduction with only high qualtiy cells
seu.obj <- dataVisUMAP(seu.obj = seu.obj.sub, outDir = "../output/s3/", outName = "tumor_QCfiltered_2_3000", final.dims = 40, final.res = 0.5, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.5, n.neighbors = 50, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )


####################################################
#######   Begin tumor/fibroblast analysis   ########
####################################################

#load in data
seu.obj <- readRDS(file = "../output/s3/tumor_QCfiltered_2_3000_res0.5_dims40_dist0.5_neigh50_S3.rds")
sorted_labels <- paste(sort(as.integer(levels(seu.obj$clusterID_sub))))
seu.obj$clusterID_sub <- factor(seu.obj$clusterID_sub, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "orig.ident_2", metaAdd = "name")
sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "name", metaAdd = "colz")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "name", metaAdd = "subtype")
outName <- "tumor_naive6"

#rename cell types
Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "Osteoblast_1", "1" = "Osteoblast_2", 
                                   "2" = "Osteoblast_3", "3" = "Osteoblast_cycling_1",
                                   "4" = "Hypoxic_osteoblast","5" = "Osteoblast_cycling_2",
                                   "6" = "Fibroblast", "7" = "Osteoblast_cycling_3", 
                                   "8" = "Osteoblast_cycling_4", "9" = "IFN-osteoblast")
                       )
seu.obj$majorID_sub <- Idents(seu.obj)

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "Osteoblast_1 (c0)", "1" = "Osteoblast_2 (c1)", 
                                   "2" = "Osteoblast_3 (c2)", "3" = "Osteoblast_cycling_1 (c3)",
                                   "4" = "Hypoxic_osteoblast (c4)","5" = "Osteoblast_cycling_2 (c5)",
                                   "6" = "Fibroblast (c6)", "7" = "Osteoblast_cycling_3 (c7)", 
                                   "8" = "Osteoblast_cycling_4 (c8)", "9" = "IFN-osteoblast (c9)")
                       )
seu.obj$majorID_subWclus <- Idents(seu.obj)


### Fig extra: Create violin plots for cell ID
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_sub", numOfFeats = 24, outName = "tumor_QCfiltered_3000", outDir = "../output/viln/tumor/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T
                     )


### Fig extra: check QC parameters
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("../output/", outName, "/", outName, "_QC_feats.png", sep = ""), width = 9, height = 3)


#create df with cell type colors for pretty plotting
majorColors.df <- as.data.frame(levels(seu.obj$majorID_sub))
colnames(majorColors.df) <- "ClusterID"
majorColors.df$colz <- c("#6D0026","#EB2C31","#EB2C31","#FFAA93",
                         "#F17D00", "#F6B31E", "#DABC9C", 
                         "#F49900", "#FFD6C6", "#AC0535")


majorColors.df$title <- "Cell Type"
majorColors.df$labCol <- c("white","black","black","black",
                           "black","black","black",
                           "black","black","white")

majorColors.df$labz <- levels(seu.obj$clusterID_sub)

### Fig 2a/b: Create labels to be cropped onto UMAP
leg <- cusLeg(legend = majorColors.df, clusLabel = "labz",legLabel = "ClusterID", colorz = "colz",labCol = "labCol",colz = 1, rowz = NULL, groupLabel = "title", dotSize = 6, groupBy = "title",sortBy = "labz", compress_x = 1.25, topBuffer = 1.1, ymin = 0, compress_y = 4)

ggsave(paste("../output/", outName, "/", outName, "_leg_forUMAP_labels.png", sep = ""), width = 3, height = 7)


### Fig 2a: Create raw UMAP
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_sub",
              cols = majorColors.df$colz,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              shuffle = TRUE
) + NoLegend()
pi <- cusLabels(plot = pi, shape = 21, size = 10, alpha = 0.8, labCol = majorColors.df$labCol, textSize = 5, smallAxes = TRUE) 
ggsave(paste("../output/", outName, "/", outName, "_rawUMAP.png", sep = ""), width = 7, height = 7)


### Fig extra: Create raw UMAP with orig clusID
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID",
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              shuffle = TRUE
)
pi <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8) + NoLegend()
ggsave(paste("../output/", outName, "/", outName, "_rawUMAP_clusterID.png", sep = ""), width = 7, height = 7)


### Fig extra: Create UMAP by cell cycle score
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "Phase",
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              shuffle = TRUE
)
pi <- formatUMAP(plot = pi) + NoLegend()
ggsave(paste("../output/", outName, "/", outName, "_UMAP_phase.png", sep = ""), width = 7, height = 7)


### Fig extra: Create UMAP by sample
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "subtype",
              split.by = "orig.ident_2",
              cols = levels(seu.obj$colz),
              pt.size = 0.25,
              ncol = 3,
              label = F,
              label.box = F,
              shuffle = TRUE
)

pi <- formatUMAP(plot = pi) + NoLegend()
ggsave(paste("../output/", outName, "/", outName, "_UMAPbySample.png", sep = ""), width = 10.5, height = 7)


### Fig extra: Make stacked bar graph
p <- stackedBar(seu.obj = seu.obj, downSampleBy = "name", groupBy = "name", 
                clusters = "clusterID_sub") + scale_fill_manual(labels = levels(seu.obj$name), 
                                                                values = levels(seu.obj$colz)
                                                               ) + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1, byrow =T)
                                                                                                             )
ggsave(file = paste("../output/", outName, "/", outName, "_stackedBar.png", sep = ""), width = 8, height = 12)


### Fig supp 2a: Search for fibs -- gene list from Azimuth_Cell_Types_2021
modulez <- list("FIBROBLAST_SIG" = c("FBLN1", "ADH1B", "DCN", "LUM", "C1S", "C1R", "SLIT2", "C3", "CFD", "ABLIM1"))

seu.obj <- AddModuleScore(seu.obj, features = modulez, name = "_score")
names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)

features <- names(modulez)
ecScores <- majorDot(seu.obj = seu.obj, groupBy = "clusterID_sub",
                     features = features
                    ) + theme(axis.title = element_blank(),
                              legend.position = "right",
                              legend.direction = "vertical",
                             ) + guides(color = guide_colorbar(title = 'Scaled\nenrichment\nscore'))

ggsave(paste("../output/", outName, "/", outName, "_fib_ecScore.png", sep = ""), width = 3.25, height=7)


### Fig extra: plot CopyKAT output on the subset data
#load in CopyKat classifications
cellCounts <- read.csv("../output/copyKat/cnvStat.csv")
cellCounts$X <- NULL
cnts <- cellCounts$cnvStat
names(cnts) <- cellCounts$code

seu.obj <- AddMetaData(
      object = seu.obj,
      metadata = cnts,
      col.name = 'copyKatPred'
      )

#exculde low quality sample from visual
seu.obj.sub <- subset(seu.obj,invert = T,
                  subset = 
                  orig.ident ==  "run_count_tumor_no_tx_7")

#plot the data
pi <- DimPlot(seu.obj.sub, 
                  reduction = "umap", 
                  group.by = "copyKatPred",
                  pt.size = 0.25,
                  label = FALSE,
                  label.box = FALSE,
               shuffle = T
                 )
pi <- formatUMAP(pi) + theme(legend.position = c(0.01, 0.95)) #+ NoLegend()
ggsave(paste("../output/", outName, "/", outName, "_uMAP_by_ploidy.png", sep = ""),width = 7,height=7)


### Fig 2b: Create heatmap of defining feats
#load in defining features determined using FinDAllMarkers
all.markers <- read.csv("../output/viln/tumor/tumor_QCfiltered_3000_gene_list.csv")
key.genes <- all.markers[!grepl("^ENSCAFG", all.markers$gene),] 
key.genes.sortedByPval = key.genes[order(key.genes$p_val),]

#sort and remove duplicate feats
features <- key.genes.sortedByPval %>%  group_by(cluster) %>% do(head(., n=5))
features <- as.data.frame(features[!duplicated(features$gene),])

#load in metadata and count data for hierarchical clustering
seu.obj$type <- seu.obj$majorID_sub
metadata <- seu.obj@meta.data
expression <- as.data.frame(t(seu.obj@assays$RNA@data)) #use log noralized count
expression$anno_merge <- seu.obj@meta.data[rownames(expression),]$type
expression <- expression[,colnames(expression) %in% c(features$gene,"anno_merge") ] #try subsetting instead if converting to a function

#get average expression by cell group and move cell types to rownames
clusAvg_expression <- expression %>% group_by(anno_merge) %>% summarise(across(where(is.numeric), mean)) %>% as.data.frame()
rownames(clusAvg_expression) <- clusAvg_expression$anno_merge
clusAvg_expression$anno_merge <- NULL
                                     
#complete hc
M <- (1- cor(t(clusAvg_expression),method="pearson"))/2
hc <- hclust(as.dist(M),method="complete")
mat <- scale(as.matrix(clusAvg_expression))
mat <- mat[,features$gene]

#plot on heatmap
outfile <- paste("../output/", outName, "/", outName, "_heatMap.png", sep = "")
png(file = outfile, width=4000, height=3000, res=400)
par(mfcol=c(1,1))                    
ht <- Heatmap(mat,
        name = "Scaled expression     ",
        rect_gp = gpar(col = "white", lwd = 2),
        border_gp = gpar(col = "black"),
        show_column_dend = F,
        cluster_rows = F,
        cluster_columns = F,
        row_dend_reorder = F,
              show_row_names = F,
        column_names_rot = 45, 
        column_names_side = "top",
              heatmap_legend_param = list(legend_direction = "horizontal", title_position = "lefttop", just = c("right", "top")))

draw(ht, padding = unit(c(2, 7, 2, 7), "mm"), heatmap_legend_side = "top")
dev.off()


### Fig 2c: Create featplot of key feats
features <- c("ALPL","COL13A1","COL3A1","FBLN1","VEGFA","ACTA2" )

p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 6, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 16)
ggsave(paste("../output/", outName, "/", outName, "_key_feats.png", sep = ""), width = 15, height = 3)


### Fig supp: umap split by sample
seu.obj$orig.subtype <- paste0(seu.obj$name,"_",seu.obj$subtype)
### Create UMAP by sample
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "orig.subtype",
              split.by = "orig.subtype",
              pt.size = 0.25,
              ncol = 2,
              #label = T,
              #label.box = T,
              shuffle = TRUE
)

pi <- formatUMAP(plot = pi) + NoLegend() + theme(axis.title = element_blank())
ggsave(paste("../output/", outName, "/", outName, "_UMAPbySample.png", sep = ""), width = 7, height = 10.5)


### Fig 2d: Complete gsea of tumor cell pops using singleseqgset

#get gene sets from msigdb
can_gene_sets <- msigdbr(species = "dog", category = "H")
can_gene_lists <- unique(can_gene_sets$gs_name)
can.sets <- vector("list",length=length(can_gene_lists))
names(can.sets) <- can_gene_lists

for (i in names(can.sets)) {
    can.sets[[i]] <- pull(can_gene_sets[can_gene_sets$gs_name == i, "gene_symbol"])
}

#run singleseqgset
logfc.data <- logFC(cluster.ids = seu.obj@meta.data$clusterID_sub,
                    expr.mat = seu.obj@assays$RNA@data)

gse.res <- wmw_gsea(expr.mat = seu.obj@assays$RNA@data,
                    cluster.cells = logfc.data[[1]],
                    log.fc.cluster = logfc.data[[2]],
                    gene.sets=can.sets)

res.stats <- gse.res[["GSEA_statistics"]]
res.pvals <- gse.res[["GSEA_p_values"]]

res.pvals <- apply(res.pvals, 2, p.adjust, method="fdr") #Correct for multiple comparisons

mat <- scale(as.matrix(res.stats))
rownames(mat) <- gsub("HALLMARK_", "", rownames(mat))


colz <- majorColors.df$colz
names(colz) <- as.character(majorColors.df$labz)
ha = HeatmapAnnotation(
    Cluster = colnames(mat),
    border = TRUE,
    col = list(Cluster = colz),
    show_legend = c(Cluster = FALSE),
    annotation_name_side = "none"
)

#plot the results in heatmap
outfile <- paste("../output/", outName, "/", outName, "_heatMap_hallmarks.png", sep = "")
png(file = outfile, width=4500, height=6000, res=400)
par(mfcol=c(1,1))                    
ht <- Heatmap(mat,
        name = "HALLMARK pathway activity",
        rect_gp = gpar(col = "white", lwd = 2),
        border_gp = gpar(col = "black"),
        show_column_dend = F,
        cluster_rows = T,
        cluster_columns = T,
        row_dend_reorder = F,
        column_names_rot = 0, 
        column_names_side = "top",
              column_names_gp = grid::gpar(fontsize = 16),
              top_annotation = ha,
              col = colorRamp2(c(-2, 0, 2), rev(brewer.pal(n=3, name="RdYlBu"))),
              heatmap_legend_param = list(title_position = "leftcenter-rot", legend_height = unit(6, "cm"))
             )
draw(ht, padding = unit(c(10, 10, 2, 30), "mm"), heatmap_legend_side = "left")
dev.off()


### Fig 2e: Comare gene expression between fibs and osteoblasts
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "clusterID_sub", idents.1 = "6", idents.2 = c('0',"1","2"), bioRep = "name",padj_cutoff = 0.05, lfcCut = 0.58, 
                        minCells = 25, outDir = paste0("../output/", outName, "/"), title = "c6_vs_c012", idents.1_NAME = "FIBROBLAST", idents.2_NAME = "OSTEOBLAST", returnVolc = T, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F, setSeed = 24
                    )

p  <- prettyVolc(plot = p_volc[[1]], rightLab = "Up in c6", leftLab = "Up in c0,1,2") + labs(x = "log2(FC) Fibroblasts vs Osteoblasts") + NoLegend()
ggsave(paste("../output/", outName, "/", outName, "_c6vc012_volcPlot.png", sep = ""), width = 7, height = 7)


### Fig extra: gsea of the DGE results
p <- plotGSEA(pwdTOgeneList = "../output/tumor_naive6/FIBROBLAST_vs_OSTEOBLAST_all_genes.csv", category = "C5", subcategory = NULL)
ggsave(paste("../output/", outName, "/", outName, "_enriched_terms_c6.png", sep = ""), width = 9, height =7)


### Fig 2e: Comare gene expression between hypoxic and non-hypoxic osteoblasts
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "clusterID_sub", idents.1 = "4", idents.2 = c("0","1","2"), bioRep = "name",padj_cutoff = 0.05, lfcCut = 0.58, 
                        minCells = 25, outDir = paste0("../output/", outName, "/"), title = "c4_vs_c012", idents.1_NAME = "HYPOXIC_OSTEOBLAST", idents.2_NAME = "OSTEOBLAST", returnVolc = T, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F, setSeed = 24
                    )

p  <- prettyVolc(plot = p_volc[[1]], rightLab = "Up in c4", leftLab = "Up in c0,1,2") + labs(x = "log2(FC) Hypoxic osteoblasts vs Osteoblasts") + NoLegend()
ggsave(paste("../output/", outName, "/", outName, "_c4vc012_volcPlot.png", sep = ""), width = 7, height = 7)


### Fig extra: gsea of the DGE results
p <- plotGSEA(pwdTOgeneList = "../output/tumor_naive6/HYPOXIC_OSTEOBLAST_vs_OSTEOBLAST_all_genes.csv", category = "C5", subcategory = NULL)
ggsave(paste("../output/", outName, "/", outName, "_enriched_terms_c4.png", sep = ""), width = 9.5, height =7)


### Fig extra: Comare gene expression between c4 and c1
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "clusterID_sub", idents.1 = "4", idents.2 = "1", bioRep = "name",padj_cutoff = 0.05, lfcCut = 0.58, 
                        minCells = 25, outDir = paste0("../output/", outName, "/"), title = "c4_vs_c1", idents.1_NAME = "c4", idents.2_NAME = "c1", returnVolc = T, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F, setSeed = 24
                    )

p  <- prettyVolc(plot = p_volc[[1]], rightLab = "Up in c4", leftLab = "Up in c1") + labs(x = "log2(FC) c4 vs c1")
ggsave(paste("../output/", outName, "/", outName, "_c4vc1_volcPlot.png", sep = ""), width = 7, height = 7)


### Fig extra: Comare gene expression between c0 and c1
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "clusterID_sub", idents.1 = "0", idents.2 = "1", bioRep = "name",padj_cutoff = 0.05, lfcCut = 0.58, 
                        minCells = 25, outDir = paste0("../output/", outName, "/"), title = "c0_vs_c1", idents.1_NAME = "c0", idents.2_NAME = "c1", returnVolc = T, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F, setSeed = 24
                    )
p  <- prettyVolc(plot = p_volc[[1]], rightLab = "Up in c0", leftLab = "Up in c1") + labs(x = "log2(FC) c0 vs c1")
ggsave(paste("../output/", outName, "/", outName, "_c0vc1_volcPlot.png", sep = ""), width = 7, height = 7)


### Fig extra: Comare gene expression between c2 and c1
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "clusterID_sub", idents.1 = "2", idents.2 = "1", bioRep = "name",padj_cutoff = 0.05, lfcCut = 0.58, 
                        minCells = 25, outDir = paste0("../output/", outName, "/"), title = "c2_vs_c1", idents.1_NAME = "c2", idents.2_NAME = "c1", returnVolc = T, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F, setSeed = 24
                    )
p  <- prettyVolc(plot = p_volc[[1]], rightLab = "Up in c2", leftLab = "Up in c1") + labs(x = "log2(FC) c2 vs c1")
ggsave(paste("../output/", outName, "/", outName, "_c2vc1_volcPlot.png", sep = ""), width = 7, height = 7)


##################################################
#######   End tumor/fibroblast analysis   ########
##################################################
