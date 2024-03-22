#!/usr/bin/Rscript

#load custom functions & packages
source("./customFunctions.R")
library(singleseqgset)

#######################################################################
#######   subset and recluster on macropahge and osteoclasts   ########
#######################################################################

#load in all cell object and add metadata
seu.obj <- readRDS(file = "../output/s3/naive6_QCfilter_2000Feats_res0.8_dims45_dist0.35_neigh40_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "majorID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "freqID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "id")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "tumorO")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "orig.ident_2", metaAdd = "name")

sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "name", metaAdd = "colz")

#subset on macs as ID'd on the total myeloid subset
seu.obj <- subset(seu.obj,
                  subset = 
                  majorID ==  "mac" | majorID ==  "oc") 


seu.obj <- indReClus(seu.obj = seu.obj, outDir = "../output/s2/", subName = "macOC_QCfiltered_2500", preSub = T, nfeatures = 2500,
                      vars.to.regress = c("percent.mt")
                       )


seu.obj <- readRDS(file = "../output/s2/macOC_QCfiltered_2500_S2.rds")
clusTree(seu.obj = seu.obj, dout = "../output/clustree/", outName = "macOC_QCfiltered_2500", test_dims = c(45,40,35), algorithm = 3, prefix = "integrated_snn_res.")



seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = "macOC_QCfiltered_2500", final.dims = 40, final.res = 0.5, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.25, n.neighbors = 40, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA"))


seu.obj.sub <- subset(seu.obj, invert = T,
                  subset = 
                  clusterID_sub ==  "7") 


seu.obj <- dataVisUMAP(seu.obj = seu.obj.sub, outDir = "../output/s3/", outName = "macOC_QCfiltered_2500", final.dims = 40, final.res = 0.6, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.25, n.neighbors = 40, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA"))

#######################################
#######   Begin MAC analysis   ########
#######################################

### Resume the journey
seu.obj <- readRDS(file = "../output/s3/macOC_QCfiltered_2500_res0.6_dims40_dist0.25_neigh40_S3.rds")
sorted_labels <- paste(sort(as.integer(levels(seu.obj$clusterID_sub))))
seu.obj$clusterID_sub <- factor(seu.obj$clusterID_sub, levels = sorted_labels)
sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "name", metaAdd = "colz")
outName <- "mac_naive6"


### Create violin plots for cell ID
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID_sub", numOfFeats = 24, outName = "macOC_QCfiltered_2500", outDir = "../output/viln/myeloid/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T
                     )

#stash new IDs
Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "TAM_ACT", "1" = "TAM_INT", 
                                   "2" = "LA-TAM_SPP2_hi", "3" = "LA-TAM_C1QC_hi",
                                   "4" = "CD4-_TIM","5" = "Cycling_OC",
                                   "6" = "Mature_OC", "7" = "ANGIO_TAM", 
                                   "8" = "Cycling_OC", "9" = "CD320_OC", 
                                   "10" = "IFN-TAM", "11" = "CD4+_TIM")
                       )


seu.obj$majorID <- Idents(seu.obj)

seu.obj$majorID <- factor(seu.obj$majorID, levels = c("CD4-_TIM (c4)","CD4+_TIM (c11)","Angio-TAM (c2)","LA-TAM_TREM2hi (c0)","LA-TAM_C1QChi (c5)",
                                                 "IFN-TAM (c8)","Cycling-OC (c4)","Reactive-OC (c7)","Mature-OC (c6)","ECM-OC (c3)"))

#used later
ct.l3 <- c(ct.l3,seu.obj$majorID)


### Fig extra: QC check
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("../output/", outName, "/", outName, "_QC_feats.png", sep = ""), width = 9, height = 3)


#set colors
colz <- gg_color_hue(12)
colz[c(6,7,9,10)] <- c("grey60","grey70","grey80","grey75")
colz[c(1,2,3,4,5,8,11,12)] <- c("#8C94BF","#00366C","#023FA5","#9BBFDD","#5295D4","#3267AD","#0066A5","#D0E4FF")


### Fig 5a: Create raw UMAP
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_sub",
              pt.size = 0.25,
              label = TRUE,
              cols = colz,
              label.box = TRUE,
              shuffle = TRUE
)
pi <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.85, 
                labCol = c("black", "white","white","black","black","black","black","white","black","black","white","black")) + NoLegend() + theme(axis.title = element_blank(),
                                                                                   panel.border = element_blank())
ggsave(paste("../output/", outName, "/", outName, "_rawUMAP_clusterID_sub_mac.png", sep = ""), width = 7, height = 7)


### Fig extra: Create raw UMAP with original clusters
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID",
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              shuffle = TRUE
)
pi <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.85) + NoLegend() + theme(axis.title = element_blank(),
                                                                                   panel.border = element_blank())
ggsave(paste("../output/", outName, "/", outName, "_rawUMAP_clusterID_mac.png", sep = ""), width = 7, height = 7)


### Fig 5b-1: Plot key feats
features <- c("DLA-DRA", "DLA-DMA", "CD80",
              "SDC2","ALDOC","RBP4",
              "TREM2","APOC1","SPP2",
              "APOE","C1QC","IGF1")
              

p <- prettyFeats(seu.obj = seu.obj, nrow = 4, ncol =  3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 24, noLegend = T)
ggsave(paste("../output/", outName, "/", outName, "_key_feats_1.png", sep = ""), width = 6, height = 8, scale = 2)


### Fig 5b-2: Plot key feats
features <- c("SELL", "IL1B","VCAN",
              "VEGFC","PMEPA1","IL18BP",
              "RSAD2","ISG15","CD40",
              "CD4","IL18","LTF")
              

p <- prettyFeats(seu.obj = seu.obj, nrow = 4, ncol =  3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 24, noLegend = T)
ggsave(paste("../output/", outName, "/", outName, "_key_feats_2.png", sep = ""), width = 6, height = 8, scale = 2)


### Fig 5c feats that define on heatmap
#stash obj b/c need to subset on just the macs
seu.obj.backup <- seu.obj
Idents(seu.obj) <- "clusterID_sub"

#subset and re-run normalization
seu.obj <- subset(seu.obj,invert = T,
                 idents = c(5,6,8,9))
seu.obj <- NormalizeData(seu.obj)

#set feats
feats_forHeat <- c("CXCL10", "CD40", "CCL19","CD5L","DLA-DRA",
                   "GCAT", "IFIT", "IL7R","KCNH7","NEURL3",
                   "GPNMB", "PTGR1", "CD36","SCIN","SPP2","SPP1","TREM2",
                   "C1QC", "IGF1", "NIBAN3","SERPING1","CLDN1",
                   "SELL", "ITGB7", "PYGL","NOD2","OSM",
                   "VEGFA", "PRDX6", "TIMP1", "PTGES","VEGFC",
                   "RSAD2", "OAS3", "IFIT3","ISG15","TNFSF10",
                   "LTF", "PADI3", "PTGS2","CD4","THBS1"
                  )

#stash new metadata slot for ease
seu.obj$type <- seu.obj$clusterID_sub

#extract data
metadata <- seu.obj@meta.data
expression <- as.data.frame(t(seu.obj@assays$RNA@data)) #use log noralized count
expression$anno_merge <- seu.obj@meta.data[rownames(expression),]$type

expression <- expression[,colnames(expression) %in% c(feats_forHeat,"anno_merge") ] #try subsetting instead if converting to a function

#get average expression for each cluster
clusAvg_expression <- expression %>% group_by(anno_merge) %>% summarise(across(where(is.numeric), mean)) %>% as.data.frame()
rownames(clusAvg_expression) <- clusAvg_expression$anno_merge
clusAvg_expression$anno_merge <- NULL

#run hierarchical clustering
M <- (1- cor(t(clusAvg_expression),method="pearson"))/2
hc <- hclust(as.dist(M),method="complete")

#extract matrix
mat <- scale(as.matrix(clusAvg_expression))

#plot the data
outfile <- paste("../output/", outName, "/", outName, "_heatMap.png", sep = "")
png(file = outfile, width=6000, height=2000, res=400)
par(mfcol=c(1,1))                    
ht <- Heatmap(mat,
        name = "Scaled\nexpression         ",
        rect_gp = gpar(col = "white", lwd = 2),
        border_gp = gpar(col = "black"),
        show_column_dend = F,
        cluster_rows = hc,
        row_dend_reorder = F,
              column_names_gp = grid::gpar(fontsize = 12),
              row_names_gp = grid::gpar(fontsize = 16),
        column_names_rot = 45, column_names_side = "top",
                      heatmap_legend_param = list(legend_direction = "vertical", title_position = "topleft",  title_gp = gpar(fontsize = 16), 
                                              labels_gp = gpar(fontsize = 12), legend_width = unit(6, "cm"))

       )
draw(ht, padding = unit(c(2, 5, 2, 5), "mm"), heatmap_legend_side = "right")
dev.off()


### Fig 5d: Feats that define LA-TAM subsets
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "clusterID_sub", idents.1 = "3", idents.2 = "2", bioRep = "name",padj_cutoff = 0.05, lfcCut = 0.58, 
                        minCells = 25, outDir = paste0("../output/", outName, "/"), title = "c3_vs_c2", idents.1_NAME = "LA_TAM_C1QC_hi", idents.2_NAME = "LA_TAM_SPP1_hi", returnVolc = T, doLinDEG = F, paired = T, addLabs = NULL,lowFilter = T, dwnSam = F
                    )

p  <- prettyVolc(plot = p_volc[[1]], rightLab = "Up in c3", leftLab = "Up in c2") + labs(x = "log2(FC) c3 vs c2")

ggsave(paste("../output/", outName, "/", outName, "_c3vc2_volcPlot.png", sep = ""), width = 7, height = 7)


### Fig 5e & Supp: M1 vs M2 signatures
#set modules
modulez <- list("Pro-inflammatory" = c("AZIN1", "CD38","CD86","CXCL10","FPR2","GPR18","IL12B","IL18","IRF5","NFKBIZ","NOS2","PTGS2","TLR4","TNF"),
                "Anti-inflammatory" = c("ALOX15", "ARG1", "CHIL3", "CHIL4","EGR2", "IL10","IRF4","KLF4","MRC1","MYC","SOCS2","TGM2")
                )

seu.obj <- AddModuleScore(seu.obj,
                          features = modulez,
                         name = "_score")
names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)

features <- names(modulez)
ecScores <- majorDot(seu.obj = seu.obj, groupBy = "clusterID_sub",
                     features = rev(features)
                    ) + coord_flip() + theme(plot.margin = margin(7, 7, 7, 7, "pt"),
                                                                   axis.title = element_blank(),
                                                                   legend.position = "right",
                                                                   legend.direction = "vertical",
                                                                   axis.text.x = element_text(angle=0, hjust = 0.5)
                                            ) + scale_y_discrete(position = "right") + scale_colour_continuous(name="Enrichment score", type = "viridis") + NoLegend()

ggsave(paste("../output/", outName, "/", outName, "_dotPlot_ecScores.png", sep = ""), width = 6,height = 1)


### Fig extra: plot each indiviudal feature
modulez <- c(list("Enrichment score" = names(modulez)), modulez)

plots <- lapply(modulez, function(x){
    majorDot(seu.obj = seu.obj, groupBy = "clusterID_sub",
                     features = rev(unname(unlist(x)))
                    ) + coord_flip() + theme(axis.text.x = element_blank(),
                                                          axis.ticks.x = element_blank(),
                                                          legend.position = "right",
                                                                   legend.direction = "vertical",
                                                          axis.title = element_blank(),
                                                          plot.margin = margin(3, 0, 3, 0, "pt")
                                                         ) + scale_colour_viridis(option="magma", name='Average\nexpression', limits = c(-2.5,2.5))
    
})

    patch <- area()
    nrow <- length(modulez)
    ncol <- 1
    counter=0
    for (i in 1:nrow) {
        for (x in 1:ncol) {
                patch <- append(patch, area(t = i, l = x, b = i, r = x))
        }
    }


plots$`Enrichment score` <- plots$`Enrichment score` + theme(axis.text.x = element_text(angle=0, hjust = 0.5)
                            ) + scale_y_discrete(position = "right") + scale_colour_viridis() + guides(color = guide_colorbar(title = 'Module\nscore'), limits = c(-2.5,2.5))

p <- Reduce( `+`, plots ) +  plot_layout(guides = "collect", design = patch, 
                                                                             height = unname(unlist(lapply(modulez, length)))/sum(unname(unlist(lapply(modulez, length)))))

ggsave(paste("../output/", outName, "/", outName, "_dotPlot_ecScores_2.png", sep = ""), width = 6,height=7)


### Figure 5f: genes that define M1 vs M2
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "clusterID_sub", idents.1 = "11", idents.2 = "3", bioRep = "name",padj_cutoff = 0.05, lfcCut = 0.58, 
                        minCells = 5, outDir = paste0("../output/", outName, "/"), title = "c11_vs_c3", idents.1_NAME = "CD4pos_TIM", idents.2_NAME = "LA_TAM_C1QC_hi", returnVolc = T, doLinDEG = F, paired = T, addLabs = NULL,lowFilter = T, dwnSam = F
                    )

p  <- prettyVolc(plot = p_volc[[1]], rightLab = "Up in c11", leftLab = "Up in c3") + labs(x = "log2(FC) c11 vs c3")

ggsave(paste("../output/", outName, "/", outName, "_CD4pos_TIMvLA_TAM_C1QC_hi_volcPlot.png", sep = ""), width = 7, height = 7)


### Fig 6g: make heatmap M1 vs M2
geneLists <- read.csv("../output/mac-oc_naive6/CD4pos_TIM_vs_LA_TAM_C1QC_hi_all_genes.csv")

#extract top 20 feats for each direction of the conntrast
geneListUp <- geneLists %>% arrange(padj) %>% filter(log2FoldChange > 0) %>% .$gene
geneListDwn <- geneLists %>% arrange(padj) %>% filter(log2FoldChange < 0) %>% .$gene
feats_forHeat <- c(head(geneListUp,20), head(geneListDwn,20))

#transfer meta slot for ease
seu.obj$type <- seu.obj$clusterID_sub

#extract counts and meta
metadata <- seu.obj@meta.data
expression <- as.data.frame(t(seu.obj@assays$RNA@data)) #use log noralized count
expression$anno_merge <- seu.obj@meta.data[rownames(expression),]$type
expression <- expression[,colnames(expression) %in% c(feats_forHeat,"anno_merge") ] #try subsetting instead if converting to a function

#calc avg expressionfor each cluster
clusAvg_expression <- expression %>% group_by(anno_merge) %>% summarise(across(where(is.numeric), mean)) %>% as.data.frame()
rownames(clusAvg_expression) <- clusAvg_expression$anno_merge
clusAvg_expression$anno_merge <- NULL

#complete herichcial cluster and spin clades to get in visually pleasing presnetion
M <- (1- cor(t(clusAvg_expression),method="pearson"))/2
hc <- hclust(as.dist(M),method="complete")
hc <- dendextend::rotate(hc, c(8,5,6,2,1,7,3,4))
hc <- dendextend::rotate(hc, c("11","4","7","1","0","10","2","3"))

#get the data
mat <- scale(as.matrix(clusAvg_expression))
mat <- mat[,feats_forHeat]

#plot the data
outfile <- paste("../output/", outName, "/", outName, "_m1vm2_heatMap.png", sep = "")
png(file = outfile, width=4000, height=4000, res=400)
par(mfcol=c(1,1))                    
ht <- Heatmap(t(mat),
        name = "Scaled expression     ",
        rect_gp = gpar(col = "white", lwd = 2),
        border_gp = gpar(col = "black"),
        show_column_dend = T,
                      show_row_dend = F,

        cluster_rows = T,
        cluster_columns = hc,
        row_dend_reorder = F,
              column_names_gp = grid::gpar(fontsize =18),
              row_names_gp = grid::gpar(fontsize = 14),
        column_names_rot = 45, 
        column_names_side = "bottom",
              heatmap_legend_param = list(legend_direction = "horizontal", title_position = "topleft",title_gp = gpar(fontsize = 16), 
                                              labels_gp = gpar(fontsize = 12),legend_width = unit(6, "cm")))
draw(ht, padding = unit(c(2, 5, 2, 5), "mm"), heatmap_legend_side = "top")
dev.off()


### Fig 5h: make heatmap GSEA HALLMARK terms using singleseqgset
seu.obj$clusterID_sub <- droplevels(seu.obj$clusterID_sub)

can_gene_sets <- msigdbr(species = "dog", category = "H")

can_gene_lists <- unique(can_gene_sets$gs_name)
can.sets <- vector("list",length=length(can_gene_lists))
names(can.sets) <- can_gene_lists

for (i in names(can.sets)) {
    can.sets[[i]] <- pull(can_gene_sets[can_gene_sets$gs_name == i, "gene_symbol"])
}

logfc.data <- logFC(cluster.ids = seu.obj@meta.data$clusterID_sub,
                    expr.mat = seu.obj@assays$RNA@data)

gse.res <- wmw_gsea(expr.mat = seu.obj@assays$RNA@data,
                    cluster.cells = logfc.data[[1]],
                    log.fc.cluster = logfc.data[[2]],
                    gene.sets=can.sets)

res.stats <- gse.res[["GSEA_statistics"]]
res.pvals <- gse.res[["GSEA_p_values"]]

res.pvals <- apply(res.pvals, 2, p.adjust, method="fdr") #Correct for multiple comparisons
res.pvals <- res.pvals[rowMins(res.pvals) < 0.05,]


mat <- scale(as.matrix(res.stats))
rownames(mat) <- gsub("HALLMARK_", "", rownames(mat))

top10Terms <- do.call(c,lapply(1:ncol(mat),function(x){names(tail(sort(mat[,x]),10))}))

mat <- mat[rownames(mat) %in% top10Terms,]

library("circlize")
library("RColorBrewer")

#plot the results
outfile <- paste("../output/", outName, "/", outName, "_hallmarks_heatMap.png", sep = "")
png(file = outfile, width=6000, height=6000, res=400)
par(mfcol=c(1,1))                    
ht <- Heatmap(mat,
        name = "HALLMARK pathway activity",
        rect_gp = gpar(col = "white", lwd = 2),
        border_gp = gpar(col = "black"),
        show_column_dend = F,
        cluster_rows = T,
        row_dend_reorder = F,
              column_names_gp = grid::gpar(fontsize = 16),
              row_names_gp = grid::gpar(fontsize = 12),
           #  width = ncol(mat)*unit(5, "mm"),
            #  height = nrow(mat)*unit(5, "mm"),
        column_names_rot = 90, 
              column_names_side = "top",
              col = colorRamp2(c(-2, 0, 2), rev(brewer.pal(n=3, name="RdYlBu"))),
                      heatmap_legend_param = list(legend_direction = "horizontal", title_position = "topleft",  title_gp = gpar(fontsize = 16), 
                                              labels_gp = gpar(fontsize = 8), legend_width = unit(6, "cm"))

       )
draw(ht, padding = unit(c(0, 0, 0, 250), "mm"), heatmap_legend_side = "top")

dev.off()


##########################################
#######   END mac cell analysis   ########
##########################################

###########################################
#######   BEGIN oc cell analysis   ########
###########################################

### OC analysis
### Resume the journey
seu.obj <- readRDS(file = "../output/s3/macOC_QCfiltered_2500_res0.6_dims40_dist0.25_neigh40_S3.rds")
sorted_labels <- paste(sort(as.integer(levels(seu.obj$clusterID_sub))))
seu.obj$clusterID_sub <- factor(seu.obj$clusterID_sub, levels = sorted_labels)
sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "name", metaAdd = "colz")
outName <- "oc_naive6"

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "TAM_ACT", "1" = "TAM_INT", 
                                   "2" = "LA-TAM_SPP2_hi", "3" = "LA-TAM_C1QC_hi",
                                   "4" = "CD4-_TIM","5" = "Cycling_OC_1",
                                   "6" = "Mature_OC", "7" = "ANGIO_TAM", 
                                   "8" = "Cycling_OC_2", "9" = "CD320_OC", 
                                   "10" = "IFN-TAM", "11" = "CD4+_TIM")
                       )

seu.obj$majorID <- Idents(seu.obj)
ct.l3 <- c(ct.l3,seu.obj$majorID)

colz <- gg_color_hue(12)
colz[c(1,2,3,4,5,8,11,12)] <- c("grey60","grey75","grey85","grey90","grey70","grey80","grey88","grey79")
colz[c(6,7,9,10)] <- c("#645A9F","#AE8FD1","#EDCEF7","#CFADE5")

### Fig 6a: Create raw UMAP
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_sub",
              pt.size = 0.25,
              label = TRUE,
              cols = colz,
              label.box = TRUE,
              shuffle = TRUE
)
pi <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.85, labCol = c("black","black","black","black","black","white","black","black","black","black","black","black")) + NoLegend() + theme(axis.title = element_blank(),
                                                                                   panel.border = element_blank())
ggsave(paste("../output/", outName, "/", outName, "_rawUMAP_clusterID_sub_oc.png", sep = ""), width = 7, height = 7)


### Fig 6b: Plot key feats
features <- c("STMN1","ATP6V0D2", "TNIP3", 
              "TOP2A","CTSK", "CD320",
              "MYBL2","HYAL1","SRM")
p <- prettyFeats(seu.obj = seu.obj, nrow = 3, ncol =  3, features = features, 
                 color = "black", order = F, pt.size = 0.01, title.size = 20, titles = features, noLegend = T)
ggsave(paste("../output/", outName, "/", outName, "_key_feats.png", sep = ""), width = 6, height = 6, scale = 2)



### Fig 6c: complete heirchical clustering
colz <- gg_color_hue(12)
colz[c(1,2,3,4,5,8,11,12)] <- c("#8C94BF","#00366C","#023FA5","#9BBFDD","#5295D4","#3267AD","#0066A5","#D0E4FF")
colz[c(6,7,9,10)] <- c("#645A9F","#AE8FD1","#EDCEF7","#CFADE5")

seu.obj$type <- seu.obj$majorID

metadata <- seu.obj@meta.data
expression <- as.data.frame(t(seu.obj@assays$RNA@data)) #use log noralized count
expression$anno_merge <- seu.obj@meta.data[rownames(expression),]$type

clusAvg_expression <- expression %>% group_by(anno_merge) %>% summarise(across(where(is.numeric), mean)) %>% as.data.frame()
rownames(clusAvg_expression) <- clusAvg_expression$anno_merge
clusAvg_expression$anno_merge <- NULL                                                                       

M <- (1- cor(t(clusAvg_expression),method="pearson"))/2
hc <- hclust(as.dist(M),method="complete")

ggtree(as.phylo(hc)) + geom_tiplab(offset = 0.006) + xlim(NA,.10) + geom_tippoint(shape = 21,size = 8,alpha = 1, colour="black", fill = colz) + geom_tiplab(aes(label = as.character(seq(1:23)-1)), color = c("black","white","white","black","black","white","black","black","black","black","white","black"), hjust = 0.5)

ggsave(paste("../output/", outName, "/", outName, "_oc_hc.png", sep = ""), width = 4, height = 7)


### Fig 6d: c6 vs mac
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "clusterID_sub", idents.1 = "6", idents.2 = c("0","1","2","3"), bioRep = "name",
                      padj_cutoff = 0.05, lfcCut = 0.58, minCells = 25, outDir = paste0("../output/", outName, "/"), 
                      title = "oc_c6_vs_mac", idents.1_NAME = "MATURE_OC", idents.2_NAME = "MAC_c0-c3", 
                      returnVolc = T, doLinDEG = F, paired = T, addLabs = c("CTSK","ACP5"),lowFilter = T, dwnSam = F
                     )

p  <- prettyVolc(plot = p_volc[[1]], rightLab = "Up in c6", leftLab = "Up in mac") + labs(x = "log2(FC) oc_c3 vs mac") + NoLegend()

ggsave(paste("../output/", outName, "/", outName, "_c6vcmac_volcPlot.png", sep = ""), width = 7, height = 7)


### Fig supp: c9 vs c6
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "clusterID_sub", idents.1 = "9", idents.2 = "6", bioRep = "name",
                      padj_cutoff = 0.05, lfcCut = 0.58, minCells = 5, outDir = paste0("../output/", outName, "/"), 
                      title = "oc_c9_vs_6", idents.1_NAME = "CD320_OC", idents.2_NAME = "MATURE_OC", 
                      returnVolc = T, doLinDEG = F, paired = T, addLabs = "",lowFilter = T, dwnSam = F
                     )

p  <- prettyVolc(plot = p_volc[[1]], rightLab = "Up in c9", leftLab = "Up in c6") + labs(x = "log2(FC) oc_c9 vs c6") + NoLegend()

ggsave(paste("../output/", outName, "/", outName, "_c9vc6_volcPlot.png", sep = ""), width = 7, height = 7)

### Fig supp: c9 vs mac
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "clusterID_sub", idents.1 = "9", idents.2 = c("0","1","2","3"), bioRep = "name",
                      padj_cutoff = 0.05, lfcCut = 0.58, minCells = 5, outDir = paste0("../output/", outName, "/"), 
                      title = "oc_c9_vs_mac", idents.1_NAME = "CD320_OC", idents.2_NAME = "MAC_c0-c3", 
                      returnVolc = T, doLinDEG = F, paired = T, addLabs = "",lowFilter = T, dwnSam = F
                     )

p  <- prettyVolc(plot = p_volc[[1]], rightLab = "Up in c9", leftLab = "Up in mac") + labs(x = "log2(FC) oc_c9 vs mac") + NoLegend()

ggsave(paste("../output/", outName, "/", outName, "_c9vmac_volcPlot.png", sep = ""), width = 7, height = 7)


p <- plotGSEA(pwdTOgeneList = "../output/oc_naive6/MATURE_OC_vs_MAC_c0-c3_all_genes.csv", termsTOplot = 15, upOnly = T, category = "C5", subcategory = "GO:BP", dwnCol = "red")
ggsave(paste("../output/", outName, "/", outName, "_enriched_terms_c6vMac.png", sep = ""), width = 9, height = 7)


### Fig 6e: intersection of CD320+ OC DEGs
# NOTE: this code was abopted from the human - canine comparision code.
# Althought the variables and values indicate two species are being used,
# that is not the case. All analysis in this case is canine and the "human"
# associated varaibles are comparing CD320_OC to macropahge, while the
# "canine" associated varaibles are comparing CD320_OC to mature osteoclasts

geneLists.hu <- read.csv("../output/oc_naive6/CD320_OC_vs_MAC_c0-c3_all_genes.csv")
geneLists.can <- read.csv("../output/oc_naive6/CD320_OC_vs_MATURE_OC_all_genes.csv")
cONtrast <- c("c9vMAC", "c9vc6")

geneLists.hu <- geneLists.hu[geneLists.hu$log2FoldChange > 0,]
geneLists.can <- geneLists.can[geneLists.can$log2FoldChange > 0,]

nlab <- 20
axisCnt <- 10
overlapGenes <- rownames(seu.obj)

geneLists.hu <- geneLists.hu %>% mutate(species = "human",
                                        signed_pVal_hu = ifelse(log2FoldChange > 0, -log10(padj), log10(padj))
                                       ) %>% filter(gene %in% overlapGenes)# %>% top_n(n = ds, wt = signed_pVal_hu)

geneLists.can <- geneLists.can %>% mutate(species = "canine",
                                          signed_pVal_can = ifelse(log2FoldChange > 0, -log10(padj), log10(padj))
                                         ) %>% filter(gene %in% overlapGenes)# %>% top_n(n = ds, wt = signed_pVal_can)

geneList <- unique(c(geneLists.hu$gene, geneLists.can$gene)) %>% as.data.frame()
colnames(geneList) <- "gene"

colUp <- "red"
colDwn <- "blue"
geneList <- geneList %>% left_join(geneLists.hu[,c("gene","species","signed_pVal_hu")], 
                                   by = "gene") %>% left_join(geneLists.can[,c("gene","signed_pVal_can")], 
                                                              by = "gene") %>% replace(is.na(.), 0) %>% 
mutate(sig = ifelse(signed_pVal_can*signed_pVal_hu == 0,signed_pVal_can+signed_pVal_hu,signed_pVal_can*signed_pVal_hu) ,
       direction = case_when(signed_pVal_can*signed_pVal_hu != 0 & signed_pVal_hu > 0 & signed_pVal_can > 0 ~ "up",
                             signed_pVal_can*signed_pVal_hu != 0 & signed_pVal_hu < 0 &  signed_pVal_can < 0 ~ "dwn",
                             signed_pVal_can*signed_pVal_hu == 0 & signed_pVal_hu > 0 ~ "axis1",
                             signed_pVal_can*signed_pVal_hu == 0 & signed_pVal_can > 0 ~ "axis2",
                             signed_pVal_can*signed_pVal_hu == 0 & signed_pVal_hu < 0 ~ "axis3",
                             signed_pVal_can*signed_pVal_hu == 0 & signed_pVal_can < 0 ~ "axis4",                             
                             signed_pVal_can*signed_pVal_hu != 0 & signed_pVal_hu > 0 & signed_pVal_can < 0 ~ "conflict1",
                             signed_pVal_can*signed_pVal_hu != 0 & signed_pVal_hu < 0 & signed_pVal_can > 0 ~ "conflict2",
                            ),
       
       label = ifelse(sig != 0, gene, NA)
      )  %>% arrange(-abs(sig)) %>% group_by(direction) %>% 
            mutate(lab=ifelse(row_number() <= nlab
                              & signed_pVal_can*signed_pVal_hu != 0 , gene, ifelse(signed_pVal_can*signed_pVal_hu == 0 & row_number() <= axisCnt, gene, NA)), 
                   lab_col=case_when(direction == "up" ~ colUp,
                                     direction == "dwn" ~ colDwn,
                                     direction == "axis1" ~ "black",
                                     direction == "axis2" ~ "black",
                                     direction == "axis3" ~ "black",
                                     direction == "axis4" ~ "black",
                                     direction == "conflict1" ~"hotpink",
                                     direction == "conflict2" ~"hotpink")
                  ) 

annotationz <- geneList %>% group_by(direction) %>% summarize(cntz = n()) %>% mutate(xpos = c(Inf,0,Inf),
                                                                                     ypos =  c(0,Inf,Inf),
                                                                                     hjustvar = c(1,0.5,1),
                                                                                     vjustvar = c(1,1,1),
                                                                                     colz = c("black","black",colUp)
                                                                                    )



ggplot(geneList, aes(x=signed_pVal_hu,y=signed_pVal_can)) + geom_hline(yintercept=0,linetype=2) +
geom_vline(xintercept=0,linetype=2) + 
geom_point() + 
geom_label_repel(max.overlaps = Inf, size=1.5, label = geneList$lab, color = geneList$lab_col, show.legend = F,seed = 777) + 
theme_classic() + 
labs(x = "c9 vs mac DE signed log10(p.adj)",
     y = "c9 vs c6 DE signed log10(p.adj)",
     title = "Intersection of DEGs c9 vs mac/c6"
    ) + coord_cartesian(clip = 'off') + theme(panel.border = element_rect(color = "black",fill = NA,size = 1),
                                              plot.title = element_text(size = 12),
                                             axis.line = element_blank())


ggsave(paste("../output/", outName, "/", outName, "_pSigned_c9.png", sep = ""), width = 4,height=4)

#########################################
#######   END oc cell analysis   ########
#########################################
