#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")


#################################################################
#######   Integrate samples using SCTransform workflow   ########
#################################################################

#filter, integrate, and visulaize naive tumor data
load10x(din = "./input/", dout = "./output/", outName = "final_notx6_12_5mt", testQC = F,
       nFeature_RNA_high = 5500, nFeature_RNA_low = 200, percent.mt_high = 12.5, nCount_RNA_high = 75000, nCount_RNA_low = 100)

sctIntegrate(din = "./output/s1_naive/", outName = "naive_n6", vars.to.regress = c("percent.mt"), nfeatures = 2000)

seu.obj <- readRDS(file = "/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/naive_n6_seu.integrated.obj_S2.rds")
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "naive_n6_12-5", test_dims = c(50,45,40,35,30,25), algorithm = 3, prefix = "integrated_snn_res.")

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/", outName = "naive_n6", final.dims = 45, final.res = 0.8, stashID = "clusterID", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = F,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

#rename idents to get technical reps into one sample name
Idents(seu.obj) <- "orig.ident"
seu.obj <- RenameIdents(seu.obj, 
                        'run_count_tumor_no_tx_1_1' = 'Naive_1', 'run_count_tumor_no_tx_1_2' = 'Naive_1', 
                        'run_count_tumor_no_tx_2_1' = 'Naive_2', 'run_count_tumor_no_tx_2_2' = 'Naive_2', 
                        'run_count_tumor_no_tx_4' = 'Naive_3', 'run_count_tumor_no_tx_4' = 'Naive_3',
                        'run_count_tumor_no_tx_5' = 'Naive_4','run_count_no_tx_tumor_6' = 'Naive_5', 
                        'run_count_tumor_no_tx_7' = 'Naive_6')

seu.obj$orig.ident_2 <- seu.obj@active.ident

#save the processed file
saveRDS(seu.obj, file = "./output/s3/naive_n6_res0.8_dims45_dist0.3_neigh30_S3.rds")


### NOTE: during deeper analysis susspected low quality clusters were ID'd, so the following steps go back and remove them
seu.obj <- readRDS(file = "./output/s3/naive_n6_res0.8_dims45_dist0.3_neigh30_S3.rds")


### Fig extra: check QC parameters -- clus 1, 11, and 28 look suss
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("./output/", outName, "/", outName, "_QC_feats.png", sep = ""), width = 9, height = 3)


### Fig extra: see what happens when modifying filtering parameters
#set mod parameters
percent.mt_high = 10
nCount_RNA_high = 75000
nFeature_RNA_high = 5500 
nFeature_RNA_low = 500
nCount_RNA_low = 500

#subset based on parameters
seu.obj.sub <- subset(seu.obj,
                          subset = nFeature_RNA < nFeature_RNA_high & nFeature_RNA > nFeature_RNA_low &
                          percent.mt < percent.mt_high & nCount_RNA < nCount_RNA_high & nCount_RNA > nCount_RNA_low
                         )


#create plot
pi <- DimPlot(seu.obj.sub, 
              reduction = "umap", 
              group.by = "freqID",
              pt.size = 0.25,
              #cols = majorColors.df$colz,
              label = F,
              label.box = F,
              shuffle = TRUE
)
pi <- formatUMAP(plot = pi) + NoLegend() + theme(axis.title = element_blank(),
                                                 panel.border = element_blank())
ggsave(paste("./output/", outName, "/", outName, "_UMAPbyMajorID_loQual.png", sep = ""), width = 7, height = 7)


#CONCLUSION: clus 1, 11, and 28 are reduced with increase in nFeature_RNA_low filter, but so do neuts etc, so just remove the sus clusters

#remove the clus and repeat integration & vis
seu.obj.sub <- subset(seu.obj,invert = T,
                  subset = 
                  clusterID ==  "11" | clusterID ==  "1" | clusterID ==  "28")

#re-run integration
seu.obj <- indReClus(seu.obj = seu.obj.sub, outDir = "./output/s2/", subName = "naive6_QCfilter_2000Feats", preSub = T, nfeatures = 2000,
                      vars.to.regress = "percent.mt"
                       )

#check ideal clustering parameters
seu.obj <- readRDS(file = "/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/s2/naive6_QCfilter_2000Feats_S2.rds")
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "naive_n6_12-5", test_dims = c(50,45,40,35,30,25), algorithm = 3, prefix = "integrated_snn_res.")

#visualize the data
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "naive6_QCfilter_2000Feats", final.dims = 45, final.res = 0.8, stashID = "clusterID", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.35, n.neighbors = 40, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )


###########################################################
#######   Begin high-level analysis of all cells   ########
###########################################################

#load in processed data and metadata
seu.obj <- readRDS(file = "./output/s3/naive6_QCfilter_2000Feats_res0.8_dims45_dist0.35_neigh40_S3.rds")
outName <- "allCells"
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "majorID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "freqID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "id")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "tumorO")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident_2", metaAdd = "name")

sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "name", metaAdd = "colz")


### Fig extra: check QC parameters
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("./output/", outName, "/", outName, "_QC_feats.png", sep = ""), width = 9, height = 3)


### Fig extra: generate violin plots showing expression of defining markers
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID", numOfFeats = 24, outName = "naive_6", outDir = "./output/viln/allCells/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T
                     )


### Fig extra: run singleR classifier
singleR(seu.obj = seu.obj, outName = "naive6_QCfilter_200", clusters = "clusterID", outDir = "./output/singleR/")


### Fig extra: create raw UMAP colorized by cluster
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID",
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = TRUE
)
pi <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8)
ggsave(paste("./output/", outName, "/", outName, "_raw_UMAP.png", sep = ""), width = 10, height = 7)


#load in color data for major cell types
majorColors.df <- as.data.frame(levels(seu.obj$freqID))
colnames(majorColors.df) <- "ClusterID"
majorColors.df$colz <- c("#DD3900","#C89074","#99BFEF","#9CD6BA",
                         "#F17D00", "#3267AD", "#EEA4BE", "#023FA5",
                         "#9A8FC4", "#DABC9C")


### Fig 1a: Create UMAP by majorID
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "freqID",
              pt.size = 0.5,
              cols = majorColors.df$colz,
              label = F,
              label.box = F,
              shuffle = TRUE
)
pi <- formatUMAP(plot = pi) + NoLegend() + theme(axis.title = element_blank(),
                                                 panel.border = element_blank())
ggsave(paste("./output/", outName, "/", outName, "_UMAPbyMajorID.png", sep = ""), width = 7, height = 7)


### Fig 1a: create labels to be cropped onto UMAP
majorColors.df$title <- "All cells"
majorColors.df$labCol <- "black"
majorColors.df$labz <- ""
leg <- cusLeg(legend = majorColors.df, clusLabel = "labz",legLabel = "ClusterID", colorz = "colz",labCol = "labCol",colz = 1, rowz = NULL, groupLabel = "title", dotSize = 8, groupBy = "title",sortBy = "labz", compress_x = 0.9)

ggsave(paste("./output/", outName, "/", outName, "_leg_forUMAP_labels.png", sep = ""), width = 2.75, height = 7)


### Fig 1b: create faceted UMAP by sample
#downsample to get equal cell numbers for each biological replicate
set.seed(12)
Idents(seu.obj) <- "orig.ident_2"
seu.obj.sub <- subset(x = seu.obj, downsample = min(table(seu.obj$orig.ident_2)))

#create the plot
pi <- DimPlot(seu.obj.sub, 
              reduction = "umap", 
              group.by = "orig.ident_2",
              split.by = "orig.ident_2",
              pt.size = 0.25,
              cols = levels(seu.obj$colz),
              ncol = 3,
              label = F,
              label.box = F,
              shuffle = TRUE
)

pi <- formatUMAP(plot = pi) + NoLegend() + theme(axis.title = element_blank())
ggsave(paste("./output/", outName, "/", outName, "_UMAPbySample.png", sep = ""), width = 10.5, height = 7)

seu.obj$orig.ident_2 <- as.factor(seu.obj$orig.ident_2)
seu.obj$freqID <- factor(seu.obj$freqID, levels = rev(levels(cluster_freq.table$ClusterID)))

stackEM <- stackedBar(seu.obj = seu.obj, downSampleBy = "orig.ident_2", groupBy = "orig.ident_2", clusters = "freqID")
stackEM <- stackEM + scale_fill_manual(name="Cell source",
                                       labels = levels(seu.obj$orig.ident_2),
                                       values = levels(seu.obj$colz)
                                      ) + theme(axis.title.y = element_blank(),
                                                legend.text=element_text(size=16),
                                                legend.spacing = unit(1.0, 'cm'),
                                                legend.key.size = unit(0.8, 'cm')
                                               )
ggsave(paste("./output/", outName, "/", outName, "_stacked.png", sep = ""), width = 10.5, height = 7)


### Fig 1c: make pie chart
#get cell names from seu metadata
clusterList <- seu.obj$freqID

#convert into df with # of cells ("Count") by cell type ("ClusterID")
cluster_freq.table <- as.data.frame(table(clusterList)) %>% melt()
cluster_freq.table <- cluster_freq.table[,-2]
colnames(cluster_freq.table) <- c("ClusterID", "Count")

#convert raw number to pct
cluster_freq.table <- cluster_freq.table %>% mutate(pct = round(Count/sum(Count),3),
                                                    pos = 1-(cumsum(pct) - (0.5 * pct))) %>% arrange(desc(pct))

cluster_freq.table$pos <- rev(cluster_freq.table$pos)

#get color data incoperated
cluster_freq.table <- cluster_freq.table %>% left_join(majorColors.df, by = "ClusterID")
cluster_freq.table$ClusterID <- as.factor(cluster_freq.table$ClusterID )

# #reorder
cluster_freq.table <- cluster_freq.table[c(1,4,2,6,7,5,10,3,9,8),]
cluster_freq.table$ClusterID <- factor(cluster_freq.table$ClusterID, levels = cluster_freq.table$ClusterID)

#do calc to reposition the labels
df2 <- cluster_freq.table %>% 
  mutate(csum = rev(cumsum(rev(pct))), 
         pos = pct/2 + lead(csum, 1),
         pos = if_else(is.na(pos), pct/2, pos))

#create the plot
p <- ggplot(cluster_freq.table, aes(x = "" , y = pct, fill = ClusterID)) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
scale_fill_manual(labels = cluster_freq.table$ClusterID, 
               values = cluster_freq.table$colz) +
  geom_text_repel(data = df2,
                   aes(y = pos, label = paste0(pct*100, "%")),
                   size = 7, nudge_x = c(1,1,1,1,1,0.3,1.5,1,1,1),show.legend = FALSE
                   #colour = c("black","black","black","white","white", "white","black","black","black","black")
                  ) +
  guides(fill = guide_legend(title = "Cell type", nrow = 10)) +
  theme_void() + theme(legend.text=element_text(size=20),
                       legend.title=element_text(size=24),
                       legend.position="left",
                       legend.spacing = unit(1.0, 'cm'),
                       legend.key.size = unit(1.0, 'cm'),
                       legend.box.margin=margin(10,10,10,14))

ggsave(paste("./output/", outName, "/", outName, "_pieChart.png", sep = ""), width = 8, height = 6)


### Fig supp 1a: Module scoring to support classification - Liu 2021 (https://doi.org/10.3389/fonc.2021.709210)
#stash modules
modulez <- list("OS" = c("ALPL", "RUNX2", "IBSP"),
                "MYELOID" = c("LYZ", "CD68"),
                "OC" = c("ACP5", "CTSK"),
                "CAF" = c("COL1A1", "FAP", "VIM"),
                "NK-T" = c("CD2", "CD3D", "CD3E", "CD3G", "GNLY", "NKG7", "KLRD1", "KLRB1"),
                "BC"= c("MS4A1", "PAX5"),
                "PLASMA"= c("JCHAIN", "MZB1"),
                "EC" = c("EGFL7", "PLVAP", "ESAM")
                )

names(modulez) <- paste0(names(modulez),"_SIG")

#run Seurat's module score function
seu.obj <- AddModuleScore(seu.obj,
                          features = modulez,
                          name = "_score")
names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)

features <- names(modulez)

#create the plot
ecScores <- majorDot(seu.obj = seu.obj, groupBy = "freqID",
                     features = features, yAxis = c("Tumor/Fibroblast","Cycling Tumor","Neutrophil","TIM/TAM","Dendritic cell","Osteoclast","T cell","B cell","Mast cell" ,"Endothelial")
                    ) + theme(plot.margin = margin(3, 0, 3, 0, "pt"),
                                             axis.text.y=element_text(size=10),
                              axis.title = element_blank(),
                                                                   #axis.title = element_blank(),
                                                                   legend.position = "right",
                                                                   legend.direction = "vertical"#,
                                                                   #axis.text.x = element_text(angle=0, hjust = 0.5)
                                            )  + scale_colour_continuous(type = "viridis") + labs(x = "Human enrichment term", y = "Cell classification") + guides(color = guide_colorbar(title = 'Module\nscore'),size=guide_legend(nrow=3,byrow=TRUE))
ggsave(paste("./output/", outName, "/", outName, "_dotPlot_ecScores.png", sep = ""), width = 6,height=6)


### Fig 1d: Plot key feats
features <- c("PTPRC", "CD3E", "CD4",
              "GZMB", "DLA-DRA", "FLT3", 
              "CD68", "CTSK", "S100A12",
              "ALPL", "COL1A1", "FBLN1", 
              "TOP2A","MS4A1", "ESAM")
titles <- c("PTPRC (CD45)", "CD3E", "CD4",
            "GZMB", "DLA-DRA (MHCII)", "FLT3", 
            "CD68", "CTSK", "S100A12", 
            "ALPL","COL1A1", "FBLN1", 
            "TOP2A", "MS4A1 (CD20)", "ESAM")
colz <- c("black", "#9CD6BA", "#9CD6BA",
          "#9CD6BA", "black", "#023FA5",
          "#99BFEF", "#3267AD", "#C89074",
          "#DD3900", "#DD3900", "#DD3900",
          "black", "#9A8FC4", "#EEA4BE")
p <- prettyFeats(seu.obj = seu.obj, nrow = 5, ncol = 3, features = features, 
                 color = colz, order = F, pt.size = 0.0000001, title.size = 14, titles = titles, noLegend = T)
ggsave(paste("./output/", outName, "/", outName, "_key_feats.png", sep = ""), width = 9, height = 12)


### Fig supp 1b: plot key fibroblast feats
features <- c("FAP", "ACTA", 
              "DCN", "LUM"
             )
titles <- c("FAP", "ACTA", 
              "DCN", "LUM")
colz <- c("black", "black",
          "black","black")
p <- prettyFeats(seu.obj = seu.obj, nrow = 2, ncol = 2, features = features, showAxis = F,smallAxis=F,
                 color = colz, order = F, pt.size = 0.0000001, title.size = 14, titles = titles, noLegend = T)
ggsave(paste("./output/", outName, "/", outName, "_key_feats_poster.png", sep = ""), width = 6, height = 6)


features <- c("AOC1","HNMT","MAOB","ALDH7A1"
             )

p <- prettyFeats(seu.obj = seu.obj, nrow = 2, ncol = 2, features = features, showAxis = F,smallAxis=F, order = F, pt.size = 0.0000001, title.size = 14, noLegend = T)
ggsave(paste("./output/", outName, "/", outName, "_key_feats_poster.png", sep = ""), width = 6, height = 6)




### Run copyKat for CNV analysis ####

# SKIP ! ### <<<<<<< results provided in cnvStat.csv, but can run code if desired
# seu.obj <- readRDS(file = "./output/s3/naive_n6_res0.8_dims45_dist0.3_neigh30_S3.rds")
# dout <- "./output/copyKat/"

# seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./tumor_naive_6.csv", groupBy = "clusterID", metaAdd = "majorID")

# seu.sub.list <- SplitObject(seu.obj, split.by = "orig.ident")
# 
# df.list <- list()
# cntr <- 0
# lapply(seu.sub.list, function(seu.obj) {
    
#     rm(copykat_run)
#     gc()
    
#     #check the distribution of cells that have more than 2000 UMIs
#     table(seu.tumor@meta.data$nCount_RNA>2000,seu.tumor$clusterID)
    
#     #select best cells
#     seu.tumor.clean <- seu.tumor[,seu.tumor@meta.data$nCount_RNA>2000]
    
#     #try with "normal" cells called
#     norm_cells <- colnames(subset(seu.obj, majorID == "myeloid" | majorID == "tcell"))
#     cntr <- cntr+1
#     rawData <- as.matrix(seu.obj@assays$RNA@counts)

#     outName <- unique(seu.obj$orig.ident)
    
#     #run copyKAT
#     copykat_run <- copykat(rawmat=rawData, 
#                            id.type="S", 
#                            ngene.chr=5, 
#                            win.size=25, 
#                            KS.cut=0.1, 
#                            sam.name=outName, 
#                            distance="euclidean", 
#                            norm.cell.names=norm_cells, 
#                            n.cores=12,
#                            output.seg="FALSE")

#     pred.score <- data.frame(copykat_run$prediction)
#     pred.score$cell.names <- NULL 
#     CNA.val <- data.frame(copykat_run$CNAmat)

#     seu.obj <- AddMetaData(
#       object = seu.obj,
#       metadata = pred.score,
#       col.name = 'copyKatPred'
#       )
    
#     #df <- as.data.frame(seu.obj@meta.data)
#     #df.list[[seu.obj@meta.data]] <- df
    
#     pi <- DimPlot(seu.obj, 
#                   reduction = "umap", 
#                   group.by = "copyKatPred",
#                   pt.size = 0.5,
#                   label = FALSE,
#                   label.box = FALSE
#                  )
    
#     df <- data.frame(matrix(ncol = 0, nrow = dim(seu.obj)[2]))
    
#     df[["code"]] <- rownames(seu.obj@meta.data)
#     df[["cnvStat"]] <- seu.obj@meta.data$copyKatPred
    
#     df.list[[cntr]] <- df
    
#     outfile <- paste(dout, outName, "_copyKatPred_seu.rds", sep = "")
#     saveRDS(seu.obj, file = outfile)

#     outfile <- paste(dout, outName, "_uMAP_by_ploidy.png", sep = "")
#     #save the final plot!
#     ggsave(plot = pi, outfile, width = 10, height = 10)
    
# })
       
# cellCounts <- do.call(rbind, df.list)
# outfile <- paste(dout,"cnvStatz.csv", sep = "")
# write.csv(cellCounts, file = outfile)

# files <- list.files(path = "./output/copyKat_2/", pattern="copyKatPred_seu.rds", all.files=FALSE,
#                         full.names=T)

# create.seu.call <- function(x) {
#         readRDS(x)
#     }

# #load all seurat objects and store in a large list
# katOuts.list <- mapply(create.seu.call, files)

# df.list <- list()

# for (obj in seq(1:length(katOuts.list))) {
#     seu <- katOuts.list[[obj]]
#     df <- data.frame(matrix(ncol = 0, nrow = dim(seu)[2]))
# 
#     df[["code"]] <- rownames(seu@meta.data)
#     df[["cnvStat"]] <- seu@meta.data$copyKatPred
#     
#     df.list[[obj]] <- df
# }
# 
# cellCounts <- do.call(rbind, df.list)
# outfile <- paste("./output/copyKat_2/cnvStat.csv", sep = "")
# write.csv(cellCounts, file = outfile)


### Fig 1e: plot CopyKAT output
cellCounts <- read.csv("./output/copyKat_2/cnvStat.csv")
cellCounts$X <- NULL
cnts <- cellCounts$cnvStat
names(cnts) <- cellCounts$code

seu.obj <- AddMetaData(
      object = seu.obj,
      metadata = cnts,
      col.name = 'copyKatPred'
      )

#exculde low quality sample
seu.obj.sub <- subset(seu.obj,invert = T,
                  subset = 
                  orig.ident ==  "run_count_tumor_no_tx_7")


#create the plot
pi <- DimPlot(seu.obj.sub, 
                  reduction = "umap", 
                  group.by = "copyKatPred",
                  pt.size = 0.25,
                  label = FALSE,
                  label.box = FALSE,
               shuffle = T
                 )
pi <- formatUMAP(pi) + theme(legend.position = c(0.01, 0.95)) + theme(axis.title = element_blank(),
                                                 panel.border = element_blank())
ggsave(paste("./output/", outName, "/", outName, "_uMAP_by_ploidy_noAxis.png", sep = ""),width = 7,height=7)


axes <- ggplot() + labs(x = "UMAP1", y = "UMAP2") + 
theme(axis.line = element_line(colour = "black", 
                               arrow = arrow(angle = 30, length = unit(0.1, "inches"),
                                             ends = "last", type = "closed"),
                              ),
      axis.title.y = element_text(colour = "black", size = 20),
      axis.title.x = element_text(colour = "black", size = 20),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()
     )

p <- pi + inset_element(axes,left= 0,
  bottom = 0,
  right = 0.25,
  top = 0.25,
                       align_to = "full")
ggsave(paste0("./output/", outName, "/", outName, "_uMAP_by_ploidy.png"), width = 7, height = 7)


###### INFER CNV ANALYSIS AND PLOTTING IN RESPONSE TO REVIEWERS ######

fileNames <- list.files(path = "./output_cnvIndv/", pattern="output_noNeuts_cleanData_run_count", all.files=FALSE,
                    full.names=F)

files <- paste0( "./output_cnvIndv/", fileNames, "/infercnv.observations.txt")

df.list <- lapply(files, read.table, sep = " ")

files <- paste0( "./output_cnvIndv/", fileNames, "/infercnv.references.txt")
df.list2 <- c(df.list, lapply(files, read.table, sep = " "))

# cnvStatus <- unlist(lapply(df.list, function(df){df %>% filter(rownames(.) == "MYC") %>% -1}))
cnvStatus <- unlist(lapply(df.list2, function(df){df <- colSums(abs(1-df))}))
# cnvStatus <- setNames(names(cnvStatus), cnvStatus)
names(cnvStatus) <- gsub("[.]","-",names(cnvStatus))

df <- as.data.frame(cnvStatus)
df <- df %>% mutate(qunatile = case_when(cnvStatus < quantile(cnvStatus)[2] ~ "q1 (low CNV)",
                                         cnvStatus > quantile(cnvStatus)[2] & cnvStatus < quantile(cnvStatus)[3] ~ "q2",
                                         cnvStatus > quantile(cnvStatus)[3] & cnvStatus < quantile(cnvStatus)[4] ~ "q3",
                                         cnvStatus > quantile(cnvStatus)[4] ~ "q4 (high CNV)"),
                    samID = substr(rownames(df),20,20)
                   ) %>% group_by(samID) %>% mutate(qunatile_grpd = case_when(cnvStatus < quantile(cnvStatus)[2] ~ "q1 (low CNV)",
                                                                              cnvStatus > quantile(cnvStatus)[2] & cnvStatus < quantile(cnvStatus)[3] ~ "q2",
                                                                              cnvStatus > quantile(cnvStatus)[3] & cnvStatus < quantile(cnvStatus)[4] ~ "q3",
                                                                              cnvStatus > quantile(cnvStatus)[4] ~ "q4 (high CNV)"),
                                                    cnvCall = ifelse(qunatile_grpd == "q4 (high CNV)", "aneuploid", "diploid")
                   )



seu.obj <- readRDS(file = "./output/s3/canine_naive_n6_annotated.rds")

cnvStatus1 <- df$cnvStatus
names(cnvStatus1) <- names(cnvStatus)
seu.obj <- AddMetaData(seu.obj, metadata = cnvStatus1, col.name = "cnvStatus")

cnvStatus2 <- df$qunatile
names(cnvStatus2) <- names(cnvStatus)
seu.obj <- AddMetaData(seu.obj, metadata = cnvStatus2, col.name = "quant")

cnvStatus3 <- df$cnvCall
names(cnvStatus3) <- names(cnvStatus)
seu.obj <- AddMetaData(seu.obj, metadata = cnvStatus3, col.name = "cnvCall")

cnvStatus4 <- df$qunatile_grpd
names(cnvStatus4) <- names(cnvStatus)
seu.obj <- AddMetaData(seu.obj, metadata = cnvStatus4, col.name = "qunatile_grpd")


p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 1, features = "cnvStatus", showAxis = F, smallAxis=F,
                 color = "black", order = F, pt.size = 0.0000001, title.size = 14, titles = "CNV Status", noLegend = T)
ggsave(paste("./", outName, "/", outName, "_inferCNV_status.png", sep = ""), width = 6, height = 6)

cellz <- names(seu.obj$quant[!is.na(seu.obj$quant)])
cellz <- names(seu.obj$qunatile_grpd[!is.na(seu.obj$qunatile_grpd)])
seu.obj.sub <- subset(seu.obj, cells = cellz)
pi <- DimPlot(seu.obj.sub, 
              reduction = "umap", 
              group.by = "qunatile_grpd",
              pt.size = 0.1,
              split.by = "name",
              ncol = 3,
              label = F,
              label.box = F,
              shuffle = TRUE
) 
p <- formatUMAP(pi)
ggsave(paste("./output/", outName, "/", outName, "test_quant.png", sep = ""), width = 10.5, height = 7)

levs <- seu.obj@meta.data[ ,c("celltype.l1", "cnvStatus")] %>% group_by(celltype.l1) %>% summarize(mea = mean(cnvStatus)) %>% arrange(desc(mea)) %>% pull(celltype.l1)
seu.obj$celltype.l1 <- factor(seu.obj$celltype.l1, levels = levs)

levs <- seu.obj@meta.data[ ,c("celltype.l2", "cnvStatus")] %>% group_by(celltype.l2) %>% summarize(mea = mean(cnvStatus)) %>% arrange(desc(mea)) %>% pull(celltype.l2)
seu.obj$celltype.l2 <- factor(seu.obj$celltype.l2, levels = levs)

p <- stackedBar(seu.obj = seu.obj, downSampleBy = "cellSource", groupBy = "qunatile_grpd", clusters = "celltype.l2") +
 theme(legend.position = "top") + guides(fill = guide_legend(nrow = 1, byrow =T))
ggsave(paste("./output/", outName, "/", outName, "_stackedBar.png", sep = ""), width = 8, height = 12)


pi <- VlnPlot(object = seu.obj,
    pt.size = 0,
    same.y.lims = F,
    group.by = "celltype.l1",
    combine = F,
#     cols = c("#FF755F", "#FFD6C6", "#AC0535", "#EB2C31", "tomato"),
    stack = F,
    fill.by = "ident",
    flip = T,
    features = "cnvStatus"
        ) & NoLegend() & theme(axis.ticks = element_blank(),
                               axis.text.y = element_blank(),
                               axis.title.x = element_blank())

#plot <- prettyViln(plot = pi, colorData = NULL, nrow = 2, ncol = 4)
ggsave(paste("./", outName, "/", outName, "_selectViln.png", sep = ""), width = 7, height =5)



seu.obj <- readRDS(file = "./output/s3/tumor_QCfiltered_2_3000_res0.5_dims40_dist0.5_neigh50_S3.rds")

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

seu.obj <- AddMetaData(seu.obj, metadata = cnvStatus1, col.name = "cnvStatus")
seu.obj <- AddMetaData(seu.obj, metadata = cnvStatus2, col.name = "quant")
seu.obj <- AddMetaData(seu.obj, metadata = cnvStatus3, col.name = "cnvCall")
seu.obj <- AddMetaData(seu.obj, metadata = cnvStatus4, col.name = "qunatile_grpd")




cellz <- names(seu.obj$quant[!is.na(seu.obj$quant)])
cellz <- names(seu.obj$qunatile_grpd[!is.na(seu.obj$qunatile_grpd)])
seu.obj.sub <- subset(seu.obj, cells = cellz)
pi <- DimPlot(seu.obj.sub, 
              reduction = "umap", 
              group.by = "qunatile_grpd",
              pt.size = 0.2,
              split.by = "name",
              ncol = 3,
              label = F,
              label.box = F,
              shuffle = TRUE
) 
p <- formatUMAP(pi)
ggsave(paste("./output/", outName, "/", outName, "test_tumor_quant.png", sep = ""), width = 10.5, height = 7)



levs <- seu.obj@meta.data[ ,c("majorID_subWclus", "cnvStatus")] %>% group_by(majorID_subWclus) %>% summarize(mea = median(cnvStatus)) %>% arrange(desc(mea)) %>% pull(majorID_subWclus)
seu.obj$majorID_subWclus <- factor(seu.obj$majorID_subWclus, levels = levs)

pi <- VlnPlot(object = seu.obj,
    pt.size = 0,
    same.y.lims = F,
    group.by = "majorID_subWclus",
    combine = F,
#      cols = unique(seu.obj$colz),
    stack = F,
    fill.by = "ident",
#               split.by = "name",
    flip = T,
    features = "cnvStatus"
        ) & NoLegend() & theme(axis.ticks = element_blank(),
                               axis.text.y = element_blank(),
                               axis.title.x = element_blank())

#plot <- prettyViln(plot = pi, colorData = NULL, nrow = 2, ncol = 4)
ggsave(paste("./output/", outName, "/", outName, "_selectViln_tumor.png", sep = ""), width = 8, height =5)


#########################################################
#######   End high-level analysis of all cells   ########
#########################################################
