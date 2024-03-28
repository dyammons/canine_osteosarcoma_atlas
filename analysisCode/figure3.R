#!/usr/bin/Rscript

#load custom functions & packages
source("./customFunctions.R")

####################################################
#######   subset and recluster on t cells   ########
####################################################

#load in metadata for all cells
seu.obj <- readRDS(file = "../output/s3/naive6_QCfilter_2000Feats_res0.8_dims45_dist0.35_neigh40_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "majorID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "freqID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "id")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "tumorO")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "orig.ident_2", metaAdd = "name")

sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "name", metaAdd = "colz")

#subset on t cells
seu.obj <- subset(seu.obj,
                  subset = 
                  majorID ==  "tcell")

#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj, outDir = "../output/s2/", subName = "tcell_2500Feats_qcFiltered", preSub = T, nfeatures = 2500,
                      vars.to.regress = "percent.mt"
                       )

# seu.obj <- readRDS(file = "../output/s2/tcell_2500Feats_qcFiltered_S2.rds")
clusTree(seu.obj = seu.obj, dout = "../output/clustree/", outName = "tcell_2500Feats_qcFiltered", test_dims = 40, algorithm = 3, prefix = "integrated_snn_res.")

#complete initial dim reduction and vis
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = "tcell_2500Feats_qcFiltered", final.dims = 40, final.res = 0.5, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 50, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

#remove a low quality cell population
seu.obj.sub <- subset(seu.obj,invert = T,
                  subset = 
                  clusterID_sub ==  "7")

clusTree(seu.obj = seu.obj.sub, dout = "../output/clustree/", outName = "tcell_2500Feats_qcFiltered", test_dims = 40, algorithm = 3, prefix = "integrated_snn_res.")

#finalize dim reduction
seu.obj <- dataVisUMAP(seu.obj = seu.obj.sub, outDir = "../output/s3/", outName = "tcell_2500Feats_qcFiltered", final.dims = 40, final.res = 0.6, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 50, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )


### Code extra
# #try appraoch to determine where the cycling cells are comming from
# seu.obj.regCYCLE <- indReClus(seu.obj = seu.obj, outDir = "../output/s2/", subName = "tcell_2500Feats_qcFiltered", preSub = T, nfeatures = 2500,
#                       vars.to.regress = c("S.Score", "G2M.Score","percent.mt")
#                        )

# seu.obj.regCYCLE <- dataVisUMAP(seu.obj = seu.obj.regCYCLE, outDir = "../output/s3/", outName = "tcell_2500Feats_qcFiltered", final.dims = 40, final.res = 0.6, stashID = "clusterID_sub2", 
#                         algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 50, assay = "integrated", saveRDS = F,
#                         features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
#                                      "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
#                                      "CD4", "MS4A1", "PPBP","HBM")
#                        )



# pi <- DimPlot(seu.obj.regCYCLE, 
#               reduction = "umap", 
#               group.by = "clusterID_sub2",
#               cols = majorColors.df$colz,
#               pt.size = 0.25,
#               label = TRUE,
#               label.box = TRUE,
#               shuffle = TRUE
# )
# pi <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = majorColors.df$labCol) + NoLegend() + theme(axis.title = element_blank(),
#                                                  panel.border = element_blank())
# ggsave(paste("../output/", outName, "/", outName, "_rawUMAP.png", sep = ""), width = 7, height = 7)


# pi <- DimPlot(seu.obj.regCYCLE, 
#               reduction = "umap", 
#               group.by = "clusterID_sub",
#               cols = majorColors.df$colz,
#               pt.size = 0.25,
#               label = TRUE,
#               label.box = TRUE,
#               shuffle = TRUE
# )
# pi <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = majorColors.df$labCol) + NoLegend() + theme(axis.title = element_blank(),
#                                                  panel.border = element_blank())
# ggsave(paste("../output/", outName, "/", outName, "_rawUMAP.png", sep = ""), width = 7, height = 7)

# Idents(seu.obj.regCYCLE) <- "clusterID_sub"
# highlight <- colnames(seu.obj.regCYCLE)[colnames(seu.obj.regCYCLE) %in% WhichCells(seu.obj.regCYCLE, ident = "4")]

# pi <- DimPlot(seu.obj.regCYCLE,#seu.integrated.obj, 
#         reduction = "umap", 
#         cols = "grey",
#         pt.size = 0.5,
#         label = F,
#         cells.highlight= highlight,
#         cols.highlight = "orchid",
#         order = F
#         #label.box = TRUE
#  )

# p <- formatUMAP(pi) + NoLegend()
# ggsave(paste("../output/", outName, "/", outName, "_supp_highlight_trdcPos.png", sep = ""), width =7, height = 7)


# features <- c("CD3E","GZMA", "CD40LG", "FOXP3","CXCL13","IL7R",
#               "CD4","GZMB", "DLA-DRA", "CTLA4","LAG3","SELL",
#                "CD8A", "GZMK", "CXCR4", "TNFRSF4","PDCD1","TOP2A")
# p <- prettyFeats(seu.obj = seu.obj.regCYCLE, nrow = 3, ncol =  6, features = features, 
#                  color = "black", order = T, pt.size = 0.5, title.size = 20, noLegend = T)
# ggsave(paste("../output/", outName, "/", outName, "_key_feats.png", sep = ""), width = 14, height = 8, scale = 2)




# vilnPlots(seu.obj = seu.obj.regCYCLE, groupBy = "clusterID_sub2", numOfFeats = 24, outName = "tcell_2500Feats_seu.obj.regCYCLE_qcFiltered", outDir = "../output/viln/tcells/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS"), assay = "RNA", 
#                       min.pct = 0.25, only.pos = T
#                      )


##########################################
#######   begin t cell analysis   ########
##########################################

#load metadata -- see line 168 if loading in data from Zenodo
seu.obj <- readRDS(file = "../output/s3/tcell_2500Feats_qcFiltered_res0.6_dims40_dist0.3_neigh50_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "orig.ident_2", metaAdd = "name")
sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./metaData/refColz.csv", groupBy = "name", metaAdd = "colz")
sorted_labels <- paste(sort(as.integer(levels(seu.obj$clusterID_sub))))
seu.obj$clusterID_sub <- factor(seu.obj$clusterID_sub, levels = sorted_labels)
outName <- "tCells_naive6"
Idents(seu.obj) <- "clusterID_sub"

#set idents based on inital analysis
Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "CD8_ex", "1" = "CD8_eff", 
                                   "2" = "CD4_act", "3" = "CD4_reg",
                                   "4" = "T_cycling","5" = "CD4_fh",
                                   "6" = "CD8_SPP1_hi", "7" = "CD4_naive", 
                                   "8" = "T_IFN", "9" = "NK")
                       )


seu.obj$majorID <- Idents(seu.obj)

seu.obj$majorID <- factor(seu.obj$majorID, levels = c("CD4_naive","CD4_act","CD4_reg","CD4_fh",
                                                "CD8_eff", "CD8_ex", "CD8_SPP1_hi",
                                                "T_IFN","T_cycling","NK"))

#export the annotated dataset for Zenodo - no need to run
# saveRDS(seu.obj, "../output/s3/tcell_subset_annotated.rds")

### If loading from Zenodo repository, can start here
# seu.obj <- readRDS("../output/s3/tumor_subset_annotated.rds")
# outName <- "tCells_naive6"


### Fig extra = QC check
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("../output/", outName, "/", outName, "_QC_feats.png", sep = ""), width = 9, height = 3)


### Create violin plots for cell ID
vilnPlots(seu.obj = seu.obj, groupBy = "majorID", numOfFeats = 24, 
          outName = "tcell_2500Feats_qcFiltered", outDir = "../output/viln/tcells/", 
          outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA",
          min.pct = 0.25, only.pos = T
         )


### Fig 3a - Create raw UMAP

#create df with cell type colors for pretty plotting
majorColors.df <- as.data.frame(levels(seu.obj$clusterID_sub))
colnames(majorColors.df) <- "ClusterID"
majorColors.df$colz <- c("#0B4151","#009DA5","#13B9AD","#47D2AF",
                         "#82E6AE", "#B7F1B2", "#236474", 
                         "#6CAFA3", "#259E98", "#00656B")

majorColors.df$title <- "Cell Type"
majorColors.df$labCol <- c("white","black","black","black",
                           "black","black","white",
                           "black","black","white")

#plot the data
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_sub",
              cols = majorColors.df$colz,
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              shuffle = TRUE
)
pi <- cusLabels(plot = pi, shape = 21, size = 8, alpha = 0.8, labCol = majorColors.df$labCol) + NoLegend() + theme(axis.title = element_blank(),
                                                 panel.border = element_blank())
ggsave(paste("../output/", outName, "/", outName, "_rawUMAP.png", sep = ""), width = 7, height = 7)


### Fig extra - Create raw UMAP with orig clusID
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


### Fig extra - Create UMAP by sample
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "orig.ident_2",
              split.by = "orig.ident_2",
              pt.size = 0.25,
              ncol = 2,
              label = F,
              label.box = F,
              shuffle = TRUE
)
pi <- formatUMAP(plot = pi) + NoLegend()
ggsave(paste("../output/", outName, "/", outName, "_UMAPbySubtype.png", sep = ""), width = 7, height = 10.5)


### Fig extra - Make stacked bar graph
p <- stackedBar(seu.obj = seu.obj, downSampleBy = "name", groupBy = "name", 
                clusters = "clusterID_sub") + scale_fill_manual(labels = levels(seu.obj$name), 
                                                                values = levels(seu.obj$colz)
                                                               ) + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 1, byrow =T)
                                                                                                             )
ggsave(file = paste("../output/", outName, "/", outName, "_stackedBar.png", sep = ""), width = 8, height = 12)


### Fig 3b - key feature plots
features <- c("CD3E","GZMA", "CD40LG", "FOXP3","CXCL13","IL7R",
              "CD4","GZMB", "DLA-DRA", "CTLA4","IL21","SELL",
               "CD8A", "GZMK", "CXCR4", "TNFRSF4","CD70","TOP2A")
p <- prettyFeats(seu.obj = seu.obj, nrow = 3, ncol =  6, features = features, 
                 color = "black", order = T, pt.size = 0.5, title.size = 20, noLegend = T)
ggsave(paste("../output/", outName, "/", outName, "_key_feats.png", sep = ""), width = 14, height = 8, scale = 2)


### Fig 3c - Use enrichment scoring help ID cells
#load in gene lists as a named list
modulez <- list("Naïve" = c("CCR7", "LEF1", "SELL", "TCF7"),
                "Cytotoxic" = c("CST7", "GZMA", "GZMB", "IFNG", "PRF1", "TNFSF10"),
                "Regulatory" = c("IL2RA", "IL4R", "IL7", "TGFB1", "TGFB3", "TGFBI", "TGFBR1"),
                "Exhausted" = c("BTLA", "CTLA4", "HAVCR2", "LAG3", "PDCD1", "TIGIT"),
                "Costimulatory" = c("ICOS", "CD226", "SLAMF1", "TNFRSF14", "TNFRSF25", "TNFRSF9"),
                "NK cell" = c("KLRF1", "STMN2", "NCR3", "F2RL3", "CD96", "IL2RB"),
                "Cycling T cell" = c("TOP2A", "MKI67", "RRM2", "H1-5", "DIAPH3", "TK1"),
                "IFN signature" = c("CXCL10", "IFI44", "OAS1", "ISG15", "IFI44L", "IFGGB2")
                )

#run module score
seu.obj <- AddModuleScore(seu.obj,
                          features = modulez,
                         name = "_score")
names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)

#plot the results of enrichment scores
features <- names(modulez)
ecScores <- majorDot(seu.obj = seu.obj, groupBy = "clusterID_sub",
                     features = rev(features)
                    ) + coord_flip() + theme(plot.margin = margin(3, 0, 3, 0, "pt"),
                                             axis.text.y=element_text(size=10),
                                                                   axis.title = element_blank(),
                                                                   legend.position = "right",
                                                                   legend.direction = "vertical",
                                                                   axis.text.x = element_text(angle=0, hjust = 0.5)
                                            ) + scale_y_discrete(position = "right") + scale_colour_continuous(name="Enrichment score", type = "viridis")

ggsave(paste("../output/", outName, "/", outName, "_dotPlot_ecScores.png", sep = ""), width = 6,height=4)

#plot indivdual members of each term
modulez <- c(list("Enrichment score" = names(modulez)), modulez)
labelz <- as.data.frame(names(modulez))
colnames(labelz) <- "labz"
labelz$modLen <- unname(unlist(lapply(modulez, length)))
cntr <- 0
plots <- lapply(modulez, function(x){
    cntr <<- cntr+1
    labz.df <- labelz[cntr,]

    majorDot(seu.obj = seu.obj, groupBy = "clusterID_sub",
                     features = rev(unname(unlist(x)))
                    ) + theme(axis.text.x = element_blank(),
                                                          axis.ticks = element_blank(),
                                                          legend.position = "right",
                                                                   legend.direction = "vertical",
                                                          axis.title = element_blank(),
                                                          plot.margin = margin(3, 0, 3, 0, "pt")
                                                         ) + scale_colour_distiller(palette = "RdYlBu", name='Average\nexpression', limits = c(-2.5,2.5)) + 
    geom_text(data = labz.df, aes(label = labz, y = 10.85, x = (modLen+1)/2),
             angle = 270, vjust = 0.5, hjust=0.5, size = 12*0.36) + coord_flip(ylim = c(1,10.75), clip = "off") + annotate("segment", x = -Inf, y = 10.5, xend = Inf, yend = 10.5, lineend = "round", linejoin = "bevel", linetype ="solid", colour = "grey70", alpha = 0.7,size = 0.5)
    
})

#get all the plots together
    patch <- area()
    nrow <- length(modulez)
    ncol <- 1
    counter=0
    for (i in 1:nrow) {
        for (x in 1:ncol) {
                patch <- append(patch, area(t = i, l = x, b = i, r = x))
        }
    }

#change to the color of the module scores for visual distinction & plot final
plots$`Enrichment score` <- plots$`Enrichment score` + theme(axis.text.x = element_text(angle=0, hjust = 0.5)
                            ) + scale_y_discrete(position = "right") + scale_colour_viridis() + guides(color = guide_colorbar(title = 'Module\nscore'), limits = c(-2.5,2.5))

p <- Reduce( `+`, plots ) +  plot_layout(guides = "collect", design = patch, 
                                                                             height = unname(unlist(lapply(modulez, length)))/sum(unname(unlist(lapply(modulez, length)))))

ggsave(paste("../output/", outName, "/", outName, "_dotPlot_ecScores_2.png", sep = ""), width = 3.5,height=7.25, scale = 2)


### Figure 3d - genes that define tregs
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "clusterID_sub", idents.1 = "3", idents.2 = "2", bioRep = "name",padj_cutoff = 0.05, lfcCut = 0.58, 
                        minCells = 25, outDir = paste0("../output/", outName, "/"), title = "treg_vs_tact", idents.1_NAME = "T_REG", idents.2_NAME = "T_ACT", returnVolc = T, doLinDEG = F, paired = T, addLabs = NULL,lowFilter = T, dwnSam = F
                    )

p  <- prettyVolc(plot = p_volc[[1]], rightLab = "Up in c3 (Treg)", leftLab = "Up in c2 (Tact)") + labs(x = "log2(FC) T_regulatory vs T_activated") + NoLegend()

ggsave(paste("../output/", outName, "/", outName, "_c3vc2_volcPlot.png", sep = ""), width = 7, height = 7)


### Figure 3e - genes that define tfh
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "clusterID_sub", idents.1 = "5", idents.2 = "2", bioRep = "name",padj_cutoff = 0.05, lfcCut = 0.58, 
                        minCells = 25, outDir = paste0("../output/", outName, "/"), title = "tfh_vs_tact", idents.1_NAME = "T_FH", idents.2_NAME = "T_ACT", returnVolc = T, doLinDEG = F, paired = T, addLabs = NULL,lowFilter = T, dwnSam = F
                    )

p  <- prettyVolc(plot = p_volc[[1]], rightLab = "Up in c5 (Tfh)", leftLab = "Up in c2 (Tact)") + labs(x = "log2(FC) T_follicular helper vs T_activated") + NoLegend()

ggsave(paste("../output/", outName, "/", outName, "_c5vc2_volcPlot.png", sep = ""), width = 7, height = 7)


### Figure extra - genes that define CD8_ex
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "clusterID_sub", idents.1 = "0", idents.2 = "1", bioRep = "name",
                      padj_cutoff = 0.05, lfcCut = 0.58, minCells = 25, lowFilter = T, dwnSam = F,
                      outDir = paste0("../output/", outName, "/"), 
                      title = "tc_c0_vs_c1", idents.1_NAME = "tc_c0", idents.2_NAME = "tc_c1", 
                      returnVolc = T, doLinDEG = F, paired = T, addLabs = NULL
                    )

p  <- prettyVolc(plot = p_volc[[1]], rightLab = "Up in c0", leftLab = "Up in c1") + labs(x = "log2(FC) c0 vs c1") + NoLegend()

ggsave(paste("../output/", outName, "/", outName, "_c0vc1_volcPlot.png", sep = ""), width = 7, height = 7)


### Fig 3f - calc % immune and total cells
#load in annotated object with all cells
# seu.obj.all <- readRDS(file = "../output/s3/naive6_QCfilter_2000Feats_res0.8_dims45_dist0.35_neigh40_S3.rds")
# seu.obj.all <- loadMeta(seu.obj = seu.obj.all, metaFile = "./metaData/tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "majorID")
# seu.obj.all <- loadMeta(seu.obj = seu.obj.all, metaFile = "./metaData/tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "freqID")
# seu.obj.all <- loadMeta(seu.obj = seu.obj.all, metaFile = "./metaData/tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "id")
# seu.obj.all <- loadMeta(seu.obj = seu.obj.all, metaFile = "./metaData/tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "tumorO")
# seu.obj.all <- loadMeta(seu.obj = seu.obj.all, metaFile = "./metaData/refColz.csv", groupBy = "orig.ident_2", metaAdd = "name")
# seu.obj.all <- AddMetaData(seu.obj.all, metadata = seu.obj$majorID, col.name = "clusterID_sub")

#load in final annotated dataset for most accurate analysis
seu.obj.all <- readRDS(file = "../output/s3/canine_naive_n6_annotated.rds")
seu.obj.all <- AddMetaData(seu.obj.all, metadata = seu.obj$majorID, col.name = "clusterID_sub")

#load in sample colors
colz.df <- read.csv("./metaData/refColz.csv")

#calculate percentage of all cells
pct.df <- table(seu.obj.all$clusterID_sub, seu.obj.all$name) %>% melt()
pct.df$Var.1 <- as.factor(pct.df$Var.1)

parentFreq <- table(seu.obj.all$name) %>% melt()
pct.df <- pct.df %>% left_join(parentFreq, by = c("Var.2" = "Var.1")) %>% left_join(colz.df, by = c("Var.2" = "name"))
pct.df$pct <- pct.df$value.x/pct.df$value.y*100

#reoder vars
pct.df$Var.1 <- factor(pct.df$Var.1, levels = c("CD4_naive","CD4_act","CD4_reg","CD4_fh",
                                                "CD8_eff", "CD8_ex", "CD8_SPP1_hi",
                                                "T_IFN","T_cycling","NK"))

p1 <- ggplot(pct.df, aes(x = Var.1, y = pct)) +
stat_summary(fun = mean, geom = "bar", fill = "grey", width = 0.7) +
geom_jitter(aes(colour = Var.2),width = 0.1) +
scale_y_continuous(limits = c(0, 10), expand = expansion(mult = c(0, 0))) + 
theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    text = element_text(colour = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.line = element_line(colour="black"),
    legend.key = element_rect(colour = NA, fill = NA),
    axis.ticks = element_line(colour="black"),
    axis.title = element_text(face = "bold"),
    axis.title.y = element_text(angle = 90, vjust = 2),
    axis.title.x = element_blank(),
    axis.text = element_text(face = "bold"),
      legend.background = element_rect(colour = "transparent", fill = NA),
    legend.box.background = element_rect(colour = "transparent", fill = NA),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.background = element_blank()
  ) + labs(y = "% all cells") + scale_colour_manual(labels = unique(pct.df$Var.2),
                                                             values = unique(pct.df$colz),
                                                             name = "Cell source")

ggsave(paste("../output/", outName, "/", outName, "_barchart.png", sep = ""), width = 6, height = 3)

pct.df.all <- pct.df

#subset on immune cells and re-calc percentage out of immune cells this time
seu.obj.all$immune <- ifelse(seu.obj.all$majorID == "endo" | seu.obj.all$majorID == "tumor" | seu.obj.all$majorID== "cyclingTumor", "tumor", "immune")

parentFreq <- table(seu.obj.all$name, seu.obj.all$immune) 
parentFreq <- parentFreq[,-2] %>% as.data.frame() %>% rownames_to_column()
colnames(parentFreq)[2] <- "immuneCnt"
pct.df <- pct.df %>% left_join(parentFreq, by = c("Var.2" = "rowname")) #%>% left_join(colz.df, by = c("Var.2" = "name"))
pct.df$pct <- pct.df$value.x/pct.df$immuneCnt*100

### Figue xx - pct of immune cells
p2 <- ggplot(pct.df, aes(x = Var.1, y = pct)) +
stat_summary(fun = mean, geom = "bar", fill = "grey", width = 0.7) +
# stat_summary(fun.data = mean_se, 
#              geom = "errorbar", width = 0) +
geom_jitter(aes(colour = Var.2),width = 0.1) +
scale_y_continuous(limits = c(0, 20), expand = expansion(mult = c(0, 0))) + 
theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    text = element_text(colour = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.line = element_line(colour="black"),
    axis.ticks = element_line(colour="black"),
    axis.title = element_text(face = "bold"),
    axis.title.y = element_text(angle = 90, vjust = 2),
    axis.title.x = element_blank(),
    axis.text = element_text(face = "bold"),
    legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = "transparent", fill = NA),
    legend.box.background = element_rect(colour = "transparent", fill = NA),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.background = element_blank()
  ) + labs(y = "% immune cells") + scale_colour_manual(labels = unique(pct.df$Var.2),
                                                             values = unique(pct.df$colz)) + NoLegend()

ggsave(paste("../output/", outName, "/", outName, "_barchart_immune.png", sep = ""), width = 6,height = 3)

pct.df.immune <- pct.df

p <- p2 + p1 + plot_layout(guides = 'collect')
ggsave(paste("../output/", outName, "/", outName, "_barchart_combined.png", sep = ""), width = 8,height = 3)


### Supp table - extract the summary data used in % bar charts
pct.df.immune <- pct.df.immune %>% group_by(Var.1) %>% summarize(Average = round(mean(pct),2),
                                                Range = paste0(round(min(pct),2),"-",round(max(pct),2)),
                                                Parent = "Immune cells")

pct.df.all <- pct.df.all %>% group_by(Var.1) %>% summarize(Average = round(mean(pct),2),
                                                Range = paste0(round(min(pct),2),"-",round(max(pct),2)),
                                                Parent = "All cells")

pct.df.merge <- rbind(pct.df.immune,pct.df.all)
colnames(pct.df.merge)[1] <- "Cell type"
write.csv(pct.df.merge, "../output/supplmental_data/tcell_pct_table.csv", row.names = F)


### Fig 3g - Surface marker screen
#load in required data
surface.markers <- read.csv("./metaData/surface_master.csv") %>% filter(Surfaceome.Label == "surface")
cluster.markers <- read.csv("../output/viln/tcells/tcell_2500Feats_qcFiltered_gene_list.csv")

### Supp data - surface markers
write.csv(cluster.markers %>% left_join(read.csv("./metaData/surface_master.csv")[c("UniProt.gene", "UniProt.description", "Surfaceome.Label", "Surfaceome.Label.Source")], by = c("gene" = "UniProt.gene")),
          file = "../output/supplmental_data/tcell_surf.csv", row.names = F)

#filter
select.surface <- cluster.markers %>% filter(gene %in% surface.markers$UniProt.gene) %>% group_by(cluster) %>% top_n(wt = avg_log2FC,n = 5)
feats_forHeat <- unique(select.surface$gene)
seu.obj$type <- seu.obj$majorID

#extract data for heatmap
metadata <- seu.obj@meta.data
expression <- as.data.frame(t(seu.obj@assays$RNA@data)) #use log noralized count
expression$anno_merge <- seu.obj@meta.data[rownames(expression),]$type
expression <- expression[,colnames(expression) %in% c(feats_forHeat,"anno_merge") ] #try subsetting instead if converting to a function

#calc avg expression by cluster
clusAvg_expression <- expression %>% group_by(anno_merge) %>% summarise(across(where(is.numeric), mean)) %>% as.data.frame()
rownames(clusAvg_expression) <- clusAvg_expression$anno_merge
clusAvg_expression$anno_merge <- NULL                                                                       

#scale and order matrix
mat <- scale(as.matrix(clusAvg_expression))
mat <- mat[,feats_forHeat]

#plot the data
outfile <- paste("../output/", outName, "/", outName, "_heatMap.png", sep = "")
png(file = outfile, width=6000, height=4000, res=400)
par(mfcol=c(1,1))                    
ht <- Heatmap(mat,
        name = "Scaled expression     ",
        rect_gp = gpar(col = "white", lwd = 2),
        border_gp = gpar(col = "black"),
              col = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"),
        show_column_dend = F,
              row_names_gp = gpar(fontsize = 18),
              column_names_gp = gpar(fontsize = 18),
        cluster_rows = F,
        cluster_columns = F,
        row_dend_reorder = F,
        column_names_rot = 45, 
        column_names_side = "top",
              heatmap_legend_param = list(legend_direction = "horizontal", title_position = "topleft",  title_gp = gpar(fontsize = 16), 
                                          labels_gp = gpar(fontsize = 12), legend_width = unit(6, "cm"))
             )

draw(ht, padding = unit(c(10, 2, 2, 15), "mm"), heatmap_legend_side = "top")
dev.off()


### Fig supp - Surface marker screen
p <- prettyFeats(seu.obj = seu.obj, nrow = 8, ncol =  5, features = feats_forHeat, 
                 color = "black", order = T, pt.size = 0.5, title.size = 20, noLegend = T)
ggsave(paste("../output/", outName, "/", outName, "_surface_feats.png", sep = ""), width = 8, height = 14, scale = 2)

########################################
#######   END t cell analysis   ########
########################################
