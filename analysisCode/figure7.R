#!/usr/bin/Rscript

#load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")
               
###################################################################
#######   subset and recluster on non-neutrophil myeloid   ########
###################################################################

##### Load in processed data and complete analysis on all cells
seu.obj <- readRDS(file = "./output/s3/naive6_QCfilter_2000Feats_res0.8_dims45_dist0.35_neigh40_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "majorID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "freqID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "id")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "tumorO")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident_2", metaAdd = "name")

sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "name", metaAdd = "colz")

#subset on non-granuloatue myeloid cells cells
seu.obj <- subset(seu.obj,
                  subset = 
                  majorID ==  "oc" | majorID ==  "mac" | majorID == "dc")

#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj, outDir = "./output/s2/", subName = "myeloid_QCfiltered_3000", preSub = T, nfeatures = 3000,
                      vars.to.regress = "percent.mt"
                       )

seu.obj <- readRDS(file = "./output/s2/myeloid_QCfiltered_3000_S2.rds")
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "myeloid_QCfiltered_3000", test_dims = c(45,40,35), algorithm = 3, prefix = "integrated_snn_res.")

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "myeloid_QCfiltered_3000", final.dims = 40, final.res = 0.6, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.2, n.neighbors = 40, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )


################################################
#######   subset and recluster on non-neutrophil myeloid   ########
################################################

#extract cell barcode annoations from DC file
seu.obj <- readRDS(file = "./output/s3/dc_QCfiltered_2000_res0.3_dims35_dist0.3_neigh50_S3.rds")

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "cDC2", "1" = "mregDC", 
                                   "2" = "cDC1", "3" = "pDC",
                                   "4" = "preDC")
                       )
ct.l3 <- Idents(seu.obj)


#extract cell barcode annoations from mac-oc file
seu.obj <- readRDS(file = "./output/s3/macOC_QCfiltered_2500_res0.6_dims40_dist0.25_neigh40_S3.rds")

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "TAM_ACT", "1" = "TAM_INT", 
                                   "2" = "LA-TAM_SPP2_hi", "3" = "LA-TAM_C1QC_hi",
                                   "4" = "CD4-_TIM","5" = "Cycling_OC",
                                   "6" = "Mature_OC", "7" = "ANGIO_TAM", 
                                   "8" = "Cycling_OC", "9" = "CD320_OC", 
                                   "10" = "IFN-TAM", "11" = "CD4+_TIM")
                       )
ct.l3 <- c(ct.l3,Idents(seu.obj))


#transfer annoations to all myeloid cell object
seu.obj <- readRDS(file = "./output/s3/myeloid_QCfiltered_3000_res0.6_dims40_dist0.2_neigh40_S3.rds")
sorted_labels <- paste(sort(as.integer(levels(seu.obj$clusterID_sub))))
seu.obj$clusterID_sub <- factor(seu.obj$clusterID_sub, levels = sorted_labels)
sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "name", metaAdd = "colz")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./myeloid.csv", groupBy = "clusterID_sub", metaAdd = "major")
outName <- "myeloid_naive6"

seu.obj$majorID_sub <- Idents(seu.obj)
seu.obj <- AddMetaData(seu.obj, metadata =ct.l3, col.name = "celltype.l3")
seu.obj$celltype.l3 <- ifelse(is.na(seu.obj$celltype.l3),"NA",seu.obj$celltype.l3)
seu.obj <- subset(seu.obj, invert = T, 
                  subset = celltype.l3 == "NA")
seu.obj <- AddMetaData(seu.obj, metadata =ct.l3, col.name = "celltype.l3")


sorted_labels <- c("CD4-_TIM", "CD4+_TIM", "ANGIO_TAM", "TAM_ACT", "TAM_INT", "LA-TAM_SPP2_hi", "LA-TAM_C1QC_hi","IFN-TAM", "CD320_OC", "Cycling_OC", "Mature_OC", "cDC2", "cDC1", "mregDC", "pDC", "preDC")

seu.obj$celltype.l3 <- factor(seu.obj$celltype.l3, levels = sorted_labels)



### Fig extra: Create raw UMAP with orig clusID
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "celltype.l3",
              pt.size = 0.25,
              label = TRUE,
              label.box = TRUE,
              shuffle = TRUE
)
pi <- formatUMAP(plot = pi)
ggsave(paste("./output/", outName, "/", outName, "_myeloid_celltype_l3.png", sep = ""), width = 10, height = 7)


vilnPlots(seu.obj = seu.obj, groupBy = "celltype.l3", numOfFeats = 24, outName = "macOC_QCfiltered_final_2500", outDir = "./output/viln/myeloid/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T
                     )

### Fig 7a: IHC marker screen

features <- c("CD68","CD163","AIF1","MSR1")


color <- "black"
pi <- VlnPlot(
    object = seu.obj,
    pt.size = 0,
    same.y.lims = T,
    group.by = "celltype.l3",
   combine = F,
    stack = T,
    cols = c("#5295D4","#D0E4FF","#3267AD","#8C94BF","#00366C","#023FA5","#9BBFDD","#0066A5","#CFADE5","#645A9F","#AE8FD1","#FF755F", "#FFD6C6", "#AC0535", "#EB2C31", "tomato"),
    fill.by = "ident",
    flip = T,
    features = features
        ) + NoLegend()+ theme(axis.text.x = element_text(angle = 45, hjust = 1),
                              axis.title.x = element_blank(),
                              plot.margin = margin(t = 7, r = 7, b =7, l = 21, unit = "pt")
                           ) 

pi <- pi + theme(axis.title.y = element_blank(),
        strip.text.y.left = element_text(angle=0, size = 1),
        strip.placement = "outside")

ggsave(paste("./output/", outName, "/", outName, "_viln1.png", sep = ""), height = 4)




### Fig 7b: Surface marker screen
surface.markers <- read.csv("./surface_master.csv") %>% filter(Surfaceome.Label == "surface")
cluster.markers <- read.csv("./output/viln/myeloid/macOC_QCfiltered_final_2500_gene_list.csv")

### Make supplemental data
write.csv(cluster.markers %>% left_join(read.csv("./surface_master.csv")[c("UniProt.gene", "UniProt.description", "Surfaceome.Label", "Surfaceome.Label.Source")], by = c("gene" = "UniProt.gene")),
          file = "./output/supplmental_data/myeloid_surf.csv", row.names = F)

select.surface <- cluster.markers %>% filter(gene %in% surface.markers$UniProt.gene) %>% group_by(cluster) %>% top_n(wt = avg_log2FC,n = 5)

feats_forHeat <- unique(select.surface$gene)

seu.obj$type <- seu.obj$celltype.l3

metadata <- seu.obj@meta.data
expression <- as.data.frame(t(seu.obj@assays$RNA@data)) #use log noralized count
expression$anno_merge <- seu.obj@meta.data[rownames(expression),]$type
expression <- expression[,colnames(expression) %in% c(feats_forHeat,"anno_merge") ] #try subsetting instead if converting to a function

clusAvg_expression <- expression %>% group_by(anno_merge) %>% summarise(across(where(is.numeric), mean)) %>% as.data.frame()
rownames(clusAvg_expression) <- clusAvg_expression$anno_merge
clusAvg_expression$anno_merge <- NULL                                                                       

M <- (1- cor(t(clusAvg_expression),method="pearson"))/2
hc <- hclust(as.dist(M),method="complete")

mat <- scale(as.matrix(clusAvg_expression))
mat <- mat[,feats_forHeat]


outfile <- paste("./output/", outName, "/", outName, "_heatMap.png", sep = "")
png(file = outfile, width=6000, height=3000, res=400)
par(mfcol=c(1,1))                    
ht <- Heatmap(mat,
        name = "Scaled expression     ",
        rect_gp = gpar(col = "white", lwd = 2),
        border_gp = gpar(col = "black"),
              col = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"),
#                col = c("-4" = "#3B9AB2", "0" = "#EBCC2A", "4" = "#F21A00"),
        show_column_dend = F,
        cluster_rows = T,
        cluster_columns = T,
        row_dend_reorder = F,
        column_names_rot = 45, 
        column_names_side = "top",
              heatmap_legend_param = list(legend_direction = "horizontal", title_position = "lefttop",  title_gp = gpar(fontsize = 16), 
                                              labels_gp = gpar(fontsize = 8), legend_width = unit(6, "cm")))

draw(ht, padding = unit(c(10, 2, 2, 15), "mm"), heatmap_legend_side = "top")
dev.off()

### Fig supp: Surface marker screen

p <- prettyFeats(seu.obj = seu.obj, nrow = 9, ncol =  6, features = feats_forHeat, 
                 color = "black", order = F, pt.size = 0.5, title.size = 16, noLegend = T)
ggsave(paste("./output/", outName, "/", outName, "_surface_feats.png", sep = ""), width = 8.5, height = 11, scale = 2)



### Fig: calc % immune and total cells
seu.obj.all <- readRDS(file = "./output/s3/naive6_QCfilter_2000Feats_res0.8_dims45_dist0.35_neigh40_S3.rds")
seu.obj.all <- loadMeta(seu.obj = seu.obj.all, metaFile = "./tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "majorID")
seu.obj.all <- loadMeta(seu.obj = seu.obj.all, metaFile = "./tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "freqID")
seu.obj.all <- loadMeta(seu.obj = seu.obj.all, metaFile = "./tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "id")
seu.obj.all <- loadMeta(seu.obj = seu.obj.all, metaFile = "./tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "tumorO")
seu.obj.all <- loadMeta(seu.obj = seu.obj.all, metaFile = "./refColz.csv", groupBy = "orig.ident_2", metaAdd = "name")

seu.obj.all <- AddMetaData(seu.obj.all, metadata = seu.obj$celltype.l3, col.name = "clusterID_sub")

colz.df <- read.csv("./refColz.csv",)

pct.df <- table(seu.obj.all$clusterID_sub, seu.obj.all$name) %>% melt()
pct.df$Var.1 <- as.factor(pct.df$Var.1)

parentFreq <- table(seu.obj.all$name) %>% melt()
pct.df <- pct.df %>% left_join(parentFreq, by = c("Var.2" = "Var.1")) %>% left_join(colz.df, by = c("Var.2" = "name"))
pct.df$pct <- pct.df$value.x/pct.df$value.y*100

### Figue xx: pct of all cells
p1 <- ggplot(pct.df, aes(x = Var.1, y = pct)) +
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
      legend.background = element_rect(colour = "transparent", fill = NA),
    legend.box.background = element_rect(colour = "transparent", fill = NA),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.background = element_blank()
  ) + labs(y = "% all cells") + scale_colour_manual(labels = unique(pct.df$Var.2),
                                                             values = unique(pct.df$colz),
                                                             name = "Cell source")

ggsave(paste("./output/", outName, "/", outName, "_barchart.png", sep = ""), width = 6, height = 3)

pct.df.all <- pct.df


seu.obj.all$immune <- ifelse(seu.obj.all$majorID == "endo" | seu.obj.all$majorID == "tumor" | seu.obj.all$majorID== "cyclingTumor", "tumor", "immune")

parentFreq <- table(seu.obj.all$name, seu.obj.all$immune) 
parentFreq <- parentFreq[,-2] %>% as.data.frame() %>% rownames_to_column()
colnames(parentFreq)[2] <- "immuneCnt"
pct.df <- pct.df %>% left_join(parentFreq, by = c("Var.2" = "rowname")) #%>% left_join(colz.df, by = c("Var.2" = "name"))
pct.df$pct <- pct.df$value.x/pct.df$immuneCnt*100

### Figue xx: pct of immune cells
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
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.background = element_blank()
  ) + labs(y = "% immune cells") + scale_colour_manual(labels = unique(pct.df$Var.2),
                                                             values = unique(pct.df$colz)) + NoLegend()

ggsave(paste("./output/", outName, "/", outName, "_barchart_immune.png", sep = ""), width = 6,height = 3)

pct.df.immune <- pct.df

p <- p2 + p1 + plot_layout(guides = 'collect') & theme(legend.key = element_rect(colour = "transparent", fill = "white"))
ggsave(paste("./output/", outName, "/", outName, "_barchart_combined.png", sep = ""), width = 8,height = 3)


pct.df.immune$Parent <- "Immune cells"
pct.df.immune$immuneCnt <- NULL
pct.df.all$Parent <- "All cells"

long.df <- rbind(pct.df.immune,pct.df.all)
long.df <- long.df[ ,c(1,2,8,9)]
colnames(long.df) <- c("Cell_type", "Dog_ID", "Percentage", "Parent")
write.csv(long.df, paste0("./output/", outName, "/", outName, "myeloid_pct_long_table.csv"), row.names = F)


### Supp table: extract the summary data used in % bar charts
pct.df.immune <- pct.df.immune %>% group_by(Var.1) %>% summarize(Average = round(mean(pct),2),
                                                                 Range = paste0(round(min(pct),2),"-",round(max(pct),2)),
                                                                 Parent = "Immune cells"
                                                                )

pct.df.all <- pct.df.all %>% group_by(Var.1) %>% summarize(Average = round(mean(pct),2),
                                                           Range = paste0(round(min(pct),2),"-",round(max(pct),2)),
                                                           Parent = "All cells"
                                                          )

pct.df.merge <- rbind(pct.df.immune,pct.df.all)
colnames(pct.df.merge)[1] <- "Cell type"
write.csv(pct.df.merge, paste0("./output/", outName, "/", outName, "myeloid_pct_table.csv"), row.names = F)
