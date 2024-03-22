#!/usr/bin/Rscript

#load custom functions & packages
source("./customFunctions.R")
library(CellChat)

#load tumor class
seu.obj <- readRDS(file = "./output/s3/tumor_QCfiltered_2_3000_res0.5_dims40_dist0.5_neigh50_S3.rds")

#load cell types
Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "Osteoblast_1", "1" = "Osteoblast_2", 
                                   "2" = "Osteoblast_3", "3" = "Osteoblast_cycling_1",
                                   "4" = "Hypoxic_osteoblast","5" = "Osteoblast_cycling_2",
                                   "6" = "Fibroblast", "7" = "Osteoblast_cycling_3", 
                                   "8" = "Osteoblast_cycling_4", "9" = "IFN-osteoblast")
                       )
seu.obj$celltype.l3 <- Idents(seu.obj)

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "Osteoblast_1", "1" = "Osteoblast_2", 
                                   "2" = "Osteoblast_3", "3" = "Osteoblast_cycling",
                                   "4" = "Hypoxic_osteoblast","5" = "Osteoblast_cycling",
                                   "6" = "Fibroblast", "7" = "Osteoblast_cycling", 
                                   "8" = "Osteoblast_cycling", "9" = "IFN-osteoblast")
                       )
seu.obj$celltype.l2 <- Idents(seu.obj)

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "Osteoblast", "1" = "Osteoblast", 
                                   "2" = "Osteoblast", "3" = "Osteoblast_cycling",
                                   "4" = "Hypoxic_osteoblast","5" = "Osteoblast_cycling",
                                   "6" = "Fibroblast", "7" = "Osteoblast_cycling", 
                                   "8" = "Osteoblast_cycling", "9" = "IFN-osteoblast")
                       )
seu.obj$celltype.l1 <- Idents(seu.obj)

ct.l1 <- as.factor(seu.obj$celltype.l1)
ct.l2 <- as.factor(seu.obj$celltype.l2)
ct.l3 <- as.factor(seu.obj$celltype.l3)


### Fig extra: generate violin plots showing expression of defining markers
vilnPlots(seu.obj = seu.obj, groupBy = "celltype.l3", returnViln = F, outName = "tumor_for_CB", outDir = "./output/viln/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T
                     )


### Export data to make cell browser
ExportToCB_cus(seu.obj = seu.obj, dataset.name = "tumor", dir = "./output/cb_input/", 
               markers = "/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/viln/tumor_for_CB_gene_list.csv", 
               reduction = "umap",  colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase",
                                                   "clusterID_sub", "name", 
                                                   "celltype.l1", "celltype.l2", "celltype.l3"), skipEXPR = T,
               test = F,
                           feats = c("ALPL","COL13A1","COL3A1",
                                     "FBLN1","VEGFA","ACTA2",
                                     "DCN", "LUM", "C1R")
                          
                          )
    



#load t classifications
seu.obj <- readRDS(file = "./output/s3/tcell_2500Feats_qcFiltered_res0.6_dims40_dist0.3_neigh50_S3.rds")

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "CD8_ex", "1" = "CD8_eff", 
                                   "2" = "CD4_act", "3" = "CD4_reg",
                                   "4" = "T_cycling","5" = "CD4_fh",
                                   "6" = "CD8_SPP1_hi", "7" = "CD4_naive", 
                                   "8" = "T_IFN", "9" = "NK")
                       )
seu.obj$celltype.l3 <- Idents(seu.obj)

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "CD8 T cell", "1" = "CD8 T cell", 
                                   "2" = "CD4 T cell", "3" = "CD4 T cell",
                                   "4" = "T_cycling","5" = "CD4 T cell",
                                   "6" = "CD8 T cell", "7" = "CD4 T cell", 
                                   "8" = "T_IFN", "9" = "NK")
                       )
seu.obj$celltype.l2 <- Idents(seu.obj)

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "T cell", "1" = "T cell", 
                                   "2" = "T cell", "3" = "T cell",
                                   "4" = "T_cycling","5" = "T cell",
                                   "6" = "T cell", "7" = "T cell", 
                                   "8" = "T_IFN", "9" = "NK")
                       )
seu.obj$celltype.l1 <- Idents(seu.obj)

ct.l1 <- c(ct.l1,seu.obj$celltype.l1)
ct.l2 <- c(ct.l2,seu.obj$celltype.l2)
ct.l3 <- c(ct.l3,seu.obj$celltype.l3)


### Fig extra: generate violin plots showing expression of defining markers
vilnPlots(seu.obj = seu.obj, groupBy = "celltype.l3", returnViln = F, outName = "tcell_for_CB", outDir = "./output/viln/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T
                     )



### Export data to make cell browser
ExportToCB_cus(seu.obj = seu.obj, dataset.name = "tcell", dir = "./output/cb_input/", 
               markers = "/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/viln/tcell_for_CB_gene_list.csv", 
               reduction = "umap",  colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase",
                                                   "clusterID_sub", "name", 
                                                   "celltype.l1", "celltype.l2", "celltype.l3"), 
               skipEXPR = F, test = F,
               
               feats = c("CD3E","GZMA", "CD40LG", "FOXP3","CXCL13","IL7R",
                         "CD4","GZMB", "DLA-DRA", "CTLA4","IL21","SELL",
                         "CD8A", "GZMK", "CXCR4", "TNFRSF4","CD70","TOP2A")
                          
                          )


#load myeloid class
seu.obj <- readRDS(file = "./output/s3/dc_QCfiltered_2000_res0.3_dims35_dist0.3_neigh50_S3.rds")

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "cDC2", "1" = "mregDC", 
                                   "2" = "cDC1", "3" = "pDC",
                                   "4" = "preDC")
                       )
seu.obj$celltype.l3 <- Idents(seu.obj)
seu.obj$celltype.l2 <- Idents(seu.obj)

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "DC", "1" = "DC", 
                                   "2" = "DC", "3" = "DC",
                                   "4" = "DC")
                       )
seu.obj$celltype.l1 <- Idents(seu.obj)

ct.l1 <- c(ct.l1,seu.obj$celltype.l1)
ct.l2 <- c(ct.l2,seu.obj$celltype.l2)
ct.l3 <- c(ct.l3,seu.obj$celltype.l3)


### Fig extra: generate violin plots showing expression of defining markers
vilnPlots(seu.obj = seu.obj, groupBy = "celltype.l3", returnViln = F, outName = "dc_for_CB", outDir = "./output/viln/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T
                     )


ExportToCB_cus(seu.obj = seu.obj, dataset.name = "dc", dir = "./output/cb_input/", 
               markers = "/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/viln/dc_for_CB_gene_list.csv", 
               reduction = "umap",  colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase",
                                                   "clusterID_sub", "name", 
                                                   "celltype.l1", "celltype.l2", "celltype.l3"), 
               skipEXPR = F, test = F,
               
               feats = c("IL3RA", "PGLYRP2", "DDR2", 
                         "IGF1","RARRES2", "IGHM","IGKC",
                         "CADM1","DNASE1L3", "CLEC1B", 
                         "IDO1","IL4I1", "CCR7",
                         "CD1C","PID1","CD300H")
                          
                          )



#load OC-mac class
seu.obj <- readRDS(file = "./output/s3/macOC_QCfiltered_2500_res0.6_dims40_dist0.25_neigh40_S3.rds")

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "TAM_ACT", "1" = "TAM_INT", 
                                   "2" = "LA-TAM_SPP2_hi", "3" = "LA-TAM_C1QC_hi",
                                   "4" = "CD4-_TIM","5" = "Cycling_OC",
                                   "6" = "Mature_OC", "7" = "ANGIO_TAM", 
                                   "8" = "Cycling_OC", "9" = "CD320_OC", 
                                   "10" = "IFN-TAM", "11" = "CD4+_TIM")
                       )

seu.obj$celltype.l3 <- Idents(seu.obj)

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "TAM_ACT", "1" = "TAM_INT", 
                                   "2" = "LA-TAM_SPP2_hi", "3" = "LA-TAM_C1QC_hi",
                                   "4" = "CD4-_TIM","5" = "Cycling_OC",
                                   "6" = "Mature_OC", "7" = "ANGIO_TAM", 
                                   "8" = "Cycling_OC", "9" = "CD320_OC", 
                                   "10" = "IFN-TAM", "11" = "CD4+_TIM")
                       )
seu.obj$celltype.l2 <- Idents(seu.obj)

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "TAM", "1" = "TAM", 
                                   "2" = "TAM", "3" = "TAM",
                                   "4" = "TIM","5" = "Cycling_OC",
                                   "6" = "Mature_OC", "7" = "TAM", 
                                   "8" = "Cycling_OC", "9" = "CD320_OC", 
                                   "10" = "IFN-TAM", "11" = "TIM")
                       )
seu.obj$celltype.l1 <- Idents(seu.obj)

ct.l1 <- c(ct.l1,seu.obj$celltype.l1)
ct.l2 <- c(ct.l2,seu.obj$celltype.l2)
ct.l3 <- c(ct.l3,seu.obj$celltype.l3)



### Fig extra: generate violin plots showing expression of defining markers
vilnPlots(seu.obj = seu.obj, groupBy = "celltype.l3", returnViln = F, outName = "macOC_for_CB", outDir = "./output/viln/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T
                     )


ExportToCB_cus(seu.obj = seu.obj, dataset.name = "macOC", dir = "./output/cb_input/", 
               markers = "/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/viln/macOC_for_CB_gene_list.csv", 
               reduction = "umap",  colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase",
                                                   "clusterID_sub", "name", 
                                                   "celltype.l1", "celltype.l2", "celltype.l3"), 
               skipEXPR = F, test = F,
               
               feats = c("CXCL10", "CD40", "CCL19","CD5L","DLA-DRA",
                   "GCAT", "IFIT", "IL7R","KCNH7","NEURL3",
                   "GPNMB", "PTGR1", "CD36","SCIN","SPP2","SPP1","TREM2",
                   "C1QC", "IGF1", "NIBAN3","SERPING1","CLDN1",
                   "SELL", "ITGB7", "PYGL","NOD2","OSM",
                   "VEGFA", "PRDX6", "TIMP1", "PTGES","VEGFC",
                   "RSAD2", "OAS3", "IFIT3","ISG15","TNFSF10",
                   "LTF", "PADI3", "PTGS2","CD4","THBS1",
                        "STMN1","ATP6V0D2", "TNIP3", 
              "TOP2A","CTSK", "CD320",
              "MYBL2","HYAL1","SRM")
                          
                          )


#load in all data and extract last few classes
seu.obj <- readRDS(file = "./output/s3/naive6_QCfilter_2000Feats_res0.8_dims45_dist0.35_neigh40_S3.rds")

#select only cells annotated in this object
Idents(seu.obj) <- "clusterID"
seu.obj <- subset(seu.obj,
                  ident = c("2","16","20","21","22")
                 )

Idents(seu.obj) <- "clusterID"
seu.obj <- RenameIdents(seu.obj, c("2" = "Neutrophil", "16" = "Endothelial cell", "20" = "Plasma cell", 
                                   "21" = "B cell", "22" = "Mast cell")
                       )
seu.obj$celltype.l3 <- Idents(seu.obj)
seu.obj$celltype.l2 <- Idents(seu.obj)


Idents(seu.obj) <- "clusterID"
seu.obj <- RenameIdents(seu.obj, c("2" = "Neutrophil", "16" = "Endothelial cell", "20" = "B cell", 
                                   "21" = "B cell", "22" = "Mast cell")
                       )
seu.obj$celltype.l1 <- Idents(seu.obj)

ct.l1 <- c(ct.l1,seu.obj$celltype.l1)
ct.l2 <- c(ct.l2,seu.obj$celltype.l2)
ct.l3 <- c(ct.l3,seu.obj$celltype.l3)

#load in "all" object and transfer classes
seu.obj <- readRDS(file = "./output/s3/naive6_QCfilter_2000Feats_res0.8_dims45_dist0.35_neigh40_S3.rds")
outName <- "allCells"
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "majorID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "freqID")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./tumor_naive_6_v2.csv", groupBy = "clusterID", metaAdd = "id")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "orig.ident_2", metaAdd = "name")
sorted_labels <- sort(unique(seu.obj$name))
seu.obj$name <- factor(seu.obj$name, levels = sorted_labels)
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./refColz.csv", groupBy = "name", metaAdd = "colz")
seu.obj <- AddMetaData(seu.obj, metadata =ct.l1, col.name = "celltype.l1")
seu.obj <- AddMetaData(seu.obj, metadata =ct.l2, col.name = "celltype.l2")
seu.obj <- AddMetaData(seu.obj, metadata =ct.l3, col.name = "celltype.l3")


#select only cells annotated in this object
seu.obj$celltype.l1 <- as.factor(ifelse(is.na(seu.obj$celltype.l1),"NA",as.character(seu.obj$celltype.l1)))
seu.obj$celltype.l2 <- as.factor(ifelse(is.na(seu.obj$celltype.l2),"NA",as.character(seu.obj$celltype.l2)))
seu.obj$celltype.l3 <- as.factor(ifelse(is.na(seu.obj$celltype.l3),"NA",as.character(seu.obj$celltype.l3)))

Idents(seu.obj) <- "celltype.l3"
seu.obj <- subset(seu.obj,invert = TRUE,
                  ident = "NA"
                 )

seu.obj$celltype.l3 <- droplevels(seu.obj$celltype.l3)

#generate a cluster ID number in order of cluster size
clusterID_final <- table(seu.obj$celltype.l3) %>% as.data.frame() %>% arrange(desc(Freq)) %>%
mutate(clusterID_final=row_number()-1) %>% arrange(clusterID_final) 

#stash the numerical ID
newID <- clusterID_final$clusterID_final
names(newID) <- clusterID_final$Var1
Idents(seu.obj) <- "celltype.l3"
seu.obj <- RenameIdents(seu.obj, newID)
table(Idents(seu.obj))
seu.obj$clusterID_final <- Idents(seu.obj)

#export the data if desired
# group.df <- as.data.frame(table(seu.obj$celltype.l1,seu.obj$celltype.l2))
# group.df <- group.df[group.df$Freq > 0,]
# colnames(group.df) <- c("celltype.l1","celltype.l2","count")
# write.csv(group.df, "./celltype_df.csv", row.names = F)

#save the object
# saveRDS(seu.obj, file = "./output/s3/canine_naive_n6_annotated.rds")


#load in the final annotated dataset
seu.obj <- readRDS(file = "./output/s3/canine_naive_n6_annotated.rds")
outName <- "allCells"

# #generate violin plots for cell type level 3
vilnPlots(seu.obj = seu.obj, inFile = NULL, groupBy = "celltype.l3", numOfFeats = 24, outName = "ct.l3",
                      outDir = "./output/viln/allCells/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T, resume = F, resumeFile = NULL
                     )

cluster.markers <- read.csv("./output/viln/allCells/ct.l3_gene_list.csv")
write.csv(cluster.markers %>% left_join(read.csv("./surface_master.csv")[c("UniProt.gene", "UniProt.description", "Surfaceome.Label", "Surfaceome.Label.Source")], by = c("gene" = "UniProt.gene")),
          file = "./output/supplmental_data/all_surf.csv", row.names = F)


### Fig extra: generate violin plots showing expression of defining markers

ExportToCB_cus(seu.obj = seu.obj, dataset.name = "allCells", dir = "./output/cb_input/", 
               markers = "/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/viln/allCells/ct.l3_gene_list.csv", 
               reduction = "umap",  colsTOkeep = c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "Phase",
                                                   "clusterID", "name", 
                                                   "celltype.l1", "celltype.l2", "celltype.l3"), 
               skipEXPR = F, test = F,
               
               feats = c("PTPRC", "CD3E", "CD4",
                         "GZMB", "DLA-DRA", "FLT3", 
                         "CD68", "CTSK", "S100A12",
                         "ALPL", "COL1A1", "FBLN1", 
                         "TOP2A","MS4A1", "ESAM")
                          
                          )

#load in color data for pretry plotting
colz <- c(c("#6D0026","#EB2C31","#EB2C31","#FFAA93", "#F17D00", "#F6B31E", "#DABC9C", "#F49900", "#FFD6C6", "#AC0535"),
c("#0B4151","#009DA5","#13B9AD","#47D2AF","#82E6AE", "#B7F1B2", "#236474", "#6CAFA3", "#259E98", "#00656B"),
c("#5295D4","#D0E4FF","#3267AD","#8C94BF","#00366C","#023FA5","#9BBFDD","#0066A5","#CFADE5","#645A9F","#AE8FD1","#FF755F", "#FFD6C6", "#AC0535", "#EB2C31", "tomato"),
         c("#DABC9C","#8A72BD","#EDCEF7","#C89074","#D87FAB"))

names(colz) <- c(c("Osteoblast_1","Osteoblast_2", "Osteoblast_3","Osteoblast_cycling_1","Hypoxic_osteoblast","Osteoblast_cycling_2","Fibroblast","Osteoblast_cycling_3", "Osteoblast_cycling_4","IFN-osteoblast"),
c("CD4_naive","CD4_act","CD4_reg","CD4_fh","CD8_eff", "CD8_ex", "CD8_SPP1_hi","T_IFN","T_cycling","NK"),
c("CD4-_TIM", "CD4+_TIM", "ANGIO_TAM", "TAM_ACT", "TAM_INT", "LA-TAM_SPP2_hi", "LA-TAM_C1QC_hi","IFN-TAM", "CD320_OC", "Cycling_OC", "Mature_OC", "cDC2", "cDC1", "mregDC", "pDC", "preDC"),
                c("Mast cell","B cell","Plasma cell","Neutrophil","Endothelial cell"))

names(colz)[!names(colz) %in% seu.obj$celltype.l3]


### Fig extra: Create UMAP by celltype.l3
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "celltype.l3",
              pt.size = 0.25,
              cols = colz,
              label = F,
              label.box = F,
              shuffle = TRUE
)
pi <- formatUMAP(plot = pi) +theme(legend.key = element_rect(colour = "transparent", fill = "transparent"), legend.text=element_text(size=7),legend.position="bottom") + guides(colour=guide_legend(nrow=6,byrow=F, override.aes = list(size=4)))
ggsave(paste("./output/", outName, "/", outName, "_UMAPbycelltype_l3.png", sep = ""), width = 9, height = 10)


### Fig 8 legend: extract and save legend from the above UMAP
leg <- get_legend(pi)
leg_p <- as_ggplot(leg)
ggsave(leg_p, file = paste("./output/", outName, "/", outName, "_UMAP_leg.png", sep = ""), width = 9, height = 2)


#########################################################################

#run cellchat
seu.obj <- readRDS(file = "./output/s3/canine_naive_n6_annotated.rds")
outName <- "allCells"
seu.obj$immune <- ifelse(seu.obj$majorID == "endo" | seu.obj$majorID == "tumor" | seu.obj$majorID== "cyclingTumor", "tumor", "immune")


# #if running on immune only un-note
# seu.obj <- subset(seu.obj, subset = immune == "immune")

#extract data for cellchat
seu.obj$celltype.l3 <- droplevels(seu.obj$celltype.l3)
seu.obj <- NormalizeData(seu.obj)
cnts <- seu.obj@assays$RNA@data
cnts <- orthogene::convert_orthologs(gene_df = cnts,
                                        gene_input = "rownames", 
                                        gene_output = "rownames", 
                                        input_species = "dog",
                                        output_species = "human",
                                        non121_strategy = "drop_both_species") 
rownames(cnts) <- unname(rownames(cnts))

meta <- seu.obj@meta.data
cell.use <- rownames(meta)

#create cellchat object
cellchat <- createCellChat(object = cnts, meta = meta, group.by = "celltype.l3")
cellchat@idents <- factor(cellchat@idents, levels = as.character(str_sort(levels(cellchat@idents),numeric = TRUE)))
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)

#run cellchat as described in vignette
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

#save the object - if running immune unnote the top saveRDS call
# saveRDS(cellchat, "/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/allCells/allCells_cellChatobj_IMMUNE_BYfinal_230724.rds")
# saveRDS(cellchat,"/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/allCells/allCells_cellChatobjBYfinalID_ctl3_230724.rds")
saveRDS(cellchat,"/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/allCells/allCells_cellChatobjBYfinalID_ctl3_CONVERTED_231228.rds")
saveRDS(cellchat, "/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/allCells/allCells_cellChatobj_IMMUNE_BYfinal_CONVERTED_231228.rds")

#load in saved object if not already loaded in
# cellchat <- readRDS("/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/allCells/allCells_cellChatobjBYfinalID_ctl3_CONVERTED_231228.rds")
# cellchat <- readRDS("/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/allCells/allCells_cellChatobjBYfinalID_ctl3_230724.rds")

#compute cellchat netP value
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

#check relavent numbers for reporting
chat_cnts <- as.data.frame(cellchat@net$count)
source_cnts <- rowSums(chat_cnts)
target_cnts <- colSums(chat_cnts)

totalNumInts <- sum(target_cnts)
print(totalNumInts)


### Supp Fig 8a: stacked bar graphs of interactions
#get the data
source_cnts.df <- as.data.frame(source_cnts)
target_cnts.df <- as.data.frame(target_cnts)

#clean the data
cnt.df <- cbind(source_cnts.df,target_cnts.df)
cnt.df$source_pct <- cnt.df$source_cnts/(cnt.df$source_cnts+cnt.df$target_cnts)
cnt.df$target_pct <- cnt.df$target_cnts/(cnt.df$source_cnts+cnt.df$target_cnts)
cnt.df.2 <- cnt.df[,c(1,2)] %>% rownames_to_column() %>% melt()

#calculat total number of interactions
cnt.df$total <- cnt.df$source_cnts + cnt.df$target_cnts

cnt.df.2$rowname <- factor(cnt.df.2$rowname, levels=rownames(cnt.df[order(cnt.df$total, decreasing=TRUE),]))

#make the first part of plot
p1 <- ggplot(cnt.df.2, aes(x = rowname, y = value, fill = factor(variable))) +
    geom_bar(stat = "identity", position = "stack", width = 1, colour="white") +
    theme_classic() + theme(axis.text.x = element_blank(),
                            axis.title.x = element_blank(),
                            legend.position = "top",
                            legend.direction = "horizontal"
                           ) + 
    ylab(label = "Count sender/reciver") + scale_fill_discrete(name = "Interaction type", labels = c("Sender", "Receiver")) + guides(fill = guide_legend(nrow = 1))
ggsave(paste("./output/", outName, "/", outName, "_stacked_barchart.png", sep = ""), width = 7,height = 4)

#get other data
cnt.df.2 <- cnt.df[,c(3,4)] %>% rownames_to_column() %>% melt()
cnt.df.2$rowname <- factor(cnt.df.2$rowname, levels=rownames(cnt.df[order(cnt.df$total, decreasing=TRUE),]))

#make second part of plot
p2 <- ggplot(cnt.df.2, aes(x = rowname, y = value, fill = factor(variable))) +
    geom_bar(stat = "identity", position = "fill", width = 1, colour="white") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                            axis.title.x = element_blank(),
                            legend.position = "top",
                            legend.direction = "horizontal"
                           ) + 
    ylab(label = "Proportion sender/reciver") + NoLegend()

p <- p1/p2
ggsave(paste0("./output/", outName, "/", outName, "_stacked_barchart_2.png"), width = 7,height = 7)



#check valid pathways
pathways <- cellchat@netP$pathways
length(pathways)

#set pathway type groupings
# old
# immune_path <- c("SPP1","CD45","CCL","SELPLG","CD86","SELL","IL1","CD80","PARs","IL16","FASLG","FLT3","PD-L1",
#                  "CD96","CD23","OSM","CD22","CADM","MPZ","BAFF","THY1","ADGRE5","SEMA4", "CD40", "CSF")
# immuneRelated_path <- c("GALECTIN","FN1","THBS","LAMININ","SEMA7","CLEC","BAG","RESISTIN","VISFATIN","NECTIN","ICAM","NRG","NOTCH", "HSPG")
# nonImmune_path <- c("COLLAGEN","APP","PTN","VEGF","BSP","PECAM1","IGF","FGF","CDH","TENASCIN","APRIL","JAM","CDH1","ESAM","VWF","PERIOSTIN","CDH5","PDGF","NCAM","CALCR","ANGPT", "GAS")


#set pathway type groupings
immune_path <- c("SPP1","CD45","CCL","SELPLG","CD86","SELL","IL1","CD80","PARs","IL16","FASLG","FLT3","PD-L1",
                 "CD96","CD23","OSM","CADM","MPZ","BAFF","THY1","ADGRE5","SEMA4", "CD40", "CSF", "MHC-II")
immuneRelated_path <- c("FN1","THBS","LAMININ","SEMA7","CLEC","BAG","RESISTIN","VISFATIN","NECTIN","ICAM","NRG","NOTCH", "HSPG")
nonImmune_path <- c("COLLAGEN","APP","PTN","VEGF","BSP","PECAM1","IGF","CDH","TENASCIN","APRIL","JAM","CDH1","ESAM","VWF","PERIOSTIN","CDH5","PDGF","NCAM","CALCR","ANGPT", "GAS")


length(unique(c(immune_path,immuneRelated_path,nonImmune_path)))

#write the pathway classifications for supplemental table
path.df <- as.data.frame(c(immune_path,immuneRelated_path,nonImmune_path))
colnames(path.df) <- "Pathway"
path.df$Path_type <- c(rep("immune_path",length(immune_path)),rep("immuneRelated_path",length(immuneRelated_path)),rep("nonImmune_path",length(nonImmune_path)))
write.csv(path.df, file = "/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/supplmental_data/pathTypes.csv", row.names = F)


### Fig 8a/b: plot the interaction numbers and strength
#make plot using cellchat funciton
gg1 <- netAnalysis_signalingRole_scatter(cellchat, color.use = colz, title = "All interactions")
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = immune_path, color.use = colz, title = "Immune")
gg3 <- netAnalysis_signalingRole_scatter(cellchat, signaling = immuneRelated_path, color.use = colz, title = "Immune-related")
gg4 <- netAnalysis_signalingRole_scatter(cellchat, signaling = nonImmune_path, color.use = colz, title = "Non-immune")

#extract data and customize plots
gg1.df <- gg1$data
gg1.df$data_type <- "All interactions"

gg2.df <- gg2$data
gg2.df$data_type <- "Immune"

gg3.df <- gg3$data
gg3.df$data_type <- "Immune-related"

gg4.df <- gg4$data
gg4.df$data_type <- "Non-immune"

gg.df <- rbind(gg1.df,gg2.df,gg3.df,gg4.df)

colz.df <- as.data.frame(colz)
colz.df$labels <- rownames(colz.df)

gg.df <- gg.df %>% mutate(strength = x*y) %>% left_join(colz.df, by = "labels")

#set cell types to label
gg.df <- gg.df %>% group_by(data_type) %>% arrange(desc(strength)) %>% mutate(lab=ifelse(row_number() <= 5, labels, NA)) %>% ungroup()

#manually add in labels for 3 more immune related populations
gg.df <- gg.df %>% mutate(lab=ifelse(data_type != "Immune-related", lab,
                           ifelse(!is.na(lab), lab,
                                 ifelse(labels %in% c("NK", "CD8_eff","Mature_OC"),labels,NA))))


gg.df.sub <- gg.df %>% filter(data_type == "All interactions")

#make fig 4a
pi1 <- ggplot(data=gg.df.sub, aes(x = x, y = y, size=Count, color = labels, label=lab)) + 
            ggtitle("All interactions") +
            geom_point(color = gg.df.sub$colz) + 
            labs(x = "Outgoing interaction strength", y = "Incoming interaction strength") +
            geom_text_repel(max.overlaps = Inf, size=3, color = "black") + 
            theme_classic() + 
    theme(axis.title = element_text(size= 10),
          axis.text = element_text(size= 8),
          title = element_text(size= 11),
          legend.title=element_text(size=10), 
          legend.text=element_text(size=8)
         )


gg.df2 <- gg.df %>% filter(data_type != "All interactions")

#make fig 4b
pis <- lapply(c("Immune","Immune-related","Non-immune"),function(z){
    gg.df.sub <- gg.df2 %>% filter(data_type == z)

    ggplot(data=gg.df.sub, aes(x = x, y = y, size=Count, color = labels, label=lab)) + 
            ggtitle(z) +
            geom_point(color = gg.df.sub$colz) + 
            labs(x = "Outgoing interaction strength", y = "Incoming interaction strength") +
            geom_text_repel(max.overlaps = Inf, size=3, color = "black") + 
            theme_classic() + 
    theme(axis.title = element_text(size= 10),
          axis.text = element_text(size= 8),
          title = element_text(size= 11),
          legend.title=element_text(size=10), 
          legend.text=element_text(size=8)
                 )
})

pi2 <- Reduce( `+`, pis ) + plot_layout(ncol = 3, guides = 'collect') & scale_size_continuous(limits = c(min(gg.df2$Count), max(gg.df2$Count)))
pi <- pi1 + pi2 + plot_layout(ncol = 2, widths = c(0.22,0.78)) 
ggsave(pi, file = paste0("./output/", outName, "/", outName, "_interactionScater.png"), width = 12.5, height = 3)



#load in immune only cellchat
cellchat <- readRDS("/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/allCells/allCells_cellChatobj_IMMUNE_BYfinal_CONVERTED_231228.rds")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

#set group names for prety plotting
groupzNames <- c("ANGIO_TAM","IFN-TAM","LA-TAM_C1QC_hi","LA-TAM_SPP2_hi" ,"CD4+_TIM","TAM_ACT","TAM_INT","CD4-_TIM",
                 "B cell","Plasma cell" ,
                 "CD4_act","CD4_fh", "CD4_naive", "CD4_reg","CD8_eff","CD8_ex", "CD8_SPP1_hi","NK", "T_cycling","T_IFN",
                 "CD320_OC","Cycling_OC","Mature_OC",
                 "cDC1","cDC2","mregDC","pDC","preDC",
                 "Mast cell","Neutrophil")
                 
groupzNames[c("ANGIO_TAM","IFN-TAM","LA-TAM_C1QC_hi","LA-TAM_SPP2_hi" ,"CD4+_TIM","TAM_ACT","TAM_INT","CD4-_TIM")] <- "Mac_Mono"
groupzNames[c("Neutrophil","Mast cell")] <- "Other immune"
groupzNames[c("B cell","Plasma cell")] <- "B cell"
groupzNames[c("CD320_OC","Cycling_OC","Mature_OC")] <- "OC"
groupzNames[c("CD4_act","CD4_fh", "CD4_naive", "CD4_reg","CD8_eff","CD8_ex", "CD8_SPP1_hi","NK", "T_cycling","T_IFN")] <- "T cell"
groupzNames[c("cDC1","cDC2","mregDC","pDC","preDC")] <- "DC"
groupzNames <- groupzNames[31:60]

#tag pathways of insterest
pathways <- cellchat@netP$pathways
immuneSuppressive <- c("PD-L1","CD80","CD86")


### Fig 8c: plot immune suppresive networks
#extract data to customize plot
gg1 <- netAnalysis_signalingRole_scatter(cellchat, signaling = immuneSuppressive, color.use = colz, title = "Immune suppressive")
gg1.df <- gg1$data

#bring color data over
colz.df <- as.data.frame(colz)
colz.df$labels <- rownames(colz.df)
gg.df <- gg1.df %>% mutate(strength = x*y) %>% left_join(colz.df, by = "labels")

#select which datapoints to label
gg.df <- gg.df %>% arrange(desc(x)) %>% mutate(labx=ifelse(row_number() <= 3, labels, "NA")) %>% arrange(desc(y)) %>% mutate(laby=ifelse(row_number() <= 3, labels, "NA"),
                                                                                                                          lab=ifelse(labx != "NA", labx, ifelse(laby != "NA", laby,NA)))

#make the plot
pi <- ggplot(data=gg.df, aes(x = x, y = y, size=Count, color = labels, label=lab)) + 
            ggtitle("Immune suppressive") +
            geom_point(color = gg.df$colz) + 
            labs(x = "Outgoing interaction strength", y = "Incoming interaction strength") +
            geom_text_repel(max.overlaps = Inf,force = 5, max.iter = 100000, ylim = c(0.12, Inf), size=3, color = "black", seed = 77) + 
            theme_classic() + 
    theme(axis.title = element_text(size= 10),
          axis.text = element_text(size= 8),
          title = element_text(size= 11),
          legend.title=element_text(size=10), 
          legend.text=element_text(size=8)
                 )

ggsave(pi, file = paste0("./output/", outName, "/", outName, "_immSup_ints.png"), width = 3.25, height = 3)


### Fig 8d/8e + Fig supp 8b-d: circos plots & violin plots
pathwayz <- immuneSuppressive
subName <- "cellCom"
lapply(pathwayz, function(pathway){
    
    #extract required plotting data
    lrData <- as.data.frame(cellchat@LR)
    net <- cellchat@net
    prob <- net$prob
    prob <- prob[,,rownames(lrData[lrData$LRsig.pathway_name == pathway,])]
    prob.sum <- apply(prob, c(1,2), sum)

    #identify which cell types are active in pathway
    idx1 <- which(Matrix::rowSums(prob.sum) == 0)
    idx2 <- which(Matrix::colSums(prob.sum) == 0)
    idx <- intersect(idx1, idx2)
    net <- prob.sum[-idx, ]
    net <- net[, -idx]
    cellTypeOFinterest <- rownames(net)

    #grey out cell types not involved
    colz2 <- colz
    colz2[names(colz2)[!names(colz2) %in% cellTypeOFinterest]] <- "grey"

    #save the plot
    outfile <- paste("./output/", outName, "/", subName, "/", pathway ,"_cell_cell_3.png", sep = "")
    png(file = outfile, width=2500, height=2500, res=400)
    par(mfrow=c(1,1), mar = c(0,0,0,0))
    gg7 <- netVisual_aggregate(cellchat, layout = "chord", signaling = pathway, group = groupzNames, color.use = colz2, remove.isolate = F, big.gap = 5)
    dev.off()

    #get the active features in the pathway
    genez <- lapply(pathways, function(x){extractEnrichedLR(cellchat, signaling = x, geneLR.return = TRUE, enriched.only = T)[["geneLR"]]})
    names(genez) <- pathways

    Idents(seu.obj) <- "celltype.l3"
    seu.obj.sub <- subset(seu.obj, idents = cellTypeOFinterest)
    
    #plot expression using Seurat function
    pi <- VlnPlot(
        object = seu.obj.sub,
        pt.size = 0,
        same.y.lims = T,
        flip = T,
        group.by = "celltype.l3",
        fill.by = "ident",
        cols = colz2,
        stack = TRUE,
        combine = FALSE,
        features = unlist(genez[pathway])
    ) + NoLegend() + theme(axis.title.x = element_blank(),
                           axis.title.y.right = element_blank())
    
    ggsave(paste("./output/", outName, "/", subName, "/", pathway ,"_viln.png", sep = ""), height = 7, width = 7)
    
})



### Fig 4g: mregDC Tcell interactions
p <- netVisual_bubble(cellchat, sources.use = 21, targets.use = c(3,4,5,6,9,10,11,23,28), remove.isolate = FALSE)

df <- p$data
df <- df[df$pval > 1,]

df$source.target <- factor(df$source.target, levels = c("mregDC -> CD4_naive","mregDC -> CD4_act","mregDC -> CD4_reg","mregDC -> CD4_fh","mregDC -> CD8_eff","mregDC -> CD8_ex","mregDC -> CD8_SPP1_hi","mregDC -> T_IFN","mregDC -> NK"))


plot <- ggplot(data = df, mapping = aes(x = source.target, y = interaction_name_2)) +
    geom_point(mapping = aes(color = prob), size = 4) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.background = element_rect(fill = "white",colour = "gray"),
          plot.background = element_rect(fill = "white",colour = NA),
          legend.background = element_rect(fill = "transparent",colour = NA),
          legend.key = element_rect(fill = "transparent",colour = NA),
          panel.grid.major = element_line(color = "gray"), 
          panel.grid.minor = element_line(color = "gray"),
          legend.spacing.y = unit(0.5, 'cm')
          ) +
    scale_colour_viridis(name='Interaction\nprobability', limits = c(min(df$prob), max(df$prob)), labels = c("low","high"), breaks = c(min(df$prob), max(df$prob))) +
    guides(size=guide_legend(title="P value", override.aes=list(fill=NA))) +
    coord_cartesian(expand = TRUE, clip = "off")

#check path is correct
outfile <- paste("./output/", outName, "/", outName, "_netVisual_bubble_mregDC-Tcell.png", sep = "")
ggsave(file = outfile, plot, height = 4.5)

##########################################
#######   End CellChat analysis   ########
##########################################
