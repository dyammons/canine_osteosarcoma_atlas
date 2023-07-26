#!/usr/bin/Rscript

### Load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")


##############################
### Human homolgy analysis ###
##############################

#####
#process count maricies from Liu et. al, 2021 "Single-Cell Transcriptomics Reveals the Complexity of the Tumor Microenvironment of Treatment-Naive Osteosarcoma"
#https://doi.org/10.3389/fonc.2021.709210
load10x(din = "./human_input/", dout = "./output/", outName = "human_naive6", testQC = F,
       nFeature_RNA_high = 5500, nFeature_RNA_low = 200, percent.mt_high = 12.5, nCount_RNA_high = 75000, nCount_RNA_low = 100)


sctIntegrate(din = "./output/s1/", outName = "human_naive_n6", vars.to.regress = c("percent.mt"), nfeatures = 2000)

seu.obj <- readRDS(file = "/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/s2/human_naive_n6_seu.integrated.obj_S2.rds")
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "naive_n6_12-5", test_dims = c(50,45,40,35,30,25), algorithm = 3, prefix = "integrated_snn_res.")

seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "human_naive_n6", final.dims = 45, final.res = 0.8, stashID = "clusterID", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

# #process nature paper atlas -- not used
# load10x(din = "./1_raw_data/", dout = "./2_output/s1/", outName = "human_chen_naive6", testQC = F,
#        nFeature_RNA_high = 7000, nFeature_RNA_low = 200, percent.mt_high = 12.5, nCount_RNA_high = 75000, nCount_RNA_low = 100)

# sctIntegrate(din = "./2_output/s1/", outName = "human_chen_naive6", vars.to.regress = c("percent.mt"), nfeatures = 2000)

# seu.obj <- readRDS(file = "/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/human_chen_naive6_seu.integrated.obj_S2.rds")
# clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "naive_n6_12-5", test_dims = c(50,45,40,35,30,25), algorithm = 3, prefix = "integrated_snn_res.")

# seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "human_chen_naive_n6", final.dims = 45, final.res = 0.8, stashID = "clusterID", 
#                         algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = F,
#                         features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
#                                      "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
#                                      "CD4", "MS4A1", "PPBP","HBM")
#                        )


#load in processed data
seu.obj <- readRDS(file = "./human_output/s3/human_naive_n6_res0.8_dims45_dist0.3_neigh30_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./all_celltype_human.csv", groupBy = "clusterID", metaAdd = "celltype.l1")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./all_celltype_human.csv", groupBy = "clusterID", metaAdd = "celltype.l2")
ct.l1_human <- seu.obj$celltype.l1
ct.l2_human <- seu.obj$celltype.l2
outName <- "human_naive"

### Check QC parameters
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("./output/", outName, "/", outName, "_QC_feats.png", sep = ""), width = 9, height = 3)


### Generate violin plots
vilnPlots(seu.obj = seu.obj, groupBy = "clusterID", numOfFeats = 24, outName = "human_naive_6", outDir = "./output/viln/allCells/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T
                     )


### Fig extra: plot key feats
features <- c("PTPRC", "CD3E", "GZMB", "HLA-DRA", 
              "FLT3", "CD68", "CTSK", "S100A12",
              "COL1A1", "FBLN1", "TOP2A", "ESAM")

p <- prettyFeats(seu.obj = seu.obj, nrow = 3, ncol = 4, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("./output/", outName, "/", outName, "_key_feats.png", sep = ""), width = 15, height = 9)


### Fig extra: Create raw UMAP
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


### Fig extra: Create UMAP by celltype.l2
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "celltype.l2",
              pt.size = 0.25,
              #cols = majorColors.df$colz,
              label = F,
              label.box = F,
              shuffle = TRUE
)
pi <- formatUMAP(plot = pi) #+ NoLegend() + theme(axis.title = element_blank(),
                            #                     panel.border = element_blank())
ggsave(paste("./output/", outName, "/", outName, "_UMAPbycelltype_l2.png", sep = ""), width = 10, height = 7)


### Fig extra: Create UMAP by sample
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "orig.ident",
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = TRUE
)

pi <- formatUMAP(plot = pi)
ggsave(paste("./output/", outName, "/", outName, "_UMAPbySample.png", sep = ""), width = 10, height = 7)


### Fig extra: Complete module scoring as in original article
modulez <- list("OS" = c("ALPL", "RUNX2", "IBSP"),
                "Myeloid" = c("LYZ", "CD68"),
                "OC" = c("ACP5", "CTSK"),
                "CAF" = c("COL1A1", "FAP", "VIM"),
                "NK-T cell" = c("CD2", "CD3D", "CD3E", "CD3G", "GNLY", "NKG7", "KLRD1", "KLRB1"),
                "Endothelial" = c("EGFL7", "PLVAP", "ESAM"),
                "B cell"= c("MS4A1", "PAX5"),
                "Plasma cell"= c("JCHAIN", "MZB1")
                )

seu.obj <- AddModuleScore(seu.obj,
                          features = modulez,
                         name = "_score")
names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)

features <- names(modulez)

ecScores <- majorDot(seu.obj = seu.obj, groupBy = "clusterID",
                     features = features
                    ) + theme(legend.position = "bottom",
                              axis.title.y = element_blank(),
                              plot.margin = margin(7, 7, 0, 24, "pt"))
ggsave(paste("./output/", outName, "/", outName, "_dotPlot_ecScores.png", sep = ""), width = 6,height=6)

p <- prettyFeats(seu.obj = seu.obj, nrow = 3, ncol = 4, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("./output/", outName, "/", outName, "_modulez.png", sep = ""), width = 15, height = 9)




### Data for Fig 9b-d: Comare gene expression between clusters
btwnClusDEG(seu.obj = seu.obj, groupBy = "celltype.l2", idents.1 = "Fibroblast", idents.2 = "Endothelial", bioRep = "orig.ident",padj_cutoff = 0.05, lfcCut = 0.58, 
                        minCells = 25, outDir = paste0("./output/", outName, "/"), title = "hu_fib_vs_endo", idents.1_NAME = "hu_fib", idents.2_NAME = "hu_endo", returnVolc = F, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F
                    )

btwnClusDEG(seu.obj = seu.obj, groupBy = "celltype.l2", idents.1 = "Mature-OC", idents.2 = "CD14_monocyte", bioRep = "orig.ident",padj_cutoff = 0.05, lfcCut = 0.58, 
                        minCells = 25, outDir = paste0("./output/", outName, "/"), title = "hu_matOC_vs_cd14mono", idents.1_NAME = "hu_matOC", idents.2_NAME = "hu_cd14mono", returnVolc = F, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F
                    )

btwnClusDEG(seu.obj = seu.obj, groupBy = "celltype.l2", idents.1 = "pDC", idents.2 = "cDC2", bioRep = "orig.ident",padj_cutoff = 0.05, lfcCut = 0.58, 
                        minCells = 5, outDir = paste0("./output/", outName, "/"), title = "hu_pDC_vs_cDC2", idents.1_NAME = "hu_pDC", idents.2_NAME = "hu_cDC2", returnVolc = F, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F
                    )

#extra contrasts:
# btwnClusDEG(seu.obj = seu.obj, groupBy = "celltype.l2", idents.1 = "IFN-TAM", idents.2 = "FABP5_Macrophage", bioRep = "orig.ident",padj_cutoff = 0.05, lfcCut = 0.58, 
#                         minCells = 25, outDir = paste0("./output/", outName, "/"), title = "hu_ifnTAMC_vs_fabp5_mac", idents.1_NAME = "hu_ifnTAMC", idents.2_NAME = "hu_fabp5_mac", returnVolc = F, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F
#                     )

# btwnClusDEG(seu.obj = seu.obj, groupBy = "celltype.l2", idents.1 = "NR4A3_Macrophage", idents.2 = "TXNIP_Macrophage", bioRep = "orig.ident",padj_cutoff = 0.05, lfcCut = 0.58, 
#                         minCells = 5, outDir = paste0("./output/", outName, "/"), title = "hu_NR4A3_Macrophage_vs_TXNIP_Macrophage", idents.1_NAME = "hu_NR4A3_Macrophage", idents.2_NAME = "hu_TXNIP_Macrophage", returnVolc = F, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F
#                     )


# btwnClusDEG(seu.obj = seu.obj, groupBy = "celltype.l2", idents.1 = "CD14_monocyte", idents.2 = "TXNIP_Macrophage", bioRep = "orig.ident",padj_cutoff = 0.05, lfcCut = 0.58, 
#                         minCells = 5, outDir = paste0("./output/", outName, "/"), title = "hu_CD14_monocyte_vs_TXNIP_Macrophage", idents.1_NAME = "hu_CD14_monocyte", idents.2_NAME = "hu_TXNIP_Macrophage", returnVolc = F, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F
#                     )

# btwnClusDEG(seu.obj = seu.obj, groupBy = "celltype.l2", idents.1 = "Mature-OC", idents.2 = "Mast", bioRep = "orig.ident",padj_cutoff = 0.05, lfcCut = 0.58, 
#                         minCells = 5, outDir = paste0("./output/", outName, "/"), title = "hu_matureOC_vs_mast", idents.1_NAME = "hu_mature_OC", idents.2_NAME = "hu_mast", returnVolc = F, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F
#                     )

# btwnClusDEG(seu.obj = seu.obj, groupBy = "celltype.l2", idents.1 = "CD14_monocyte", idents.2 = "Mast", bioRep = "orig.ident",padj_cutoff = 0.05, lfcCut = 0.58, 
#                         minCells = 5, outDir = paste0("./output/", outName, "/"), title = "hu_CD14_monocyte_vs_mast", idents.1_NAME = "hu_CD14_monocyte", idents.2_NAME = "hu_mast", returnVolc = F, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F
#                     )

# btwnClusDEG(seu.obj = seu.obj, groupBy = "celltype.l2", idents.1 = "Mast", idents.2 = NULL, bioRep = "orig.ident",padj_cutoff = 0.05, lfcCut = 0.58, 
#                         minCells = 5, outDir = paste0("./output/", outName, "/"), title = "hu_Mast_v_all", idents.1_NAME = "hu_mast", idents.2_NAME = "hu_all", returnVolc = F, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F
#                     )


### Fig supp 6: IHC marker screen in human data
seu.obj <- subset(seu.obj,
                  subset = 
                  celltype.l1 ==  "Macrophage" | celltype.l1 ==  "Osteoclast" | celltype.l1 ==  "DC")

seu.obj$celltype.l2 <- droplevels(seu.obj$celltype.l2)
seu.obj$celltype.l2 <- factor(seu.obj$celltype.l2, levels = c("CD14_monocyte", "NR4A3_Macrophage","TXNIP_Macrophage","FABP5_Macrophage","IFN-TAM","Pre-OC","Mature-OC","cDC2","pDC"))

features <- c("CD68","CD163","AIF1","MSR1","KIT")


color <- "black"
pi <- VlnPlot(
    object = seu.obj,
    pt.size = 0,
    same.y.lims = T,
    group.by = "celltype.l2",
   combine = F,
    stack = T,
    cols = c("#5295D4","#D0E4FF","#3267AD","#8C94BF","#00366C","#023FA5","#9BBFDD","#0066A5","#CFADE5","#645A9F","#AE8FD1","#FF755F", "#FFD6C6", "#AC0535", "#EB2C31", "tomato"),
    fill.by = "ident",
    flip = T,
    features = features
        ) + NoLegend()+ theme(axis.text.x = element_text(angle = 45, hjust = 1),
                              axis.title.x = element_blank(),
                              plot.margin = margin(t = 7, r = 7, b =7, l = 34, unit = "pt")
                           ) 

pi <- pi + theme(axis.title.y = element_blank(),
        strip.text.y.left = element_text(angle=0, size = 1),
        strip.placement = "outside")

ggsave(paste("./output/", outName, "/", outName, "_viln_IHC.png", sep = ""), height = 4)


################################
### Integrate human with dog ###
################################


#read in processed objects
seu.obj.human <- readRDS(file = "./human_output/s3/human_naive_n6_res0.8_dims45_dist0.3_neigh30_S3.rds")
seu.obj.dog <- readRDS(file = "./output/s3/canine_naive_n6_annotated.rds")

#merge then integrate
seu.obj.crossSpecies <- merge(seu.obj.dog,seu.obj.human)
indReClus(seu.obj = seu.obj.crossSpecies, outDir = "./output/s2/", subName = "cross_Species_3000", preSub = T, nfeatures = 3000,
                      vars.to.regress = "percent.mt"
                       )

#check ideal clustering parameteres
seu.obj <- readRDS(file = "./output/s2/cross_Species_3000_S2.rds")
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "cross_Species_3000", test_dims = c(50,45,40,35,30,25), algorithm = 3, prefix = "integrated_snn_res.")

#complete dim reduction and visulaization
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "cross_Species_3000", final.dims = 45, final.res = 0.6, stashID = "clusterID_2", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )

#load in integrated dataset
seu.obj <- readRDS(file = "./human_output/s3/cross_Species_2000_res0.6_dims45_dist0.3_neigh30_S3.rds")
seu.obj$cellSource <- ifelse(grepl("no_tx", seu.obj$orig.ident), "Canine", "Human")
seu.obj$name <- ifelse(grepl("no_tx", seu.obj$orig.ident), seu.obj$name, seu.obj$orig.ident)
outName <- "human_dog_naive"


#bring annotations over to ensure they are up-to-date
seu.obj.hu <- readRDS(file = "./human_output/s3/human_naive_n6_res0.8_dims45_dist0.3_neigh30_S3.rds")
seu.obj.hu <- loadMeta(seu.obj = seu.obj.hu, metaFile = "./all_celltype_human.csv", groupBy = "clusterID", metaAdd = "celltype.l1")
seu.obj.hu <- loadMeta(seu.obj = seu.obj.hu, metaFile = "./all_celltype_human.csv", groupBy = "clusterID", metaAdd = "celltype.l2")

seu.obj.hu$celltype.l2 <- droplevels(seu.obj.hu$celltype.l2)
# #calculate defining features of each human cluster
# vilnPlots(seu.obj = seu.obj.hu, inFile = NULL, groupBy = "celltype.l2", numOfFeats = 24, outName = "human",
#                       outDir = "./output/viln/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS"), assay = "RNA", 
#                       min.pct = 0.25, only.pos = T, resume = F, resumeFile = NULL
#                      )

seu.obj.k9 <- readRDS(file = "./output/s3/canine_naive_n6_annotated.rds")

ct.l1_human <- seu.obj.hu$celltype.l1[!is.na(seu.obj.hu$celltype.l1)]
ct.l1_can <- seu.obj.k9$celltype.l1[!is.na(seu.obj.k9$celltype.l1)]
names(ct.l1_can) <- paste0(names(seu.obj.k9$celltype.l1), "_1")
names(ct.l1_human) <- paste0(names(seu.obj.hu$celltype.l1), "_2")
seu.obj <- AddMetaData(seu.obj, metadata = c(ct.l1_can,ct.l1_human), col.name = "celltype.l1_merged")


ct.l2_human <- seu.obj.hu$celltype.l2[!is.na(seu.obj.hu$celltype.l2)]
ct.l2_can <- seu.obj.k9$celltype.l2[!is.na(seu.obj.k9$celltype.l2)]
names(ct.l2_can) <- paste0(names(seu.obj.k9$celltype.l2), "_1")
names(ct.l2_human) <- paste0(names(seu.obj.hu$celltype.l2), "_2")
seu.obj <- AddMetaData(seu.obj, metadata = c(ct.l2_can,ct.l2_human), col.name = "celltype.l2_merged")

seu.obj$type <- ifelse(grepl("Canine", seu.obj$cellSource), paste0("can_",seu.obj$celltype.l2_merged), paste0("hu_",seu.obj$celltype.l2_merged))    

#ensure it looks correct
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "type",
              pt.size = 0.25,
              label = F,
              label.box = F,
              shuffle = TRUE,
              na.value = NA
)
pi <- formatUMAP(plot = pi) + theme(legend.position = "top")
ggsave(paste("./output/", outName, "/", outName, "_raw_UMAP_test.png", sep = ""), width = 10, height = 10)

#remove NAs from canine dataset
seu.obj <- subset(seu.obj, subset = type == "can_NA", invert = T)
seu.obj$type <- as.factor(seu.obj$type)
levels(seu.obj$type) <- gsub(" ","_",levels(seu.obj$type))
levels(seu.obj$type)

### Fig 9a: complete hierarchical clustering
#extract data
metadata <- seu.obj@meta.data
expression <- as.data.frame(t(seu.obj@assays$integrated@data)) #use sct count
expression$anno_merge <- seu.obj@meta.data[rownames(expression),]$type

#calc average expression
clusAvg_expression <- expression %>% group_by(anno_merge) %>% summarise(across(where(is.numeric), mean)) %>% as.data.frame()
rownames(clusAvg_expression) <- clusAvg_expression$anno_merge
clusAvg_expression$anno_merge <- NULL
                  
#do clustering
M <- (1- cor(t(clusAvg_expression),method="pearson"))/2
hc <- hclust(as.dist(M),method="complete")
mat <- scale(as.matrix(clusAvg_expression))

#create plot
node.df <- as.data.frame(list("node" = c(70,71,72,61,60,66,68,69),
                              "colour" = rep(c("lightblue","lightgrey"),4),
                              "cellType" = c("T & NK cells","B cell & pDCs","Osteoclasts","Osteoblasts & stroma", "Cycling", "TAM","TIM & Neuts", "DCs")))

ggtree(as.phylo(hc)) + geom_tiplab(size = 6) + xlim(NA, 4)+ geom_hilight(data=node.df, mapping=aes(node=node, fill=colour)
                                                   ) + scale_fill_manual(values=c("lightblue","lightgrey")
                                                                        ) + geom_cladelab(data=node.df, mapping=aes(node=node, label=cellType), align=TRUE, angle=270, offset = 0.75, hjust = 0.5, vjust = 0,fontsize = 5) +  NoLegend()

ggsave(paste0("./output/", outName, "/", outName, "_ggTree.png"), width = 9, height = 8, scale = 2) 


### Fig 9b-d: compelte cross species dge
p <- crossSpeciesDEG(pwdTOspecies1 = "/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/human_naive/hu_fib_vs_hu_endo_all_genes.csv", 
                     pwdTOspecies2 = "/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/human_dog_naive/can_fib_vs_can_endo_all_genes.csv", 
                     species1 = "Human", species2 = "Canine", cONtrast = c("fibroblast","endothelial"),
                     nlab = 10, nlab_axis = 2, colUp = "red", colDwn = "blue", seed = 12,
                     overlapGenes = overlapGenes
                    )
ggsave(paste("./output/", outName, "/", outName, "_pSigned_cfibVend.png", sep = ""), width = 7,height=5)


p <- crossSpeciesDEG(pwdTOspecies1 = "/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/human_naive/hu_matOC_vs_hu_cd14mono_all_genes.csv", 
                     pwdTOspecies2 = "/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/human_dog_naive/can_matOC_vs_can_tim_all_genes.csv", 
                     species1 = "Human", species2 = "Canine", cONtrast = c("mature OC", "TIM"),
                     nlab = 10, nlab_axis = 2, colUp = "red", colDwn = "blue", seed = 12,
                     overlapGenes = overlapGenes
                    )
ggsave(paste("./output/", outName, "/", outName, "_pSigned_cOCvTIM.png", sep = ""), width = 7,height=5)



p <- crossSpeciesDEG(pwdTOspecies1 = "/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/human_naive/hu_pDC_vs_hu_cDC2_all_genes.csv", 
                     pwdTOspecies2 = "/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/human_dog_naive/can_pDC_vs_can_cDC2_all_genes.csv", 
                     species1 = "Human", species2 = "Canine", cONtrast = c("pDC", "cDC2"), seed = 666,
                     nlab = 10, nlab_axis = 2, colUp = "red", colDwn = "blue",
                     overlapGenes = overlapGenes,
                     hjustvar = c(1,0.85,0,0.5,1,0,1,0),
                     vjustvar = c(0.75,1,0.5,0,0,1,1,0)
                    )
ggsave(paste("./output/", outName, "/", outName, "_pSigned_cpDCvcDC2.png", sep = ""), width = 7,height=5)
