#!/usr/bin/Rscript

### Load custom functions & packages
source("/pl/active/dow_lab/dylan/repos/K9-PBMC-scRNAseq/analysisCode/customFunctions.R")


###################################
### Human homolgy preprocessing ###
###################################

### Process count maricies from Liu et. al, 2021 "Single-Cell Transcriptomics Reveals the Complexity of the Tumor Microenvironment of Treatment-Naive Osteosarcoma"
### https://doi.org/10.3389/fonc.2021.709210

#load in their processed count matricies
load10x(din = "./human_input/", dout = "./output/", outName = "human_naive6", testQC = F,
       nFeature_RNA_high = 5500, nFeature_RNA_low = 200, percent.mt_high = 12.5, nCount_RNA_high = 75000, nCount_RNA_low = 100)

#integrate data
seu.obj <- sctIntegrate(din = "./output/s1/", outName = "human_naive_n6", vars.to.regress = c("percent.mt"), nfeatures = 2000)

#check clustering resolution
clusTree(seu.obj = seu.obj, dout = "./output/clustree/", outName = "naive_n6_12-5", test_dims = c(50,45,40,35,30,25), algorithm = 3, prefix = "integrated_snn_res.")

#dim reduction and visualization
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "./output/s3/", outName = "human_naive_n6", final.dims = 45, final.res = 0.8, stashID = "clusterID", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 30, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       )


##############################
### Human homolgy analysis ###
##############################


#load in processed data
seu.obj <- readRDS(file = "./human_output/s3/human_naive_n6_res0.8_dims45_dist0.3_neigh30_S3.rds")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./all_celltype_human.csv", groupBy = "clusterID", metaAdd = "celltype.l1")
seu.obj <- loadMeta(seu.obj = seu.obj, metaFile = "./all_celltype_human.csv", groupBy = "clusterID", metaAdd = "celltype.l2")
ct.l1_human <- seu.obj$celltype.l1
ct.l2_human <- seu.obj$celltype.l2
outName <- "human_naive"

#check QC parameters
features <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- prettyFeats(seu.obj = seu.obj, nrow = 1, ncol = 3, features = features, 
                 color = "black", order = F, pt.size = 0.0000001, title.size = 18)
ggsave(paste("./output/", outName, "/", outName, "_QC_feats.png", sep = ""), width = 9, height = 3)

#generate violin plots
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

cnts <- seu.obj.dog@assays$RNA@counts
cnts <- orthogene::convert_orthologs(gene_df = cnts,
                                        gene_input = "rownames", 
                                        gene_output = "rownames", 
                                        input_species = "dog",
                                        output_species = "human",
                                        non121_strategy = "drop_both_species") 
rownames(cnts) <- unname(rownames(cnts))

seu.obj.dog <- CreateSeuratObject(cnts, project = "humanConvert", assay = "RNA",
                                  min.cells = 0, min.features = 0, names.field = 1,
                                  names.delim = "_", meta.data = seu.obj.dog@meta.data)


#merge then integrate
seu.obj.crossSpecies <- merge(seu.obj.dog,seu.obj.human)
seu.obj <- indReClus(seu.obj = seu.obj.crossSpecies, outDir = "./output/s2/", subName = "cross_Species_3000", preSub = T, nfeatures = 3000,
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


###############################
### Analyze integrated data ###
###############################


#load in integrated dataset
seu.obj <- readRDS(file = "./output/s3/cross_Species_3000_res0.6_dims45_dist0.3_neigh30_S3.rds")
seu.obj$cellSource <- ifelse(grepl("no_tx", seu.obj$orig.ident), "Canine", "Human")
seu.obj$name <- ifelse(grepl("no_tx", seu.obj$orig.ident), seu.obj$name, seu.obj$orig.ident)
outName <- "human_dog_naive"

#bring annotations over to ensure they are up-to-date
seu.obj.hu <- readRDS(file = "./human_output/s3/human_naive_n6_res0.8_dims45_dist0.3_neigh30_S3.rds")
seu.obj.hu <- loadMeta(seu.obj = seu.obj.hu, metaFile = "./all_celltype_human.csv", groupBy = "clusterID", metaAdd = "celltype.l1")
seu.obj.hu <- loadMeta(seu.obj = seu.obj.hu, metaFile = "./all_celltype_human.csv", groupBy = "clusterID", metaAdd = "celltype.l2")

#extract cell type anotations by barcode for human
ct.l1_human <- seu.obj.hu$celltype.l1[!is.na(seu.obj.hu$celltype.l1)]
names(ct.l1_human) <- paste0(names(seu.obj.hu$celltype.l1), "_2")
seu.obj.hu$celltype.l2 <- droplevels(seu.obj.hu$celltype.l2)
ct.l2_human <- seu.obj.hu$celltype.l2[!is.na(seu.obj.hu$celltype.l2)]
names(ct.l2_human) <- paste0(names(seu.obj.hu$celltype.l2), "_2")

#clean env
rm(seu.obj.hu)
gc()

#load in annotated canine data to extract values
seu.obj.k9 <- readRDS(file = "./output/s3/canine_naive_n6_annotated.rds")

#extract cell type anotations by barcode for canine
ct.l1_can <- seu.obj.k9$celltype.l1[!is.na(seu.obj.k9$celltype.l1)]
names(ct.l1_can) <- paste0(names(seu.obj.k9$celltype.l1), "_1")
seu.obj <- AddMetaData(seu.obj, metadata = c(ct.l1_can,ct.l1_human), col.name = "celltype.l1_merged")
ct.l2_can <- seu.obj.k9$celltype.l2[!is.na(seu.obj.k9$celltype.l2)]
names(ct.l2_can) <- paste0(names(seu.obj.k9$celltype.l2), "_1")
seu.obj <- AddMetaData(seu.obj, metadata = c(ct.l2_can,ct.l2_human), col.name = "celltype.l2_merged")

#clean env
rm(seu.obj.k9)
gc()

#stash new metadata with cell type name prefixed with the specices
seu.obj$type <- ifelse(grepl("Canine", seu.obj$cellSource), paste0("can_",seu.obj$celltype.l2_merged), paste0("hu_",seu.obj$celltype.l2_merged))    
seu.obj$type <- as.factor(seu.obj$type)
levels(seu.obj$type) <- gsub(" ","_",levels(seu.obj$type))
levels(seu.obj$type)

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
node.df <- as.data.frame(list("node" = c(66,67,71,70,59,65,64),
                              "colour" = rep(c("lightblue","lightgrey"),4)[1:7],
                              "cellType" = c("T & NK","B cell","OC","Mast","Osteoblast & stromal", "TAM, TIM, & Neutrophil", "DC")))

# ggtree(as.phylo(hc))  + geom_nodelab(aes(label = as.character(seq(1:107))))+ geom_tiplab(size = 6) + xlim(NA, 4) 
ggtree(as.phylo(hc)) + geom_tiplab(size = 6) + xlim(NA, 4) + geom_hilight(data=node.df, mapping=aes(node=node, fill=colour)
                                                   ) + scale_fill_manual(values=c("lightblue","lightgrey")
                                                                        ) + geom_cladelab(data=node.df, mapping=aes(node=node, label=cellType), align=TRUE, angle=270, offset = 0.75, hjust = 0.5, vjust = 0,fontsize = 5) +  NoLegend()

ggsave(paste0("./output/", outName, "/", outName, "_ggTree.png"), width = 9, height = 8, scale = 2) 


### Fig 9b-d: Comare gene expression between clusters

#subset on the human cells only and extract required data
seu.obj.hu <- subset(seu.obj, subset = cellSource == "Human")

#calculate defining features of each human cluster
seu.obj.hu$celltype.l2_merged <- droplevels(seu.obj.hu$celltype.l2_merged)
vilnPlots(seu.obj = seu.obj.hu, inFile = NULL, groupBy = "celltype.l2_merged", numOfFeats = 24, outName = "human",
                      outDir = "./output/viln/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T, resume = F, resumeFile = NULL, returnViln = F
                     )

#complete DGE analysis between clusters
btwnClusDEG(seu.obj = seu.obj, groupBy = "celltype.l2", idents.1 = "Fibroblast", idents.2 = "Endothelial", bioRep = "orig.ident",padj_cutoff = 0.05, lfcCut = 0.58, 
                        minCells = 25, outDir = paste0("./output/", outName, "/"), title = "hu_fib_vs_endo", idents.1_NAME = "hu_fib", idents.2_NAME = "hu_endo", returnVolc = F, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F
                    )

btwnClusDEG(seu.obj = seu.obj, groupBy = "celltype.l2", idents.1 = "Mature-OC", idents.2 = "CD14_monocyte", bioRep = "orig.ident",padj_cutoff = 0.05, lfcCut = 0.58, 
                        minCells = 25, outDir = paste0("./output/", outName, "/"), title = "hu_matOC_vs_cd14mono", idents.1_NAME = "hu_matOC", idents.2_NAME = "hu_cd14mono", returnVolc = F, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F
                    )

btwnClusDEG(seu.obj = seu.obj, groupBy = "celltype.l2", idents.1 = "pDC", idents.2 = "cDC2", bioRep = "orig.ident",padj_cutoff = 0.05, lfcCut = 0.58, 
                        minCells = 5, outDir = paste0("./output/", outName, "/"), title = "hu_pDC_vs_cDC2", idents.1_NAME = "hu_pDC", idents.2_NAME = "hu_cDC2", returnVolc = F, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F
                    )

#subset on the canine cells only and extract required data
seu.obj.k9 <- subset(seu.obj, subset = cellSource == "Canine")

#calculate defining features of each canine cluster
vilnPlots(seu.obj = seu.obj.k9, inFile = NULL, groupBy = "celltype.l2", numOfFeats = 24, outName = "canine",
                      outDir = "./output/viln/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T, resume = F, resumeFile = NULL
                     )

btwnClusDEG(seu.obj = seu.obj.k9, groupBy = "celltype.l2", idents.1 = "Fibroblast", idents.2 = "Endothelial cell", bioRep = "orig.ident",padj_cutoff = 0.05, lfcCut = 0.58, 
                        minCells = 25, outDir = paste0("./output/", outName, "/"), title = "can_fib_vs_endo", idents.1_NAME = "can_fib", idents.2_NAME = "can_endo", returnVolc = F, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F
                    )

btwnClusDEG(seu.obj = seu.obj.k9, groupBy = "celltype.l2", idents.1 = "Mature_OC", idents.2 = "TIM", bioRep = "orig.ident",padj_cutoff = 0.05, lfcCut = 0.58, 
                        minCells = 25, outDir = paste0("./output/", outName, "/"), title = "can_matOC_vs_tim", idents.1_NAME = "can_matOC", idents.2_NAME = "can_tim", returnVolc = F, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F
                    )

btwnClusDEG(seu.obj = seu.obj.k9, groupBy = "celltype.l2", idents.1 = "pDC", idents.2 = "cDC2", bioRep = "orig.ident",padj_cutoff = 0.05, lfcCut = 0.58, 
                        minCells = 5, outDir = paste0("./output/", outName, "/"), title = "can_pDC_vs_cDC2", idents.1_NAME = "can_pDC", idents.2_NAME = "can_cDC2", returnVolc = F, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = T, dwnSam = F
                    )

#plot Fig 9b-d using the output of the above contrasts
p <- crossSpeciesDEG(pwdTOspecies1 = "/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/human_naive/hu_fib_vs_hu_endo_all_genes.csv", 
                     pwdTOspecies2 = "/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/human_dog_naive/can_fib_vs_can_endo_all_genes.csv", 
                     species1 = "Human", species2 = "Canine", cONtrast = c("fibroblast","endothelial"),
                     nlab = 10, nlab_axis = 2, colUp = "red", colDwn = "blue", seed = 12,
                     overlapGenes = overlapGenes, saveGeneList = T
                    )
ggsave(paste("./output/", outName, "/", outName, "_pSigned_cfibVend.png", sep = ""), width = 7,height=5)

p <- crossSpeciesDEG(pwdTOspecies1 = "/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/human_naive/hu_matOC_vs_hu_cd14mono_all_genes.csv", 
                     pwdTOspecies2 = "/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/human_dog_naive/can_matOC_vs_can_tim_all_genes.csv", 
                     species1 = "Human", species2 = "Canine", cONtrast = c("mature OC", "TIM"),
                     nlab = 10, nlab_axis = 2, colUp = "red", colDwn = "blue", seed = 12,
                     overlapGenes = overlapGenes, saveGeneList = T 
                    )
ggsave(paste("./output/", outName, "/", outName, "_pSigned_cOCvTIM.png", sep = ""), width = 7,height=5)

p <- crossSpeciesDEG(pwdTOspecies1 = "/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/human_naive/hu_pDC_vs_hu_cDC2_all_genes.csv", 
                     pwdTOspecies2 = "/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/human_dog_naive/can_pDC_vs_can_cDC2_all_genes.csv", 
                     species1 = "Human", species2 = "Canine", cONtrast = c("pDC", "cDC2"), seed = 666,
                     nlab = 10, nlab_axis = 2, colUp = "red", colDwn = "blue",
                     overlapGenes = overlapGenes,saveGeneList = T,
                     hjustvar = c(1,0.85,0,0.5,1,0,1,0),
                     vjustvar = c(0.75,1,0.5,0,0,1,1,0)
                    )
ggsave(paste("./output/", outName, "/", outName, "_pSigned_cpDCvcDC2.png", sep = ""), width = 7,height=5)


### Supplemental fig - jaccard similarity

#load in cell type gene signatures
dog.df <- read.csv("/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/viln/canine_gene_list.csv")
human.df <- read.csv("/pl/active/dow_lab/dylan/k9_OS_tumor_scRNA/analysis/output/viln/human_gene_list.csv")

#check the data quality
dog.df %>% group_by(cluster) %>% summarize(nn = n()) %>% summarize(min = min(nn))
human.df %>% group_by(cluster) %>% summarize(nn = n()) %>% summarize(min = min(nn))

#calc the jaccard similarity index
res <- lapply(unique(dog.df$cluster), function(x){
    dog.list <- dog.df[dog.df$cluster == x, ] %>% .$gene
    
    res_pre <- lapply(unique(human.df$cluster), function(y){
        human.list <- human.df[human.df$cluster == y, ] %>% .$gene
        
        interSect <- length(intersect(human.list, dog.list)) 
        uni <- length(human.list) + length(dog.list) - interSect 
        JaccardIndex <- interSect/uni
        names(JaccardIndex) <- x
        
        return(JaccardIndex)
    })
    
    res_pre <- do.call(rbind, res_pre)
    rownames(res_pre) <- unique(human.df$cluster)
    return(res_pre)
    
})
    
res1 <- do.call(cbind, res)

#order the celltypes
rowTarg <- c("Endothelial","Fibroblast" ,"Osteoblast","Cycling osteoblast",
       "CD14_monocyte","NR4A3_Macrophage","TXNIP_Macrophage","FABP5_Macrophage",
       "IFN-TAM","Pre-OC","Mature-OC","pDC","cDC2","Mast",
       "CD8 T cell","CD4 T cell","IFN-sig T cell","NK cell",
       "Plasma cell","B cell")
colTarg <- c("Endothelial cell","Fibroblast","Osteoblast_1","Osteoblast_2","Osteoblast_3","Hypoxic_osteoblast",
         "IFN-osteoblast","Osteoblast_cycling","CD4+_TIM","CD4-_TIM","ANGIO_TAM","TAM_INT","TAM_ACT","LA-TAM_SPP2_hi",
         "LA-TAM_C1QC_hi","IFN-TAM","Cycling_OC","CD320_OC","Mature_OC","pDC","cDC1","cDC2","mregDC",
         "Mast cell","Neutrophil","CD8 T cell","CD4 T cell","T_IFN","T_cycling","NK","Plasma cell","B cell")
res1 <- res1[match(rowTarg, rownames(res1)),]        
res1 <- res1[ ,match(colTarg, colnames(res1))]   
       
#plot the data
png(file = paste0("./output/", outName, "/", outName, "_jaccard.png"), width=4000, height=4000, res=400)
par(mfcol=c(1,1))         
ht <- Heatmap(t(res1), #name = "mat", #col = col_fun,
              name = "Jaccard similarity index",
              cluster_rows = F,
              row_title = "Canine cell types",
              row_title_gp = gpar(fontsize = 24),
              col=viridis(option = "magma",100),
              cluster_columns = F,
              column_title = "Human cell types",
              column_title_gp = gpar(fontsize = 24),
              column_title_side = "bottom",
              column_names_rot = 45,
              heatmap_legend_param = list(legend_direction = "horizontal", title_position = "topleft",  title_gp = gpar(fontsize = 16), 
                                          labels_gp = gpar(fontsize = 8), legend_width = unit(6, "cm")),
             )
draw(ht, padding = unit(c(2, 12, 2, 5), "mm"),heatmap_legend_side = "top")
dev.off()

####################################
### End integrated data analysis ###
####################################
