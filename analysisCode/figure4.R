#!/usr/bin/Rscript

#load custom functions & packages
source("./customFunctions.R")
               
################################################
#######   subset and recluster on dcs   ########
################################################

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

#subset on dcs
seu.obj <- subset(seu.obj,
                  subset = 
                  majorID ==  "dc")

#complete independent reclustering
seu.obj <- indReClus(seu.obj = seu.obj, outDir = "../output/s2/", subName = "dc_QCfiltered_2000", preSub = T, nfeatures = 2000,
                      vars.to.regress = "percent.mt"
                       )

# seu.obj <- readRDS(file = "../output/s2/dc_QCfiltered_2000_S2.rds")
clusTree(seu.obj = seu.obj, dout = "../output/clustree/", outName = "dc_QCfiltered_2000", test_dims = c(45,40,35), algorithm = 3, prefix = "integrated_snn_res.")

#complete dim reduction and vis
seu.obj <- dataVisUMAP(seu.obj = seu.obj, outDir = "../output/s3/", outName = "dc_QCfiltered_2000", final.dims = 35, final.res = 0.3, stashID = "clusterID_sub", 
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.3, n.neighbors = 50, assay = "integrated", saveRDS = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA"))

######################################
#######   Begin DC analysis   ########
######################################

#load in processed object -- see line 70 if loading in data from Zenodo
seu.obj <- readRDS(file = "../output/s3/dc_QCfiltered_2000_res0.3_dims35_dist0.3_neigh50_S3.rds")

#stash new idents
Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "cDC2", "1" = "mregDC", 
                                   "2" = "cDC1", "3" = "pDC",
                                   "4" = "preDC")
                       )
seu.obj$majorID_sub <- Idents(seu.obj)

Idents(seu.obj) <- "clusterID_sub"
seu.obj <- RenameIdents(seu.obj, c("0" = "cDC2 (c0)", "1" = "mregDC (c1)", 
                                   "2" = "cDC1 (c2)", "3" = "pDC (c3)",
                                   "4" = "preDC (c4)")
                       )
seu.obj$majorID_sub2 <- Idents(seu.obj)

#set output path
outName <- "dc_naive6"

#export the annotated dataset for Zenodo - no need to run
# saveRDS(seu.obj, "../output/s3/dc_subset_annotated.rds")

### If loading from Zenodo repository, can start here
# seu.obj <- readRDS("../output/s3/dc_subset_annotated.rds")

#used later
ct.l3 <- seu.obj$majorID_sub

### Fig extra: create violin plots for cell ID
vilnPlots(seu.obj = seu.obj, groupBy = "majorID_sub2", numOfFeats = 24, outName = "dc_QCfiltered_2000", outDir = "../output/viln/myeloid/", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^ENSCAF", "^RPS"), assay = "RNA", 
                      min.pct = 0.25, only.pos = T
                     )


### Fig 4a: Create raw UMAP
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "clusterID_sub",
              pt.size = 0.5,
              label = F,
              cols = c("#FF755F", "#FFD6C6", "#AC0535", "#EB2C31", "tomato"),
              label.box = F,
              shuffle = TRUE
)
pi <- formatUMAP(plot = pi) + NoLegend() + theme(axis.title = element_blank(), panel.border = element_blank())
ggsave(paste("../output/", outName, "/", outName, "_rawUMAP.png", sep = ""), width = 7, height = 7)


### Create labels to be cropped onto UMAP
df <- as.data.frame(c("0" = "cDC2", "1" = "mregDC", 
                                   "2" = "cDC1", "3" = "pDC",
                                   "4" = "preDC"))
colnames(df) <- "cellType"
df$cluster <- rownames(df)
df$colorz <- c("#FF755F", "#FFD6C6", "#AC0535", "#EB2C31", "tomato")
df$title <- "DCs"
df$labCol <- c("black","black","white","white","black")
leg <- cusLeg(legend = df, clusLabel = "cluster",colorz = "colorz", legLabel = "cellType",labCol = "labCol",colz = 1, rowz = NULL, groupLabel = "title", groupBy = "title",sortBy = "cluster",
             compress_x = 0.7)

ggsave(paste("../output/", outName, "/", outName, "_leg_forUMAP_labels.png", sep = ""), width = 1.25, height = 7)


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


### Fig extra: Create UMAP by sample
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
ggsave(paste("../output/", outName, "/", outName, "_UMAPbySample.png", sep = ""), width = 7, height = 10.5)


### Fig extra: ref map with canine pbmc dataset from Ammons 2023 (https://doi.org/10.3389/fimmu.2023.1162700)
#set the path to the location in which the reference file is saved
#downloaded using `wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225599/suppl/GSE225599_final_dataSet_HvO.rds.gz`
reference <- readRDS(file = "../../../k9_PBMC_scRNA/analysis/output/s3/final_dataSet_HvO.rds") 

#prepare the reference
reference[['integrated']] <- as(object = reference[['integrated']] , Class = "SCTAssay")
DefaultAssay(reference) <- "integrated"

#find conserved anchors with query and reference
anchors <- FindTransferAnchors(
    reference = reference,
    query = seu.obj,
    normalization.method = "SCT",
    reference.reduction = "pca",
    dims= 1:50
)

#select meta.data slot to use for label transfer -- change refdata value to use alternate labels (i.e., refdata = reference$celltype.l1)
predictions <- TransferData(anchorset = anchors, refdata = reference$celltype.l3,
    dims = 1:50)
seu.obj <- AddMetaData(seu.obj, metadata = predictions)

#generate and save the image
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "predicted.id",
             # split.by = "predicted.id",
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = F#,
             # ncol = 3
)
ggsave(paste("../output/", outName, "/", outName, "_referenceMap.png", sep = ""), width = 10, height = 7)


### Fig 4b: Create violin plots for key feats
features = c("IL3RA", "PGLYRP2", "DDR2", 
             "IGF1","RARRES2", "IGHM","IGKC",
             "CADM1","DNASE1L3", "CLEC1B", 
             "IL21R","IL4I1", "IDO1",
             "CD1C","PID1","CD300H"
            )

pi <- VlnPlot(
    object = seu.obj,
    pt.size = 0,
    same.y.lims = F,
    group.by = "majorID_sub2",
    combine = T,
    cols = c("#FF755F", "#FFD6C6", "#AC0535", "#EB2C31", "tomato"),
    stack = T,
    fill.by = "ident",
    flip = T,
    features = rev(features)
        ) + NoLegend() + theme(axis.ticks = element_blank(),
                               axis.text.y = element_blank(),
                               axis.title.x = element_blank())

#plot <- prettyViln(plot = pi, colorData = NULL, nrow = 2, ncol = 4)
ggsave(paste("../output/", outName, "/", outName, "_selectViln.png", sep = ""), width = 5, height =6)


### Fig 4c: complete heirchical clustering
#set id of interest to type for ease
seu.obj$type <- seu.obj$majorID_sub

#extract metadata and data
metadata <- seu.obj@meta.data
expression <- as.data.frame(t(seu.obj@assays$RNA@data)) #use log noralized count
expression$anno_merge <- seu.obj@meta.data[rownames(expression),]$type

#get cell type expression averages
clusAvg_expression <- expression %>% group_by(anno_merge) %>% summarise(across(where(is.numeric), mean)) %>% as.data.frame()
rownames(clusAvg_expression) <- clusAvg_expression$anno_merge
clusAvg_expression$anno_merge <- NULL                                                                       

#complete clustering
M <- (1- cor(t(clusAvg_expression),method="pearson"))/2
hc <- hclust(as.dist(M),method="complete")

#plot the results
ggtree(as.phylo(hc)) + geom_tiplab(offset = 0.003) + xlim(NA,.05) + geom_tippoint(shape = 21,size = 8,alpha = 1, colour="black", fill = c("#FF755F", "#FFD6C6", "#AC0535", "#EB2C31", "tomato")) + geom_tiplab(aes(label = c("0","1","2","3","4","5","6","7","8")),colour=c("black","black","white","white","black"),offset = -0.001)
ggsave(paste("../output/", outName, "/", outName, "_dc_hc.png", sep = ""), width = 3, height = 7)


### Fig supp 4a: key feature plots
features <- c("MS4A1","JCHAIN", "FLT3", "DLA-DRA")
p <- prettyFeats(seu.obj = seu.obj, nrow = 2, ncol =  2, features = features, 
                 color = "black", order = T, pt.size = 0.5, title.size = 20, noLegend = T)
ggsave(paste("../output/", outName, "/", outName, "_key_feats.png", sep = ""), width = 6, height = 6)


#load in human genessets to define DCs
modulez <- list("TLRs" = c("MYD88","MAVS","TLR1","TLR2","TLR3","TLR4","TLR5","TLR6",
                          "TLR7","TLR8","TLR9"),
                "Maturation" = c("CD40", "CD80","CD86","RELB","CD83"),
                "Regulatory" = c("CD274", "PDCD1LG2", "CD200", "FAS","ALDH1A2", "SOCS1","SOCS2"),
                "Migration" = c("CCR7","MYO1G","CXCL16","ADAM8","ICAM1","FSCN1","MARCKS","MARCKSL1")
                )

seu.obj <- AddModuleScore(seu.obj,
                          features = modulez,
                         name = "_score")
names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)


### Fig extra: module scoring of  the terms
features <- names(modulez)
ecScores <- majorDot(seu.obj = seu.obj, groupBy = "majorID_sub",
                     features = rev(features)
                    ) + coord_flip() + theme(plot.margin = margin(7, 7, 7, 7, "pt"),
                                                                   axis.title = element_blank(),
                                                                   legend.position = "right",
                                                                   legend.direction = "vertical",
                                                                   axis.text.x = element_text(angle=0, hjust = 0.5)
                                            ) + scale_y_discrete(position = "right") + scale_colour_continuous(name="Enrichment score", type = "viridis") + NoLegend()

ggsave(paste("../output/", outName, "/", outName, "_dotPlot_ecScores.png", sep = ""), width = 6,height = 2)


### Fig 4d/e: use human genes to define DCs
#need to manually mod output to
labelz <- as.data.frame(names(modulez))
colnames(labelz) <- "labz"
labelz$modLen <- unname(unlist(lapply(modulez, length)))

cntr <- 0
plots <- lapply(modulez, function(x){
    cntr <<- cntr+1
    labz.df <- labelz[cntr,]

    majorDot(seu.obj = seu.obj, groupBy = "majorID_sub",
                     features = unname(unlist(x))
                    ) + theme(axis.ticks = element_blank(),
                              
                              plot.title = element_text(size = 12),
                                                          legend.position = "right",
                                                                   legend.direction = "vertical",
                                                          axis.title = element_blank(),
                                                          plot.margin = margin(3, 0, 3, 0, "pt")
                                                         ) + scale_colour_distiller(palette = "RdYlBu", name='Average\nexpression', limits = c(-2.5,2.5)) + ggtitle(labz.df$labz)+ NoLegend() + {if(cntr > 2){theme(axis.text.y = element_blank())}}
    
})

modulez2 <- list("TLRs" = c("MYD88","MAVS","TLR1","TLR2","TLR3","TLR4","TLR5","TLR6",
                          "TLR7","TLR8","TLR9"),
                "Maturation" = c("CD40", "CD80","CD86","RELB","CD83"),
                "Regulatory" = c("CD274", "PDCD1LG2", "CD200", "FAS","ALDH1A2", "SOCS1","SOCS2"),
                "Migration" = c("CCR7","MYO1G","CXCL16","ADAM8","ICAM1","FSCN1","MARCKS","MARCKSL1")
                )

    patch <- area()
    ncol <- length(modulez)
    nrow <- 1
    counter=0
    for (i in 1:nrow) {
        for (x in 1:ncol) {
                patch <- append(patch, area(t = i, l = x, b = i, r = x))
        }
    }

p <- Reduce( `+`, plots ) +  plot_layout(guides = "collect", design = patch, 
                                                                             width = unname(unlist(lapply(modulez, length)))/sum(unname(unlist(lapply(modulez, length)))))

ggsave(paste("../output/", outName, "/", outName, "_dotPlot_ecScores_2.png", sep = ""), width = 10,height=3)


### Fig 4f: genes that define
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "clusterID_sub", idents.1 = "1", idents.2 = "0", bioRep = "name",padj_cutoff = 0.05, lfcCut = 0.58, 
                      minCells = 25, outDir = paste0("../output/", outName, "/"), title = "mregDC_vs_cDC2", idents.1_NAME = "mregDC", idents.2_NAME = "cDC2",
                      returnVolc = T, doLinDEG = F, paired = T, addLabs = NULL,lowFilter = T, dwnSam = F
                     )

p  <- prettyVolc(plot = p_volc[[1]], rightLab = "Up in c1 (mregDC)", leftLab = "Up in c0 (cDC2)") + labs(x = "log2(FC) mregDC vs cDC2")

ggsave(paste("../output/", outName, "/", outName, "_c1vc0_volcPlot.png", sep = ""), width = 7, height = 7)


p <- plotGSEA(pwdTOgeneList = "../output/dc_naive6/mregDC_vs_cDC2_all_genes.csv", category = "C2", upOnly = T, termsTOplot = 15)
ggsave(paste("../output/", outName, "/", outName, "_enriched_terms.png", sep = ""), width = 9, height =7)

p <- plotGSEA(pwdTOgeneList = "../output/dc_naive6/mregDC_vs_cDC2_all_genes.csv", category = "C2",subcategory = "CP:REACTOME")
ggsave(paste("../output/", outName, "/", outName, "_enriched_terms_2.png", sep = ""), width = 9, height =7)


### Fig extra: genes that define cDC1
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "clusterID_sub", idents.1 = "2", idents.2 = "0", bioRep = "name",padj_cutoff = 0.05, lfcCut = 0.58, 
                        minCells = 5, outDir = paste0("../output/", outName, "/"), title = "cDC1_vs_cDC2", idents.1_NAME = "c2", idents.2_NAME = "c0", returnVolc = T, doLinDEG = F, paired = T, addLabs = NULL,lowFilter = T, dwnSam = F
                    )

p  <- prettyVolc(plot = p_volc[[1]], rightLab = "Up in c1 (cDC1)", leftLab = "Up in c0 (cDC2)") + labs(x = "log2(FC) cDC1 vs cDC2")

ggsave(paste("../output/", outName, "/", outName, "_c2vc0_volcPlot.png", sep = ""), width = 7, height = 7)


### Fig extra: genes that define pDC1
p_volc <- btwnClusDEG(seu.obj = seu.obj, groupBy = "clusterID_sub", idents.1 = "3", idents.2 = "0", bioRep = "name",padj_cutoff = 0.05, lfcCut = 0.58, 
                        minCells = 5, outDir = paste0("../output/", outName, "/"), title = "pDC_vs_cDC2", idents.1_NAME = "c3", idents.2_NAME = "c0", returnVolc = T, doLinDEG = F, paired = T, addLabs = NULL,lowFilter = T, dwnSam = F
                    )

p  <- prettyVolc(plot = p_volc[[1]], rightLab = "Up in c3 (pDC)", leftLab = "Up in c0 (cDC2)") + labs(x = "log2(FC) pDC vs cDC2")

ggsave(paste("../output/", outName, "/", outName, "_c3vc0_volcPlot.png", sep = ""), width = 7, height = 7)


################################################
#######   END dendritic cell analysis   ########
################################################
