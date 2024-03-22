#!/usr/bin/Rscript

library(Seurat)
library(tidyverse)
library(clustree)
library(stringr)
library(DoubletFinder)
library(patchwork)
library(scales)
library(cowplot)
library(ggrepel)
library(colorspace)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(SeuratDisk)
library(SingleR)
library(viridis)
library(reshape)
library(lemon)
library(ggsankey)
library(msigdbr)
library(clusterProfiler)
library(slingshot)
library(ggpubr)
library(scRNAseq)
library(scuttle)
library(ape)
library(ggtree)
library(ComplexHeatmap)

############ gg_color_hue ############

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
############ hsv2rgb ############
#taken from: https://datascienceconfidential.github.io/r/graphics/2017/11/23/soothing-pastel-colours-in-r.html

hsv2rgb <- function(x){  
    # convert an hsv colour to rgb  
    # input:  a 3 x 1 matrix (same as output of rgb2hsv() function)  
    # output: vector of length 3 with values in [0,1]    
        
    # recover h, s, v values  
    h <- x[1,1]  
    s <- x[2,1]  
    v <- x[3,1]    
        
    # follow the algorithm from Wikipedia  
    C <- s*v   
        
    # in R, h takes values in [0,1] rather than [0, 360], so dividing by  
    # 60 degrees is the same as multiplying by six  
    hdash <- h*6  
    X <- C * (1 - abs(hdash %% 2 -1))
        
    if (0 <= hdash & hdash <=1) RGB1 <- c(C, X, 0)  
    if (1 <= hdash & hdash <=2) RGB1 <- c(X, C, 0)  
    if (2 <= hdash & hdash <=3) RGB1 <- c(0, C, X)  
    if (3 <= hdash & hdash <=4) RGB1 <- c(0, X, C)  
    if (4 <= hdash & hdash <=5) RGB1 <- c(X, 0, C)  
    if (5 <= hdash & hdash <=6) RGB1 <- c(C, 0, X) 
        
    # the output is a vector of length 3. This is the most convenient  
    # format for using as the col argument in an R plotting function 
    RGB1 + (v-C)
}

############ pastellize ############
#taken from: https://datascienceconfidential.github.io/r/graphics/2017/11/23/soothing-pastel-colours-in-r.html

pastellize <- function(x, p){
    # x is a colour
    # p is a number in [0,1]
    # p = 1 will give no pastellization
  
    # convert hex or letter names to rgb
    if (is.character(x)) x <- col2rgb(x)/255
  
    # convert vector to rgb
    if (is.numeric(x)) x <- matrix(x, nr=3)
  
    col <- rgb2hsv(x, maxColorValue=1)
    col[2,1] <- col[2,1]*p
    col <- hsv2rgb(col)
  
    # return in convenient format for plots
    rgb(col[1], col[2], col[3])
}

############ load10x ############

load10x <- function(din = NULL, dout = NULL, outName = NULL, testQC = FALSE,
                     nFeature_RNA_high = 4500, nFeature_RNA_low = 200, nCount_RNA_high = 20000, nCount_RNA_low = 100, percent.mt_high = 10, mt_pattern = "^MT-",
                     removeDubs = TRUE, removeRBC_pal = TRUE, pal_feats = NULL,
                     featPlots = c("PTPRC", "CD3E", "CD8A", "GZMA", "IL7R", "ANPEP", "FLT3", "DLA-DRA", "CD4", "MS4A1", "PPBP","HBM"), isolatePalRBC = FALSE, readCnts = F){
    
    if (testQC == FALSE){
        print(paste0("The QC parameters are: nFeature_RNA < ", nFeature_RNA_high, " & nFeature_RNA > ", nFeature_RNA_low, " & percent.mt < ", percent.mt_high, " & nCount_RNA < ", nCount_RNA_high," & nCount_RNA > ", nCount_RNA_low, sep = ""))
    }

    fpath <-  paste("./", din,"/", sep = "") 

    files <- list.files(path = fpath, pattern=NULL, all.files=FALSE,
                        full.names=F)

    df.list <- list()
    #foreach (infile = files, .combine = 'c') %dopar% {
    for (infile in files) {
        #set up df for export
        if (removeRBC_pal == TRUE){
            df <- data.frame(matrix(ncol = 4, nrow = 0))
            colnames(df) <- c("Initial","Filtered","Platelet_rbc_rm","Singlets")
        } else {
            df <- data.frame(matrix(ncol = 3, nrow = 0))
            colnames(df) <- c("Initial","Filtered","Singlets")
        }
    
        #set import path
        pwd <- paste("./", din,"/", infile, sep = "") 
        
        if(readCnts){
            cntMat.pwd <- list.files(path = paste0(pwd,"/"), pattern=NULL, all.files=F,
                        full.names=T)
            
            cnts.df <- read.table(cntMat.pwd)
            
            seu_obj <- CreateSeuratObject(counts = cnts.df,
                                          project = infile,
                                          min.cells = 3,
                                          min.features = 200)
        }else{
            
            #read 10X data
            indata <- Read10X(pwd)
            
            #create seurat object
            
            seu_obj <- CreateSeuratObject(counts = indata,
                                      project = infile,
                                      min.cells = 3,
                                      min.features = 200)
        }

        #Add mitochondrial QC data to seurat metadata
        seu_obj[["percent.mt"]] <- PercentageFeatureSet(seu_obj, pattern = mt_pattern)
        seu_obj[["percent.hbm"]] <- PercentageFeatureSet(seu_obj, pattern = "HBM")
        seu_obj[["percent.ppbp"]] <- PercentageFeatureSet(seu_obj, pattern = "PPBP")
        if(!is.null(pal_feats)){            
            data <- lapply(pal_feats, function(x){PercentageFeatureSet(seu_obj, pattern = x)})
            data <- as.data.frame(bind_cols(data))
            seu_obj[["percent.pal"]] <- rowSums(data)
        }
    
        #visualize QC metrics as a violin plot
        outfile <- paste("./",dout,"/", infile,"_QC_S1.png", sep = "") 
        p <- VlnPlot(seu_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
        ggsave(outfile)

        if (testQC == TRUE){
            next
        }
    
        df[1,1] <- dim(seu_obj)[2]

        #set QC cutoffs based on above plots ### MODIFY VALUES based on your data!! ###
        seu_obj <- subset(seu_obj,
                          subset = nFeature_RNA < nFeature_RNA_high & nFeature_RNA > nFeature_RNA_low &
                          percent.mt < percent.mt_high & nCount_RNA < nCount_RNA_high & nCount_RNA > nCount_RNA_low
                         )
  
        #next steps normalize, scale, and run UMAP
        seu_obj <- NormalizeData(seu_obj,
                                 normalization.method = "LogNormalize",
                                 Scale.factor = 10000) ###change method to scran???

        seu_obj <- FindVariableFeatures(seu_obj,
                                        selection.method = "vst", 
                                        nfeatures = 2500) #can change number of feats used
  
        all.genes <- rownames(seu_obj)
        seu_obj <- ScaleData(seu_obj, features = all.genes)
        seu_obj <- RunPCA(seu_obj, features = VariableFeatures(object = seu_obj))
        seu_obj <- FindNeighbors(seu_obj,
                                 dims = 1:10
                                ) #can change dims
        seu_obj <- FindClusters(seu_obj,
                                resolution = 0.1
                               ) #can change resolution
        seu_obj <- RunUMAP(seu_obj, 
                           dims = 1:15
                          ) #can change dims

        #visulize UMAP featPlots
        outfile <- paste("./",dout,"/", infile,"_featPlotDefault_S1.png", sep = "")
        p <- FeaturePlot(seu_obj,features = featPlots)
        ggsave(outfile, width = 12, height = 8)
        
        #visulize UMAP
        outfile <- paste("./",dout,"/", infile,"_uMAP_S1.png", sep = "") 
        p <- DimPlot(seu_obj, reduction = "umap")
        ggsave(outfile)
        
        #stash initial cluster IDs
        seu_obj[["int.clusID"]] <- Idents(object = seu_obj)
        
        df[1,2] <- dim(seu_obj)[2]
        
        if (removeRBC_pal == TRUE){
            #find the rbc cluster
            if(length(AverageExpression(seu_obj, features = "HBM")) != 0){
                rbc.df <- as.data.frame(AverageExpression(seu_obj, features = "HBM"), header = TRUE)
                rbc.clus <- str_split(as.character(colnames(rbc.df)[max.col(rbc.df,ties.method="first")]),"[.]")[[1]][2]
                
                sus.rbc.size <- as.numeric(length(WhichCells(seu_obj, expression = HBM > 0)))
                sus.rbc.clus.size <- as.numeric(length(WhichCells(seu_obj, idents = rbc.clus )))
                
                sus.rbc.comp.pct <- sus.rbc.size/sus.rbc.clus.size
                rbc.clus <- ifelse(sus.rbc.comp.pct > 0.8, rbc.clus, "NULL")
                print(rbc.clus)
            } else{rbc.clus <- "NULL"}
            
            #find the palelet cluster
            if(length(AverageExpression(seu_obj, features = "PPBP")) != 0){
                pal.df <- as.data.frame(AverageExpression(seu_obj, features = "PPBP"), header = TRUE)
                pal.clus <- str_split(as.character(colnames(pal.df)[max.col(pal.df,ties.method="first")]),"[.]")[[1]][2]
                
                sus.pal.size <- as.numeric(length(WhichCells(seu_obj, expression = PPBP > 0)))
                sus.pal.clus.size <- as.numeric(length(WhichCells(seu_obj, idents = pal.clus )))
                
                sus.pal.comp.pct <- sus.pal.size/sus.pal.clus.size
                pal.clus <- ifelse(sus.pal.comp.pct > 0.8, pal.clus, "NULL")
                print(pal.clus)
            } else{pal.clus <- "NULL"}
            
            if(rbc.clus != "NULL" | pal.clus != "NULL"){
                if(isolatePalRBC == F){
                    seu_obj <- subset(seu_obj,
                                      subset = int.clusID != pal.clus & int.clusID != rbc.clus)
                }
                else{
                    seu_obj <- subset(seu_obj,
                                      subset = int.clusID == pal.clus | int.clusID == rbc.clus)
                }
            } else{print("No platelets of rbcs detected in sample!")}
            
            df[1,3] <- dim(seu_obj)[2]
        }
        
        #next steps complete doublet identification using DoubletFinder
        sweep.res.list_seu_obj <- paramSweep_v3(seu_obj, PCs = 1:10, sct = FALSE)
        sweep.stats_seu_obj <- summarizeSweep(sweep.res.list_seu_obj, GT = FALSE)
        bcmvn_seu_obj <- find.pK(sweep.stats_seu_obj)
        pk_val <- as.numeric(as.character(bcmvn_seu_obj$pK[bcmvn_seu_obj$BCmetric == max(bcmvn_seu_obj$BCmetric)])) 
        annotations <- seu_obj@meta.data$seurat_clusters
        homotypic.prop <- modelHomotypic(annotations)
        
        #determine expected doublet rate -- this formula assumes 0.5% doublet per 1000 cells (lower than 10x recommendation to account for homoypic doublets and cells removed during QC)
        expceted_dub_rate <- dim(seu_obj)[2]/1000*0.5/100
        
        nExp_poi <- round(expceted_dub_rate*length(seu_obj@meta.data$orig.ident))  ### MODIFY VALUE; should be percentage expected to be doublets based on 10X expected values ###
        nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
        
        pN_value <- 0.25 ### MODIFY VALUE as needed ###
        
        pANN_value <- paste0("pANN_",pN_value,"_",pk_val,'_',nExp_poi, sep = "")
        seu_obj <- doubletFinder_v3(seu_obj, PCs = 1:10, pN = pN_value, pK = pk_val, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
        dfClass <- paste0("DF.classifications_",pN_value,"_",pk_val,'_',nExp_poi, sep = "")
        dfmeta <- paste0("seu_obj@meta.data$",dfClass, sep = "")
        seu_obj[["doublet"]] <- eval(parse(text = dfmeta))
        
        #export UMAP highlighting doublet vs singlet cells
        outfile <- paste("./",dout,"/", infile,"_DF_S1.png", sep = "") 
        DimPlot(seu_obj,pt.size = 1,label=TRUE, label.size = 5,reduction = "umap",group.by = dfClass )
        ggsave(outfile)
        
        if (removeDubs == TRUE){
            #remove putative doublets
            seu_obj <- subset(seu_obj,
                              subset = doublet == "Singlet"
                             )
        }    
        
        #update user
        if (removeRBC_pal == TRUE){
            df[1,4] <- dim(seu_obj)[2]
        } else {df[1,3] <- dim(seu_obj)[2]
               }
        
        rownames(df) <- infile
        df.list[[which(infile == files)]] <- df
        
        #recluster the data with all unwanted cells removed
        seu_obj <- RunPCA(seu_obj, features = VariableFeatures(object = seu_obj))
        seu_obj <- FindNeighbors(seu_obj, 
                                 dims = 1:10
                                ) #can change dims
        seu_obj <- FindClusters(seu_obj, 
                                resolution = 0.1
                               ) #can change resolution
        seu_obj <- RunUMAP(seu_obj, 
                           dims = 1:15
                          ) #can change dims
        
        outfile <- paste("./",dout,"/", infile,"_featPlotDefault_postRBC_pal_rm_S1.png", sep = "")
        p <- FeaturePlot(seu_obj,features = featPlots)
        ggsave(outfile, width = 12, height = 8)
        
        #export final umap
        outfile <- paste("./",dout,"/", infile,"_uMAP_postRBC_pal_rm_S1.png", sep = "") 
        DimPlot(seu_obj, reduction = "umap")
        ggsave(outfile)
        
        #identify cycling cells
        seu_obj <- CellCycleScoring(
            object = seu_obj,
            s.features = cc.genes.updated.2019$s.genes,
            g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = FALSE
        )
        
        #export UMAP with cycling data
        outfile <- paste("./",dout,"/", infile,"_uMAP_cellCylce_postRBC_pal_rm_S1.png", sep = "") 
        DimPlot(seu_obj, reduction = "umap", group.by = "Phase")
        ggsave(outfile)
        
        #store cell cycle state
        seu_obj[["clusters"]] <- Idents(object = seu_obj)
        
        #export processed seurat object as an .RDS file
        outfile <- paste("./",dout,"/", infile,"_S1.rds", sep = "") 
        saveRDS(seu_obj, file = outfile)
    }
    
    cellCounts <- do.call(rbind, df.list)
    outfile <- paste("./",dout,"/", outName, "_cell_counts_S1.csv", sep = "")
    write.csv(cellCounts, file = outfile)
}

############ sctIntegrate ############

sctIntegrate <- function(din = "", dout = "./output/", outName = "", vars.to.regress = c("percent.mt"), nfeatures = 2000,
                         featTOexclude = NULL, pattern = "S1.rds", returnObj = T,
                         varFeatsToUse = NULL
                        ) {
    #get seurat objects to process and integrate
    fpath <- paste("./", din,"/", sep = "") 
    files <- list.files(path = fpath, pattern = pattern, all.files=FALSE,
                        full.names=TRUE)

    create.seu.call <- function(x) {
        readRDS(x)
    }

    #load all seurat objects and store in a large list
    seu.obj <- mapply(create.seu.call, files)

    #use SCT inegration to merge the samples
    seu.obj <- lapply(seu.obj, 
                      SCTransform, 
                      vars.to.regress = vars.to.regress,
                      verbose = TRUE,
                      conserve.memory=TRUE)
    
    if(is.null(varFeatsToUse)){
        SelectedFeatures <- SelectIntegrationFeatures(object.list = seu.obj,
                                                      nfeatures = nfeatures) 
    
        if(!is.null(featTOexclude)){
            SelectedFeatures <- SelectedFeatures[!SelectedFeatures %in% featTOexclude]
            if(nfeatures != length(SelectedFeatures)){
                message <- paste("NOTE: ", featTOexclude, " was/were excluded from the variable features used in integration!",sep = "")
                print(message)
                SelectedFeatures <- SelectIntegrationFeatures(object.list = seu.obj,
                                                              nfeatures = nfeatures+(nfeatures-length(SelectedFeatures))
                                                             )
                SelectedFeatures <- SelectedFeatures[!SelectedFeatures %in% featTOexclude]

            }else{
                message <- paste("NOTE: The features to exclude (", featTOexclude, ") was/were not included in the variable features used in integration, so the option was not used.",sep = "")
                print(message)
            }
        }
    }else{
        SelectedFeatures <- varFeatsToUse
    }
    
    seu.integrated <- PrepSCTIntegration(
        object.list = seu.obj,
        anchor.features = SelectedFeatures,
        verbose = TRUE
    )

    gc()

    seu.integrated.anchors <- FindIntegrationAnchors(
        object.list = seu.integrated,
        normalization.method = "SCT",
        anchor.features = SelectedFeatures,
        verbose = TRUE
    )

    #clean up environment a bit
    rm(seu.integrated)
    gc()

    #integrate data and keep full gene set - still might not be retaining all genes
    seu.integrated.obj <- IntegrateData(
        anchorset = seu.integrated.anchors,
        normalization.method = "SCT",
        verbose = TRUE
    )

    #clean up environment a bit
    rm(seu.integrated.anchors)
    gc()
    
    seu.integrated.obj <- RunPCA(seu.integrated.obj)

    p <- ElbowPlot(seu.integrated.obj, ndims = 50)
    outfile <- paste(dout, outName, "_seu.integrated.obj_S2_elbow.png", sep = "")
    ggsave(outfile)

    DefaultAssay(seu.integrated.obj) <- "integrated"

    outfile <- paste(dout, outName, "_seu.integrated.obj_S2.rds", sep = "")
    saveRDS(seu.integrated.obj, file = outfile)
    
    if(returnObj){
        return(seu.integrated.obj)
    }
}

############ clusTree ############
#uses too much ram if paralellized

clusTree <- function(seu.obj = NULL, dout = "./output/", outName = "", test_dims = c(50,45,40,35), algorithm = 3, resolution = c(0.01, 0.05, 0.1, seq(0.2, 2, 0.1)), 
                     prefix = "integrated_snn_res."
                    ) {
    
    for (dimz in test_dims){
        seu.test <- FindNeighbors(object = seu.obj, dims = 1:dimz)
        seu.test <- FindClusters(object = seu.test, algorithm = algorithm, resolution = resolution)
        
        p <- clustree::clustree(seu.test, prefix = prefix) + ggtitle(paste0("The number of dims used:", dimz, sep = ""))
        
        outfile <- paste(dout,dimz ,"_clustree_", outName,".png", sep = "") 
        png(file = outfile, height = 1100, width = 2000)
        print(p)
        dev.off()
        
    }
}

############ dataVisUMAP ############

dataVisUMAP <- function(file = NULL, seu.obj = NULL, outDir = "", outName = "", final.dims = NULL, final.res = NULL, stashID = "clusterID", returnFeats = T,
                        algorithm = 3, prefix = "integrated_snn_res.", min.dist = 0.6, n.neighbors = 75, assay = "integrated", saveRDS = T, return_obj = T,
                        features = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                       ) {
    
    if(!is.null(file)){
        try(seu.integrated.obj <- readRDS(file), silent = T)
    }else if(!is.null(seu.obj)){
        seu.integrated.obj <- seu.obj
    }

    #perform clustering: Neighbor finding and resolution parameters
    DefaultAssay(seu.integrated.obj ) <- assay
    seu.integrated.obj <- FindNeighbors(object = seu.integrated.obj, dims = 1:final.dims)
    seu.integrated.obj <- FindClusters(object = seu.integrated.obj, algorithm = algorithm, resolution = final.res)

    #choose appropriate clustering resolution
    res <- paste(prefix, final.res,sep = "") 
    seu.integrated.obj <- SetIdent(object = seu.integrated.obj, value = res)
    
    # Run UMAP embedding
    seu.integrated.obj <- RunUMAP(seu.integrated.obj,
                                  dims = 1:final.dims,
                                  min.dist = min.dist,
                                  n.neighbors = n.neighbors
                                 )

    #save cluster number as "clusterID"
    seu.integrated.obj[[stashID]] <- seu.integrated.obj@active.ident

    #build hierarchical clustering based on clustering results
    DefaultAssay(seu.integrated.obj) <- "RNA"
    seu.integrated.obj <- NormalizeData(seu.integrated.obj)
    seu.integrated.obj <- BuildClusterTree(seu.integrated.obj, assay = "RNA", dims = 1:final.dims)

    #export UMAP colored by sample ID and cluster
    outfile <- paste(outDir, outName, "_res", final.res, "_dims", final.dims, "_dist",min.dist, "_neigh",n.neighbors, "_cluster_S3.png", sep = "") 
    p <- DimPlot(seu.integrated.obj, label = TRUE, reduction = "umap", group.by = stashID)
    ggsave(outfile)
    
    outfile <- paste(outDir, outName, "_res", final.res, "_dims", final.dims, "_dist",min.dist, "_neigh",n.neighbors, "_sample_S3.png", sep = "") 
    p <- DimPlot(seu.integrated.obj, reduction = "umap", group.by = "orig.ident")
    ggsave(outfile)

    if(returnFeats == T){
     
        #visulize UMAP featPlots
        outfile <- paste(outDir, outName, "_res", final.res, "_dims", final.dims, "_dist",min.dist, "_neigh",n.neighbors, "_featPlotDefault_S3.png", sep = "")
        p <- FeaturePlot(seu.integrated.obj,features = features)
        ggsave(outfile, width = 12, height = 8)
    }

    #stash identy by healthy vs OS & store as cellSource - add feature
    #Idents(seu.integrated.obj) <- "orig.ident" #need to do smthg about this
    #seu.integrated.obj$cellSource <- ifelse(grepl("ealthy", seu.integrated.obj@meta.data$orig.ident), "Healthy", "Osteosarcoma")

    #save seurat object as rds
    if(saveRDS == T){
        outfile <- paste(outDir, outName, "_res", final.res, "_dims", final.dims, "_dist",min.dist, "_neigh",n.neighbors,"_S3.rds", sep = "")
        saveRDS(seu.integrated.obj, file = outfile)
    }
    
    if(return_obj == T){
        return(seu.integrated.obj)
    }
}



############ indReClus ############
#adopteed from final post on Dec 4, 2021  at https://github.com/satijalab/seurat/issues/1883
indReClus <- function(seu.obj = NULL, group.by = NULL, sub = NULL, outDir = "", subName = "", preSub = F, seu.list = NULL, featTOexclude = NULL, varFeatsToUse = NULL, nfeatures = 2000, k = NULL, saveRDS = T, returnObj = T,ndims = 50,
                      vars.to.regress = "percent.mt", z = 30
                       ) {
    
    #print(paste0("The seurat object loaded: ", seu.obj, sep = ""))
    #print(paste0("The seurat object will be subset on: ", sub, sep = ""))
    print(paste0("The subclustered object will be output as: ", outDir, subName,"_S2.rds", sep = ""))
    
    #read in data
    if(preSub == F){
        Idents(seu.obj) <- group.by
        seu.sub <- subset(seu.obj, idents = sub)
    }else{
        seu.sub <- seu.obj
    }
      
    options(future.globals.maxSize = 100000 * 1024^2)

    if(is.null(seu.list)){
        seu.sub.list <- SplitObject(seu.sub, split.by = "orig.ident")
        z <- table(seu.sub@meta.data$orig.ident)
        print(z)
        k <- ifelse(is.null(k), ifelse(min(z)>100,100, min(z)), k)
        print(k)
    }else{
        seu.sub.list <- seu.list
        k <- ifelse(is.null(k),100,k)
        print(k)
    }
    
    if(is.null(varFeatsToUse)){
        seu.obj <- lapply(seu.sub.list,
                          SCTransform, 
                          vars.to.regress = vars.to.regress,
                          verbose = FALSE,
                          conserve.memory=TRUE)
    
        SelectedFeatures <- SelectIntegrationFeatures(object.list = seu.obj,
                                                      nfeatures = nfeatures) 
    
        if(!is.null(featTOexclude)){
            SelectedFeatures <- SelectedFeatures[!SelectedFeatures %in% featTOexclude]
            if(nfeatures != length(SelectedFeatures)){
                message <- paste("NOTE: ", featTOexclude, " was/were excluded from the variable features used in integration!",sep = "")
                print(message)
                SelectedFeatures <- SelectIntegrationFeatures(object.list = seu.obj,
                                                              nfeatures = nfeatures+(nfeatures-length(SelectedFeatures))
                                                             )
                SelectedFeatures <- SelectedFeatures[!SelectedFeatures %in% featTOexclude]

            }else{
                message <- paste("NOTE: The features to exclude (", featTOexclude, ") was/were not included in the variable features used in integration, so the option was not used.",sep = "")
                print(message)
            }
        }
        
        seu.integrated <- PrepSCTIntegration(object.list = seu.obj,
                                            anchor.features = SelectedFeatures,
                                            verbose = FALSE
                                            )
    }else{
        
        seu.obj <- lapply(seu.sub.list,
                          SCTransform, 
                          vars.to.regress = vars.to.regress,
                          verbose = FALSE,
                          variable.features.n = 3000,
                          return.only.var.genes = FALSE)
        
        seu.obj <- lapply(seu.obj, function(obj){
            varKeep <- varFeatsToUse[varFeatsToUse %in% rownames(obj@assays$SCT@scale.data)]
            obj@assays$SCT@scale.data = obj@assays$SCT@scale.data[varKeep, ]
            obj@assays$SCT@var.features = varFeatsToUse
            return(obj)
        })
        
         varFeats <- lapply(seu.obj, function(obj){
             varKeep <- varFeatsToUse[varFeatsToUse %in% rownames(obj@assays$SCT@scale.data)]
             return(varKeep)
         })
        
        SelectedFeatures <- Reduce(intersect, varFeats)
        
        #may want to skip this step
        seu.integrated <- PrepSCTIntegration(object.list = seu.obj,
                                             anchor.features = SelectedFeatures,
                                             verbose = FALSE
                                             )
        
    }



    gc()

    seu.integrated.anchors <- FindIntegrationAnchors(
        object.list = seu.integrated,
        normalization.method = "SCT",
        anchor.features = SelectedFeatures,
        dims = 1:ifelse(min(z)>30, 30, min(z)-1),
        k.filter = ifelse(min(z)>200, 200, min(z)-1),
        k.score = ifelse(min(z)>30, 30, min(z)-1)
    )

    #clean up environment a bit
    rm(seu.sub)
    rm(seu.sub.list)
    rm(seu.obj)
    rm(seu.integrated)
    gc()

    #integrate data and keep full gene set - still might not be retaining all genes
    seu.integrated.obj <- IntegrateData(
        anchorset = seu.integrated.anchors,
        normalization.method = "SCT",
        k.weight = k,
        verbose = FALSE
    )

    #clean up environment a bit
    rm(seu.integrated.anchors)
    gc()

    seu.integrated.obj <- RunPCA(seu.integrated.obj)


    outfile <- paste(outDir, subName,"_S2_elbow.png", sep = "")
    p <- ElbowPlot(seu.integrated.obj, ndims = ndims)
    ggsave(outfile)

    DefaultAssay(seu.integrated.obj) <- "integrated"

    if(saveRDS){
        outfile <- paste(outDir, subName,"_S2.rds", sep = "")
        saveRDS(seu.integrated.obj, file = outfile)
    }
    
    if(returnObj){
        return(seu.integrated.obj)
    }
}


############ prettyFeats ############

prettyFeats <- function(seu.obj = NULL, nrow = 3, ncol = NULL, features = "", color = "black", order = FALSE, titles = NULL, noLegend = F, bottomLeg = F, min.cutoff = NA, pt.size = NULL, title.size = 18, legJust = "bottom",showAxis = F,smallAxis=F, legInLine = F, returnPlots = F
                       ) {
    
    DefaultAssay(seu.obj) <- "RNA"
    features <- features[features %in% c(unlist(seu.obj@assays$RNA@counts@Dimnames[1]),unlist(colnames(seu.obj@meta.data)))]
    
    if(is.null(ncol)){
        ncol = ceiling(sqrt(length(features)))
    }
    
    if(is.null(titles)){
        titles <- features #- add if statement
    }
    
    
    #strip the plots of axis and modify titles and legend -- store as large list
    plots <- Map(function(x,y,z) FeaturePlot(seu.obj,features = x, pt.size = pt.size, order = order, min.cutoff = min.cutoff) + labs(x = "UMAP1", y = "UMAP2") +
                 theme(axis.text= element_blank(), 
                       axis.ticks = element_blank(),
                       axis.title = element_blank(), 
                       axis.line = element_blank(),
                       title = element_text(size= title.size, colour = y),
                       legend.position = "none"
                      ) + 
                 scale_color_gradient(breaks = pretty_breaks(n = 3), limits = c(NA, NA), low = "lightgrey", high = "darkblue") + 
                 ggtitle(z), x = features, y = color, z = titles) 


    asses <- ggplot() + labs(x = "UMAP1", y = "UMAP2") + 
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
    if(!noLegend){
        leg <- FeaturePlot(seu.obj,features = features[1], pt.size = 0.1) + 
        theme(legend.position = 'bottom',
              legend.direction = 'vertical',
              legend.justification = "center",

             ) + 
        scale_color_gradient(breaks = pretty_breaks(n = 1), labels = c("low", "high"), limits = c(0,1), low = "lightgrey", high = "darkblue") + 
        guides(color = guide_colourbar(barwidth = 1)) 
        
        if(bottomLeg){
            leg <- leg + theme(legend.direction = 'horizontal') + guides(color = guide_colourbar(barwidth = 8))
        }

        legg <- get_legend(leg)
    }

    #nrow <- ceiling(length(plots)/ncol) - add if statement
    patch <- area()

    counter=0
    for (i in 1:nrow) {
        for (x in 1:ncol) {
            counter = counter+1
            if (counter <= length(plots)) {
                patch <- append(patch, area(t = i, l = x, b = i, r = x))
            }
        }
    }
    
    if(!noLegend){
        if(!bottomLeg & !legInLine){
            legPos <- ifelse(legJust == "bottom",ceiling(length(features)/ncol),1)
            patch <- append(patch, area(t = legPos, l = ncol+1, b = legPos, r = ncol+1))
        }else if(legInLine){
            patch <- append(patch, area(t = nrow, l = ncol, b = nrow, r = ncol))
        }else{
            patch <- append(patch, area(t = ceiling(length(features)/ncol)+1, l = ncol, b = ceiling(length(features)/ncol)+1, r = ncol))
        }
    }else{
        if(smallAxis){
            patch <- append(patch, area(t = nrow, l = 1, b = nrow, r = 1))
        }else{
            patch <- append(patch, area(t = 1, l = 1, b = nrow, r = ncol))
        }
    }

    
    if(returnPlots){
        return(plots)
    }else{
        p <- Reduce( `+`, plots ) +  {if(noLegend & showAxis){asses}} +
        {if(!noLegend){legg}} + plot_layout(guides = "collect") +
        {if(!noLegend & bottomLeg){plot_layout(design = patch, heights = c(rep.int(1, nrow),0.2))}else if(!noLegend & !bottomLeg){plot_layout(design = patch, widths = c(rep.int(1, ncol),0.2))}else{plot_layout(design = patch, widths = rep.int(1, ncol))}}
        
        return(p)
    }
}

############ linDEG ############

linDEG <- function(seu.obj = NULL, threshold = 1, thresLine = F, groupBy = "clusterID", 
                   comparision = "cellSource", outDir = "./output/", outName = "", cluster = NULL, 
                   labCutoff = 20, noTitle = F, contrast = NULL, colUp = "red", colDwn = "blue", 
                   subtitle = T, returnUpList = F, returnDwnList = F, forceReturn = F, useLineThreshold = F, 
                   pValCutoff = 0.01, logfc.threshold = 0.58,
                   saveGeneList = F, addLabs = "", labsHide = c("^MT-","^RPL","^ENS","^RPS"), returnPlots = F
                  ) {
    
    #set active.ident
    Idents(seu.obj) <- groupBy    
    
    if(is.null(contrast)){
        message(paste0("It is now reccomened that a contrast be provided. This should be two levels present in the comparision metadata slot. Please provide this information as c(",  "'comp1'", ",",  "'comp2'",") and the analysis will run comp1 VS comp2."))
        contrast <- levels(seu.obj@meta.data[[comparision]])[c(1:2)]
        message(paste0("Contrasting - (1) ", contrast[1]," - VS - (2) ", contrast[2]))
    }
    
    #loop through active.ident to make plots for each group
    lapply(if(is.null(cluster)){levels(seu.obj)}else{cluster}, function(x) {
        
        #subset on population of cells to analyze -- if all cells are included this will effectively not do anything
        seu.sub <- subset(seu.obj, idents = x)
        seu.sub@meta.data[[groupBy]] <- droplevels(seu.sub@meta.data[[groupBy]])
               
        #use FindMarkers buried in the the custom vilnSplitComp function to complete stats to ID DEGs
        geneList <- vilnSplitComp(seu.obj = seu.sub, groupBy = groupBy, refVal = comparision, 
                                  outDir = outDir, outName = outName, 
                                  saveOut = F, saveGeneList = saveGeneList, returnGeneList = T, 
                                  contrast = contrast, pValCutoff = pValCutoff, logfc.threshold = logfc.threshold
                                 ) 
        geneList <- geneList[[1]]
        
        #filter FindMarkers output to only retain features that pass the p_val threshold
        geneList <- geneList[geneList$p_val_adj < pValCutoff,]
        geneList$gene <- rownames(geneList)
        
        #extract average expression values for plotting (FYI: AverageExpression() is applied to non-logged data!!)
        Idents(seu.sub) <- comparision
        avg.seu.sub <- log1p(AverageExpression(seu.sub, verbose = FALSE)$RNA)
        avg.seu.sub <- as.data.frame(avg.seu.sub)
        
        #convert groups to "X" and "Y" for plotting where contrast[2] is on x-axis and contrast[1] is on y-axis
        avg.seu.sub <- avg.seu.sub[ ,contrast[c(2,1)]]
        colnames(avg.seu.sub) <- c("X","Y")

        #do some manupulation to label select data points
        #this is recomened approach
        if(!useLineThreshold){
            #join the coordinates df with the DEG results then run filtering to select data points to label
            #this labels the top n selected by the user, ranked by abs(X*Y) to find the data points furthest from n.s.

            avg.seu.sub <- avg.seu.sub %>% mutate(gene=rownames(avg.seu.sub)) %>% 
            left_join(geneList, by = "gene") %>% 
            mutate(direction=case_when(avg_log2FC > 0 ~ "up",
                                       avg_log2FC < 0 ~ "dwn",
                                       is.na(avg_log2FC) ~ "NA"),
                   direction=ifelse(gene %in% rownames(geneList)[grepl(paste(labsHide,collapse="|"), row.names(geneList))], "exclude", direction),
                   residual=case_when(direction == "up" ~ abs(Y-X),
                                      direction == "dwn" ~ abs(Y-X),
                                      direction == "NA" ~ 0,
                                      direction == "exclude" ~ 0)) %>% arrange(desc(residual)) %>% group_by(direction) %>% 
            mutate(lab=ifelse(row_number() <= labCutoff
                              & direction != "NA" 
                              & direction != "exclude"
                              & gene %in% c(rownames(geneList),addLabs), gene, ifelse(gene %in% addLabs, gene, NA)),
                   lab_col=case_when(direction == "up" ~ colUp,
                                     direction == "dwn" ~ colDwn,
                                     direction == "NA" ~ "black",
                                     direction == "exclude" & avg_log2FC > 0 ~ colUp,
                                     direction == "exclude" & avg_log2FC < 0 ~ colDwn)
                  )
        }else{
            avg.seu.sub <- avg.seu.sub %>% mutate(gene=rownames(avg.seu.sub),
                                                  up=threshold+1*X, 
                                                  dn=-threshold+1*X,
                                                  direction=case_when(Y > up ~ "up",
                                                                      Y < dn ~ "dwn",
                                                                      Y < up && Y > dn ~ "NA"),
                                                  residual=case_when(direction == "up" ~ Y-up,
                                                                     direction == "dwn" ~ dn-Y,
                                                                     direction == "NA" ~ 0),
                                                 ) %>% group_by(direction) %>% 
            arrange(desc(residual)) %>% 
            mutate(lab=ifelse(row_number() <= labCutoff 
                              & direction != "NA" 
                              & gene %in% c(rownames(geneList),addLabs), gene, ifelse(gene %in% addLabs, gene, NA)),  
                   lab_col=case_when(direction == "up" ~ colUp,
                                     direction == "dwn" ~ colDwn,
                                     direction == "NA" ~ "black")
                  )
        }
        
        #check if there was at least 1 DGE or force return == T before saving the plot
        if(length(na.omit(avg.seu.sub$lab)) > 0 | forceReturn == T){
            outfile <- paste0(outDir, outName, "_", gsub(" ", "_", x),"_linear_deg_by_", gsub(" ", "_", groupBy),".png")
            
            #make the plot
            p <- ggplot(data=avg.seu.sub, aes(x = X, y = Y, label=lab)) + 
            ggtitle(x, 
                    if(subtitle) {subtitle = paste("Average gene expression (", contrast[1]," vs ", contrast[2],")", sep = "")}
                   ) +
            geom_point(color = avg.seu.sub$lab_col) + 
            labs(x = contrast[2], y = contrast[1]) +
            {if(thresLine)geom_abline(intercept = threshold, slope = 1)} +
            {if(thresLine)geom_abline(intercept = -threshold, slope = 1)} + 
            geom_label_repel(max.overlaps = Inf, size=5, color = avg.seu.sub$lab_col, max.time = 2, max.iter = 10000000) + 
            theme_classic() + 
            {if(noTitle)theme(axis.title = element_text(size= 20),
                  axis.text = element_text(size= 14),
                  title = element_blank(),
                  plot.subtitle = element_text(size= 14)
                 )
             } + 
            {if(!noTitle)theme(axis.title = element_text(size= 20),
                  axis.text = element_text(size= 14),
                  title = element_text(size= 24),
                  plot.subtitle = element_text(size= 14)
                 )
                
            }
        
            #save the plot
            ggsave(outfile)
            
            #if request to return, then do so
            if(returnPlots){
                return(p)
                
            }
            
            #if request to return, then do so -- this returns features upregulated in contrast[1]
            if(returnUpList){
                up <- avg.seu.sub[avg.seu.sub$direction == "up",]
                upList <- up$gene
                return(upList)
                
            }
            
            #if request to return, then do so -- this returns features upregulated in contrast[2]
            if(returnDwnList){
                dwn <- avg.seu.sub[avg.seu.sub$direction == "dwn",]
                dwnList <- dwn$gene
                return(dwnList)
                
            }
            
            }
    })
}

############ formatUMAP ############
formatUMAP <- function(plot = NULL, smallAxes = F) {
    
    plot <- plot + labs(x = "UMAP1", y = "UMAP2") +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(),
          axis.title = element_text(size= 20),
          plot.title = element_blank(),
          title = element_text(size= 20),
          axis.line = element_blank(),
          panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = 2)
          )
    
    if(smallAxes){
       
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
        
        plot <- plot + theme(axis.title = element_blank(),
                             panel.border = element_blank(),
                             plot.margin = unit(c(-7, -7, -7, -7), "pt")
                            )

        plot <- plot + inset_element(axes,left= 0,
                                bottom = 0,
                                right = 0.25,
                                top = 0.25,
                                align_to = "full"
                               )
    }
    
    return(plot)
}


############ cusUMAP ############
#this function requires a UMAP plot gerneated using DimPlot with label = T, label.box = T
cusLabels <- function(plot = NULL, shape = 21, labCol = "black", size = 8, alpha = 1, rm.na = T, nudge_x = NULL, nudge_y = NULL, textSize = 4, smallAxes = FALSE
                  ) {
    
    #extract label coords and colors
    g <- ggplot_build(plot)
    #pointData <- as.data.frame(g$data[[1]][!duplicated(g$data[[1]][,"colour"]),])
    #rownames(pointData) <- NULL
    #pointData$group <- pointData$group-1
    #pointData$group <- as.factor(pointData$group)
    
    labCords <- as.data.frame(g$data[2]) #add error if labels are not present
    labCordz <- labCords[,c(1:4,7)]

    colnames(labCordz) <- c("colour", "UMAP1", "UMAP2", "clusterID", "labCol")
#     labCordz$clusterID <- as.character(labCordz$clusterID)
#     labCordz$clusterID <- labCordz$clusterID-1
    
    #labCordz <- left_join(labCordz, pointData[ , c("group", "colour")], by = c("clusterID" = "group"), keep = T)
    labCordz <- labCordz[order(labCordz$clusterID),]
    labCordz$labCol <- labCol
    #rownames(labCordz) = seq(length=nrow(labCordz))
    #labCordz$clusterID <- as.factor(sort(as.numeric(as.character(labCordz$clusterID))))
    #labCordz <- labCordz[order(as.numeric(as.character(labCordz$clusterID))),]
    
    if(rm.na == T){
        labCordz <- na.omit(labCordz)
    }
    
    if(!is.null(nudge_x)){
        labCordz$UMAP1 <- labCordz$UMAP1 + nudge_x
    }
    
    if(!is.null(nudge_y)){
        labCordz$UMAP2 <- labCordz$UMAP2 + nudge_y
    }

    #rownames(labCordz) <- NULL
    #remove old labels
    plot$layers[2] <- NULL

    #add labels to the stripped plot to create final image
    pi <- plot + geom_point(data = labCordz, aes(x = UMAP1, y = UMAP2),
                            shape=shape,
                            size=size,
                            fill=labCordz$colour,
                            stroke=1,
                            alpha=alpha,
                            colour="black") +
    geom_text(data = labCordz, size = textSize, mapping = aes(x = UMAP1, y = UMAP2), label = labCordz$clusterID, color = labCordz$labCol)
    
    plot <- formatUMAP(pi, smallAxes = smallAxes)

    
    return(plot)
}

############ freqPlots ############

freqPlots <- function(seu.obj = NULL, groupBy = "clusterID", refVal = "orig.ident", comp = "cellSource", colz = NULL, namez = NULL, method = 1, nrow = 3, title = F, legTitle = NULL, no_legend = F, showPval = T
                       ) {
    
    if(method == 1){
    fq <- prop.table(table(seu.obj@meta.data[[groupBy]], seu.obj@meta.data[,refVal]), 2) *100
    df <- reshape2::melt(fq, value.name = "freq", varnames = c(groupBy, 
                                                               refVal)) 
        } else if(method == 2){
        
        Idents(seu.obj) <- refVal
        set.seed(12)
        seu.sub <- subset(x = seu.obj, downsample = min(table(seu.obj@meta.data[[refVal]])))

        fq <- prop.table(table(seu.sub@meta.data[[refVal]], seu.sub@meta.data[,groupBy]), 2) *100       
        df <- reshape2::melt(fq, value.name = "freq", varnames = c(refVal,
                                                                   groupBy))
    
    } else{print("Method is not recognized. Options are method = 1 (by sample) OR method = 2 (by cluster).")
          break()
          }
    
    df <- df %>% left_join(seu.obj@meta.data[ ,c(refVal,comp)][!duplicated(seu.obj@meta.data[ ,c(refVal,comp)][,1]),], by = refVal)
    #df <- df[grepl("\\b9\\b|12|28|33|41",df$clusterID),] #12 # add grep for only certain groupBy vals
    #colnames(df) <- c("clusterID","orig.ident","freq","cellSource")
    
    if(!is.null(colz)){
        if(length(colz) == 1){
            df <- df %>% left_join(seu.obj@meta.data[ ,c(refVal,colz)][!duplicated(seu.obj@meta.data[ ,c(refVal,colz)][,1]),], by = refVal)
            warn1 <- F
        }else{print("Warning: It is ideal if colors are stored in meta.data slot of the Seurat object. Colors may be mismatched.")
              warn1 <- T
             }
    }else{warn1 <- F}
    
    if(is.null(namez)){namez <- refVal}
    
    if(ifelse(length(namez) > 1, "", namez) != refVal){
         if(length(namez) == 1){
             df <- df %>% left_join(seu.obj@meta.data[ ,c(refVal,namez)][!duplicated(seu.obj@meta.data[ ,c(refVal,namez)][,1]),], by = refVal)
             warn2 <- F
         }else{print("Warning: It is ideal if names are stored in meta.data slot of the Seurat object. Names may be mismatched.")
               warn2 <- T
               }
    }else{warn2 <- F}
    
    p <- ggplot(df, aes_string(y = "freq", x = comp)) + 
    labs(x = NULL, y = "Proportion (%)") + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(), 
          strip.background = element_rect(fill = NA, color = NA), 
          strip.text = element_text(face = "bold"), 
          axis.ticks.x = element_blank(), 
          axis.text = element_text(color = "black")         )
    
    if(is.null(legTitle)){
        legTitle <- refVal
    }
    
    
    #note these are not corrected p-values
    pi <- p + facet_wrap(groupBy, scales = "free_y", nrow = nrow) + 
    guides(fill = "none") + 
    geom_boxplot(aes_string(x = comp), alpha = 0.25, outlier.color = NA) + 
    geom_point(size = 2, position = position_jitter(width = 0.25),
               aes_string(x = comp, y = "freq", color = refVal)) +
    labs(color = legTitle) +
    {if(showPval){ggpubr::stat_compare_means(aes(label = paste0("p = ", ..p.format..)), label.x.npc = "left", label.y.npc = 1,vjust = -1, size = 3)}} + 
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    theme(panel.grid.major = element_line(color = "grey", size = 0.25),
          #legend.position = "none",
          text = element_text(size = 12)
          ) +                    
    {if(!is.null(colz)){
        scale_color_manual(labels = if(warn1){
            namez
        }else{
            levels(df[[namez]])
        },
                           values = if(warn2){
                               colz
                           }else{
                               levels(df[[colz]])
                           })}} + {if(no_legend == T){NoLegend()}}
                                           
    return(pi)

}
                                           
############ vilnSplitComp ############

vilnSplitComp <- function(seu.obj = NULL, groupBy = "clusterID", refVal = "cellSource", outName = "", nPlots = 9, saveOut = T, saveGeneList = F, returnGeneList = F, contrast = NULL, filterOutFeats = c("^MT-","^RPL","^RPS"), logfc.threshold = 0.5, pValCutoff = 0.01,
                          outDir = "./output/"
                       ) {
    
    #set the groups to contrast
    DefaultAssay(seu.obj) <- "RNA"
    seu.obj$cluster.condition <- paste(seu.obj@meta.data[[groupBy]], seu.obj@meta.data[[refVal]], sep = "_")
    Idents(seu.obj) <- "cluster.condition"
    
    p <- lapply(levels(seu.obj@meta.data[[groupBy]]), function(x) {

        #split the groups by condition
        comp1 <- paste(x, gsub(" ", "_", contrast[1]), sep = "_") 
        comp2 <- paste(x, gsub(" ", "_", contrast[2]), sep = "_")
        
        #use FindMarkers to complete DGE analysis -- default parameters used
        message(paste0("Contrasting - (1) ", comp1," - VS - (2) ", comp2))
        cluster.markers <- FindMarkers(seu.obj, ident.1 = comp1,  ident.2 = comp2, min.pct = 0, logfc.threshold = logfc.threshold) 
        message(paste0("Number of DEGs prior to p_val_adj filter:"))
        print(cluster.markers %>% mutate(direction = ifelse(avg_log2FC > 0, "Up", "Down")) %>% group_by(direction) %>% summarize(nRow = n()))
        
        #filter the results based on the user specified adjusted p-value
        cluster.markers <- cluster.markers %>% filter(p_val_adj < pValCutoff)
        message(paste0("Number of DEGs following p_val_adj filter:"))
        print(cluster.markers %>% mutate(direction = ifelse(avg_log2FC > 0, "Up", "Down")) %>% group_by(direction) %>% summarize(nRow = n()))

        #save the output DGE results if requested
        if(saveGeneList){
            cluster.markers$cellType <- x
            write.csv(cluster.markers, file = paste0(outDir, outName, "_", gsub(" ", "_", x) ,"_geneList.csv"))
        }

        #return the gene list if requested
        if(returnGeneList){
            return(cluster.markers)
        }
        
        #filter out requested features using regex
        cluster.markers <- cluster.markers[!grepl(paste(filterOutFeats,collapse="|"), row.names(cluster.markers)), ]
        
        #plot the desired number of DEGs by subsetting the DEGs
        cluster.markers.toPlot <- cluster.markers %>% top_n(n = -nPlots, wt = p_val)
        if(length(rownames(cluster.markers.toPlot)) > 0){
            seurat.vlnplot <- VlnPlot(object = seu.obj,
                                      idents = c(comp1,comp2),
                                      features = rownames(cluster.markers.toPlot),
                                      split.by = refVal
            )
            
            if(saveOut){
                png(file = paste0(outDir, outName,"_", x,"_top",nPlots,"_comp_vlnPlot.png"), width=2520, height=1460)
                print(seurat.vlnplot)
                dev.off()
            }
        }
    })
}
############ loadMeta ############

loadMeta <- function(seu.obj = NULL, seu.file = NULL, metaFile = "", groupBy = "clusterID", metaAdd = "majorID", 
                     save = FALSE, outName = "", header = TRUE
                    ){
    
    #make this an lapply incase you have a vector of meta data to add
    metaData <- read.csv(file = metaFile, header = header)
    new.cluster.ids <- metaData[[metaAdd]]
    seu.obj <- SetIdent(seu.obj, value = groupBy)
    
    names(new.cluster.ids) <- metaData[[groupBy]]
    seu.obj <- RenameIdents(seu.obj, new.cluster.ids)
    table(Idents(seu.obj))
    seu.obj@meta.data[[metaAdd]] <- Idents(seu.obj)
    
    return(seu.obj)
    
    #make it so you can save the updated rds

}
                       
############ vilnSplitCompxGene ############
    
vilnSplitCompxGene <- function(seu.obj = NULL, groupBy = "clusterID_sub", comp = "cellSource", metaAdd = "majorID", features = "", labelData = NULL,
                               cols = c("mediumseagreen","mediumpurple1"), save = FALSE, outName = "", outDir = "",
                               height = 4, width = 8
                              ){
    
    DefaultAssay(seu.obj) <- "RNA"
    Idents(seu.obj) <- groupBy

    rectangles <- data.frame(
        xmin = seq(1,length(levels(seu.obj)), by = 2)+0.99,
        xmax = seq(1,length(levels(seu.obj)), by = 2)+1.01,
        ymin = 6,
        ymax = 0
    )

    
    p <- VlnPlot(
        features = features,
        object = seu.obj,
        split.by = comp,
        cols = cols,
        stack = TRUE,
        flip = TRUE
    ) + theme(legend.position = "top",
                 axis.title = element_blank(),
                 axis.text = element_blank(),
                 axis.ticks.y = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.text.x = element_blank(),
                  plot.margin = margin(0, 0, 0, 3, "pt"),
                  legend.text = element_text(margin = margin(r = 10, unit = "pt"))
                 ) & theme(axis.text.x = element_blank())
        
        if(save == T){
            outfile <- paste(outDir, outName,"_", x,"_comp_vlnPlot.png", sep = "") 
            ggsave(plot = p, outfile, height = height, width = width) 
        }
        
    
#         #plot lablels for bottom of plot
#         g <- ggplot(colArray, aes(x = clusterID_sub, y = 1, fill = clusterID_sub, label = clusterID_sub)) + 
#     geom_tile(fill = "transparent") +
#     geom_point(shape = 21, stroke = 0.75, size = 7) +
#     geom_text(fontface = "bold", size = 3, color = colArray$labCol) + theme_bw(base_size = 12) +
#     scale_fill_manual(values = colArray$colour) + scale_y_discrete(expand = c(0, 0)) +
#     theme(legend.position = "none", panel.spacing = unit(0, "lines"),
#               panel.background = element_blank(), 
#               panel.border = element_blank(),
#               plot.background = element_blank(), 
#               plot.margin = margin(0, 0, 0, 3, "pt"),
#               axis.title.y = element_blank(),
#               axis.ticks.y = element_blank(),
#               axis.ticks.x = element_blank(),
#               axis.text = element_blank(),
#                   panel.grid.major = element_blank(), 
#                   panel.grid.minor = element_blank()) + xlab("Cluster") 
#      pi <- plot_grid(p, g, ncol = 1, rel_heights = c(0.85, 0.15), align = "v", axis = "lr") 
    
    return(p)
}

############ getPB ############

#define the function
#where bioRep == annotation vector cell & orig.ident data (with levels for each patient)
getPb <- function(mat.sparse, bioRep) {
   mat.summary <- do.call(cbind, lapply(levels(bioRep$cellSource), function(rep) {
     cells <- row.names(bioRep)[bioRep$cellSource==rep]
     pseudobulk <- Matrix::rowSums(mat.sparse[,cells])
     return(pseudobulk)
   }))
    
   colnames(mat.summary) <- levels(bioRep$cellSource)
   return(mat.summary)
}
                       
############ createPB ############

createPB <- function(seu.obj = NULL, groupBy = "clusterID_sub", comp = "cellSource", biologicalRep = "orig.ident",
                     clusters = NULL, outDir = "", min.cell = 25, lowFilter = F, dwnSam =T, featsTOexclude = NULL, cnts = T,
                     grepTerm = NULL, grepLabel = NULL #improve - fix this so it is more functional
                    ){


    Idents(seu.obj) <- groupBy
    
    #extract number of biological replicates in the grouping system
    ztest <- dim(table(seu.obj@meta.data[[biologicalRep]]))
    
    groupz <- ifelse(is.null(clusters),levels(seu.obj),list(clusters)) #improve - fix this so you can choose not to run the conversion on all samples...
    #loop though all the groups
    test <- lapply(levels(seu.obj), function(x) {
        
        #subset the data by group
        seu.sub <- subset(seu.obj, idents = x)
        
        #determine value to downsample by and extract metadata - uses bioRep with lowest number of cells > "25"
        z <- table(seu.sub@meta.data[[biologicalRep]])
        z <- z[z>min.cell]
        zKeep <- names(z)

        if(length(zKeep) >= 0.5*ztest){
        
            #remove samples that have insufficent cell numbers then normalize the data
            Idents(seu.sub) <- biologicalRep
            seu.sub.clean <- subset(seu.sub, idents = zKeep)
            seu.sub.clean <- NormalizeData(seu.sub.clean)
    
            
            seu.sub.clean@meta.data[[biologicalRep]] <- as.factor(seu.sub.clean@meta.data[[biologicalRep]])
            seu.sub.clean@meta.data[[biologicalRep]] <- droplevels(seu.sub.clean@meta.data[[biologicalRep]])
            
            
            z <- table(seu.sub.clean@meta.data[[biologicalRep]])
            
            #extract and save metadata in a data.frame
            z <- as.data.frame(z)
            colnames(z) <- c("sampleID", "nCell")
            z$clusterID <- toString(x) #add clusterID value
            z$groupID <- ifelse(grepl(grepTerm, z$sampleID), grepLabel[1], grepLabel[2]) #fix this - hardcoded and will only work for PBMCs improve this
    
            ds <- min(z$nCell[z$nCell > min.cell])

            if(dwnSam){
            message <- paste("Downsampling cluster: ", x," at a level of ", ds, " cells per replicate", sep = "")
            print(message)
        
            #randomly downsample the subset data
            Idents(seu.sub.clean) <- biologicalRep
            set.seed(12)
            seu.sub.clean <- subset(x = seu.sub.clean, downsample = ds) #works, but stops working if the min is "0"
            }
            

            #extract required data for pseudobulk conversion
            if(cnts){
                mat <- seu.sub.clean@assays$RNA@counts
            }else{
                mat <- seu.sub.clean@assays$RNA@data
            }
            
            
            if(lowFilter){
                mat <- mat[rowSums(mat > 1) >= 10, ]
            }
            
            bioRep <- as.data.frame(seu.sub.clean@meta.data[[biologicalRep]])
            colnames(bioRep) <- "cellSource"
            row.names(bioRep) <- colnames(seu.sub.clean)
            bioRep$cellSource <- as.factor(bioRep$cellSource)

            #use custom function to convert to pseudobulk
            pbj <- getPb(mat, bioRep)
            
            if(!is.null(featsTOexclude)){
                pbj <- pbj[!rownames(pbj) %in% featsTOexclude,]
            }
        
            #remove genes with 50% "0" values
            #need to add code here
        
            #log the number of reps included
            if(length(colnames(pbj))-1 != ztest){
                message <- paste("The following replicates were used for psudobluk conversion: ",as.list(colnames(pbj)), sep = "")
                print(message)
            } else {
                message <- "All replicates were used for psudobluk conversion"
                print(message)
            }
            
            outfile <- paste(outDir,x,"_pb_matrix.csv", sep = "")
            write.csv(pbj, file = outfile)
        } else {
            message <- paste("Unable to downsample cluster: ",x," due to insufficent cell numbers", sep = "")
            print(message)
        } 
        return(z) 
    })
    
    df <- do.call(rbind, test)
    
    csvOut <- paste(outDir,groupBy ,"_deg_metaData.csv", sep = "")
    write.csv(df, file = csvOut)
}
                       
############ pseudoDEG ############
# contrast will be idents.1_NAME vs idents.2_NAME !!!
pseudoDEG <- function(metaPWD = "", metaPass = F, padj_cutoff = 0.1, lfcCut = 0.58, outDir = "", outName = "", idents.1_NAME = NULL, idents.2_NAME = NULL, returnDDS = F,
                     inDir = "", title = "", fromFile = T, meta = NULL, pbj = NULL, returnVolc = F, paired = F, pairBy = "", minimalOuts = F, saveSigRes = T, topn=c(20,20), degFormula = NULL,saveFullRes = F,
                     filterTerm = "^ENSCAF", addLabs = NULL, mkDir = F, dwnCol = "blue", stblCol = "grey",upCol = "red", labSize = 3
                     ){
    if(fromFile){
        files <- list.files(path = inDir, pattern="pb_matrix.csv", all.files=FALSE,full.names=FALSE)
        clusters <- unname(sapply(files, function(x) {unlist(strsplit(x, split = "_pb_"))[1]}))
        outfileBase <- paste(outDir, outName, "_cluster_", sep = "")
    }else{
        clusters <- outName
        outfileBase <- outDir
    }
        

    lapply(clusters, function(x) {
        if(fromFile){
            inFile <- paste(inDir, x,"_pb_matrix.csv", sep = "")
            pbj <- read.csv(file = inFile, row.names = 1) #will likely want to remove alll rbc/platlet related genes
            pbj <- pbj[!apply(pbj==0, 1, all),]
            
            if(!metaPass){
                meta <- read.csv(file = metaPWD, row.names = 1)
                meta[ ,colnames(meta)] <- lapply(meta[ ,colnames(meta)] , factor)
            }
            
            meta <- meta[meta$clusterID == x, ]
        }
        
        if(mkDir){
            outDir <- paste(outDir, "/", x, "/", sep = "")
            dir.create(outDir)
            outfileBase <- paste(outDir, outName, "_cluster_", sep = "")
            }
        
        if(paired){
            dds <- DESeqDataSetFromMatrix(round(pbj), 
                                          colData = meta, ### add pt meta data here
                                          design = formula(paste("~ groupID + ",noquote(pairBy), sep = "")))
        } else if(!is.null(degFormula)){
            dds <- DESeqDataSetFromMatrix(round(pbj), 
                                          colData = meta, ### add pt meta data here
                                          design = formula(noquote(degFormula))
                                         )
        }
        else{
            dds <- DESeqDataSetFromMatrix(round(pbj), 
                                          colData = meta, ### add pt meta data here
                                          design = ~ groupID)
        }
        
        if(returnDDS){
            return(dds)
        }
        
        #transforma and plot the data with PCA
        rld <- varianceStabilizingTransformation(dds, blind=TRUE)
        if(!minimalOuts){
            outfile <- paste(outfileBase, x,"_pca.png", sep = "")
            print(outfile)
            p <- DESeq2::plotPCA(rld, intgroup = "groupID")
            ggsave(outfile, width = 7, height = 7)
            
            outfile <- paste(outfileBase, x,"_pca2.png", sep = "")
            p <- DESeq2::plotPCA(rld, intgroup = "sampleID")
            ggsave(outfile, width = 7, height = 7)
            
        }
        
        #set up QC to evaluate treatment seperation by heatmap & plot
        rld <- rlog(dds, blind=FALSE)
        rld_mat <- assay(rld)
        rld_cor <- cor(rld_mat)

        if(!minimalOuts){
            outfile <- paste(outfileBase, x,"_pheatmap.png", sep = "")
            p <- pheatmap::pheatmap(rld_cor)
            ggsave(p, file = outfile)
        }

        #run DESeq2 then plot dispersion estimates
        dds <- DESeq(dds)
        
        #extract Wald test results

        if(!is.null(idents.1_NAME) | !is.null(idents.2_NAME)){
            contrast <- c("groupID", idents.1_NAME, idents.2_NAME)
        }else{
            contrast <- c("groupID", unique(meta$groupID)[1], unique(meta$groupID)[2])
            print("Variables idents.1_NAME and/or idents.2_NAME not specify, this may impact directionality of contrast - confirm results are as expect and/or specify ident names.")
        }
        

        #perform logFoldChange and shrinkage
        res <- results(dds, 
                       contrast = contrast,
                       alpha = padj_cutoff,
                       lfcThreshold = lfcCut,
                       cooksCutoff=FALSE)
        summary(res)
        
        res <- lfcShrink(dds, 
                         contrast = contrast,
                         res=res,
                         type="normal") #would prefer to use something other than normal
        
#         res <- lfcShrink(dds, 
#                          coef = resultsNames(res)[1],
#                          type="apeglm") #would prefer to use something other than normal
        
        #extract the results and save as a .csv
        res_tbl <- res %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()
        
        if(saveFullRes){
            write.csv(res_tbl,
                      file = paste(outfileBase, x,"_full_all_genes.csv", sep = ""),
                      quote = FALSE,
                      row.names = FALSE)
        }
        
        sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
        dplyr::arrange(padj)

        if(saveSigRes){
            sig_res$gs_base <- toupper(x)
            
            write.csv(sig_res,
                      file = paste(outfileBase, x,"_all_genes.csv", sep = ""),
                      quote = FALSE,
                      row.names = FALSE)
        }
        
        #get nomarlized counts and plot top 20 DEGs
        normalized_counts <- counts(dds, 
                                    normalized = TRUE)

        top20_sig_genes <- sig_res %>%
        dplyr::arrange(padj) %>%
        dplyr::pull(gene) %>%
        head(n=20)

        top20_sig_norm <- data.frame(normalized_counts) %>%
        rownames_to_column(var = "gene") %>%
        dplyr::filter(gene %in% top20_sig_genes)

        gathered_top20_sig <- top20_sig_norm %>%
        gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")

        gathered_top20_sig <- meta %>% inner_join(gathered_top20_sig, by = c("sampleID" = "samplename")) #need more metadata
    
        if(dim(gathered_top20_sig)[1] > 0){ 
            
            #plot with ggplot2
            if(!minimalOuts){
                outfile <- paste(outfileBase, x,"_genePlot.png", sep = "")
                p <- ggplot(gathered_top20_sig) +
                geom_point(aes(x = gene, 
                               y = normalized_counts, 
                               color = groupID), #need to fix samplename & change to groupID
                           position=position_jitter(w=0.1,h=0)) +
                scale_y_log10() +
                xlab("Genes") +
                ylab("log10 Normalized Counts") +
                ggtitle("Top 20 Significant DE Genes") +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                theme(plot.title = element_text(hjust = 0.5))
                ggsave(p, file = outfile)
            }
            
            #extract sig DEGs based on normalized data and plot heat map
            sig_norm <- data.frame(normalized_counts) %>%
            rownames_to_column(var = "gene") %>%
            dplyr::filter(gene %in% sig_res$gene)

            rownames(sig_norm) <- sig_norm$gene
            colAnn <- as.data.frame(meta[,"groupID"], colname = "groupID")
            rownames(colAnn) <- meta[,"sampleID"]
            colnames(colAnn) <- "groupID"

            heat_colors <- brewer.pal(6, "YlOrRd")

#             if(!minimalOuts){
#                 outfile <- paste(outfileBase, x,"_pheatmapComp.png", sep = "")
#                 p <- pheatmap(sig_norm[ , 2:length(colnames(sig_norm))], 
#                               color = heat_colors, 
#                               cluster_rows = FALSE, 
#                               show_rownames = TRUE,
#                               annotation = colAnn, 
#                               #border_color = "grey", 
#                               fontsize = 10, 
#                               scale = "row", 
#                               fontsize_row = 10, 
#                               height = 20)    
#                 ggsave(p, file = outfile)
#             }

            #set threshold and flag data points then plot with ggplot
            
            res_table_thres <- res_tbl %>% mutate(threshold = ifelse(padj < padj_cutoff & abs(log2FoldChange) >= lfcCut, 
                                            ifelse(log2FoldChange > lfcCut ,'Up','Down'),'Stable')
                                                 )
            
            res_table_thres <- res_table_thres[!is.na(res_table_thres$padj),]


            res_table_thres.sortedByPval = res_table_thres[order(res_table_thres$padj),]
            res_table_thres.sortedByPval <- res_table_thres.sortedByPval[!grepl(filterTerm, res_table_thres.sortedByPval$gene),]
            top20_up <- res_table_thres.sortedByPval[res_table_thres.sortedByPval$threshold == "Up",] %>%  do(head(., n=topn[1]))
            top20_down <- res_table_thres.sortedByPval[res_table_thres.sortedByPval$threshold == "Down",] %>%  do(head(., n=topn[2]))

            res_table_thres <- res_table_thres %>% mutate(label = ifelse(gene %in% top20_up$gene | gene %in% top20_down$gene | gene %in% addLabs, gene, NA)
                         )
        
            cntUp <- nrow(res[which(res$log2FoldChange > lfcCut & res$padj < padj_cutoff),])
            cntDwn <- nrow(res[which(res$log2FoldChange < -lfcCut & res$padj < padj_cutoff),])

            outfile <- paste(outfileBase, x,"_volcano.png", sep = "")
            
            if(fromFile){
                if(is.null(title)){
                    title <- paste(idents.1_NAME, " vs ",idents.2_NAME, "within", x, sep="")
                }
            }
            
            p <- ggplot(data = res_table_thres,
                        aes(x = log2FoldChange, 
                            y = -log10(padj), 
                            colour = threshold)) +
            geom_vline(xintercept = c(-lfcCut,lfcCut), lty = 4, col="black", lwd = 0.8) +
            geom_hline(yintercept = -log10(padj_cutoff), lty = 4, col="black", lwd = 0.8) +
            geom_point(alpha = 0.4, size = 3.5) +
            geom_label_repel(max.overlaps = Inf, size=labSize, label = res_table_thres$label, show.legend = FALSE) +
            labs(x="log2(fold change)",
                 y="-log10(padj)",
                 title=title) + 
            scale_color_manual(values=c("Down" = dwnCol, "Stable" = stblCol,"Up" = upCol), labels=c(paste("Down (", cntDwn,")", sep = ""), "Stable", paste("Up (", cntUp,")", sep = ""))) +
            theme_bw() +
            theme(plot.title = element_text(size = 20, hjust=0.5), 
                  legend.position = "right", 
                  legend.title = element_blank()
                 )
            ggsave(p, file = outfile)
            
            if(returnVolc){
                return(p)
            }
            

        }
    })
}
############ btwnClusDEG ############
#work in progress             - need to to fix doLinDEG option ### NOTE: cannot have special char in ident name
btwnClusDEG <- function(seu.obj = NULL,groupBy = "majorID_sub", idents.1 = NULL, idents.2 = NULL, bioRep = "orig.ident",padj_cutoff = 0.01, lfcCut = 0.5, topn=c(20,20),
                        minCells = 25, outDir = "", title = NULL, idents.1_NAME = "", idents.2_NAME = "", returnVolc = F, doLinDEG = F, paired = T, addLabs = NULL, lowFilter = F, dwnSam = T, setSeed = 12, dwnCol = "blue", stblCol = "grey",upCol = "red", labSize = 3
                    ){
    
    seu.integrated.obj <- seu.obj
    DefaultAssay(seu.integrated.obj) <- "RNA"

    Idents(seu.integrated.obj) <- groupBy
    if(is.null(idents.2)){
        idents.2 <- levels(seu.integrated.obj)[!levels(seu.integrated.obj) %in% idents.1]
    }
    
    seu.sub <- subset(seu.integrated.obj, idents = as.list(append(idents.1,idents.2)))
    
    if(all(idents.1 %in% levels(seu.integrated.obj))){
        if(length(idents.1) > 1){
            grepTerm <- paste(idents.1,collapse="|")
        }else{
            grepTerm <- idents.1
        }
    }else{
        print("Idents not found in seurat object. Please check that they are correct and try again.")
        break()
    }

    if(is.null(title)){
        title <- paste(gsub(" ", "_", idents.1_NAME), "_vs_",gsub(" ", "_",idents.2_NAME))
    }
    seu.sub$compare <- ifelse(grepl(grepTerm, seu.sub@meta.data[[groupBy]]), "idents.1", "idents.2")
    
    Idents(seu.sub) <- "compare"

    #extract number of biological replicates in the grouping system
    ztest <- dim(table(seu.sub@meta.data[[bioRep]]))

    seu.sub@meta.data$conditional <- paste(seu.sub@meta.data[[bioRep]], seu.sub@meta.data$compare, sep = "_")

    z <- table(seu.sub@meta.data$conditional)
    z <- z[z>minCells]
    zKeep <- names(z)

    Idents(seu.sub) <- "conditional"
    seu.sub.clean <- subset(seu.sub, idents = zKeep)
    #seu.sub.clean <- NormalizeData(seu.sub.clean)
    
    z <- table(seu.sub.clean@meta.data$conditional)
    
    #extract and save metadata in a data.frame
    z <- as.data.frame(z)
    colnames(z) <- c("sampleID", "nCell")
    z$groupID <- ifelse(grepl("idents.1", z$sampleID), idents.1_NAME, idents.2_NAME)
    
    ds <- min(z$nCell)

    if(dwnSam){
    message <- paste("Downsampling at a level of ", ds, " cells per replicate", sep = "")
    print(message)
        
    #randomly downsample the subset data
        Idents(seu.sub.clean) <- "conditional"
        set.seed(setSeed)
        seu.sub.clean <- subset(x = seu.sub.clean, downsample = ds)
    }

    #extract required data for pseudobulk conversion
    mat <- seu.sub.clean@assays$RNA@counts
    if(lowFilter){
        mat <- mat[rowSums(mat > 1) >= 10, ]
    }
    
    bioRep <- as.data.frame(seu.sub.clean@meta.data$conditional)
    colnames(bioRep) <- "cellSource"
    row.names(bioRep) <- colnames(seu.sub.clean)
    bioRep$cellSource <- as.factor(bioRep$cellSource)

    #use custom function to convert to pseudobulk
    pbj <- getPb(mat, bioRep)
    pbj <- pbj[rowSums(pbj[])>0,] # remove any rows that are all zeros

    ### pseudobulk DEG using DEseq2 backbone ###
    meta <- z
    
    meta$bioRepPair <- str_split_fixed(meta$sampleID, "_ident.", 2)[,1]
    if(doLinDEG == T){
        seu.sub.clean$compareLinDEG <- ifelse(grepl(grepTerm, seu.sub.clean@meta.data$clusterID), "idents.1", "idents.2")
        seu.sub.clean$groupz <- "cellz"
        linDEG(seu.obj = seu.sub.clean, threshold = 1, thresLine = T, groupBy = "groupz", comparision = "compareLinDEG", outDir = outDir, outName = paste(gsub(" ", "_", idents.1_NAME), "_vs_",gsub(" ", "_",idents.2_NAME), sep = ""), colUp = "red", colDwn = "blue", subtitle = T, returnUpList = F, returnDwnList = F, forceReturn = T
              )
    }
    print(meta)
    p <- pseudoDEG(padj_cutoff = padj_cutoff, lfcCut = lfcCut, outName = paste(gsub(" ", "_", idents.1_NAME), "_vs_",gsub(" ", "_",idents.2_NAME), sep = ""), 
              outDir = outDir, title = title, fromFile = F, meta = meta, pbj = pbj, returnVolc = returnVolc, paired = paired, pairBy = "bioRepPair",
                   idents.1_NAME = idents.1_NAME, idents.2_NAME = idents.2_NAME, minimalOuts = T, saveSigRes = T, addLabs = addLabs, topn = topn, dwnCol = dwnCol, stblCol = stblCol,upCol = upCol, labSize = labSize
                     )    
    return(p)
}
############ volcFromFM ############
volcFromFM <- function(seu.obj = NULL, padj_cutoff = 0.01, lfcCut = 0.58, title = "",
                     clusters = NULL, outDir = "", preSub = F
                    ){    
    
    if(preSub == F & !is.null(clusters)){
        seu.obj <- subset(seu.obj,
                          subset = 
                          clusterID_sub == clusters) #this aint right fix
        }

    clusMarkers <- FindMarkers(seu.obj, ident.1 = "OS", group.by = "cellSource")

    clusMarkers <- tibble::rownames_to_column(clusMarkers, "gene")

    #set threshold and flag data points then plot with ggplot
    res_table_thres <- clusMarkers %>% mutate(threshold = ifelse(p_val_adj < padj_cutoff & abs(avg_log2FC) >= lfcCut, 
                                                                 ifelse(avg_log2FC > lfcCut ,'Up','Down'),'Stable')
                                             )

    res_table_thres <- res_table_thres[!is.na(res_table_thres$p_val_adj),]

    res_table_thres.sortedByPval = res_table_thres[order(res_table_thres$p_val_adj),]
    cntUp <- nrow(res_table_thres[res_table_thres$threshold == "Up",])
    cntDwn <- nrow(res_table_thres[res_table_thres$threshold == "Down",])
    res_table_thres.sortedByPval <- res_table_thres.sortedByPval[!grepl("^ENSCAF|^MT-|^RPS|^RPL", res_table_thres.sortedByPval$gene),] 
    top20_up <- res_table_thres.sortedByPval[res_table_thres.sortedByPval$threshold == "Up",] %>% do(head(., n=20))
    top20_down <- res_table_thres.sortedByPval[res_table_thres.sortedByPval$threshold == "Down",] %>%  do(head(., n=20))


    res_table_thres <- res_table_thres %>% mutate(label = ifelse(gene %in% top20_up$gene | gene %in% top20_down$gene, gene, NA),
                                                  shape = ifelse(p_val_adj == 0, 17, 19)
                                                 )
    
    
    outfile <- paste(outDir, title,"_volcano.png", sep = "")

    p <- ggplot(data = res_table_thres, 
                aes(x = avg_log2FC, 
                y = -log10(p_val_adj), 
                colour = threshold)) +
    geom_vline(xintercept = c(-0.58,0.58), lty = 4, col="black", lwd = 0.8) +
    geom_hline(yintercept = -log10(padj_cutoff), lty = 4, col="black", lwd = 0.8) +
    geom_point(alpha = 0.4, size = 3.5,
                shape = res_table_thres$shape) +
    geom_label_repel(max.overlaps = Inf, size=3, label = res_table_thres$label, show.legend = FALSE) +
    #xlim(c(-2.5, 2.5)) +
    #ylim(c(0, 3)) +
    labs(x="log2(fold change)",
         y="-log10(padj)",
         title=title) + #paste("Differential expression (",idents.1_NAME ," vs ",idents.2_NAME ,")", sep = "")
    #scale_colour_discrete(labels=c(paste("Down (", cntDwn,")", sep = ""), "Stable", paste("Up (", cntUp,")", sep = ""))) + 
    scale_color_manual(values=c("Down" = "blue", "Stable" = "grey","Up" = "red"), labels=c(paste("Down (", cntDwn,")", sep = ""), "Stable", paste("Up (", cntUp,")", sep = ""))) +
    theme_bw() +
    theme(plot.title = element_text(size = 20), 
          legend.position = "right", 
          legend.title = element_blank()
         )
    
    ggsave(plot = p, file = outfile)
    
}

############ vilnPlots ############

vilnPlots <- function(seu.obj = NULL, inFile = NULL, groupBy = "clusterID", numOfFeats = 24, outName = "",
                      outDir = "", outputGeneList = T, filterOutFeats = c("^MT-", "^RPL", "^RPS"), assay = "RNA", returnViln = T, features = NULL,
                      min.pct = 0.25, only.pos = T, resume = F, resumeFile = NULL
                     ){ 
    
    if(!resume){
        if(!is.null(inFile)){
            try(seu.obj <- readRDS(file), silent = T)
        }else if(!is.null(seu.obj)){
            seu.obj <- seu.obj
        }

        DefaultAssay(seu.obj) <- assay

        Idents(seu.obj) <- groupBy
        ident.level <- levels(seu.obj@active.ident)

        cluster.markers <- FindAllMarkers(seu.obj, only.pos = only.pos, min.pct = min.pct, features = features)
        cluster.markers <- cluster.markers[cluster.markers$p_val_adj < 0.01, ]
        
        if(length(filterOutFeats) > 0){
            cluster.markers <- cluster.markers[!grepl(paste(filterOutFeats,collapse="|"), rownames(cluster.markers)), ]
        }     
        
        if (outputGeneList == TRUE){
            outfile <- paste(outDir, outName, "_gene_list.csv", sep = "")
            write.csv(cluster.markers, file = outfile)
        }  
        
    }else{
        cluster.markers <- read.csv(file = resumeFile, header = T)
        Idents(seu.obj) <- groupBy
    }
    
    if(returnViln){
        
        lapply(unique(cluster.markers$cluster), function(x) {

            geneList <- cluster.markers[cluster.markers$cluster == x, ]
            geneList <- geneList$gene

            plotList <- head(geneList, n = numOfFeats)

            seurat.vlnplot <- VlnPlot(
                object = seu.obj,
                features = rev(plotList)
            )

            outfile <- paste(outDir, outName, "_", x, "_top", numOfFeats, "_vlnPlot.png", sep = "") 
            png(file = outfile, width=2520, height=1460)

            print(seurat.vlnplot)
            dev.off()
        })
    }
    
}

############ singleR ############

singleR <- function(seu.obj = NULL, outName = "", clusters = "clusterID", outDir = ""
                     ){
    
    cntData <- GetAssayData(seu.obj, slot = "data", assay = "RNA")
    Kotliarov <- scRNAseq::KotliarovPBMCData()
    refs <- list(human_ref1 = celldex::HumanPrimaryCellAtlasData(), 
                 human_ref2 = BlueprintEncodeData(),
                 human_ref3 = DatabaseImmuneCellExpressionData(),
                 human_ref4 = NovershternHematopoieticData(), 
                 human_ref5 = MonacoImmuneData()
                )
   
    
    #Perform annotation using each reference
    sapply(names(refs), FUN = function(ref_idx) {
        ref <- refs[[ref_idx]]
        ref_name <- paste("SingleR", ref_idx, sep = ".")
        
        #Annotation is performed on cluster-level rather than default single-cell level
        rst <- SingleR(test = cntData, ref = ref, clusters = seu.obj@meta.data[[clusters]], labels = ref$label.fine)
        
        # Assign predicted cell labels to seurat object
        # SingleR assigns labels for all clusters, even for those not in the reference
        # So use pruned.labels to remove uncertain labels 
        seu.obj@meta.data[[ref_name]] <- rst$pruned.labels[match(seu.obj@meta.data[[clusters]], rownames(rst))]
        
        #Visualize cell labels
        p <- DimPlot(seu.obj, reduction = "umap", group.by = ref_name, label = TRUE)
        ggsave(filename = paste(outDir,outName,"_",ref_name, ".png", sep = ""), plot = p, width = 10, height = 7)
        
        return(p)
    })
}
    
############ cusLeg ############

cusLeg <- function(legend = NULL, colz = 2, rowz = NULL, clusLabel = "clusterID", legLabel = NULL, groupLabel = NULL, colorz = NULL, dotSize = 8,
                   groupBy = "", sortBy = "", labCol = "", headerSize = 6, #add if len labCol == 1 then add col to the df
                   cellHeader = T, bump = 2, nudge_left = 0, nudge_right = 0, topBuffer = 1.05, ymin = 0, compress_y = 0, compress_x = 0.75, titleOrder = NULL, spaceBtwnCols = NULL, breakGroups = F, horizontalWrap = F, returnData = F, overrideData = NULL
                     ){
    
    if(is.null(colorz)){
        legend$colorz <- gg_color_hue(length(legend[[clusLabel]]))
        colorz <- "colorz"
        
    }
    
    if(cellHeader == T){
        if(!is.null(titleOrder)) {
            legend[[groupLabel]] <- factor(legend[[groupLabel]], levels = titleOrder)
        }else{
            legend[[groupLabel]] <- factor(legend[[groupLabel]])
        }
        
        header <- legend  %>% group_by(!!as.name(groupLabel)) %>% summarize(n = n()) %>% arrange(factor(!!as.name(groupLabel), 
                           levels = groupBy)) %>% mutate(pos = n+1,
                                                         cum_sum = cumsum(pos),
                                                         dfNum = cum_sum-n
                                                        )  %>% as.data.frame()
        nTerms <- max(header$cum_sum) #add this var outsdide of if
        
        if(is.null(rowz)){
            rowz <- ceiling((length(legend[[clusLabel]])+dim(header)[1])/colz)
        }
    }
    
    if(length(spaceBtwnCols) < colz-1){
        message("Invalid spaceBtwnCols entry. Ensure list is correct length!")
    }
    
    #if there is one header...
    if(rowz < max(header$n) & nrow(header) == 1){
        if(!horizontalWrap){
            leg_y <- rep(seq(rowz, 1, by = -1),colz)[1:header$n]
        }else{
            leg_y <- rep(seq(rowz, 1, by = -1), each = colz)[1:header$n]
        }
        
        leg_y <- leg_y + compress_y

        if(is.null(spaceBtwnCols)){
            if(!horizontalWrap){
                leg_x <- rep(seq(1,colz)/2, each = rowz)[1:header$n]
            }else{
                leg_x <- rep(seq(1,colz)/2, rowz)[1:header$n]
            }
        }else{
            leg_x <- tryCatch(
                {
                    if(!horizontalWrap){
                        rep(cumsum(c(0.5,spaceBtwnCols)), each = rowz)[1:header$n]
                    }else{
                        rep(cumsum(c(0.5,spaceBtwnCols)), rowz)[1:header$n]
                    }
                }, 
                error = function(e) {
                    if(!horizontalWrap){
                        return(rep(seq(1,colz)/2, each = rowz)[1:header$n])
                    }else{
                        return(rep(seq(1,colz)/2, rowz)[1:header$n])
                    }
                    message("Invalid spaceBtwnCols entry. Using defeult spacing of 0.5")
                }
            )
        }
        
        header$header_x <- leg_x[header$dfNum]
        header$header_y <- leg_y[header$dfNum]+1
        headerBump = T
        
    #else there is more than one header
    }else{
        headerGroups <- seq(1,nrow(header), by = ceiling(nrow(header)/colz))
        if(length(headerGroups) != colz){
            headerGroups <- sort(c(headerGroups,seq(1,nrow(header), by = 1)[!seq(1,nrow(header), by = 1) %in% headerGroups][1]))
        }
        headData <- header[headerGroups,] %>% mutate(diff = cum_sum - lag(cum_sum, default = 0))
        
        leg_y <- do.call(c,(lapply(headData$diff,function(x){seq(max(headData$diff),max(headData$diff)-x+1, by = -1)})))
        
        if(is.null(spaceBtwnCols)){
            leg_x <- rep(seq(1,colz)/2, headData$diff)
        }else{
            leg_x <- tryCatch(
                {
                    rep(cumsum(c(0.5,spaceBtwnCols)), headData$diff)
                }, 
                error = function(e) {
                    return(rep(seq(1,colz)/2, headData$diff))
                    message("Invalid spaceBtwnCols entry. Using default spacing of 0.5")
                }
            )
        }
           
        leg_x <- leg_x[1:(dim(header)[1]+dim(legend)[1])]
        leg_y <- leg_y + compress_y
    
        header$header_x <- leg_x[header$dfNum]
        header$header_y <- leg_y[header$dfNum]
    
        leg_x <- leg_x[-c(header$dfNum)]
        leg_y <- leg_y[-c(header$dfNum)]
        headerBump = F
    }
        
    #legend <- legend %>% arrange(factor(!!as.name(groupLabel), levels = groupBy), !!as.name(sortBy))
    legend <- legend %>% arrange(!!as.name(groupLabel), !!as.name(sortBy))
    
    legend$leg_x <- leg_x
    legend$leg_y <- leg_y
    
    if(!is.null(overrideData)){
        legend <- as.data.frame(overrideData[1])
        header <- as.data.frame(overrideData[2])
    }

    p <- ggplot() + geom_point(data = legend, aes(x = leg_x, y = leg_y),
                    shape=21,
                    size=dotSize,
                    fill=legend[[colorz]],
                    stroke=1,
                    colour="black") +
          geom_text(data = legend, size = 4, mapping = aes(x = leg_x, y = leg_y), label = legend[[clusLabel]], colour = legend[[labCol]]) +
          geom_text(data = legend, size = 6, mapping = aes(x = leg_x+0.05, y = leg_y), label = legend[[legLabel]], hjust = 0) + #NoLegend() + 
          geom_text(data = header, size = headerSize, mapping = aes(x = header_x-0.03, y = header_y), label = header[[groupLabel]], hjust = 0, fontface =2) +
      theme(axis.text = element_blank(), 
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            plot.title = element_blank(),
            title = element_blank(),
            axis.line = element_blank(),
            panel.border = element_blank(),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.margin = unit(c(0, 0, 0, 0), "cm")
    ) + coord_cartesian(ylim = c(ymin, max(ifelse(headerBump == T, leg_y+1, leg_y))*topBuffer), expand = TRUE, clip = "off") + scale_x_continuous(limits = c(0.47,colz*compress_x))

    if(!returnData){
        return(p)
    }else{
        return(list(legend, header))
    }
        
}
    
############ majorDot ############

majorDot <- function(seu.obj = NULL, groupBy = "",
                     yAxis = NULL, scale = T,
                     features = "", split.by = NULL, cols = c("lightgrey", "blue"), cluster.idents = F
                    ){
    
    t <- try(head(seu.obj@assays$RNA@scale.data),silent = T)

    if("try-error" %in% class(t)){
        seu.obj <- ScaleData(seu.obj)
    }
                       
    p <- DotPlot(seu.obj,
                 assay = "RNA",
                 features = features,
                 group.by = groupBy,
                 cols = cols,
                 scale = scale,
                 split.by = split.by#,
                 #idents = levels(seu.obj@meta.data[[groupBy]]),
                 #cluster.idents = cluster.idents
                ) +
      geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
      labs(size='Percent\nexpression')  +
      theme(axis.line = element_blank(),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position = "top",
            legend.direction = "horizontal",
            legend.justification='center',
            panel.background = element_rect(fill = "white",colour = NA),
            plot.background = element_rect(fill = "white",colour = NA),
            legend.key.size = unit(1, "line"),
            panel.border = element_rect(color = "black",
                                        fill = NA,
                                        size = 1),
            ) +
      scale_colour_viridis(option="magma", name='Average\nexpression') +
      guides(size=guide_legend(override.aes = list(shape=21, colour="black", fill="white"),
                               label.position = "bottom")) +
      scale_size(range = c(0.5, 8), limits = c(0, 100)) +
      #annotate("rect", xmin = features_cnt$startVal, xmax = features_cnt$endVal, ymin = features_cnt$cluster-0.5, ymax = features_cnt$cluster+0.5, fill = NA, colour = "mediumpurple1", size = 1) +
      {if(!is.null(yAxis)){scale_y_discrete(limits=rev(yAxis))}} +
      guides(color = guide_colorbar(title = 'Scaled\nExpression')) 
      
    
    return(p)
}
                                           
############ autoDot ############
                                           
autoDot <- function(seu.integrated.obj = NULL, inFile = NULL, groupBy = "",
                     MIN_LOGFOLD_CHANGE = 0.5, MIN_PCT_CELLS_EXPR_GENE = 0.1,
                    filterTerm = "ENSCAFG"
                    ){
    
    if(!is.null(file)){
        all.markers <- read.csv(inFile)
    }else if(!is.null(seu.obj)){
        all.markers = FindAllMarkers(seu.obj,
                                     min.pct = MIN_PCT_CELLS_EXPR_GENE,
                                     logfc.threshold = MIN_LOGFOLD_CHANGE,
                                     only.pos = TRUE)
    }
    
    key.genes <- all.markers[!grepl(filterTerm, row.names(all.markers)),] 
    key.genes.sortedByPval = key.genes[order(key.genes$p_val),]

    features <- key.genes.sortedByPval %>%  group_by(cluster) %>% do(head(., n=5))
    features <- as.data.frame(features[!duplicated(features$gene),])

    #features_cnt <- features %>% count(cluster)
    #features_cnt$n <- rev(features_cnt$n)
    features_cnt <- as.data.frame(table(features$cluster))
    features_cnt$n <- rev(features_cnt$Freq)
    
    features_cnt <- features_cnt %>% mutate(endVal = cumsum(n)+0.5, startVal = endVal-n)

    features_cnt$endVal <- rev(features_cnt$endVal)
    features_cnt$startVal <- rev(features_cnt$startVal)
#     features_cnt$cluster <- as.numeric(features_cnt$cluster)
        features_cnt$Var1 <- as.numeric(features_cnt$Var1)

    
    p <- DotPlot(seu.integrated.obj,
             assay = "RNA",
             features = rev(features$gene),
             #col.min = 0,
             group.by = groupBy,
             scale = TRUE) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  labs(size='Percent\nExpression')  +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.justification='center',
        #legend.spacing.x = unit(1, "cm"),
        legend.key.size = unit(1, "line"),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1),
       # ) +
        axis.ticks.x = element_blank()) +
  scale_colour_viridis(option="magma", name='Average\nExpression') +
  guides(size=guide_legend(override.aes = list(shape=21, colour="black", fill="white"),
                           label.position = "bottom")) +
  scale_size(range = c(0.5, 8), limits = c(0, 100)) +
  annotate("rect", xmin = features_cnt$startVal, xmax = features_cnt$endVal, ymin = features_cnt$Var1-0.5, ymax = features_cnt$Var1+0.5, fill = NA, colour = "mediumpurple1", size = 1) +
  geom_tile(aes(fill = id, x = 0), size = 1, show.legend = FALSE) + 
  geom_tile(aes(fill = id, x = as.numeric(length(unique(features$gene)))+1), size = 1, show.legend = FALSE) + #need to extract data then add colors from colArray to get this corrected
  #scale_fill_manual(values=c("0" = "#CDAD00","1" = "#FFC125", "2" = "#0288D1", "3" = "blue", "4" = "#00CDCD", "5" = "gold", "6" = "brown")) +
  coord_flip() +
  guides(color = guide_colorbar(title = 'Scaled\nExpression')) +
  scale_y_discrete(expand = c(0, 0)) +
  geom_point(aes(x = 0, y=id, size = 100), shape = 21, stroke = 0.75) +
  geom_point(aes(x = as.numeric(length(unique(features$gene)))+1, y=id, size = 100), shape = 21, stroke = 0.75) +
  geom_text(aes(x = 0, label = id), size = 4.5) +
  geom_text(aes(x = as.numeric(length(unique(features$gene)))+1, label = id), size = 4.5)
    
    return(p)
}

############ stackedBar ############
                              #work in progress             
stackedBar <- function(seu.obj = NULL, downSampleBy = "orig.ident", groupBy = "orig.ident", clusters = "clusterID"
                    ){
    
    seu.int.hiQC <- seu.obj
    DefaultAssay(seu.int.hiQC) <- "RNA"
    
    Idents(seu.int.hiQC) <- downSampleBy
    
    set.seed(12)
    seu.int.hiQC.subsampled <- subset(x = seu.int.hiQC, downsample = min(table(seu.int.hiQC@active.ident)))
    
    groupByList <- seu.int.hiQC.subsampled@meta.data[[groupBy]]
    clusterList <- seu.int.hiQC.subsampled@meta.data[[clusters]]

    cluster_freq.table <- as.data.frame(table(groupByList, clusterList)) %>% melt()
    cluster_freq.table <- cluster_freq.table[,-3]
    colnames(cluster_freq.table) <- c("Sample", "ClusterID", "Count")

    cluster_freq.table <- cluster_freq.table %>% dplyr::group_by(ClusterID) %>%
    mutate(pct = round(prop.table(Count),2)) %>% filter(Count !=  0)

    p <- ggplot(cluster_freq.table, aes(x = ClusterID, y = pct, fill = factor(Sample))) +
    geom_bar(stat = "identity", position = "fill", width = 1, colour="white") +
    theme_classic() +
    theme(title = element_text(size= 14),
          legend.title = element_blank(),
          legend.text = element_text(size= 12),
          legend.position = "top",
          legend.direction = "horizontal",
          plot.title = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.key.size = unit(1,"line"),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 16),
          plot.margin = margin(t = 0, r = 21, b = 0, l = 0, unit = "pt")
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_flip() +
    ylab(label = "Fraction of cells") +
    xlab(label = "Cluster ID") + 
    scale_x_discrete(breaks = unique(cluster_freq.table$ClusterID), expand=c(0,0)) +
    guides(fill = guide_legend(nrow = 1))
    
    return(p)

}

############ runMonocle ############
                              #work in progress             
runMonocle <- function(seu.obj = NULL
                    ){
    
    cds <- as.CellDataSet(seu.obj) # make sure non-normlized data are imported
    
    cds <- clusterCells(cds)
    disp_table <- dispersionTable(cds)
    ordering_genes <- subset(disp_table, mean_expression >= 0.1)
    cds <- setOrderingFilter(cds, ordering_genes)
    cds <- reduceDimension(cds)
    cds <- orderCells(cds)
    
    diff_test_res <- differentialGeneTest(cds, fullModelFormulaStr = "~Media")
    ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
    
    cds <- setOrderingFilter(cds, ordering_genes)
    plot_ordering_genes(cds)
    
    cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
    plot_cell_trajectory(cds, color_by = "Hours")
    
    
    
    
    seu.obj.Export <- exportCDS(cds, 'Seurat')
}

### trajectory plot from slingshot
cleanSling <- function(plot = NULL, shape = 21, labCol = "black", size = 8, alpha = 1, rm.na = T, branchData = NULL
                  ) {
    
    pi <- formatUMAP(plot)

    #extract label coords and colors
    g <- ggplot_build(pi)
    
    labCords <- as.data.frame(g$data[2]) #add error if labels are not present
    labCordz <- labCords[,c(1:4,7)]

    colnames(labCordz) <- c("colour", "UMAP1", "UMAP2", "clusterID", "labCol")

    labCordz <- labCordz[order(labCordz$clusterID),]
    labCordz$labCol <- labCol
    
    if(rm.na == T){
        labCordz <- na.omit(labCordz)
    }
    
    #remove old labels
    pi$layers[2] <- NULL

    df.list <- list()
    cnt = 0 
    for (lin in 1:length(branchData)) {
        cnt <- cnt+1
        df <- labCordz[labCordz$clusterID %in% branchData[lin][[1]],]
        #df <- df[order(match(branchData[lin][[1]], labCordz$clusterID)),]
        df <- left_join(data.frame(clusterID = branchData[lin][[1]]),  
                         df,
                         by = "clusterID")
        df$lineage <- names(branchData[lin])
        df <- df %>% mutate(UMAP1_end = lead(UMAP1),
                            UMAP2_end = lead(UMAP2)
                           )

        df.list[[cnt]] <- df
    }
    lineData <- do.call(rbind, df.list)
    
    #add lines and dots to the stripped plot to create final image
    plot <- pi + geom_segment(data = lineData, aes(x = UMAP1,  xend = UMAP1_end,
                                                   y = UMAP2,  yend = UMAP2_end),
                              size = 1) + 
    geom_point(data = lineData, aes(x = UMAP1, y = UMAP2), size = 3) + 
    geom_point(data = lineData[1,], aes(x = UMAP1, y = UMAP2), size = 1.5, colour = "green")

    
    return(plot)
}
                                           
############ createCIBERsort ############
#updated 07.14.23 to calc CPM instead of using log normalizaed data
createCIBERsort <- function(seu.obj = NULL, groupBy = NULL, downSample = F, outDir = "", outName = "", normMethod = "log"
                    ){

    Idents(seu.obj) <- groupBy

    if(downSample){
        #randomly downsample the subset data
        set.seed(12)
        seu.obj <- subset(x = seu.obj, downsample = min(table(seu.obj@meta.data[[groupBy]])))
        if(normMethod == "log"){
            seu.obj <- NormalizeData(seu.obj)
        }
        
    }
    
    #extract required data for pseudobulk conversion
    if(normMethod == "log"){
        mat <- seu.obj@assays$RNA@data
    }else{
        mat <- seu.obj@assays$RNA@counts
    }
    
    bioRep <- as.data.frame(seu.obj@meta.data[[groupBy]])
    colnames(bioRep) <- "cellSource"
    row.names(bioRep) <- colnames(seu.obj)
    bioRep$cellSource <- as.factor(bioRep$cellSource)
    
    #use custom function to convert to pseudobulk
    pbj <- getPb(mat, bioRep)

    if(normMethod == "cpm"){
        pbj <- t(t(pbj)/colSums(pbj))*1e6
    }
    
    outfile <- paste(outDir,outName,"_ciberSort_matrix.csv", sep = "")
    write.csv(pbj, file = outfile,quote=T)
    
}

############ sankeyPlot ############

sankeyPlot <- function(seu_obj = NULL, new.ident = NULL, old.ident = "clusterID", old.colorz = NULL,
                       new.colorz = NULL, old.labCol = NULL, new.labCol = NULL, flowCol = "grey"
                    ){
    
    Idents(seu_obj) <- new.ident

    #get node data
    new <- levels(seu_obj@active.ident)
    new <- paste("c", new, sep="")
    nodeNum <- length(unique(seu.obj@meta.data[[old.ident]])) + length(levels(seu_obj@active.ident)) - 1
    nodes <- data.frame(node = c(0:nodeNum), 
                        name = c(as.character(sort(unique(seu.obj@meta.data[[old.ident]]))), new))

    #exctra data for plotting
    seu_obj_data <- as.data.frame(seu_obj@meta.data)
    seu_obj_data$barcode <- rownames(seu_obj_data)
    seu_obj_data <- seu_obj_data[c("barcode", old.ident,new.ident)] 
    colnames(seu_obj_data) <- c("barcode", "Initial","SubCluster")
    seu_obj_data[c("Initial","SubCluster")] <- lapply(seu_obj_data[c("Initial","SubCluster")], as.character)

    #use ggsankey package to get data in correct format
    data <- seu_obj_data %>% make_long(Initial, SubCluster)

    #prefix sub clusters with "S"
    data$next_node <- ifelse(!is.na(data$next_node),paste("c",data$next_node,sep=""),NA)
    data <- data %>% mutate(node = ifelse(x == "SubCluster",paste("c",node,sep=""),node))

    #order the groups so they are colored appropriately
    data$node <- factor(data$node, levels = nodes$name)

    #make the plot
    p <- ggplot(data, aes(x = x,
                          next_x = next_x,
                          node = node,
                          next_node = next_node,
                          fill = factor(node),
                          label = node)) + geom_sankey() +
    geom_sankey(flow.alpha = 0.5, node.fill = c(old.colorz,new.colorz), flow.fill = flowCol) + 
    geom_sankey_label(size = 3.5, color = 1, fill = "white") +
    theme_sankey(base_size = 16) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = 1),
          plot.margin = margin(t = 7, r = 14, b = 7, l = 7, unit = "pt")) + scale_x_discrete(expand = c(0, 0))
    
    return(p)
}

############ dotPlotBY_TYPE ############

dotPlotBY_TYPE <- function(seu_obj = NULL, pwdTOvilnCSVoutput = NULL, groupBy = NULL, namedCols = NULL, database = "clfamiliaris_gene_ensembl", exlcude = "", boxColor = "black"
                          ){
    
    df <- read.csv(pwdTOvilnCSVoutput)
    
    ensembl <-  biomaRt::useMart("ensembl", dataset = database)
    
    cl_go_anno <- biomaRt::getBM(
        attributes = c("external_gene_name", "description", "go_id", "name_1006", "namespace_1003", "definition_1006"),
        filters = c("external_gene_name"), 
        values = unique(as.character(df$gene)),
        mart = ensembl) %>%
    dplyr::filter(namespace_1003 %in% c("cellular_component", "molecular_function"))
    
    surface <- sort(unique(cl_go_anno$external_gene_name[grep("component of membrane", cl_go_anno$name_1006)]))
    surface.neg <- sort(unique(cl_go_anno$external_gene_name[grep("mitochondrial|integral component of Golgi membrane|intracellular membrane-bounded organelle|nuclear membrane", cl_go_anno$name_1006)]))
    surface <- surface[!(surface %in% surface.neg)]
    
    secreted <- sort(unique(cl_go_anno$external_gene_name[
        c(grep("extracellular space", cl_go_anno$name_1006),
          grep("granzyme", cl_go_anno$description))]))
    
    tf <- sort(unique(cl_go_anno$external_gene_name[
        grep("transcription factor activity|transcription activator activity|transcription repressor activity|transcription coactivator activity|DNA binding|DNA-binding",cl_go_anno$name_1006)]))
    
    intersect(surface, secreted)
    intersect(surface, tf)
    intersect(secreted, tf)
    
    surface <- surface[!(surface %in% tf)]
    secreted <- secreted[(!(secreted %in% surface)) & (!(secreted %in% tf))]
    neg <- sort(unique(df$gene[!(df$gene %in% c(surface, secreted, tf))]))
    
    anno.df <- data.frame(
        gene = c(surface, secreted, tf, neg),
        anno = c(rep("cell surface", length(surface)),
                 rep("secreted", length(secreted)),
                 rep("transcription factor", length(tf)),
                 rep("", length(neg)))) %>%
    dplyr::filter(!(gene %in% exlcude))
    
    df.plot <- df %>% dplyr::left_join(anno.df, by = c("gene")) %>%
    dplyr::select(gene, cluster, anno, avg_log2FC) 
    
    df.plot$anno <- factor(df.plot$anno, levels = c("cell surface", "secreted", "transcription factor", ""))

    if(is.null(groupBy)){
        Idents(seu.obj) <- "clusterID_sub"
    }
    
    p <- mapply(function(i, j){
        x <- df.plot %>%
        distinct(gene,.keep_all= TRUE) %>%
        dplyr::filter(anno == i) %>%
        dplyr::group_by(cluster) %>%
        dplyr::top_n(5, avg_log2FC) %>%
        #dplyr::group_by(cluster) %>%
        dplyr::arrange(gene, .by_group = T) %>% 
        dplyr::group_by(cluster) %>%
        mutate(n = n())
        
        features_cnt <- x[!duplicated(x$cluster),]#  x %>% unique(cluster)
        features_cnt$n <- rev(features_cnt$n)
        
        features_cnt <- as.data.frame(features_cnt) %>% mutate(endVal = cumsum(n)+0.5, startVal = endVal-n)
        
        features_cnt$endVal <- rev(features_cnt$endVal)
        features_cnt$startVal <- rev(features_cnt$startVal)
        features_cnt$cluster <- as.numeric(features_cnt$cluster)+1 
        
        p <- DotPlot(seu.obj, features = rev(x$gene), group.by = groupBy) +
        geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
        labs(size='Percent\nExpression', title = j)  +
        theme(axis.text.x = element_blank(),
              axis.title = element_blank(),
              axis.line = element_blank(),
              legend.position = "top",
              legend.direction = "horizontal",
              legend.justification='center',
              legend.key.size = unit(1, "line"),
              panel.border = element_rect(color = "black",
                                          fill = NA,
                                          size = 1),
              axis.ticks.x = element_blank()) +
        scale_colour_continuous_diverging(palette = 'Blue-Red', mid = 0) + 
        guides(size=guide_legend(override.aes = list(shape=21, 
                                                     colour="black", fill="white"),
                                 label.position = "bottom")) +
        scale_size(range = c(0.5, 8), limits = c(0, 100)) +
        annotate("rect", xmin = features_cnt$startVal, xmax = features_cnt$endVal, ymin = features_cnt$cluster-0.5, ymax = features_cnt$cluster+0.5, fill = NA, colour = boxColor, size = 1) +
        geom_tile(aes(fill = id, x = 0), size = 1, show.legend = FALSE) + 
        geom_tile(aes(fill = id, x = as.numeric(length(unique(x$gene)))+1), size = 1, show.legend = FALSE) + 
        scale_fill_manual(values=namedCols) + 
        coord_flip() +
        guides(color = guide_colorbar(title = 'Scaled\nExpression')) +
        scale_y_discrete(expand = c(0, 0)) +
        geom_text(aes(x = 0, label = id), size = 4.5) + #hardcoded
        geom_text(aes(x = as.numeric(length(unique(x$gene)))+1, label = id), size = 4.5)
        return(p)
    }, i = c("cell surface", "secreted", "transcription factor", ""),
                j = c("Cell\nsurface","Secreted","Transcription\nfactor","Remaining\ngenes"),
                SIMPLIFY = F,
                USE.NAMES = F)
    
    p[[1]] <- p[[1]] + NoLegend()
    p[[2]] <- p[[2]] + theme(axis.title.y = element_blank())
    p[[3]] <- p[[3]] + theme(axis.title.y = element_blank()) + NoLegend()
    p[[4]] <- p[[4]] + theme(axis.title.y = element_blank()) + NoLegend()
    
    finalPlot <- plot_grid(plotlist = p, nrow = 1, labels = "none", label_size = 8, rel_widths = c(1,1,1,1)) 
    
    return(finalPlot)
}
                                           
############ prettyViln ############

prettyViln <- function(plot = NULL, colorData = NULL, nrow = 2, ncol = NULL){
    
    if(is.null(ncol)){
        ncol = ceiling(sqrt(length(plot)))
    }
    
    
    p <- lapply(plot, function(x) x + theme(axis.title = element_blank(),
                                          axis.text = element_blank(),
                                          axis.text.x = element_text(angle = 0, vjust = 0, hjust=0.5, size =5 ),
                                          #axis.ticks.x = element_blank(),
                                          axis.ticks = element_blank(),
                                          title = element_text(size = 16),
                                          legend.position = "none",
                                          plot.margin = margin(5.5, 0, 0, 0, "pt")
                                         ) + coord_cartesian(expand = TRUE)
                )
    
    if(!is.null(colorData)){
        g <- ggplot(colArray, aes(x = clusterID_sub, y = 1, fill = clusterID_sub, label = clusterID_sub)) + 
        geom_tile(fill = "transparent") +
        geom_point(shape = 21, stroke = 0.75, size = 11) +
        geom_text(fontface = "bold", size = 6, color = colArray$labCol) + theme_bw(base_size = 12) +
        scale_fill_manual(values = colArray$colour) + scale_y_discrete(expand = c(0, 0)) +
        theme(legend.position = "none", panel.spacing = unit(0, "lines"),
              panel.background = element_blank(), 
              panel.border = element_blank(),
              plot.background = element_blank(), 
              plot.margin = margin(0, 0, 0, 0, "pt"),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank()
             )
    }
        
        patch <- area()
        counter=0
        for (i in 1:nrow) {
            for (x in 1:ncol) {
                counter = counter+1
                if (counter <= ncol*nrow) {
                    patch <- append(patch, area(t = i, l = x, b = i, r = x))
                }
            }
        }
        
    
        pi <- Reduce( `+`, p ) + plot_layout(design = patch, 
                                                             widths = c(rep.int(1, ncol)), 
                                                             heights = c(rep.int(1, nrow))
                                                            ) 
        return(pi)
}

############ prettyVolc ############                                          
prettyVolc <- function(plot = NULL, rightLab = NULL, leftLab = NULL, rightCol = "red", leftCol = "blue", arrowz = T
                    ){
    
    p <- plot + scale_x_symmetric(mid = 0) + theme(legend.position = c(0.10, 0.9),
                                                      legend.background = element_blank(),
                                                      legend.key = element_blank(),
                                                      axis.title = element_text(size = 16),
                                                   axis.text = element_text(size = 12),
                                                      panel.grid.major = element_blank(),
                                                      panel.grid.minor = element_blank(),
                                                      panel.border = element_blank(),
                                                      panel.background = element_blank(),
                                                      axis.line = element_line(color="black"),
                                                      plot.title = element_blank()
                                                     ) + 
    {if(arrowz){
        annotate("segment", x = 0.58*1.5, 
                 y = ggplot_build(plot)$layout$panel_scales_y[[1]]$range$range[2]*1.06, 
                 xend = c(max(abs(plot$data$log2FoldChange)),-max(abs(plot$data$log2FoldChange)))[1], 
                 yend = ggplot_build(plot)$layout$panel_scales_y[[1]]$range$range[2]*1.06, 
                 lineend = "round", linejoin = "bevel", linetype ="solid", colour = rightCol,
                 size = 1, arrow = arrow(length = unit(0.1, "inches"))
                ) 
        }} +
    {if(arrowz){
        annotate("segment", x = -0.58*1.5, 
                 y = ggplot_build(plot)$layout$panel_scales_y[[1]]$range$range[2]*1.06, 
                 xend = c(max(abs(plot$data$log2FoldChange)),-max(abs(plot$data$log2FoldChange)))[2],
                 yend = ggplot_build(plot)$layout$panel_scales_y[[1]]$range$range[2]*1.06, 
                 lineend = "round", linejoin = "bevel", linetype ="solid", colour = leftCol,
                 size = 1, arrow = arrow(length = unit(0.1, "inches"))
                )
        }} + 
    {if(!is.null(rightLab)){
        annotate(geom = "text", x = (max(abs(plot$data$log2FoldChange))-0.58*1.5)/2+0.58*1.5, 
                 y = ggplot_build(plot)$layout$panel_scales_y[[1]]$range$range[2]*1.09,
                 label = rightLab,
                 hjust = 0.5,
                 size = 5)
        }} + 
    {if(!is.null(leftLab)){
        annotate(geom = "text", x = -(max(abs(plot$data$log2FoldChange))-0.58*1.5)/2-0.58*1.5, 
                 y = ggplot_build(plot)$layout$panel_scales_y[[1]]$range$range[2]*1.09,
                 label = leftLab,
                 hjust = 0.5,
                 size = 5)
    }} 
    
    return(p)
}


############ plotGSEA ############
#TO DO:L fix problem with trimTerm

plotGSEA <- function(pwdTOgeneList = NULL, geneList = NULL, geneListDwn = NULL,
                     category = "C5", species = "dog", upCol = "red", dwnCol = "blue", saveRes = NULL,
                     pvalueCutoff = 0.05, subcategory = NULL, termsTOplot = 8, upOnly = F, trimTerm = T, size = 4
                    ){
    if(!is.null(pwdTOgeneList)){
        geneLists <- read.csv(pwdTOgeneList)
        
        geneListUp <- geneLists %>% arrange(padj) %>% filter(log2FoldChange > 0) %>% .$gene
        geneListDwn <- geneLists %>% arrange(padj) %>% filter(log2FoldChange < 0) %>% .$gene
        
    }else{
        geneListUp <- geneList
    }
    

    
    can_gene_sets <- as.data.frame(msigdbr(species = species, category = category, subcategory = subcategory))
    msigdbr_list <- split(x = can_gene_sets$gene_symbol, f = can_gene_sets$gs_name)
    datas <- can_gene_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
    
    enriched_up <- as.data.frame(enricher(gene = geneListUp, TERM2GENE = datas, pvalueCutoff = pvalueCutoff)
                                ) %>% arrange(p.adjust) %>% mutate(x_axis = -log10(p.adjust),
                                                                   direction = "UP")
    if(!upOnly){
        enriched_dwn <- as.data.frame(enricher(gene = geneListDwn, 
                                               TERM2GENE = datas, pvalueCutoff = pvalueCutoff)
                                     ) %>% arrange(p.adjust) %>% mutate(x_axis = log10(p.adjust),
                                                                        direction = "DOWN")
    
        enriched <- rbind(enriched_up[1:ifelse(nrow(enriched_up) > termsTOplot, termsTOplot, length(enriched_up)),], 
                          enriched_dwn[1:ifelse(nrow(enriched_dwn) > termsTOplot, termsTOplot, length(enriched_up)),])
        
    } else {
        enriched <- enriched_up[1:ifelse(nrow(enriched_up) > termsTOplot, termsTOplot, nrow(enriched_up)),]

    }
    
    if(!is.null(saveRes)){
        if(!exists("enriched_dwn")){
            write.csv(enriched_up, file = saveRes)
        } else {
             write.csv(rbind(enriched_up,enriched_dwn), file = saveRes)
        }
       
    }

    enriched <- enriched %>% arrange(desc(x_axis))
    
    if(trimTerm){
        enriched$Description <- gsub("REACTOME_", "", enriched$Description)
        enriched$Description <- gsub("_AND_", "_&_", enriched$Description)
        enriched$Description <- gsub("OTHER_", "", enriched$Description)
        enriched$Description <- gsub("GOBP_", "", enriched$Description)
        enriched$Description <- gsub("GOCC_", "", enriched$Description)
        enriched$Description <- gsub("GOMF_", "", enriched$Description)
        enriched$Description <- gsub("HP_", "", enriched$Description)
        enriched$Description <- gsub("POSITIVE_", "POS_", enriched$Description)
        enriched$Description <- gsub("ANTIGEN_PROCESSING_&_", "", enriched$Description)
        enriched$Description <- gsub("ADAPTIVE_IMMUNE_RESPONSE_BASED_ON_SOMATIC_RECOMBINATION_OF_IMMUNE_RECEPTORS_BUILT_FROM_IMMUNOGLOBULIN_SUPERFAMILY_DOMAINS", "ADAPTIVE_IMMUNE_RESPONSE_BASED_ON_SOMATIC_RECOMBINATION", enriched$Description)
        enriched$Description <- gsub("PROTON_TRANSPORTING_ATP_SYNTHASE_ACTIVITY_ROTATIONAL_MECHANISM", "PROTON_TRANSPORTING_ATP_SYNTHASE_ACTIVITY", enriched$Description)
        
#         enriched$Description <- unlist(lapply(enriched$Description, str_trunc, width = 50, side = "right", ellipsis = "*")) ### TO DO: this causes problems if trunc results in the same string -- need to ensure that is handled before adding this back in
    }
    
    #remove NAs and then plot
    enriched <- na.omit(enriched)
    p <- ggplot(data=enriched, aes(x=x_axis, y=Description, fill = direction)) +
    geom_bar(stat="identity") +  theme_classic() + scale_y_discrete(limits=rev(enriched$Description)
                                                                   ) + 
    geom_text(aes(x = ifelse(x_axis > 0, -0.05,0.05), label = Description), hjust = ifelse(enriched$x_axis > 0, 1,0), size = size) + 
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          legend.justification = "top"
         ) + coord_cartesian(clip = "off") + scale_x_symmetric(mid = 0, name = "Signed log10(padj)") + NoLegend() + scale_fill_manual(values = c(dwnCol,upCol))
    
    return(p)

    
}

############ crossSpeciesDEG ############
###TO DO: fix intercettion line
crossSpeciesDEG <- function(pwdTOspecies1 = NULL, pwdTOspecies2 = NULL, species1 = "Canine", species2 = "Human", 
                            nlab = 10, nlab_conflict = 5, nlab_axis = 2, overlapGenes = overlapGenes, colUp = "red",
                            colDwn = "blue", cONtrast = c("celltype_1","celltype_2"), seed = 12,
                            hjustvar = c(1,0.5,0,0.5,1,0,1,0),
                            vjustvar = c(0.5,1,0.5,0,0,1,1,0),
                            outDir = "./output/", saveGeneList = F
                    ){

    geneLists.spec1 <- read.csv(pwdTOspecies1)
    geneLists.spec2 <- read.csv(pwdTOspecies2)
    
#     overlapGenes <- intersect(rownames(seu.obj.k9),rownames(seu.obj.hu))
    
    geneLists.spec1 <- geneLists.spec1 %>% mutate(species = species1,
                                                  signed_pVal_spec1 = ifelse(log2FoldChange > 0, -log10(padj), log10(padj))
                                                 ) 
#     %>% filter(gene %in% overlapGenes) %>% top_n(n = ds, wt = signed_pVal_hu)
    
    geneLists.spec2 <- geneLists.spec2 %>% mutate(species = species2,
                                                  signed_pVal_spec2 = ifelse(log2FoldChange > 0, -log10(padj), log10(padj))
                                                 ) 
#     %>% filter(gene %in% overlapGenes) %>% top_n(n = ds, wt = signed_pVal_can)
    
    geneList <- unique(c(geneLists.spec1$gene, geneLists.spec2$gene)) %>% as.data.frame()
    colnames(geneList) <- "gene"
    
    #add intersect to only plot overlappling genes btwn species -- done above?
    
    
    geneList <- geneList %>% left_join(geneLists.spec1[,c("gene","species","signed_pVal_spec1")], 
                                       by = "gene") %>% left_join(geneLists.spec2[,c("gene","signed_pVal_spec2")], 
                                                                  by = "gene") %>% replace(is.na(.), 0) %>% 
    mutate(species = ifelse(species == "0", species2,species1),
           species = ifelse(signed_pVal_spec2 != "0" & signed_pVal_spec1 != "0", "Both",species1),
        sig = ifelse(signed_pVal_spec2*signed_pVal_spec1 == 0,signed_pVal_spec2+signed_pVal_spec1,signed_pVal_spec2*signed_pVal_spec1) ,
       direction = case_when(signed_pVal_spec2*signed_pVal_spec1 != 0 & signed_pVal_spec1 > 0 & signed_pVal_spec2 > 0 ~ "up",
                             signed_pVal_spec2*signed_pVal_spec1 != 0 & signed_pVal_spec1 < 0 &  signed_pVal_spec2 < 0 ~ "dwn",
                             signed_pVal_spec2*signed_pVal_spec1 == 0 & signed_pVal_spec1 > 0 ~ "axis1",
                             signed_pVal_spec2*signed_pVal_spec1 == 0 & signed_pVal_spec2 > 0 ~ "axis2",
                             signed_pVal_spec2*signed_pVal_spec1 == 0 & signed_pVal_spec1 < 0 ~ "axis3",
                             signed_pVal_spec2*signed_pVal_spec1 == 0 & signed_pVal_spec2 < 0 ~ "axis4",                             
                             signed_pVal_spec2*signed_pVal_spec1 != 0 & signed_pVal_spec1 > 0 & signed_pVal_spec2 < 0 ~ "conflict1",
                             signed_pVal_spec2*signed_pVal_spec1 != 0 & signed_pVal_spec1 < 0 & signed_pVal_spec2 > 0 ~ "conflict2",
                            ),
       
       label = ifelse(sig != 0, gene, NA)
      )  %>% arrange(-abs(sig)) %>% group_by(direction) %>% #arrange(p_val_adj) %>% 
            mutate(lab_col=case_when(direction == "up" ~ colUp,
                                     direction == "dwn" ~ colDwn,
                                     direction == "axis1" ~ "black",
                                     direction == "axis2" ~ "black",
                                     direction == "axis3" ~ "black",
                                     direction == "axis4" ~ "black",
                                     direction == "conflict1" ~"hotpink",
                                     direction == "conflict2" ~"hotpink"),
                   lab=ifelse(row_number() <= nlab & grepl("up|dwn", direction) , gene,
                              ifelse(grepl("axis", direction) & row_number() <= nlab_axis, gene, 
                                     ifelse(grepl("conflict", direction) & row_number() <= nlab_conflict, gene, NA))
                             )
                  ) 
    
    anno.df <- as.data.frame(list(
        direction = c("axis1","axis2","axis3","axis4","conflict1","conflict2","up","dwn"),
        xpos = c(Inf,0,-Inf,0,Inf,-Inf,Inf,-Inf),
        ypos =  c(0,Inf,0,-Inf,-Inf,Inf,Inf,-Inf),
        hjustvar = hjustvar,
        vjustvar = vjustvar,
        colz = c("black","black","black","black","hotpink","hotpink",colUp,colDwn)
    ))

    
    annotationz <- geneList %>% group_by(direction) %>% summarize(cntz = n()) %>% left_join(anno.df, by = "direction")

    if(saveGeneList){
        write.csv(geneList, file = paste0(outDir,species1,"_v_", species2, "_", cONtrast[1],"_",cONtrast[2],".csv"), row.names = F)
    }
    
    p <- ggplot(geneList, aes(x=signed_pVal_spec1,y=signed_pVal_spec2)) + geom_hline(yintercept=0,linetype=2) +
    geom_vline(xintercept=0,linetype=2) + 
    geom_point() + 
    geom_label_repel(max.overlaps = Inf, size=2, label = geneList$lab, color = geneList$lab_col, show.legend = F, seed = seed) + 
    theme_classic() + scale_x_symmetric(mid = 0) + scale_y_symmetric(mid = 0) + 
    geom_label(data = annotationz, aes(x = xpos,y = ypos, hjust = hjustvar, vjust = vjustvar, label = cntz), color = annotationz$colz) +
    labs(x = paste0(species1, " " ,cONtrast[1]," vs ",cONtrast[2]," DE signed log10(p.adj)"), 
                    y = paste0(species2, " ",cONtrast[1]," vs ",cONtrast[2]," DE signed log10(p.adj)"), title = paste0(species1, " and ", species2," " ,cONtrast[1]," vs ",cONtrast[2])) + coord_cartesian(clip = 'off') + 
    theme(plot.title = element_text(size = 20))

    return(p)
    
}
    

############ skewPlot ############
skewPlot <- function(seu.obj = seu.obj
                    ){
    
    pct.df <- table(seu.obj$majorID_sub, seu.obj$name) %>% melt() %>% group_by(Var.2) %>% mutate(samN = sum(value))
    pct.df$Var.1 <- as.factor(pct.df$Var.1)
    
    pct.df$pct <- pct.df$value/pct.df$samN*100
    pct.df$cellSource <- ifelse(grepl("tils",pct.df$Var.2),"TILs","Blood")

    statz <- compare_means(pct ~ cellSource, group.by = "Var.1", pct.df)

    status <- pct.df %>% group_by(Var.1,cellSource) %>% summarize(med = median(pct)) %>% spread(cellSource,med) %>% left_join(statz[c("Var.1","p.adj")], by = "Var.1") %>% mutate(lfc = abs(log2(TILs/Blood)),
    lab_col=case_when(lfc >= 3 & p.adj < 0.05 ~ "Unique",
                      p.adj < 0.05 & lfc < 3 ~ "Skewed",
                      p.adj > 0.05 ~ "Common")) %>% select(Var.1,lab_col) %>% as.data.frame()

    pct.df <- pct.df %>% left_join(status, by = "Var.1")

    pct.df$lab_col <- as.factor(pct.df$lab_col)
    
    p <- ggplot(pct.df, aes(x = Var.1, y = pct, fill = lab_col)) +
    stat_summary(fun = mean, geom = "bar", width = 0.7) +
    stat_summary(fun.data = mean_se, geom = "errorbar",width = 0) + 
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0))) + 
    labs(x = "Cell type", y = "Percent total") +
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
  ) + facet_wrap("cellSource", nrow = 2) + scale_fill_manual(labels = c("Common","Skewed","Unique"),
                                                           values = brewer.pal(n = 3, name = "Dark2"),
                                                           name = "Uniqueness\nclassification")
    return(p)
}
    
    
############ ExportToCB_cus ############

#need to update to get barcode in first column
ExportToCB_cus <- function(seu.obj = seu.obj, dataset.name = "", outDir = "./output/", markers = NULL, reduction = "umap", test = F, skipEXPR=F, colsTOkeep=NULL,
                           feats = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                                     "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                                     "CD4", "MS4A1", "PPBP","HBM")
                          
                          ){
    
    outDir <- paste0(outDir,dataset.name, "/")
    dir.create(outDir)
    
    if(test){seu.obj <- subset(x = seu.obj, downsample = 500)}

    
    meta <- seu.obj@meta.data
    meta <- rownames_to_column(meta, "barcode")
    if(!is.null(colsTOkeep)){
        meta <- meta[ ,colnames(meta) %in% c("barcode",colsTOkeep)]
    }
    write.table(meta, paste0(outDir,"meta.tsv"), quote=FALSE, sep='\t', row.names = F)
    
    quickGenes <- as.data.frame(feats[feats %in% rownames(seu.obj)])
    colnames(quickGenes) <- "symbol"
    write.csv(quickGenes, paste0(outDir,"quickGenes.csv"), quote=FALSE, row.names = F)
    
    markers <- read.csv(markers)
    markers <- markers[c("cluster","gene","avg_log2FC","p_val_adj")]
    #colnames(markers)[3,4] <- c("avg_diff","p_val")
    write.table(markers,paste0(outDir,"markers.tsv"), quote=FALSE, sep='\t', row.names = F)
    
    if(!skipEXPR){
        data.df <- as.data.frame(seu.obj@assays$RNA@data)
        data.df <- rownames_to_column(data.df, "symbol")
        write.table(data.df, paste0(outDir,"exprMatrix.tsv"), quote=FALSE, sep='\t', row.names = F)
        R.utils::gzip(paste0(outDir,"exprMatrix.tsv"), overwrite = TRUE)
    }
    
    data.df <- as.data.frame(seu.obj[[reduction]]@cell.embeddings)
    data.df <- rownames_to_column(data.df, "barcode")
    write.table(data.df,paste0(outDir,reduction,".coords.tsv"), quote=FALSE, sep='\t', row.names = F)
    
}



loadKal <- function(din = NULL, dout = NULL, seu.obj = NULL){

    fpath <-  paste0("./", din,"/") 

    files <- list.files(path = fpath, pattern=NULL, all.files=FALSE,
                        full.names=F)
    
    for (infile in files) {
        pwd <- paste0("./", din,"/", infile, "/adata.h5ad") 

        seu.kal <- MuDataSeurat::ReadH5AD(pwd)

        # gene symbols are stashed in meta.features slot, so extract them
        dups <- as.character(seu.kal@assays$RNA@meta.features[duplicated(seu.kal@assays$RNA@meta.features$gene_name), ])
        if(length(dups) > 0){
            message(paste0("Duplicate gene symbols found in the seu.kal@assays$RNA@meta.features slot. A total of ", length(dups), " duplicated gene symbols were identified and the gene ID will be retained.\n"))
            message(paste0("The impacted features are: ", paste(dups, collapse = ", ")))

            newGenez <- seu.kal@assays$RNA@meta.features %>% rownames_to_column("gene") %>% mutate(newGenez = ifelse(duplicated(seu.kal@assays$RNA@meta.features$gene_name), as.character(gene), as.character(gene_name))) %>% pull(newGenez)

        }

        message(paste0("Subseting the annotated dataset for transfer.", 
                       "Searching for ",infile," in ", paste(as.character(unique(seu.obj$orig.ident)) , collapse = ", "))
               )
        seu.obj.sub <- subset(seu.obj, cells = colnames(seu.obj)[grepl(infile, seu.obj$orig.ident)])
        goodCodez <- substr(colnames(seu.obj.sub), 1, 16) #barcodes must be nchar == 16

        message("Subseting the kallisto pseudo-aligned object to remove empty droplets and other low quality data points.")
        seu.kal <- subset(seu.kal, subset = barcode %in% goodCodez)

        cnt.mat <- as.matrix(seu.kal@assays$RNA@counts)
        colnames(cnt.mat) <- goodCodez
        rownames(cnt.mat) <- newGenez

        meta.df.orig <- seu.obj.sub@meta.data
        meta.df.orig$barcode <- goodCodez

        #test that the order matches up
        table(substr(rownames(meta.df.orig), 1, 16) == substr(colnames(seu.obj.sub), 1, 16))
        # TRUE 
        # 4745 

        meta.df <- seu.kal@meta.data %>% left_join(meta.df.orig, by = "barcode")
        rownames(meta.df) <- meta.df$barcode

        # Generate new Seurat object.
        seu.obj.kal <- CreateSeuratObject(cnt.mat,
                                      project = infile,
                                      assay = "RNA",
                                      min.cells = 0,
                                      min.features = 0,
                                      names.field = 1,
                                      names.delim = "_",
                                      meta.data = meta.df
                                     )


        # reduce the data
        seu.obj.kal <- NormalizeData(seu.obj.kal,
                                 normalization.method = "LogNormalize",
                                 Scale.factor = 10000) ###change method to scran???

        seu.obj.kal <- FindVariableFeatures(seu.obj.kal,
                                        selection.method = "vst", 
                                        nfeatures = 2500) #can change number of feats used

        all.genes <- rownames(seu.obj.kal)
        seu.obj.kal <- ScaleData(seu.obj.kal, features = all.genes)
        seu.obj.kal <- RunPCA(seu.obj.kal, features = VariableFeatures(object = seu.obj.kal))
        seu.obj.kal <- FindNeighbors(seu.obj.kal,
                                 dims = 1:10
                                ) #can change dims
        seu.obj.kal <- FindClusters(seu.obj.kal,
                                resolution = 0.1
                               ) #can change resolution
        seu.obj.kal <- RunUMAP(seu.obj.kal, 
                           dims = 1:15
                          ) #can change dims


        featPlots = c("PTPRC", "CD3E", "CD8A", "GZMA", 
                      "IL7R", "ANPEP", "FLT3", "DLA-DRA", 
                      "CD4", "MS4A1", "TRDC","FOXP3")

        ### Fig : Create UMAP by majorID_sub
        pi <- DimPlot(seu.obj.kal, 
                      reduction = "umap", 
                      group.by = "celltype.l3",
                      pt.size = 0.25,
                      label = T,
                      label.box = T,
                      shuffle = TRUE
        )

        ggsave(paste0("./",dout,"/",infile, "_rawUMAP.png"), width = 7, height = 7)

        saveRDS(seu.obj.kal, file = paste0("./",dout,"/", infile,"_kallisto_S1.rds"))
    }

}





netVisual_aggregate_mod <- function(object, signaling, signaling.name = NULL, color.use = NULL, thresh = 0.05, vertex.receiver = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL, top = 1, remove.isolate = FALSE,
                                vertex.weight = 1, vertex.weight.max = NULL, vertex.size.max = NULL,
                                weight.scale = TRUE, edge.weight.max = NULL, edge.width.max=8,
                                layout = c("circle","hierarchy","chord","spatial"),
                                pt.title = 12, title.space = 6, vertex.label.cex = 0.8,
                                alpha.image = 0.15, point.size = 1.5,
                                group = NULL,cell.order = NULL,small.gap = 1, big.gap = 10, scale = FALSE, reduce = -1, show.legend = FALSE, legend.pos.x = 20,legend.pos.y = 20,
                                ...) {
  layout <- match.arg(layout)
  if (is.null(vertex.weight)) {
    vertex.weight <- as.numeric(table(object@idents))
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    } else {
      vertex.size.max <- 15
    }
  }
  pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = T)

  if (is.null(signaling.name)) {
    signaling.name <- signaling
  }
  net <- object@net

  pairLR.use.name <- dimnames(net$prob)[[3]]
  pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
  pairLR <- pairLR[pairLR.name, ]
  prob <- net$prob
  pval <- net$pval

  prob[pval > thresh] <- 0
  if (length(pairLR.name) > 1) {
    pairLR.name.use <- pairLR.name[apply(prob[,,pairLR.name], 3, sum) != 0]
  } else {
    pairLR.name.use <- pairLR.name[sum(prob[,,pairLR.name]) != 0]
  }


  if (length(pairLR.name.use) == 0) {
    stop(paste0('There is no significant communication of ', signaling.name))
  } else {
    pairLR <- pairLR[pairLR.name.use,]
  }
  nRow <- length(pairLR.name.use)

  prob <- prob[,,pairLR.name.use]
  pval <- pval[,,pairLR.name.use]

  if (length(dim(prob)) == 2) {
    prob <- replicate(1, prob, simplify="array")
    pval <- replicate(1, pval, simplify="array")
  }
  # prob <-(prob-min(prob))/(max(prob)-min(prob))

  if (layout == "hierarchy") {
    prob.sum <- apply(prob, c(1,2), sum)
    # prob.sum <-(prob.sum-min(prob.sum))/(max(prob.sum)-min(prob.sum))
    if (is.null(edge.weight.max)) {
      edge.weight.max = max(prob.sum)
    }
    par(mfrow=c(1,2), ps = pt.title)
    netVisual_hierarchy1(prob.sum, vertex.receiver = vertex.receiver, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max, title.name = NULL, vertex.label.cex = vertex.label.cex,...)
    netVisual_hierarchy2(prob.sum, vertex.receiver = setdiff(1:nrow(prob.sum),vertex.receiver), sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max, title.name = NULL, vertex.label.cex = vertex.label.cex,...)
    graphics::mtext(paste0(signaling.name, " signaling pathway network"), side = 3, outer = TRUE, cex = 1, line = -title.space)
    # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
    # grid.echo()
    # gg <-  grid.grab()
    gg <- recordPlot()
  } else if (layout == "circle") {
    prob.sum <- apply(prob, c(1,2), sum)
    # prob.sum <-(prob.sum-min(prob.sum))/(max(prob.sum)-min(prob.sum))
    gg <- netVisual_circle(prob.sum, sources.use = sources.use, targets.use = targets.use, idents.use = idents.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max,title.name = paste0(signaling.name, " signaling pathway network"), vertex.label.cex = vertex.label.cex,...)
  }  else if (layout == "spatial") {
    prob.sum <- apply(prob, c(1,2), sum)
    if (vertex.weight == "incoming"){
      if (length(slot(object, "netP")$centr) == 0) {
        stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
      }
      vertex.weight = object@netP$centr[[signaling]]$indeg
    }
    if (vertex.weight == "outgoing"){
      if (length(slot(object, "netP")$centr) == 0) {
        stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
      }
      vertex.weight = object@netP$centr[[signaling]]$outdeg
    }
    coordinates <- object@images$coordinates
    labels <- object@idents
    gg <- netVisual_spatial(prob.sum, coordinates = coordinates, labels = labels, alpha.image = alpha.image, point.size = point.size, sources.use = sources.use, targets.use = targets.use, idents.use = idents.use, remove.isolate = remove.isolate, top = top, color.use = color.use, vertex.weight = vertex.weight, vertex.weight.max = vertex.weight.max, vertex.size.max = vertex.size.max, weight.scale = weight.scale, edge.weight.max = edge.weight.max, edge.width.max=edge.width.max,title.name = paste0(signaling.name, " signaling pathway network"), vertex.label.cex = vertex.label.cex,...)

  } else if (layout == "chord") {
    prob.sum <- apply(prob, c(1,2), sum)
    gg <- netVisual_chord_cell_internal(prob.sum, color.use = color.use, sources.use = sources.use, targets.use = targets.use, remove.isolate = remove.isolate,
                                        group = group, cell.order = cell.order,
                                        lab.cex = vertex.label.cex,small.gap = small.gap, big.gap = big.gap,
                                        scale = scale, reduce = reduce,
                                        title.name = NULL, show.legend = show.legend, legend.pos.x = legend.pos.x, legend.pos.y= legend.pos.y)
  }

  return(gg)

}