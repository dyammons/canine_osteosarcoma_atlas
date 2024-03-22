# canine_osteosarcoma_atlas

This GitHub repository contains all the analysis code used in, "Single-cell RNA sequencing reveals the cellular and molecular composition of treatment-naïve primary canine osteosarcoma tumors"

The manuscript is currently under review and the GitHub page will be finalized by the time of publication.

A pre-print is available at:
> https://www.researchsquare.com/article/rs-3232360/v1

## Repository goals: 
- provide a resource to make the data generated from this project accessible
- enable reproducible/transparent data reporting
- provide analysis code to reproduce custom figures

If you have any questions or concerns, please submit an issue, contact the corresponding author(s), and/or contact Dylan Ammons at dylan.ammons @ colostate dot edu.

## File structure:
- [:file\_folder: input](/input) contains relevant metadata files and instructions for obtaining data associated with this study
- [:file\_folder: analysis](/analysis) contains the analysis code (largely separated in scripts by figure) and source file used to complete the data analysis

## Supplemental data and potential uses:
1. [Browse the data](#1-browse-the-complete-annotated-dataset)
1. [Cell type annotations](#2-cell-type-annotations-with-defining-markers)
1. [Reference Mapping](#3-using-the-data-to-complete-reference-mapping)
1. [GSEA using dataset](#4-gene-set-enrichment-analysis)
1. [Module scoring](#5-module-scoring)
1. [Deconvoloution](#6-deconvoloution-of-bulkRNA-seq-data)

### 1. Browse the complete annotated dataset

The processed dataset will be available for browsing via the UCSC Cell Browser portal shortly.
Using the portal you can explore feature expression throughout the dataset as well as obtain the transcriptomic signatures of each cell type though an interactive webpage.

Link to the dataset: TBD

Link to UCSC Cell Browser documentation: https://cellbrowser.readthedocs.io/en/master/

### 2. Cell type annotations with defining markers

Cell markers lists will be curated for each cell type.

<details><summary>cellType.l1</summary>
<p>
  
|Cell Type         |Marker                                                             |
|------------------|-------------------------------------------------------------------|

  
</p>
</details>

<details><summary>cellType.l2</summary>
<p>

|Cell type         |                                                                    |Marker                                                              |
|------------------|--------------------------------------------------------------------|--------------------------------------------------------------------|

</p>
</details>

<details open><summary>cellType.l3 (Default)</summary>
<p>

|Cell Type     |                       |Marker                                                                   |
|--------------|-----------------------|-------------------------------------------------------------------------|


</p>
</details>

### 3. Using the data to complete reference mapping
Reference mapping is useful tool to facilitate the identification of cell types in single cell datasets. The approach described here uses Seurat functions to identify anchors between a query dataset (external/personal data) and the reference datasets generated in this study.

NOTE: this will be avalible at time of release on NCBI GEO.

Before running the reference mapping code, a Seurat object need to be preprocessed and stored as an object named `seu.obj`.
```r
#set the path to the location in which the reference file is saved
reference <- readRDS(file = "../../k9_PBMC_scRNA/analysis/output/s3/final_dataSet_HvO.rds")

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
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = F
)
ggsave("./output/referenceMap.png", width = 7, height = 7)
```

### 4. Gene set enrichment analysis

NOTE: this will be avalible at time of release on NCBI GEO.

The data generated from this work have the potential to provide supporting evidence to evaluate/confirm the cell identity of sorted bulk RNA sequencing dataset. One approach to do this is to use gene set enrichment analysis (GSEA) with the terms representing the cell type identified in our dataset.

Required input: a list of gene symbols that you wish to query. In this case the genelists are stored in a dataframe called `clus.markers`

These gene lists could be generated by simply using the features with the highest level of expression after normalizing your dataset, comparing the transcriptome of a cell population of interest (i.e., blood-derived macrophage) verses a reference (i.e., total PBMCs), or any other relevant approach to identify genes of interest.

Example data frame format:

```r
> str(clus.markers)
'data.frame':   400 obs. of  2 variables:
 $ gene   : chr  "B2M" "CD74" "DLA-64" "PPBP" ...
 $ cluster: chr  "sample_1" "sample_2" "sample_3" "sample_4" ...
```

```r
#read in the one of the supplemntal data files provided with the publication
geneLists <- read.csv(file = "./input/supplementalData_1.csv")

#clean the reference
datas <- geneLists[,c("cluster","gene")]
colnames(datas) <- c("gs_name", "gene_symbol")
datas <- datas %>% group_by(gs_name) %>% top_n(50) %>% dplyr::distinct(gene_symbol) %>% as.data.frame()

#run GSEA using clusterProfiler
clusters <- unique(clus.markers$cluster)
df.list <- list()
for (cluster in clusters) {
    clus_sub <- clus.markers[clus.markers$cluster == cluster, ]

    #run enricher
    enriched <- as.data.frame(clusterProfiler::enricher(gene = clus_sub$gene, TERM2GENE = datas, pvalueCutoff = 1))
    if(nrow(enriched) > 0){
        enriched$cluster <- cluster
        enriched <- head(enriched) #only takes the top 6 terms - can modify if desired
        df.list[[which(cluster == clusters)]] <- enriched
    }
}

cellCalls <- do.call(rbind, df.list)
outfile <- paste("./output/cell_classification.csv", sep = "")
write.csv(cellCalls, file = outfile)

#plot the data
plot <- ggplot(data = cellCalls, mapping = aes_string(x = 'cluster', y = 'ID')) +
    geom_point(mapping = aes_string(size = 'Count', color = -log2(cellCalls$p.adjust))) +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.background = element_rect(fill = "white",colour = NA),
          plot.background = element_rect(fill = "white",colour = NA),
          legend.background = element_rect(fill = "transparent",colour = NA),
          legend.key = element_rect(fill = "transparent",colour = NA),
          panel.grid.major = element_line(color = "gray"), 
          panel.grid.minor = element_line(color = "gray")
          ) + 
    scale_colour_viridis(option="magma", name='-log2(padj)') +
    guides(size=guide_legend(title="Gene count", override.aes=list(fill=NA))) +
    geom_text(aes(y = as.numeric(length(unique(ID))), label = cluster), size = 3.5,vjust = -.3, angle=45, hjust = -0.1) +
    coord_cartesian(expand = TRUE, clip = "off") +
    xlab("Sample") + ylab("GSEA term")

#check path is correct
ggsave("./output/gsea_scRNA_terms.png", width = 6, height = 4)
```

### 5. Module scoring

NOTE: this will be avalible at time of release on NCBI GEO.

Module scoring is a supplemental approach that can be applied to single cell datasets with the goal of providing further insights into cell identities. The approach described below uses the Seurat function `AddModuleScore` and the gene lists presented in Table 3 (also found in supplemental data 4) of our associated manuscript. 

The concept of the AddModuleScore() function is similar to GSEA, but also distinct in many ways. Read the [Seurat documentation](https://satijalab.org/seurat/reference/addmodulescore) and/or check out [this webpage](https://www.waltermuskovic.com/2021/04/15/seurat-s-addmodulescore-function/) for more details.

```r
#load in the reference file from supplemental data
ref.df <- read.csv("supplementalData_4.csv", header = T)

#organize the data
modulez <- split(ref.df$gene, ref.df$cellType_l2)

#complete module scoring
seu.obj <- AddModuleScore(seu.obj,
                          features = modulez,
                         name = "_score")

#correct the naming -- credit to: https://github.com/satijalab/seurat/issues/2559
names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)

#plot the results -- uses a custom function, so you will need to source the customFeunctions.R file. Alt: can also be visulized with FeaturePlot() or DotPlot()
features <- names(modulez)
ecScores <- majorDot(seu.obj = seu.obj, groupBy = "clusterID_sub", scale = T,
                     features = features
                    ) + theme(axis.title = element_blank(),
                              #axis.ticks = element_blank(),
                              #legend.justification = "left",
                              #plot.margin = margin(7, 21, 7, 7, "pt")
                              legend.direction = "vertical",
                              legend.position = "right"
                             ) + guides(color = guide_colorbar(title = 'Scaled\nenrichment\nscore')) + guides(size = guide_legend(nrow = 3, byrow = F, title = 'Percent\nenriched'))

ggsave(paste("./output/", outName, "/", outName, "_dots_celltypes.png", sep = ""),width = 10,height=6)
```

### 6. Deconvoloution of bulkRNA seq data

Under development

