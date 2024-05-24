# Canine osteosarcoma scRNA-seq atlas

[![DOI](https://zenodo.org/badge/630249079.svg)](https://zenodo.org/doi/10.5281/zenodo.10666968)

This GitHub repository contains all the analysis code used in, "Single-cell RNA sequencing reveals the cellular and molecular composition of treatment-naïve primary canine osteosarcoma tumors"

If you use our raw/processed data, extract data using the UCSC Cell Browser portal, or use portions of our code in your analysis, please cite:
> Ammons, D.T., Hopkins, L.S., Cronise, K.E., Kurihara, J., Regan, D.P. and Dow, S., 2024. Single-cell RNA sequencing reveals the cellular and molecular heterogeneity of treatment-naïve primary osteosarcoma in dogs. Communications Biology, 7(1), p.496. https://doi.org/10.1038/s42003-024-06182-w

Interested in more canine scRNA data? Check out our [healthy and OS canine PBMC atlas](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2023.1162700/full) and the [associated GitHub page](https://github.com/dyammons/Canine_Leukocyte_scRNA).

## Repository goals: 
- Provide a resource to make the data generated from this project accessible
- Enable reproducible/transparent data reporting
- Provide analysis code to reproduce custom figures

If you have any questions or concerns, please submit an issue via GitHub and/or contact the corresponding author via email: dylan.ammons @ colostate dot edu.

## File structure:
- [:file\_folder: input](/input) contains relevant metadata files and instructions for obtaining data associated with this study
- [:file\_folder: analysisCode](/analysisCode) contains the analysis code (largely separated in scripts by figure) and R source file containing custom functions used to complete data analysis
- [:file\_folder: output](/output) contains the expected output directory structure is completing a reproducible run


## Supplemental data and potential uses:
1. [Browse the data](#browse-the-complete-annotated-dataset)
1. [Cell type annotations](#cell-type-annotations-with-defining-markers)
1. [Reference Mapping to scRNA data](#using-the-data-to-complete-reference-mapping-scrna-seq)
1. [Reference Mapping to spatial data](#using-the-data-to-complete-reference-mapping-spatial-transcriptomics)
1. [GSEA using dataset](#gene-set-enrichment-analysis)
1. [Module scoring](#module-scoring)
1. [Deconvoloution](#deconvoloution-of-bulkRNA-seq-data)

<br>

### Browse the complete annotated dataset

The processed dataset is available for browsing via the UCSC Cell Browser portal.
Using the portal you can explore feature expression as well as obtain the transcriptomic signatures of each cell type though an interactive webpage.

Link to the dataset: https://cells.ucsc.edu/?ds=canine-os-atlas

Link to UCSC Cell Browser documentation: https://cellbrowser.readthedocs.io/en/master/

<br>

### Cell type annotations with defining markers

Cell markers lists were curated by evaluating the top 50 defining features (identified using `FindMarkers`) for each cell type. The specificity of the top features was evaluated using violin plots and the UCSC cell browser. Preference given to unique features only found in the top 50 of one cell type.

<details open><summary>Cell types (High-resolution)</summary>
<p>

|Cell type     |                                |Markers                                       |
|--------------|--------------------------------|----------------------------------------------|
|B cell        |                                |                                              |
|              |B cell                          |PAX5, CD22, MS4A1, FCRLA, CCR7, IGHM          |
|              |Plasma cell                     |JCHAIN, DERL3, TNFRSF17, MZB1, POU2AF1        |
|T cell        |                                |                                              |
|              |CD4 naïve                       |CCR7, CD52, LTB, LEF1, TCF7                   |
|              |CD4 activated                   |CXCR4, CD28, IL2RB, IL7R, ICOS                |
|              |CD4 follicular helper           |TNFRSF18, TNFSF8, PDCD1, CXCL13, IL4I1        |
|              |CD4 regulatory                  |IL2RB, GATA3, OCIAD2, ARID5B, IL21R           |
|              |T-IFN                           |GZMA, ISG20, CCL5, IFI44L, OAS1               |
|              |CD8 SPP1+                       |DNAJB1, HSP90AA1, FOS, HSPB1, C6H7orf50       |
|              |CD8 effector                    |GZMB, NCR3, GZMA, CD96, FASLG, IL12RB2        |
|              |CD8 exhausted                   |SEC14L1, PDCD1, GZMK, CCL5, CCL4              |
|              |Cycling T cell                  |H1-5, MKI67, CENPF, SMC2, H2AZ1               |
|Dendritic cell|                                |                                              |
|              |plasmacytoid DC                 |FCRLA, SPATS2L, IGKC, CLEC2D, RYR1, IGF1      |
|              |precursor DC                    |FCRLA, PGLYRP2, DDR2, GPHA2, TCF4             |
|              |mature regulatory DC            |CCR7, FSCN1, IL1I1, MARCKSL1, CCL19, CD274    |
|              |conventional DC subtype 1       |CPNE3, CLEC1B, BATF3, SERPINB6, SMYD3         |
|              |conventional DC subtype 2       |CD300H, CD1C, PID1, LGALS3, MAFB              |
|Tumor         |                                |                                              |
|              |Hypoxic osteoblast              |ENO1, PGF, PTGES, SFRP2, CAV1                 |
|              |Malignant osteoblast subtype 1  |MPP6, LIFR, FBLN, NPY, C1S                    |
|              |Malignant osteoblast subtype 2  |IBSP, SPARC, SMPD3, ALPL, WFDC1               |
|              |Malignant osteoblast subtype 3  |DNAJB1, HSP90AA1, FOS, HERPUD1, HSPH1         |
|              |Cycling osteoblast subtype 1    |UBE2S, DLGAP5, HMMR, TPX2, TUBA1B, CENPF      |
|              |Cycling osteoblast subtype 2    |TPX2, TK1, H1-4, RRM2, DNAJC9                 |
|              |Cycling osteoblast subtype 3    |MCM6, RAD51AP1, HELLS, CDC6, UHRF1            |
|              |Cycling osteoblast subtype 4    |CDC20, PLK1, CENPE, MIK67, DLGAP5, NUF2       |
|              |IFN-osteoblast                  |MX2, OAS1, IFI44, OAS2, IFI6                  |
|Osteoclast    |                                |                                              |
|              |Mature osteoclast               |CRYAB, ATP6V1C1, SLC4A2, CD84, NEURL3, HYAL1  |
|              |CD320 osteoclast                |VDR, HMGA1, APEX1, DDX21, RSL1D1              |
|              |Cycling osteoclast 1/2          |H2AZ1, STMN1, CENPF, CDC20, MKI67             |
|Monocyte      |                                |                                              |
|              |CD4- tumor infiltrating monocyte|CXCL8, VCAN, LYZ, PLBD1, LSP1                 |
|              |CD4+ tumor infiltrating monocyte|IL1B, PTGS2, LTF, THBS1, CXCL8, VCAN          |
|Macrophage    |                                |                                              |
|              |ANGIO-TAM                       |HBEGF, VEGFA, IL18BP, AREG, VEGFC             |
|              |Intermediate TAM                |CTSS, TPI, ENO1, LAMP2, CCL7                  |
|              |Activated TAM                   |CCL3, CD80, CCL19, CD5L, CXCL16, DLA-79       |
|              |Lipid associated TAM (C1QC high)|C1QB, C1QC, PLTP, SERPING1, DAB2, CLDN1       |
|              |Lipid associated TAM (SPP2 high)|TREM2, APOE, CD36, GPNMB, PRDX1               |
|              |IFN-TAM                         |MX2, RSAD2, CCL8, CD40, IL7R, TNFSF10         |
|Miscellaneous |                                |                                              |
|              |Neutrophil                      |SELL, SOD2, CXCL8, CD4, S100A8, PADI3         |
|              |Mast cell                       |MS4A2, IL3RA, ADORA3, CSF2RB, ACE2, HPGD, CPA3|
|              |Fibroblast                      |DCN, IGFBP7, COL3A1, COL6A3, COL12A1, COL6A1  |
|              |Endothelial cell                |CD34, PLVAP, ESM1, EGFL7, FLT1, VWF           |

</p>
</details>

<br>

### Using the data to complete reference mapping (scRNA-seq)
Reference mapping is useful tool to facilitate the identification of cell types in single cell datasets. The approach described here uses Seurat functions to identify anchors between a query dataset (external/personal data) and the reference datasets generated in this study.

Before running the reference mapping code, a query Seurat object needs to be preprocessed and stored as an varible named `seu.obj`. Additionally, the processed `reference` Seurat object needs to be downloaded and loaded into the R session. Instructions to obtain the reference dataset can be found in [:file\_folder: input](/input).

```r
### Reference mapping using Seurat
# NOTE: this was designed to be run with Seurat v4, but should run with a query dataset processed using Seurat v5.
# Please let us know if it is not working for your application.

#set path location where the downloaded reference file is saved
reference <- readRDS(file = "./final_dataSet.rds")

#prepare the reference -- NOTE: run this code block ONLY if the query data was SCT normalized, otherwise skip step
reference[['integrated']] <- as(object = reference[['integrated']] , Class = "SCTAssay")
DefaultAssay(reference) <- "integrated"

#find conserved anchors with query and reference
anchors <- FindTransferAnchors(
    reference = reference,
    query = seu.obj,
    normalization.method = "SCT", #if using log normlaized data, change to "LogNormalize"
    reference.reduction = "pca",
    dims = 1:30
)

# FYI: this is a computationally intense step
# can change refdata argument to use alternate cell type labels (i.e., refdata = reference$celltype.l1)
predictions <- TransferData(
    anchorset = anchors, 
    refdata = reference$celltype.l3,
    dims = 1:30
)

#store the metadata values
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

<br>

### Using the data to complete reference mapping (spatial transcriptomics)

In addition to transfering labels to other single-cell datasets, the data can also be integrated with canine spatial transcriptomics data. The process is very similar to the above approach, but uses the `Spatial` assay and stores the results in a new assay named `predictions`.

The code below assumes `seu.obj` is external/personal data that was preprocessed using a log normalization (+/- integration) workflow.

```r
### Integration with spatial visium data
# NOTE: this was designed to be run with Seurat v5.1.0
# Please let us know if it is not working for your application!

#import processed data
reference <- readRDS(file = "./final_dataSet.rds")

#transfer annotations from canine OS atlas
ref.anchors <- FindTransferAnchors(
    reference = reference,
    query = seu.obj,
    query.assay = "Spatial",
    dims = 1:30,
    reference.reduction = "pca", 
    features = rownames(seu.obj)[rownames(seu.obj) %in% rownames(reference)]
)
predictions.assay <- TransferData(
    anchorset = ref.anchors, 
    refdata = reference$celltype.l3, 
    prediction.assay = TRUE,
    dims = 1:30
)

#store predicted cell types in a new assay for easy plotting
seu.obj[["predictions"]] <- predictions.assay

#generate and save the image
pi <- DimPlot(seu.obj, 
              reduction = "umap", 
              group.by = "predicted.id",
              pt.size = 0.25,
              label = T,
              label.box = T,
              shuffle = F
)
ggsave("../output/referenceMap.png", width = 10, height = 7)

#plot and save the label transfer on the tissues!!!
DefaultAssay(seu.obj) <- "predictions"
lapply(rownames(seu.obj), function(x){
    p <- SpatialFeaturePlot(seu.obj, features = x, pt.size.factor = 1.6, ncol = 7, crop = TRUE) & 
    theme(legend.title = element_blank()) & patchwork::plot_annotation(title = x)
    ggsave(paste0("../output/", x, "_slide_visium.png"), height = 7, width = 12)    
})
```

<br>

### Gene set enrichment analysis

The data generated from this work have the potential to provide supporting evidence to evaluate/confirm the cell identity of sorted bulk RNA sequencing dataset. One approach to do this is to use gene set enrichment analysis (GSEA) with the terms representing the cell type identified in our dataset.

Required input: a list of gene symbols that you wish to query. In this example the genelists to be tested are stored in a dataframe called `clus.markers`

These gene lists could be generated by simply using the features with the highest level of expression after normalizing your dataset, comparing the transcriptome of a cell population of interest (i.e., blood-derived macrophage) verses a reference (i.e., total PBMCs), or any other relevant approach to identify genes of interest.

Example data frame format:

```r
> str(clus.markers)
'data.frame':   400 obs. of  2 variables:
 $ gene   : chr  "B2M" "CD74" "DLA-64" "PPBP" ...
 $ cluster: chr  "sample_1" "sample_2" "sample_3" "sample_4" ...
```

```r
library(Seurat)
library(tidyverse)

#read in the supplemental data file provided with the publication
geneLists <- read.csv(file = "./input/supplementalData_1.csv") #check file name is correct

#clean the reference
datas <- geneLists[ ,c("cluster","gene")]
colnames(datas) <- c("gs_name", "gene_symbol")

#subset on the 50 defining features (optional)
datas <- datas %>% 
    group_by(gs_name) %>% 
    top_n(50) %>% 
    distinct(gene_symbol) %>% 
    as.data.frame()

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

ggsave("gsea_scRNA_terms.png", width = 6, height = 4)
```

<br>

### Module scoring

Module scoring is a supplemental approach that can be applied to single cell datasets with the goal of providing further insights into cell identities. The approach described below uses the Seurat function `AddModuleScore` and the gene lists presented above (and in supplemental data of our associated manuscript). 

The concept of the `AddModuleScore` function is similar to GSEA, but also distinct in many ways. Read the [Seurat documentation](https://satijalab.org/seurat/reference/addmodulescore) and/or check out [this webpage](https://www.waltermuskovic.com/2021/04/15/seurat-s-addmodulescore-function/) for more details.

```r
#load in the reference file from supplemental data
ref.df <- read.csv("supplementalData_4.csv", header = T) #check file name is correct

#organize the data
modulez <- split(ref.df$gene, ref.df$cellType_l2)

#complete module scoring
seu.obj <- AddModuleScore(seu.obj,
                          features = modulez,
                         name = "_score")

#correct the naming -- credit to: https://github.com/satijalab/seurat/issues/2559
names(seu.obj@meta.data)[grep("_score", names(seu.obj@meta.data))] <- names(modulez)

#plot the results -- uses a custom function, so you will need to source the customFunctions.R file. Alt: can also be visulized with FeaturePlot() or DotPlot()
features <- names(modulez)
ecScores <- majorDot(seu.obj = seu.obj, groupBy = "clusterID_sub", scale = T,
                     features = features
                    ) + theme(axis.title = element_blank(),
                              legend.direction = "vertical",
                              legend.position = "right"
                             ) + guides(color = guide_colorbar(title = 'Scaled\nenrichment\nscore')) + guides(size = guide_legend(nrow = 3, byrow = F, title = 'Percent\nenriched'))

ggsave(paste("./output/", outName, "/", outName, "_dots_celltypes.png", sep = ""),width = 10,height=6)
```

<br>

### Deconvoloution of bulk RNA sequencing data

The data generated from this project provides the data necessary to generate a __canine-specific__ reference to deconvolute bulk RNA-seq data for canine osteosarcoma tumors.  

Currently instructions are not provided, but please reach out with questions as we can provide guidence for reference generation using CIBERSORTx, EPIC, TIMER, or other deconvolution tools.
