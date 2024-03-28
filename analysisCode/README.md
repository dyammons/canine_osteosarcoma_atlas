## Analysis code overview  

### To complete a reproducible run from raw data:
Retrieve the fastq files from SRA and then align to the canine genome. 
Instructions to download the fastq files can be found in [:file\_folder: input](/input). 
The alignment code is currently not provided (used a default Cell Ranger workflow), but can be shared if desired.

### To complete a reproducible run from Cell Ranger output count matrices:
Clone this repository (`git clone https://github.com/dyammons/canine_osteosarcoma_atlas.git`) then follow the instructions in [:file\_folder: input](/input/README.md#instructions-to-obtain-count-matrices-from-ncbi-geo) to obtain count matrices from NCBI GEO to get the count matricies into the `input` directory. Start with `figure1.R` then move through the scripts in increasing figure number.

### To reproduce/explore data using processed data:
Clone this repository (`git clone https://github.com/dyammons/canine_osteosarcoma_atlas.git`) then follow the instructions in [:file\_folder: input](/input/README.md#instructions-to-obtain-processed-seurat-objects) to obtain processed Seurat objects from Zenodo to get the `.rds` files in the `/output/s3/` directory. These objects can be loaded into scripts for figures 2-8.

### Environment to reproduce:
This work was completed on the UC boulder Alpine supercomputer on a Linux operating system. All code was run in a conda environment and I am happy to share the .yml file upon request.

Otherwise, see below for packages loaded into the system.

```r
> sessionInfo()
R version 4.1.1 (2021-08-10)
Platform: x86_64-conda-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux 8.4 (Ootpa)

Matrix products: default
BLAS/LAPACK: /projects/dyammons@colostate.edu/software/anaconda/envs/r_env/lib/libopenblasp-r0.3.18.so

locale:
 [1] LC_CTYPE=C.UTF-8    LC_NUMERIC=C        LC_TIME=C          
 [4] LC_COLLATE=C        LC_MONETARY=C       LC_MESSAGES=C      
 [7] LC_PAPER=C          LC_NAME=C           LC_ADDRESS=C       
[10] LC_TELEPHONE=C      LC_MEASUREMENT=C    LC_IDENTIFICATION=C

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] circlize_0.4.16             singleseqgset_0.1.2.9000   
 [3] Matrix_1.5-1                copykat_1.0.5              
 [5] infercnv_1.10.1             ComplexHeatmap_2.13.2      
 [7] ggtree_3.2.1                ape_5.7-1                  
 [9] scuttle_1.4.0               scRNAseq_2.8.0             
[11] ggpubr_0.4.0                slingshot_2.7.0            
[13] TrajectoryUtils_1.2.0       SingleCellExperiment_1.16.0
[15] princurve_2.1.6             clusterProfiler_4.2.2      
[17] msigdbr_7.5.1               ggsankey_0.0.99999         
[19] lemon_0.4.5                 reshape_0.8.9              
[21] viridis_0.6.2               viridisLite_0.4.1          
[23] SingleR_1.8.1               SeuratDisk_0.0.0.9019      
[25] RColorBrewer_1.1-3          pheatmap_1.0.12            
[27] DESeq2_1.34.0               SummarizedExperiment_1.24.0
[29] Biobase_2.54.0              MatrixGenerics_1.6.0       
[31] matrixStats_1.0.0           GenomicRanges_1.46.1       
[33] GenomeInfoDb_1.30.1         IRanges_2.28.0             
[35] S4Vectors_0.32.4            BiocGenerics_0.40.0        
[37] colorspace_2.0-3            ggrepel_0.9.1              
[39] cowplot_1.1.1               scales_1.2.1               
[41] patchwork_1.1.2             DoubletFinder_2.0.3        
[43] clustree_0.4.4              ggraph_2.0.5               
[45] forcats_0.5.2               stringr_1.4.1              
[47] dplyr_1.0.10                purrr_0.3.5                
[49] readr_2.1.2                 tidyr_1.2.1                
[51] tibble_3.1.8                ggplot2_3.3.6              
[53] tidyverse_1.3.1             SeuratObject_4.1.3         
[55] Seurat_4.3.0               

loaded via a namespace (and not attached):
  [1] rsvd_1.0.5                    ica_1.0-3                    
  [3] Rsamtools_2.10.0              foreach_1.5.2                
  [5] lmtest_0.9-40                 crayon_1.5.2                 
  [7] MASS_7.3-58.1                 nlme_3.1-160                 
  [9] backports_1.4.1               reprex_2.0.1                 
 [11] argparse_2.2.2                GOSemSim_2.20.0              
 [13] rlang_1.0.6                   XVector_0.34.0               
 [15] ROCR_1.0-11                   readxl_1.4.1                 
 [17] irlba_2.3.5                   limma_3.50.3                 
 [19] filelock_1.0.2                BiocParallel_1.28.3          
 [21] rjson_0.2.21                  bit64_4.0.5                  
 [23] glue_1.6.2                    sctransform_0.3.5            
 [25] parallel_4.1.1                spatstat.sparse_3.0-0        
 [27] AnnotationDbi_1.56.2          DOSE_3.20.1                  
 [29] spatstat.geom_3.0-5           haven_2.4.3                  
 [31] tidyselect_1.2.0              fitdistrplus_1.1-8           
 [33] XML_3.99-0.12                 zoo_1.8-9                    
 [35] GenomicAlignments_1.30.0      xtable_1.8-4                 
 [37] magrittr_2.0.3                phyclust_0.1-32              
 [39] cli_3.6.0                     zlibbioc_1.40.0              
 [41] rstudioapi_0.14               miniUI_0.1.1.1               
 [43] sp_1.5-1                      rjags_4-13                   
 [45] fastmatch_1.1-3               lambda.r_1.2.4               
 [47] ensembldb_2.18.3              treeio_1.18.1                
 [49] shiny_1.7.1                   BiocSingular_1.10.0          
 [51] xfun_0.39                     clue_0.3-62                  
 [53] cluster_2.1.4                 caTools_1.18.2               
 [55] tidygraph_1.2.2               KEGGREST_1.34.0              
 [57] interactiveDisplayBase_1.32.0 listenv_0.8.0                
 [59] Biostrings_2.62.0             png_0.1-7                    
 [61] future_1.29.0                 withr_2.5.0                  
 [63] bitops_1.0-7                  ggforce_0.4.1                
 [65] plyr_1.8.7                    cellranger_1.1.0             
 [67] AnnotationFilter_1.18.0       coda_0.19-4                  
 [69] pillar_1.8.1                  gplots_3.1.3                 
 [71] GlobalOptions_0.1.2           cachem_1.0.6                 
 [73] GenomicFeatures_1.46.5        multcomp_1.4-20              
 [75] fs_1.5.2                      hdf5r_1.3.8                  
 [77] GetoptLong_1.0.5              DelayedMatrixStats_1.16.0    
 [79] vctrs_0.5.1                   ellipsis_0.3.2               
 [81] generics_0.1.3                tools_4.1.1                  
 [83] munsell_0.5.0                 tweenr_2.0.2                 
 [85] fgsea_1.20.0                  DelayedArray_0.20.0          
 [87] fastmap_1.1.0                 compiler_4.1.1               
 [89] abind_1.4-5                   httpuv_1.6.5                 
 [91] rtracklayer_1.54.0            ExperimentHub_2.2.1          
 [93] plotly_4.10.1                 GenomeInfoDbData_1.2.7       
 [95] gridExtra_2.3                 edgeR_3.36.0                 
 [97] lattice_0.20-45               deldir_1.0-6                 
 [99] utf8_1.2.2                    later_1.3.0                  
[101] BiocFileCache_2.2.1           jsonlite_1.8.3               
[103] ScaledMatrix_1.2.0            tidytree_0.3.9               
[105] pbapply_1.5-0                 carData_3.0-5                
[107] sparseMatrixStats_1.6.0       genefilter_1.76.0            
[109] lazyeval_0.2.2                promises_1.2.0.1             
[111] car_3.1-0                     doParallel_1.0.17            
[113] goftest_1.2-3                 spatstat.utils_3.0-1         
[115] reticulate_1.34.0             sandwich_3.0-2               
[117] Rtsne_0.16                    downloader_0.4               
[119] uwot_0.1.14                   igraph_1.5.1                 
[121] survival_3.4-0                yaml_2.3.6                   
[123] htmltools_0.5.3               memoise_2.0.1                
[125] modeltools_0.2-23             BiocIO_1.4.0                 
[127] locfit_1.5-9.6                graphlayouts_0.8.3           
[129] digest_0.6.30                 assertthat_0.2.1             
[131] mime_0.12                     rappdirs_0.3.3               
[133] futile.options_1.0.1          RSQLite_2.2.18               
[135] yulab.utils_0.0.5             future.apply_1.10.0          
[137] data.table_1.14.4             blob_1.2.3                   
[139] futile.logger_1.4.3           splines_4.1.1                
[141] AnnotationHub_3.2.2           ProtGenerics_1.26.0          
[143] RCurl_1.98-1.12               broom_1.0.1                  
[145] hms_1.1.2                     modelr_0.1.8                 
[147] BiocManager_1.30.19           shape_1.4.6                  
[149] libcoin_1.0-9                 aplot_0.1.2                  
[151] coin_1.4-2                    Rcpp_1.0.11                  
[153] RANN_2.6.1                    mvtnorm_1.1-3                
[155] enrichplot_1.14.2             fansi_1.0.3                  
[157] tzdb_0.3.0                    parallelly_1.32.1            
[159] R6_2.5.1                      ggridges_0.5.4               
[161] lifecycle_1.0.3               formatR_1.12                 
[163] curl_4.3.3                    ggsignif_0.6.4               
[165] leiden_0.4.3                  fastcluster_1.2.3            
[167] DO.db_2.9                     qvalue_2.26.0                
[169] TH.data_1.1-1                 RcppAnnoy_0.0.20             
[171] iterators_1.0.14              spatstat.explore_3.0-5       
[173] htmlwidgets_1.5.4             beachmat_2.10.0              
[175] polyclip_1.10-4               biomaRt_2.50.3               
[177] shadowtext_0.1.1              timechange_0.1.1             
[179] gridGraphics_0.5-1            rvest_1.0.3                  
[181] globals_0.16.1                spatstat.random_3.1-3        
[183] progressr_0.11.0              codetools_0.2-18             
[185] lubridate_1.9.0               GO.db_3.14.0                 
[187] gtools_3.9.3                  prettyunits_1.1.1            
[189] dbplyr_2.1.1                  gtable_0.3.1                 
[191] DBI_1.1.3                     ggfun_0.0.8                  
[193] tensor_1.5                    httr_1.4.4                   
[195] KernSmooth_2.23-20            stringi_1.7.8                
[197] progress_1.2.2                reshape2_1.4.4               
[199] farver_2.1.1                  annotate_1.72.0              
[201] xml2_1.3.3                    BiocNeighbors_1.12.0         
[203] restfulr_0.0.15               geneplotter_1.72.0           
[205] ggplotify_0.1.0               scattermore_0.8              
[207] BiocVersion_3.14.0            bit_4.0.4                    
[209] scatterpie_0.1.7              spatstat.data_3.0-0          
[211] pkgconfig_2.0.3               babelgene_22.9               
[213] rstatix_0.7.0                 knitr_1.40    
```
