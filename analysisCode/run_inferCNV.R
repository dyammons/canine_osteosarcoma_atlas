#!/usr/bin/Rscript

### Execution:
# Rscript seurat_inferCNVindv.R -t12 -i1

#create parser object
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()

#specify our desired options 
#by default ArgumentParser will add an help option 
parser$add_argument("-t", "--ntask", type="integer", default=1,
    help="Number of cores to use [default]")
parser$add_argument("-i", "--sampleIndex", type="integer", 
    default=1, help="Index of sample to run, based on Seurat object orig.ident value")
                                        
#get command line options, if help option encountered print help and exit,
#otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

pthread <- args$ntask
index <- args$sampleIndex+1

#load da source
source("./customFunctions.R")
library(infercnv)

#read in data & load celltype metadata
seu.tumor <- readRDS(file = "../output/s3/canine_naive_n6_annotated.rds")

#split out objects to process 1 at a time
seu.sub.list <- SplitObject(seu.tumor, split.by = "orig.ident")
seu.obj <- seu.sub.list[[index]]

#remove neutrophils for this analysis
seu.obj <- subset(seu.obj, invert = T, subset = majorID == "neut")

#clean env
rm(seu.sub.list)
gc()

#DO NOT RUN -- bash cmd to make clean gene list --
#tail -n +6 CanFam3.1.104.filtered.gtf | awk '$3 == "gene" {print $0}' | sed 's/[";]//g;' | awk '{OFS="\t"; print $1, $4,$5,$14,$10}' | sort > gene_position.txt

#load clean gene list and then sort it
gene.pos <- as.data.frame(read.csv(file = "./metaDatagene_position.txt", header = FALSE, sep = "\t"))
gene.pos.sorted <- gene.pos[order( gene.pos[,1], gene.pos[,2] ),]
gene.pos.sorted$V4 <- ifelse(gene.pos.sorted$V4 == "ensembl", gene.pos.sorted$V5,gene.pos.sorted$V4)
gene.pos.sorted <- gene.pos.sorted[,c(4,1,2,3)]
gene.pos.sorted <- gene.pos.sorted[!duplicated(gene.pos.sorted$V4),]
rownames(gene.pos.sorted) <- gene.pos.sorted$V4
gene.pos.sorted <- gene.pos.sorted[,-1]

#run inferCNV
message("[", paste0(Sys.time(), "]", " INFO: initating inferCNV run on ",unique(seu.obj$orig.ident), " using ", pthread, " cores." ))
matrix <- as.matrix(seu.obj@assays$RNA@counts)
ann <- as.data.frame(seu.obj$majorID, row.names=rownames(seu.obj@meta.data))

inferobj <- CreateInfercnvObject(raw_counts_matrix=matrix,
                                 annotations_file=ann,
                                 gene_order_file=gene.pos.sorted,
                                 ref_group_names=c("endo","mac")
                                )

outDir <- paste0("../inferCNV/output_noNeuts_cleanData_", unique(seu.obj$orig.ident),"/")
dir.create(outDir)
inferobj <- infercnv::run(inferobj,
                      cutoff=0.1,
                      out_dir=outDir,
                      cluster_by_groups=FALSE,
                      plot_steps = FALSE,
                      no_prelim_plot = TRUE,
                      denoise=TRUE,
                      num_ref_groups = 2,
                      HMM=TRUE,
                          num_threads = pthread,
                      analysis_mode='subclusters',
                      png_res = 300)

message("[", paste0(Sys.time(), "]", " INFO: end of inferCNV run. This script processed ",unique(seu.obj$orig.ident), " using ", pthread, " cores." ))
