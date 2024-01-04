### Instructions to obtain processed data

From terminal (recommended): 
```sh
### may use Zenodo instead of GEO to share the processed data ###

# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE252nnn/GSE252470/suppl/GSE252470_.rds.gz 
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE252nnn/GSE252470/suppl/GSE252470_.rds.gz 
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE252nnn/GSE252470/suppl/GSE252470_.rds.gz 
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE252nnn/GSE252470/suppl/GSE252470_.rds.gz 
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE252nnn/GSE252470/suppl/GSE252470_.rds.gz 
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE252nnn/GSE252470/suppl/GSE252470_.rds.gz

# gunzip *.rds.gz
```
From R (NOTE: be sure to change `dest` to the path you want to download to!):
```r
### may use Zenodo instead of GEO to share the processed data ###

# utils::download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE252nnn/GSE252470/suppl/GSE252470_.rds.gz", dest = "/pwd/to/dir/.rds.gz")
# utils::download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE252nnn/GSE252470/suppl/GSE252470_.rds.gz", dest = "/pwd/to/dir/.rds.gz")
# utils::download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE252nnn/GSE252470/suppl/GSE252470_.rds.gz", dest = "/pwd/to/dir/.rds.gz")
# utils::download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE252nnn/GSE252470/suppl/GSE252470_.rds.gz", dest = "/pwd/to/dir/.rds.gz")
# utils::download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE252nnn/GSE252470/suppl/GSE252470_.rds.gz", dest = "/pwd/to/dir/.rds.gz")
# utils::download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE252nnn/GSE252470/suppl/GSE252470_.rds.gz", dest = "/pwd/to/dir/.rds.gz")
```
You can then use the R function `GEOquery::gunzip`, terminal `gunzip`, or prefered method to unzip the files

### Instructions to obtain count matrices from NCBI GEO
Navigate to directory in which you wish to download the data and run the following:
```sh
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE252nnn/GSE252470/suppl/GSE252470_RAW.tar
```
(Alternatively, you can download the zip folder by visiting the [GSE252470](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE252470) page)


Then upacked with:
```sh
tar -xf GSE252470_RAW.tar
```
(Or, right click and extract all)


To rearrange the file structure for easy loading into Seurat, follow the code chunks below (Alternatively, manually rearrange):

First, get the provided deCoder.tsv file in the directory you wish to house the input files then:
```sh
touch rearrange.sh
```

Then, copy the contents of the below script into the rearrange.sh file
```sh
#!/usr/bin/env bash

while read line
do
    old=$(echo "$line" | cut -f1)
    new=$(echo "$line" | cut -f2)
    
    mkdir $new
    
    barcodes="./${old}_barcodes.tsv.gz"
    feats="./${old}_features.tsv.gz"
    mtx="./${old}_matrix.mtx.gz"

    mv $barcodes ./$new/$new\_barcodes.tsv.gz
    mv $feats ./$new/$new\_features.tsv.gz
    mv $mtx ./$new/$new\_matrix.mtx.gz

done < deCoder.tsv
```

After transfering the code, run:
```sh
bash rearrange.sh
```

From there you can follow the provided analysis scripts.

### Instructions to obtain raw data from SRA
Navigate to directory of interest and run for each file you wish to pull down. Then with SRA toolkit installed run:

NOTE: all raw data files are just under 2TB of data
```sh
prefetch -v --max-size=55000000 SRR #smallest file for ex
fastq-dump --gzip --split-files SRR
```
From there you will have to modify the file names to be compatible with the cellranger input (if using cellranger). It expects `[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz` -- feel free to reach out if you have trouble

