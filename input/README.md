## Instructions to obtain processed Seurat objects  

The processed data are available on Zenodo and can be downloaded by visiting the [project repository web page](https://zenodo.org/records/10666969).
Once on the web page scroll down and select download for the file(s) of interest.

Alternatively, use `wget` in a terminal to retrieve the data:
```sh
wget https://zenodo.org/record/10891255/files/canine_naive_n6_annotated.rds  # Full dataset
wget https://zenodo.org/record/10891255/files/tumor_subset_annotated.rds     # Tumor cell dataset
wget https://zenodo.org/record/10891255/files/tcell_subset_annotated.rds     # T cell dataset
wget https://zenodo.org/record/10891255/files/dc_subset_annotated.rds        # Dendritic cell dataset
wget https://zenodo.org/record/10891255/files/macOC_subset_annotated.rds     # Macrophage and osteoclast dataset
wget https://zenodo.org/record/10891255/files/myeloid_subset_annotated.rds   # Macrophage, osteoclast, and dendritic cell dataset
```

Prefer to use tools in python or R? Check out `zenodo_get` or `inborutils` to download within the respective software. 

<details><summary>Example usage of zenodo_get </summary>
<p>

Below is the code needed to install `zendo_get` using `pip` and the command to download the repositiry specific to this project (this should be completed in an environment with python3 installed).  

Visit the [`zendo_get`](https://github.com/dvolgyes/zenodo_get) page for most up to date instructions.

```sh
#install the python tool using pip
pip3 install zenodo_get

#download the Zenodo repository
zenodo_get 10.5281/zenodo.10891255
```

</p>
</details>

## Instructions to obtain count matrices from NCBI GEO  

### Retrieve the data

To download via command line, navigate to directory in which you wish to download the data and run the following:
```sh
#pull down the data
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE252nnn/GSE252470/suppl/GSE252470_RAW.tar

#unpack the tar ball
tar -xvf GSE252470_RAW.tar
```

To download via NCBI webpage, navigate to the [GSE252470](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE252470) page and download the zip folder containing the count matrices. Once downloaded, you will likely need to extract/unpack the data before using (should be able to do by right clicking on the `GSE252470_RAW.tar` file).

### Rearrange the data (optional)

The count matrices were renamed when uploaded to NCBI. If you wish to use the analysis in it's entirety you will need to rename the files accordingly.

To rearrange the file structure for easy loading into Seurat, follow the code chunks below (alternatively, the files can be manually rearrange/renamed):

First, get the provided `deCoder.tsv` (TO DO: add deCoder.tsv !!) file in the directory you wish to house the input files then create a file called `rearrange.sh`.

```sh
#pull down deCoder.tsv
wget https://github.com/dyammons/canine_osteosarcoma_atlas/tree/main/input/deCoder.tsv

#create an empty file
touch rearrange.sh
```

Once `rearrange.sh` is created, copy the contents of the code block below into the file.

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

After transferring the code, run:
```sh
bash rearrange.sh
```

Now the files should be rearranged and renamed as needed to be identical to the nomenclature used in the associated analysis code.

## Instructions to obtain raw data from SRA
Navigate to directory of interest and run for each file you wish to pull down. Then with SRA toolkit installed run:

NOTE: all raw data files are around 1TB of data
```sh
prefetch -v --max-size=55000000 SRR #smallest file for ex
fastq-dump --gzip --split-files SRR
```
From there you will have to modify the file names to be compatible with the Cell Ranger input (if using Cell Ranger). It expects:  
`[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz`  

Feel free to reach out (email or create an issue on GitHub) if you have trouble getting the data downloaded.

