# JRC_sorter 
## A classifier for Jointly Regulated CpGs (JRCs) 

#### Benjamin Planterose Jiménez, Brontë Kolar, Manfred Kayser, Athina Vidaki
Department of Genetic Identification, Erasmus MC University Medical Center Rotterdam, Rotterdam, the Netherlands.

## Requirements

    Operating system: tested on Ubuntu 18.04.6 LTS
    R: tested on R version 4.1.2 (2021-11-01) -- "Bird Hippie" & R version 4.2.2 Patched (2022-11-10 r83330) -- "Innocent and Trusting"

For issues/questions on JRC_sorter, either [report an issue](https://github.com/BenjaminPlanterose/JRC_sorter/issues) or simply contact me via b.planterosejimenez at erasmusmc.nl

**IMPORTANT**: Test data available at [Zenodo](https://zenodo.org/record/7636817#.Y-qHlRzMJu1).

## Dependencies

#### [R](https://cran.r-project.org/)

To install R, click on the hyperlink above and follow instructions. To download dependency R-packages, do as follows. Open R by running:
```bash
R
```
And run the following R commands:

```r
# To install
install.packages('data.table')
install.packages('parallel')
install.packages('MASS')

# To verify that installation was succesful
library(data.table)
library(parallel)
library(MASS)
```

#### JRC_sorter

Obtain manual page by running:
```bash
bash <path_to_JRC_sorter>/src/JRC_sorter --h
# Usage: JRC_sorter -t FILE -i FILE [-l INT -e DOUBLE -b DOUBLE -m DOUBLE -o CHAR -c INT -f INT -a INT]
#
#   -t           A file containing target chromosomic locations (one per row) in the following format chr1:1234-3456.
#   -i           An input directory that contains a folder for each sample. Each sample folder contains files split by chromosome containing M and U counts.
#   -l           Bin length. Default value is 200 bp.
#   -e           Value for epsilon. Default value is 0.05.
#   -b           Bin likelihood tolerance. Default value is 0.15 (e.g. 15 % of the maximum log(L)).
#   -m           Model likelihood tolerance. Default value is 0.1 (e.g. 10 % of the maximum log(L)).
#   -o           output directory; default value is 'results'.
#   -c           Number of cores to employ; default value is 1.
#   -f           Flank length added left and right to each region; default value is 200 bp.
#   -a           Average count per bin and per sample below which is considered not data. Default value is 18.
```


## Test run

Download ```test_data``` from [Zenodo](https://zenodo.org/record/7636817#.Y-qHlRzMJu1) and uncompress. This is not included in the Github repository since its size exceeds the 100 MB limit. 

```test_data``` is a small part of the study of [Mordaunt et al](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-020-00785-8) (i.e. 10 samples, only chr12, chrX and chrY). The complete dataset is available under GEO accession [GSE140730](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140730).

```test_data``` includes:

* ```im_regions.txt``` - List of regions on which to run JRC_sorter
* ```PROCESSED.zip``` - Please uncompress. ```/PROCESSED/``` contains one directory per sample. Within each sub-directory, a file per chromosome is included:
	* chr12: contains the methylated/unmethylated CpG counts (paired-end WGBS experiment) on chr12. The file contains the following fields
		* Fwd_Chromosome
		* Fwd_CpG_start
		* Fwd_CpG_end
		* Fwd_M
		* Fwd_U
		* Fwd_CG, confirms they are CpG sites
		* Fwd_CpG_sequence_context
		* All upper fields but for Rv direction
	* chrX/Y: Sex is predicted from the size of these files.
* ```expected_output.zip``` - contains the expected output of JRC_sorter on this data (for debugging purposes).


To launch JRC_sorter on the example, attempt:

```bash
bash <path_to_JRC_sorter>/src/JRC_sorter -t <path_to_test_run>/im_regions.txt \ 
-i <path_to_test_run>/PROCESSED -l 200 -f 200 -b 0.15 -m 0.10 -e 0.05 -c 1 -a 18 -o results
```


## How to run on your own data?

We processed standard CpG_reports (obtained with Bismark by Mordaunt *et al*) with ```/src/opt/process_samples.sh```.

We also include some useful functions in ```/src/opt/useful_functions.R```.

