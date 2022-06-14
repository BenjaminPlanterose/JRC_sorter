# JRC_sorter 
## A classifier for Jointly Regulated CpGs (JRCs) 

#### Benjamin Planterose Jiménez<sup>1</sup>, Brontë Kolar<sup>1</sup>, Manfred Kayser<sup>1</sup>, Athina Vidaki<sup>1</sup>

<sup>1</sup> Department of Genetic Identification, Erasmus MC University Medical Center Rotterdam, Rotterdam, the Netherlands.


## Requirements

    Operating system: tested on Ubuntu 18.04.6 LTS
    R: tested on R version 4.1.2 (2021-11-01) -- "Bird Hippie"


## Installation

To begin with, make sure the right permissions are given to JRC_sorter

```bash
chmod 777 JRC_sorter
```

To install R dependencies, open R by running:
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

Finally, we need to export JRC_sorter´s path to the .bashrc file. Run in linux bash:

```bash
cd ~
vim .bashrc
```

Press i to insert text and write at the end of the file the following line:
```bash
export PATH="$PATH:/absolute_route_to_JRC_sorter/"
```
To save changes and close VIM, press *Esc* key, write *wq* and press enter. Finally, we source the .bashrc by running:

```bash
source .bashrc
```

Verify this step by executing:

```bash
# This should display the help page for binokulars
JRC_sorter --h
```

## Example

Example cord blood WGBS files can be download from GEO entry [GSE140730](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140730).

These need to be uncompressed and processed (Cs in Fwd/Rv strands need to be combined into one row). A script to do so is available in /src/process_samples.sh.


To run the example, a test file regions.txt is given under /example/. To run:

```bash
cd example
JRC_sorter -t /home/ultron/Git/JRC_sorter/example/regions.txt -i /media/ultron/2tb_disk2/0_startallover/followup_meQTLs/cord_blood/test2/PROCESSED/ -l 200 -f 100 -r 100 -e 0.05 -c 2 -o results
```

Please contact me at b.planterosejimenez at erasmusmc.nl for any questions or issues concerning the scripts.

### References and Supporting Information
B. Planterose *et al* (**2022**). Identifying jointly regulated CpGs in the human methylome. *TBD*
