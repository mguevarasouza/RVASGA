REQUIREMENTS
------------
In order to run RVASGA, the following R libraries must be installed. They can be found in CRAN. 
-GA
-WriteXLS
-hash
-SKAT

To install a library the following command can be used.
-install.packages("GA")

Installation
------------
RVASGA is an R script so no installation is needed. In runs from command line using R environment.

CONFIGURATION
-------------
The script has no menu or modifiable settings. There is no configuration. 

Preparation
-------------
The following folder structure is assumed

├── Working Directory
|   ├── 1
|   └── 2
|   └── 3
|   └── ...
|   └── 22

A function called splitFiles is included in the script. This function will create the folder
structure mentioned above. To use it, call splitFiles() from R in the directory where the .geno,
.pheno and .mapping are.

EXECUTION
-------------

First, the script should be loaded to workspace. From command line use:
-source("/path/to/script/RVASGA.r")

RVASGA is conformed by several functions but there is a main function that runs all 
the functions needed. Once the script is loaded, change to the root of the dataset 
(where the folders of the chromosomes are) and run:

runAllGA(5000000)

The parameter is used as an offset to generate a window that contains the genes to be
analyzed. The offset is applied to the left and to the right, so the size of the window
is the double of the offset.

It is also possible to run the algorithm for one chromosome at a time. In this case, 
a second parameter is needed. To run the example chromosome use:

RVASGA(3000000,17)

The results will be stored in the 17.csv file.

INPUT FILE FORMAT
-------------

RVASGA requires two text files as input. Both files should be in the folder passed as
argument. 

Phenotype file:

A phenotype file with with rows representing samples, the 1st column is sample name, the 2nd column is the quantitative or binary phenotype.

Example:

SAMP1	1
SAMP2	1
SAMP3	1
SAMP4	1


A genotype file file with rows representing variants and the columns represent sample genotypes (order of the rows matches the genotype file).

Example:

GEN1:55524145	00	00	00	00	00	00	00	00	00
GEN1:55524161	00	00	00	00	00	00	00	00	00
GEN1:55561658	00	00	00	00	00	00	00	00	00
GEN1:55561805	00	00	00	00	00	00	00	00	00
GEN1:55561861	00	00	00	00	00	00	00	00	00
GEN1:55564614	00	00	00	00	00	00	00	00	00
GEN1:55564616	00	00	00	00	00	00	00	00	00
GEN1:55570043	00	00	00	00	00	00	00	00	00
GEN1:55570067	00	00	00	00	00	00	00	00	00
GEN1:55589871	00	00	00	00	00	00	00	00	00
GEN1:55592059	00	00	00	00	00	00	00	00	00


Maintainers
------------
-Mauricio Guevara Souza (mauricio.guevara.souza@gmail.com)
