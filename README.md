# Linkage Disequilibrium Analysis

##### Package with an automatized workflow to 'rehh' package's functions. 

#### Installation
```r
# install.packages("devtools")
devtools::install_github("cmcouto-silva/LD")
```

#### Workflow instructions

Supposing you're handling a dataset with chromossomes from 1 to 22, here is your to-do steps:

##### **Initialization:**
```r
ldstart()
```

This creates five folders: 'haplotypes', 'populations', 'snps', 'rehh_in', and 'rehh_out'.

Then all you have to do (manually) is to put your haplotypes' files in the 'haplotypes' folder, the populations' files in the 'populations' folder, and the target-SNPs' files in the 'snps' folder. All haplotypes' files must be in conventional  .haps/.sample format. Population and SNP's files must be text files (.txt) with one ID per line.

For more details, see `help(ldstart)`


Following workflow steps are as simple as:
```r
# Produce required files
make.binary.haplotypes()
make.input.files()
make.scanhh()

# Load blahblahblah
load("./scanhh.RData")

# Computing Statistics and generating result files
ihs()
xpehh(scanhh.RData$population1, scanhh.RData$population2)
```

#### **At the end of this workflow you'll have:**
- iHS result files
  - Simple .csv files with original results
  - Results with SNP annotation and allele state information (.xls files)
  - Figures showing the distribution plots
- XP-EHH result files
  - Simple .csv files with original results
  - Results with SNP annotation and allele state information (.xls files)
  - Figures showing the XP-EHH-related plots
- Two R binary files
  - "haps-sample.RData" 
  - "_scanhh.RData_"  

This package has been created pursuing an easy, simple, and fast execution of a workflow for linkage disequilibrium analysis,
using algorithm's computations avaialable in 'rehh' package (which reproduce Statistics proposed by the original papers). 
Therefore, there is a reasonable kind of freedom restriction expected to users. 
Nevertheless, useful parameters are available and can be consulted in the help session for each function.

The R binary file from class "_haplohh_" named "_haps-sample.RData_" allows users to run `rehh::scanhh()` function with 
desired arguments available in rehh package's parameters, while R binary file named "_scanhh.RData_" allows users 
to perform all set of algorithms available in 'rehh' package as they want to, giving them freedom of parameters/arguments 
choice and data manipulation.

This workflow works just on Unix-based platforms (_e.g._ Linux and MAC OS).
In the next few days will be implemented an internal-function which enables Windows users to use it as well.

Vignettes and reproducible example files will be included soon.

#### Best regards,
##### Cainã Max Couto-Silva













