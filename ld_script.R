setwd("~/cmcouto.silva@usp.br/lab_files/all_datasets/HGDP_NAM/ld_analysis/")

if (!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("cmcouto-silva/LD")
library(LD)
gt::rm.all()

# Initialization
ldstart()

# Produce required files
make.binary.haplotypes()
make.input.files()
make.scanhh()

# Load 'scanhh.list'
load("./scanhh.RData")

# Computing Statistics and generating result files
ihs()
xpehh(scanhh.list$Han, scanhh.list$Nam)
rsb(scanhh.list$Han, scanhh.list$Nam)

