#### LD PACKAGE  ####
####################################################################################################!
devtools::load_all()

# Start and setting files
ldstart()

# Set populations & SNPs
# populations and SNPs were manually settled

## Main workflow

# Produce required files
make.binary.haplotypes()
make.input.files()
make.scanhh()

# Load 'scanhh.list'
load("./scanhh.RData")
names(scanhh.list) <- c('AFR', 'EAS', 'EUR', 'MAY', 'MES', 'NAM', 'OCEANIA', 'PIMA', 'SIB', 'SKOG')

# Computing Statistics and generating result files
ihs(annot = F, plot = T, plot.format = "png")