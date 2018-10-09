#' @title Adjust Ancentral/Derived Alleles According to Reference
#' @description
#' This function uses a reference file obtained from biomart R Package, setting up ancestral alleles,
#' and excluding SNPs without ancestral information from dataset.
#' @param haps Scalar character with the name of the dataset file (.haps).
#' @param ref Scalar character with the name of reference file (obtained by biomart R Package).
#' @param out Scalar character with desired name for the output file (without extension).
#' @details
#' To obtain the Ensembl annotation for ancestral/derived SNPs, please check the function avaialable at:
#' https://github.com/cmcouto-silva/scripts/blob/master/ensembl_anc_allele.R
#' @examples
#' \dontrun{
#' haps <- "dataset/mydataset_phased.haps"
#' ref <- "reference_panels/ensembl_annot.txt"
#' out <- "out_folder/newfile"
#'
#' set_anc_alleles(haps, ref, out)
#'}
#' @return Dataset with adjusted ancestral/derived alleles.
#' @import data.table
#' @importFrom magrittr set_names %>% %<>%
#' @export
#' @author Cainã Max Couto-Silva

set_anc_alleles <- function(haps, ref, out) {

  # Read refecence panel and SNPs' files
  haps <- fread(haps) %>%
    setnames(paste0("V", 1:5), c("CHR", "SNP", "POS", "A1", "A2"))

  ref <- fread(ref)[refsnp_id %in% haps[, SNP]][
    , chr_name := suppressWarnings(as.integer(chr_name))][
      complete.cases(chr_name)]

  # Equalize SNPs
  haps <- haps[SNP %in% ref[, refsnp_id]] # select only mutual alleles
  ref %<>% .[chorder(chmatch(x = .[, refsnp_id], table = haps[, SNP]))] # put ref at the same order

  # Verify multi-variant alleles (3+) flipping them
  flip <- which(haps[, A1] != ref[, allele_1] & haps[, A2] != ref[, allele_1])
  haps[flip, c('A1', 'A2') := .(gt::flip_strand(A1), gt::flip_strand(A2))]

  # Exclude SNPs with non-matches
  exclude <- which(haps[, A1] != ref[, allele_1] & haps[, A2] != ref[, allele_1])
  if (length(exclude) > 0) haps %<>% .[-exclude]

  ref %<>% .[refsnp_id %in% haps[, SNP]]
  ref %<>% .[chorder(chmatch(x = .[, refsnp_id], table = haps[, SNP]))]

  # Separating info from genotype data
  haps.info <- haps[, 1L:5L]
  haps %<>% .[, -c(1L:5L)]

  # Getting indexes to invert alleles to change SNP IDs and Genotypes
  swap_geno <- c(rbind(colnames(haps)[c(F,T)], colnames(haps)[c(T,F)]))
  snp_inv <- which(haps.info[, A1] != ref[, allele_1])
  a1 <- haps.info[A1 != ref[, allele_1], A1]
  a2 <- haps.info[A1 != ref[, allele_1], A2]

  # Invert SNP Alleles
  haps[snp_inv, colnames(haps) := .SD, .SDcols = swap_geno]
  # Invert SNP IDs
  haps.info[snp_inv, c('A1', 'A2'):= .(a2, a1)]

  # Merging info and genotypes
  cbind2 <- function(...) (setattr(do.call(c,list(...)),"class",c("data.table","data.frame")))
  haps <- cbind2(haps.info, haps)

  # Saving file
  fwrite(x = haps, file = out, quote = F, sep = " ", row.names = F, col.names = F)
}
