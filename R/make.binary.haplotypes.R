#' @title Make binary file with haplotypes
#'
#' @description This function imports all haplotypes present in the 'haplotypes' directory,
#' creating a single binary file for loading them all (as requested on subsequent functions).
#'
#' @usage make.binary.haplotypes()
#'
#' @param chrs Numeric vector with chromosome ID(s). Default set to 1:22.
#'
#' @details
#' Haplotype file names must start with "chr", followed by the chromosome ID
#' (represented as Arabic numerals), then ending with either ".haps" or ".sample"
#' extensions. Anything between the "chr[digit]" and the file extension can be written
#' for population identification.
#'
#' \preformatted{Examples: chr1<ANYTHING>.haps and chr1<ANYTHING>.sample
#'           chr2.foo-bar.haps and chr2.foo-bar.sample
#'           chr3_foobar.haps and chr3_foobar.sample
#'           chr4.haps and chr4.sample
#'           etc
#' }
#'
#' @return R binary file named 'haps-sample.RData'.
#' @seealso \code{\link{ldstart}}
#'
#' @examples
#' \dontrun{
#' make.binary.haplotypes() # analyzing all chromosomes (default = 1:22)
#' make.binary.haplotypes(chrs = 1:10) # analyzing chromosomes from 1 to 10
#' make.binary.haplotypes(chrs = c(2,9,10,15,21)) # analyzing specific chromosomes
#' }
#'
#' @export
#' @author Cainã Max Couto-Silva

make.binary.haplotypes <- function(chrs = 1:22) {
  
  # Creating empty lists for storing data
  haps <- list()
  sample <- list()
  haps.info <- list()
  sample.info <- list()
  
  # Acquiring file names
  files <- character()
  ids <- character()
  
  if (missing(chrs)) chrs <- 1:22
  
  for (i in 1:length(chrs)) {
    files[i] <- gtools::mixedsort(list.files(path = "haplotypes", pattern = paste0("^chr", chrs[i], "[^0-9+]*\\.haps$")))
    ids[i] <- gtools::mixedsort(list.files(path = "haplotypes", pattern = paste0("^chr", chrs[i], "[^0-9+]*\\.sample$")))
  }
  
  # Specifying chromosomes to be analyzed
  chrs <- as.integer(chrs)
  cat(paste("  Importing", length(chrs), "haplotypes...\n\n"))
  
  # Reading haps and sample files
  for (i in 1:length(files)) {
    cat(paste("\tReading"), files[i], "&", ids[i], "files \n")
    
    # haps (only genome)
    haps[[i]] <- data.table::fread(input = paste0("haplotypes/", files[i]), header = F, stringsAsFactors = F, drop = 1:5)
    names(haps)[[i]] <- files[i]
    # haps (only data info)
    haps.info[[i]] <- data.table::fread(input = paste0("haplotypes/", files[i]), header = F, stringsAsFactors = F, select = 1:5)
    names(haps.info)[[i]] <- files[i]
    # sample (only data)
    sample[[i]] <- data.table::fread(input = paste0("haplotypes/", ids[i]), header = F, stringsAsFactors = F)[-c(1:2), ]
    names(sample)[[i]] <- ids[i]
    # sample (header)
    sample.info[[i]] <- data.table::fread(input = paste0("haplotypes/", ids[i]), header = F, stringsAsFactors = F)[c(1:2), ]
    names(sample.info)[[i]] <- ids[i]
  }
  
  haps <- lapply(haps, as.data.frame)
  haps.info <- lapply(haps.info, as.data.frame)
  sample <- lapply(sample, as.data.frame)
  sample.info <- lapply(sample.info, as.data.frame)
  
  # Saving dataset
  cat("\n  Saving R binary file...")
  save(haps, haps.info, sample, sample.info, file = "haplotypes/haps-sample.RData")
  
  cat(paste("\n\n ", length(chrs), "haplotypes have been included\n"))
  cat("  haps-sample.RData file successfully created in \"haplotypes\" directory")
  
}
