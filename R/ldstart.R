#' @title Create required (sub)directories for LD analysis
#'
#' @description This function creates directories and subdirectories that will be used
#' for accessing input files, as well as for delivering output files from LD analysis.
#'
#' @usage ldstart()
#'
#' @details
#' Running this function generates the following folders at current directory:
#' 'haplotypes', 'populations', 'snps', 'rehh_in', and 'rehh_out'. Users must place
#' their data inside the respective directories, as explained below:
#'
#' 'haplotypes': The whole dataset of haplotypes (extensions .hap and .sample) split by chromosomes,
#' obtained by Shapeit Software (\url{http://www.shapeit.fr}).
#'
#' 'populations': Text files from each population to be analyzed (one individual per line).
#' File names must correspond to the names from populations (with or without ".txt" extension).
#'
#' 'snps': Text files with a set of SNP ID's (e.g. rs1042522), with or without ".txt" extension
#' (a single ID per line).
#'
#' @return A set of required directories.
#'
#' @examples
#' \dontrun{ldstart() }
#'
#' @export
#'

ldstart <- function() {

    dir.names <- c("haplotypes", "populations", "snps", "rehh_in", "rehh_out")

    invisible(sapply(dir.names, function(mkdir) {
        if (!dir.exists(mkdir)) {
            dir.create(paste0(mkdir), showWarnings = FALSE)
            cat(paste0("   Directory ", "\"", mkdir, "\" successfully created.\n"))
        } else {
            message(paste0("   Warning: directory ", "\"", mkdir, "\" already exists."))
        }
    }))
}
