#' @title Create R binary file with EHH-based Statistics
#'
#' @usage make.input.files(populations = "all", chrs = 1:22, haps.subset = F, output)
#'
#' @param populations Character vector with target population's file names (assumed to be in
#' 'populations' directory). If populations = 'all' (default), then all files (populations)
#' will be included on analysis.
#' @param chrs Numeric vector with chromosome ID(s).
#' @param output Character vector containing desired population names as output.
#'
#' @description
#' This function computes EHH-based Statistics which will be used in subsequent functions,
#' resulting in an R binary file containing all populations/chromosomes data analyzed.
#'
#' @details
#' The files acquired after running make.input.files() function will be used at this function.
#' Therefore, they are absolutely necessary (.inp and .thap files in 'rehh_in' directory).
#'
#' By default, all populations files in 'populations' directory are included in this function,
#' although users can also specify the populations in a vector with the exact file names
#' (without providing the path, which is supposed to be in 'populations' directory).
#'
#' Be careful when using 'output' parameter while using 'populations' parameter's default argument,
#' once you would not have control of the file names' order. If using the 'output' parameter,
#' you'd rather specify the populations' names as well. Default values to 'output' parameter
#' are automatically set according to acquired/loaded population's file names so that
#' there are no issues to worry about.
#'
#' @return R binary file named 'scanhh.RData'.
#' @seealso \code{\link{make.input.files}}
#'
#' @examples
#' \dontrun{make.input.files()}
#'
#' \dontrun{
#' make.input.files(populations = c('Afr.txt', 'EUR.txt'),
#'                  chrs = 1:10, haps.subset = T,
#'                  output = c('Africans','Europeans'))
#' }
#'
#' @export
#'

make.scanhh <- function(populations = "all", chrs = 1:22, output) {

    if (populations[1] == "all") {
        populations <- list.files(path = "populations/", all.files = TRUE)
        populations <- grep(pattern = "[[:alnum:]]", x = populations, value = TRUE)
    }

    populations <- tolower(populations)
    populations <- unlist(strsplit(x = populations, split = ".txt$"))

    if (missing(output)) {
        output <- LD::cap(populations)
    }

    if (length(populations) != length(output)) {
        stop("Output length must be equal to the number of populations analyzed.")
    }

    ## Loading files as R objects for each population ##

    scanhh.list <- as.list(replicate(length(populations), NULL))
    names(scanhh.list) <- output

    for (n in 1:length(populations)) {
        cat(paste0("\n# Population under analysis: ", output[n], "\n\n"))
        for (i in 1:length(chrs)) {
            haplohh <- rehh::data2haplohh(hap_file = paste0("./rehh_in/chr", chrs[i], "_", names(output)[n], ".rehh.thap"),
                map_file <- "./rehh_in/mapfileR.inp", chr.name = chrs[i], haplotype.in.columns = TRUE)
            results <- rehh::scan_hh(haplohh)
            if (chrs[i] == min(chrs)) {
                scanhh.list[[n]] <- results
            } else {
                scanhh.list[[n]] <- rbind(scanhh.list[[n]], results)
                scanhh.list[[n]] <- scanhh.list[[n]]
            }
        }
        comment(scanhh.list[[n]]) <- names(scanhh.list)[n]
    }
    cat(paste0("\n R binary file 'scanhh.RData' saved into \"rehh_in\" directory\n\n"))
    save(scanhh.list, file = "scanhh.RData")
}
