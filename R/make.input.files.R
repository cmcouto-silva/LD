#' @title Make input files for rehh package
#'
#' @description This function firstly extracts data according to user requirements, then convert .haps
#' in .thap files, updating their allele codification and creating a unique map (.inp) file for all
#' chromosomes and populations aimed by the user (files required as input in rehh package).
#'
#' @usage make.input.files(populations = "all", chrs = 1:22, haps.subset = F, output)
#'
#' @param populations Character vector with target population's file names (assumed to be in
#' 'populations' directory). If populations = "all" (default), then all files (populations)
#' will be included on analysis.
#' @param chrs Numeric vector with chromosome ID(s).
#' @param haps.subset Logical scalar. If TRUE, new .haps & .sample files will be generated
#' to the populations under analysis (default = FALSE).
#' @param output A character vector containing the desired population names as output.
#'
#' @details
#' R binary file acquired by make.binary.haplotypes.R() function is used at this function,
#' therefore essentially necessary ('haplotypes/haps-sample.RData').
#'
#' By default, all populations files in 'populations' directory are included in this function,
#' although users can also specify the populations in a vector with the exact file names
#' (without providing the path, which is supposed to be in 'populations' directory).
#'
#' The 'haps.subset' parameter enables the user to extract data from specific populations
#' as they are originally (.haps and .sample formats).
#'
#' Be careful when using 'output' parameter while using 'populations' parameter's default argument,
#' once you would not have control of the file names' order. If using the 'output' parameter,
#' you'd rather specify the populations' names as well. Default values to 'output' parameter
#' are automatically set according to acquired/loaded population's file names so that
#' there are no issues to worry about.
#'
#' @return Input files for rehh package: .thap, .sample, and .inp .
#' @seealso \code{\link{make.binary.haplotypes}}
#'
#' @examples
#' \dontrun{
#' make.input.files()
#' make.input.files(populations = c('Afr.txt', 'EUR.txt'), chrs = 1:10, haps.subset = T,
#' output = c('Africans','Europeans'))
#' }
#'
#' @export

make.input.files <- function(populations = "all", chrs = 1:22, haps.subset = F, output) {

    # Population's files to import
    if (populations[1] == "all") {
        populations <- list.files(path = "populations/", all.files = TRUE)
        populations <- grep(pattern = "[[:alnum:]]", x = populations, value = TRUE)
    }

    populations.data <- lapply(X = populations, FUN = function(x) readLines(paste0("populations/", x)))

    populations <- tolower(populations)

    if (any(grep(".txt$", populations))) {
        populations <- unlist(strsplit(x = populations, split = ".txt$"))
    }

    if (missing(output)) {
        output <- populations
    }

    if (length(populations) != length(output)) {
        stop("Output length must be equal to the number of populations analyzed.")
    }

    names(populations.data) <- output
    rm(populations)

    cat(" Loading haps-sample.RData file... ")
    load("./haplotypes/haps-sample.RData")
    cat("OK\n\n")
    cat(" Extracting data for specified populations... \n\n")

    for (n in 1:length(populations.data)) {

        #### Extracting populations
        phased.dataset <- mapply(function(haps, haps.info, sample, sample.info) {

            # Assessing individuals' position in `.sample` and respective columns in `.haps` files for Amazonians
            pos_sample = which(sample$V2 %in% populations.data[[n]])
            pos_in_haps_a1 = pos_sample * 2 - 1
            pos_in_haps_a2 = pos_sample * 2
            pos_in_haps_both_alleles = c(rbind(pos_in_haps_a1, pos_in_haps_a2))
            # Writing filtered haps & sample files
            new.haps = cbind(haps.info, haps[, pos_in_haps_both_alleles])
            new.sample = rbind(sample.info, sample[pos_sample, ])
            list(new.haps, new.sample)

        }, haps = haps, haps.info = haps.info, sample = sample, sample.info = sample.info)

        # splitting haps & sample files in two lists
        new.haps <- phased.dataset[c(T, F)]
        new.sample <- phased.dataset[c(F, T)]
        populations.data[[n]] <- list(new.haps, new.sample)
        cat(paste0("\tSuccessfully extracted data for \"", output[n], "\"\n"))
    }

    # Saving new haps & sample files
    if (haps.subset == T) {
        cat(paste("\n Saving new .haps & .sample files...\n\n"))
        for (n in 1:length(populations.data)) {
            for (i in 1:length(chrs)) {
                write.table(x = new.haps[[i]], file = paste0("rehh_in/chr", chrs[i], "_", output[n], ".haps"), quote = F,
                  row.names = F, col.names = F)
                write.table(x = new.sample[[i]], file = paste0("rehh_in/chr", chrs[i], "_", output[n], ".sample"), quote = F,
                  row.names = F, col.names = F)
            }
        }
        for (i in 1:length(populations.data)) cat(paste0("\tNew haps/sample files created for \"", output[i], "\" population\n"))
    }

    # Converting and saving new haps files to .thap

    cat("\n Converting and saving .thap and .sample files... \n\n")
    for (n in 1:length(populations.data)) {
        for (i in 1:length(chrs)) {
            write.table(populations.data[[n]][[1]][[i]][, -c(1:5)], file = paste0("rehh_in/chr", chrs[i], "_", output[n],
                ".rehh.thap"), quote = F, row.names = F, col.names = F)
            write.table(populations.data[[n]][[2]][[i]], file = paste0("rehh_in/chr", chrs[i], "_", output[n], ".rehh.sample"),
                quote = F, row.names = F, col.names = F)
        }
        cat(paste0("\tNew thap/sample files created for \"", output[n], "\" population\n"))
    }

    cat("\n Updating code alleles from 0/1 to 1/2 as requestes by rehh package...\n")
    OS <- Sys.info()
    if (OS["sysname"] == "Linux" | OS["sysname"] == "MAC") {
        for (n in 1:length(populations.data)) {
            for (i in 1:length(chrs)) {
                system(paste0("tr 01 12 < rehh_in/chr", chrs[i], "_", output[n], ".rehh.thap > rehh_in/chr", chrs[i],
                  ".thap && mv rehh_in/chr", chrs[i], ".thap rehh_in/chr", chrs[i], "_", output[n], ".rehh.thap; done"))
            }
        }
    }

    cat(paste0(" Conversion completed!\n\n"))

    # Making the map file
    inp <- do.call("rbind", haps.info)[, c(2, 1, 3:5)]
    write.table(inp, file = "rehh_in/mapfileR.inp", quote = F, row.names = F, col.names = F)

    cat(" Map file created,\n\n")

    cat(paste0("\t", (length(populations.data) * (length(chrs) * 2)) + 1, " files have been written into \"input\" directory:\n"))
    cat(paste0("\t - ", length(populations.data) * (length(chrs)), " files .thap \n"))
    cat(paste0("\t - ", length(populations.data) * (length(chrs)), " files .sample\n"))
    cat(paste0("\t - ", 1, " map file (.inp)\n"))

}
