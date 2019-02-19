#' @title Computes iHS Statistics (standardized IHH)
#'
#' @description This function computes iHS Statistics (standardized IHH) for all target populations and SNPs.
#'
#' @usage ihs(snp.list = "all", filter = 2, annot = T, write.xls = "both", plot = T)
#'
#' @param snp.list Character vector. It could be a scalar named "all" (default) for loading all SNPs
#' present in 'snps' directory, or a character vector with the exact names from the target-SNP files
#' (also supposed to be in 'populations' directory).
#' @param filter Numeric value for filtering significant iHS scores (default = 2).
#' @param annot Logical scalar. If TRUE (default), SNP annotation in xls files will be written with
#' correspondent SNP's genes, as well as ancestral and derived allele present in the dataset.
#' FALSE turns off this feature.
#' @param write.xls Scalar character. It could be set to "all.snps" for annotating all SNPs included
#' in this analysis, and "ss.snps" for annotating only SNPs statistically significant
#' (defined by 'filter' argument), or 'both' (default) for both annotation outputs.
#' @param plot Logical scalar. If TRUE (default), then distribution plots are generated and saved
#' as png figures.
#'
#' @details
#' All populations included in 'scanhh.RData' will be analyzed. This function creates a new directory in
#' 'rehh_out' folder, named 'ihs', with two sub-directories: 'ihs/graphics' for plots, and 'ihs/tables'
#' for csv and xls files.
#'
#' The default value for 'filter' parameter has been set to 2 accordingly the standard value described on
#' original paper. Two simple files (.csv) with iHS results will be generated for all included populations:
#' one with all SNPs included on analysis, and another with only statistically significant SNPs.
#'
#' The 'annot' parameter turns on or off annotation procedure (default is TRUE, meaning 'on'), while argument
#' passed in 'write.xls' parameter tells which SNPs must be included on analysis ("all.snps", "ss.snps"
#' or "both"). Plot argument enables the generation of two plots/figures per population.
#'
#' Except by annotation argument ('annot = T'), this function usually runs very fast. However,
#' depending on how many target-SNPs are present in the dataset, the annotation velocity may vary considerably.
#'
#' @return iHS files (csv or xls with SNP annotation) and iHS plots.
#' @seealso \code{\link{make.scanhh}}
#'
#' @examples
#' \dontrun{ihs()}
#'
#' \dontrun{
#' ihs(snp.list = 'all', filter = 2, annot = T, write.xls = 'ss.snps', plot = F)
#' }
#'
#' @export

ihs <- function(snp.list = "all", filter = 2, annot = T, write.xls = "both", plot = T, plot.format =  "png", freqbin = 0.025, minmaf = 0.05) {
  
  if(!plot %in% c(T,F))
    stop("plot must be TRUE or FALSE.")
  
  if(!annot %in% c(T,F))
    stop("annot must be TRUE or FALSE.")
  
  if(!write.xls %in% c("all.snps","ss.snps","both"))
    stop('write.xls must be "all.snps", "ss.snps", or "both".')
  
  if (!plot.format %in% c("png", "jpeg", "tiff", "pdf", "svg", "ps", "win"))
    stop('File type not supported! Supporterd formats: "png", "jpeg", "tiff", "pdf", "svg", "ps", and "win".')
  
  # Loading Target-SNPs
  if (length(snp.list) == 1 && snp.list == "all") {
    snp.list <- list.files(path = "snps/", all.files = TRUE)
    snp.list <- grep(pattern = "[[:alnum:]]", x = snp.list, value = TRUE)
    snp.list.data <- unlist(lapply(X = snp.list, FUN = function(x) readLines(paste0("snps/", x))))
  } else {
    snp.list.data <- unlist(lapply(X = snp.list, FUN = function(x) readLines(paste0("snps/", x))))
  }
  
  
  # Checking if there are any target-SNP on dataset
  dataSNPs <- row.names(scanhh.list[[1]])
  available_SNPs <- intersect(dataSNPs, snp.list.data)
  if ( length(available_SNPs) == 0L ) stop("There are no target-SNPs in the dataset.")
  
  # Loading R binary file for running annotation
  if (annot == TRUE) {
    cat(" Loading haps-sample.RData file for acquiring allele state information... \n")
    load("haplotypes/haps-sample.RData")
    haps.info.alleles <- data.table::rbindlist(haps.info)
  }
  
  # Computing iHS Statistics as described by Sabeti et al. 2007
  cat("\n Computing iHS Statistics... \n")
  ld.ihs <- lapply(scanhh.list, rehh::ihh2ihs, freqbin = freqbin, minmaf = minmaf)
  
  # All SNPs
  ihs.snp.list <- lapply(ld.ihs, function(x) {
    ihs.snp.list <- x$iHS[row.names(x$iHS) %in% snp.list.data, ]
    ihs.snp.list <- ihs.snp.list[complete.cases(ihs.snp.list), ]
  })
  
  # Statistically significant SNPs
  ihs.snp.list.stat <- lapply(ihs.snp.list, function(x) x[x$iHS >= filter | x$iHS <= -filter,])
  
  # Dealing with pop names
  for (i in 1:length(ihs.snp.list)){
    comment(ihs.snp.list[[i]]) <- names(ld.ihs)[i]
    comment(ihs.snp.list.stat[[i]]) <- names(ld.ihs)[i]
  }
  
  # Creating directories
  dir.create(path = "rehh_out/ihs", showWarnings = F)
  dir.create(path = "rehh_out/ihs/tables", showWarnings = F)
  
  # Writing iHS files
  invisible(lapply(ihs.snp.list, function(fx) {
    data.table::fwrite(x = fx, file =  paste0("./rehh_out/ihs/tables/", comment(fx), ".csv"), quote = F, row.names = T)}))
  
  invisible(lapply(ihs.snp.list.stat, function(fx) {
    data.table::fwrite(x = fx, file =  paste0("./rehh_out/ihs/tables/", comment(fx), ".stat.csv"), quote = F, row.names = T) }))
  
  # Plotting iHS distribution and qqplot
  if (plot) {
    cat(" Plotting iHS results... \n")
    dir.create(path = "rehh_out/ihs/graphics", showWarnings = FALSE)
    
    for (i in 1:length(ld.ihs)) {
      distribplot(ihs.data = ld.ihs[[i]]$iHS[, 3], main = names(ld.ihs)[i], xlab = "", 
                  file.name = paste0("rehh_out/ihs/graphics/", names(ld.ihs)[i]), 
                  file.type = plot.format)
    }
  }
  
  # ANNOTATION
  
  if (annot == TRUE) {
    
    ## Annotation for results regardless iHS score
    if (write.xls == "all.snps" || write.xls == "both") {
      cat("\n Annotation for all SNPs under analysis... \n")
      
      # Checking if there are too much SNPs for annotation
      if (any(unlist(lapply(ihs.snp.list, nrow)) > 10000)) {
        run.annot.all <- FALSE
      } else {
        run.annot.all <- TRUE
      }
      
      # Running annotation
      if(run.annot.all) {
        ihs.snp.list.ncbi.xls <- lapply(ihs.snp.list, function(fx) {
          if(nrow(fx) != 0L) {
            alleles <- haps.info.alleles[haps.info.alleles$V2 %in% rownames(fx), ]
            ihs.snp.list.ncbi.xls <- cbind(Population = comment(fx), CHR = fx$CHR, SNP = rownames(fx),
                                           Allele_A = alleles$V4, Allele_D = alleles$V5, GENE = snp.annot(rownames(fx)), fx[2:4])
          } else {
            cat("\tNo target-SNPs have been found in", comment(fx), "dataset.\n")
          }
        })
        
        # Merge data from all populations
        ihs.snp.list.ncbi.xls <- data.table::rbindlist(ihs.snp.list.ncbi.xls)
        
        if (nrow(ihs.snp.list.ncbi.xls) > 0L) {
          WriteXLS::WriteXLS(x = ihs.snp.list.ncbi.xls, ExcelFileName = "rehh_out/ihs/tables/ihs.all.snps.xls",
                             AdjWidth = TRUE, BoldHeaderRow = TRUE)
          cat("\tFile \"ihs.all.snps.xls\" successfully saved into 'rehh_out/ihs/tables/' directory \n\n")
        }
        
      } else {
        cat(" - Annotation has been skipped! Too much SNPs for annotation.\n")
      }
    }
    
    ## Annotation for results with significative iHS scores
    if (write.xls == "ss.snps" || write.xls == "both") {
      cat(paste0("\n Annotation for all SNPs with iHS score equal, less, or bigger than '", filter, "'\n"))
      
      # Checking if there are too much SNPs for annotation
      if (any(unlist(lapply(ihs.snp.list.stat, nrow)) > 10000)) {
        run.annot.stat <- FALSE
      } else {
        run.annot.stat <- TRUE
      }
      
      # Running annotation
      if(run.annot.stat) {
        ihs.snp.list.ncbi.stat.xls <- lapply(ihs.snp.list.stat, function(fx) {
          if(nrow(fx) != 0L) {
            alleles <- haps.info.alleles[haps.info.alleles$V2 %in% rownames(fx), ]
            ihs.snp.list.ncbi.stat.xls <- cbind(Population = comment(fx), CHR = fx$CHR, SNP = rownames(fx),
                                                Allele_A = alleles$V4, Allele_D = alleles$V5, GENE = snp.annot(rownames(fx)), fx[2:4])
          } else {
            cat("\tNo statistically significant SNPs have been found in", comment(fx), "dataset.\n")
          }
        })
        
        # Merge data from all populations
        ihs.snp.list.ncbi.stat.xls <- data.table::rbindlist(ihs.snp.list.ncbi.stat.xls)
        
        if (nrow(ihs.snp.list.ncbi.stat.xls) > 0L) {
          WriteXLS::WriteXLS(x = ihs.snp.list.ncbi.stat.xls, ExcelFileName = "rehh_out/ihs/tables/ihs.ss.snps.xls",
                             AdjWidth = TRUE, BoldHeaderRow = TRUE)
          cat("\tFile \"ihs.ss.snps.xls\" successfully saved into 'rehh_out/ihs/tables/' directory \n\n")
        }
        
      } else {
        cat(" - Annotation has been skipped (stat)! Too much SNPs for annotation.\n")
      }
    } 
  } 
}
