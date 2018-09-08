#' @export

snp.annot <- function(snpIDs) {

  if(!is.numeric(length(snpIDs)) || length(snpIDs) == 0) stop("Length of SNP vector must be equal or greater than 1.")
  if(!all(grepl("^rs", snpIDs))) stop("All SNPs must be codified as Reference SNP ID (starting with 'rs').")

  snp_annot_function <-  function(snp_ids) {

  url <- paste0("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?", "db=snp&id=", paste(snp_ids, collapse = ','), "&report=DocSet")
  annot <- readLines(url)

  indexes <- grep(pattern = "^rs", annot)
  indexes <- append(x = indexes, values = length(annot)+1)

  total <- length(indexes) -1
  annot.list <- as.list(replicate(total, NULL))

  for(i in 1:total) annot.list[[i]] <- annot[indexes[i]:(indexes[i+1]-1)]

  annot.list <- lapply(annot.list, function(DocSet) {
    index <- grep("GENE=", DocSet, fixed = T)
    gene <- unlist(strsplit(DocSet[index], "GENE=", fixed = T))
    gene <- gene[gene != ""]
  })

  annot.list <- lapply(annot.list, function(x){
    if(is.null(x)) {x <- ""
    } else {
      x
    }
  })

  return(do.call(c, annot.list))
  }

  i <- 0L
  annotation <- character()

  while(i < length(snpIDs)) {

    if((i+250L) > length(snpIDs)) {
      j <- i + abs(length(snpIDs) - i)
    } else {
      j <- i + 250L
    }

    annotation <- append(annotation, snp_annot_function(snpIDs[(i+1):j]) )

    i <- j

  }

  return(annotation)

}
