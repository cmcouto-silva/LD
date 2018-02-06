snp.annot <- function(snpIDs) {

  url <- paste0("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?", "db=snp&id=", paste(snpIDs, collapse = ','), "&report=DocSet")
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

  annotation <- do.call(c, annot.list)
  return(annotation)

}
