cap <- function(word) {
  sapply(word, function(x) {
    vec <- unlist(strsplit(x, split = ""))
    paste0(c(toupper(vec)[1], vec[2:length(vec)]), collapse = "")
  })
}

rm.txt <- function(file.txt) {
  if (any(grep(".txt$", file.txt))) {
    tolower(unlist(strsplit(x = file.txt, split = ".txt$")))
  } else {
    return(file)
  }
}

check.populations <- function(pop1, pop2){
  if (missing(pop1)) stop("You must specify the populations in 'scanhh.list' file.\nExample:\n> load(scanhh.list)\n> pop1 = scanhh_list$population1")
  if (missing(pop2)) stop("You must specify the populations in 'scanhh.list' file.\nExample:\n> load(scanhh.list)\n> pop2 = scanhh_list$population2")
  if(!is.data.frame(pop1) | !is.data.frame(pop2)) stop('pop1 and pop2 must be data frames from scanhh.RData')
}

check.popnames <- function(popname1, popname2){
  if (missing(popname1)) popname1 <- comment(pop1)
  if (missing(popname2)) popname2 <- comment(pop2)
}

check.snplist <- function(snp.list) {

  if (length(snp.list) == 1 && snp.list == "all") {
    snp.list <- list.files(path = "snps", all.files = TRUE)
    snp.list <- grep(pattern = "[[:alnum:]]", x = snp.list, value = TRUE)
    snp.list.data <- unlist(lapply(X = snp.list, FUN = function(x) readLines(paste0("snps/", x))))
  } else {
    snp.list.data <- unlist(lapply(X = snp.list, FUN = function(x) readLines(paste0("snps/", x))))
  }
  return(snp.list.data)
}

check.method <- function(method) {
  if (!method %in% c("bilateral", "unilateral", "both")){
    stop(paste0('"', method,'" is not recognizable. Only the values "bilateral", "unilateral" or "both" are accepted in method\'s parameter.'))
  }
  return(method)
}
