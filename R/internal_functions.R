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

