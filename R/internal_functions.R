cap <- function(word) {
  sapply(word, function(x) {
    vec <- unlist(strsplit(x, split = ""))
    paste0(c(toupper(vec)[1], vec[2:length(vec)]), collapse = "")
  })
}
