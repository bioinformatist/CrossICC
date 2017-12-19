cor.cc <- function(x, y, z, cc, k){
  platforms <- list(x = x, y = y, z = z)
  lapply(names(k), function(p) lapply(1:k[[p]], function(q) apply(platforms[[p]][, names(which(cc[[p]][[k[[p]]]]$consensusClass == q)), with = FALSE], 1, median)))
  shit <- as.data.table(unlist(fuck, recursive = FALSE))
}
