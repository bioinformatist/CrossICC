#' Title get centroids of clusters by a set of gene signature and their cluster information
#'
#' @param x expression matrix for calculating centroid
#' @param gene.signature a vector of genenames representing a certain siganture of pathway or biological function.
#' @param cluster a named cluster numbers with gene name as the cluster name and the cluster number as the value
#' @export
#'

cluster.centroid <- function(x, gene.signature, cluster) {
    x.sample.name <- colnames(x)
    # print(x.sample.name)
    cluster <- cluster[x.sample.name]
    # print(cluster) gene.signature <- intersect(rownames(x),gene.signature)
    x.scale <- t(scale(t(x[gene.signature, ])))
    centroids <- c()
    for (i in unique(cluster)) {
        sub.mat <- x.scale[gene.signature, names(which(cluster == i))]
        if (class(sub.mat) == "numeric") {
            centroids <- cbind(centroids, sub.mat)
        } else {
            centroids <- cbind(centroids, apply(sub.mat, 1, median, na.rm = FALSE))
        }
        
    }
    colnames(centroids) <- unique(cluster)
    return(centroids)
}
