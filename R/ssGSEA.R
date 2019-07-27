#' Title To get GSEA-like ranked matrix from CrossICC result.
#'
#' @param x a eSet object or eSet-like matrix.
#' @param cluster only used when plot.file is provided. Single super cluster vector of list from CrossICC result list.
#' @param gene.signature gene signatures calculated by CrossICC.
#' @param color to be deprecated.
#' @param geneset2gene a matrix contains geneset (cluster name) mapping to gene.
#'
#' @return a matrix with samples' eigenvalue in different super clusters.
#' @export
#'
#' @examples
#' \donttest{
#' CrossICC.object <- CrossICC(demo.platforms, skip.mfs = TRUE, max.iter = 1)
#' ssGSEA(demo.platforms[[1]], CrossICC.object$gene.signature, CrossICC.object$unioned.genesets)
#' }
ssGSEA <- function(x, gene.signature, geneset2gene) {
    # Our up-stream matrix is already with gene symbols, so provide fake 'genewprobe' for runFAIME. But still need geneset2gene: An
    # one-to-one mapping matrix with two columns, the 1st column is geneset ID/name, and the 2nd is its gene members
    if (is.null(gene.signature)) {
        gene.signature <- geneset2gene[, 2]
    }
    genewprobe <- gene.signature
    names(genewprobe) <- gene.signature
    fs <- runFAIME(x, genewprobe, geneset2gene, weightRank = FALSE)
    fs.scale <- t(scale(t(fs)))
    fs.scale <- replace(fs.scale, fs.scale < -2, -2)
    fs.scale <- replace(fs.scale, fs.scale > 2, 2)
    sn <- gene.signature

    fs.scale <- fs.scale[, names(sort(cluster))]
    # col.col <- c()
    # cluster.member <- unique(sort(cluster))
    # for (i in 1:length(cluster)) {
    #     cluster.names <- names(which(cluster == cluster.member[i]))
    #     col.col <- c(col.col, rep(color[i], length(cluster.names)))
    # }

    return(list(fs, fs.scale))
}
