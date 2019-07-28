#' To get GSEA-like ranked matrix from CrossICC result.
#'
#' @param x a eSet object or eSet-like matrix.
#' @param gene.signature gene signatures calculated by CrossICC.
#' @param geneset2gene a matrix contains geneset (cluster name) mapping to gene.
#' @param cluster CrossICC returned clusters. Note: Must mapping to x!
#'
#' @return a matrix with samples' eigenvalue in different super clusters.
#' @export
#'
#' @examples
#' CrossICC.object <- CrossICC(demo.platforms, skip.mfs = TRUE, max.iter = 1)
#' CrossICC.ssGSEA <- ssGSEA(x = demo.platforms[[1]], gene.signature = CrossICC.object$gene.signature, geneset2gene = CrossICC.object$unioned.genesets, cluster = CrossICC.object$clusters$clusters[[1]])
ssGSEA <- function(x, gene.signature, geneset2gene, cluster) {
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
    fs.scale <- fs.scale[, names(sort(cluster))]

    return(list(fs, fs.scale))
}
