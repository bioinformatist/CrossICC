#' Title To get GSEA-like ranked matrix from CrossICC result.
#'
#' @param x a eSet object or eSet-like matrix.
#' @param cluster only used when plot.file is provided. Single super cluster vector of list from CrossICC result list.
#' @param gene.signature gene signatures calculated by CrossICC.
#' @param plot.file to be deprecated.
#' @param color to be deprecated.
#' @param geneset2gene a matrix contains geneset (cluster name) mapping to gene.
#'
#' @return a matrix with samples' eigenvalue in different super clusters.
#' @export
#'
#' @examples
#' \donttest{
#' CrossICC.object <- CrossICC(example.matrices, max.iter = 20)
#' ssGSEA(example.matrices$x, CrossICC.object[[1]]$gene.signature, CrossICC.object[[1]]$unioned.genesets)
#' }
ssGSEA <- function(x, gene.signature, geneset2gene, plot.file = FALSE, color = c("purple", "#910202",
  "#035005", "#1273d6", "#ff8207"), cluster) {
  # Our up-stream matrix is already with gene symbols, so provide fake
  # 'genewprobe' for runFAIME. But still need geneset2gene: An one-to-one mapping matrix with two columns,
  # the 1st column is geneset ID/name, and the 2nd is its gene members
  # For testing:
  # gs <- data.frame(geneset = c(rep('set1', 8), rep('set2', 8)), gene = fuck[[1]]$gene.signature)
  # gs <- as.matrix(gs)
  genewprobe <- gene.signature
  names(genewprobe) <- gene.signature
  fs <- runFAIME(x, genewprobe, geneset2gene, weightRank = FALSE)
  fs.scale <- t(scale(t(fs)))
  fs.scale <- replace(fs.scale, fs.scale < -2, -2)
  fs.scale <- replace(fs.scale, fs.scale > 2, 2)
  sn <- gene.signature
  if (plot.file) {
    fs.scale <- fs.scale[, names(sort(cluster))]
    col.col <- c()
    cluster.member <- unique(sort(cluster))
    for (i in 1:length(cluster)) {
      cluster.names <- names(which(cluster == cluster.member[i]))
      col.col <- c(col.col, rep(color[i], length(cluster.names)))
    }
    png(plot.file, width = 1000, height = 800, pointsize = 18)
    if (require("gplots")) {
      heatmap.2(fs.scale, col = greenred(75), ColSideColors = col.col, sepcol = "white",
        sepwidth = c(0.001, 0.001), rowsep = 1:nrow(fs), scale = "none",
        trace = "none", density.info = "none", labRow = NULL, labCol = NULL)
    }
    dev.off()
  }
  return(fs)
}
