#' Title
#'
#' @param x
#' @param cluster
#' @param gene.signature
#' @param plot.file
#' @param color
#'
#' @return
#' @export
#'
#' @examples
ssGSEA<- function(x, cluster, gene.signature, plot.file, color=c("purple", "#910202", "#035005","#1273d6","#ff8207")){
  # Our up-stream matrix is already with gene symbols, so provide fake "genewprobe" and "geneset2gene" for runFAIME.
  genewprobe <- gene.signature
  names(genewprobe) <- gene.signature
  geneset2gene <- cbind(gene.signature, gene.signature)
  fs <- runFAIME(x, genewprobe, geneset2gene, weightRank = FALSE)
  fs.scale <- t(scale(t(fs)))
  fs.scale <- replace(fs.scale, fs.scale < -2, -2)
  fs.scale <- replace(fs.scale, fs.scale > 2, 2)
  sn <- gene.signature
  if(exists("plot.file")){
    fs.scale=fs.scale[,names(sort(cluster))]
    col.col=c()
    cluster.member=unique(sort(cluster))
    for(i in 1:length(cluster)){
      cluster.names <- names(which(cluster==cluster.member[i]))
      col.col <- c(col.col,rep(color[i],length(cluster.names)))
    }
    png(plot.file,width=1000,height=800,pointsize=18)
    if(require("gplots")){
      heatmap.2(fs.scale, col=greenred(75),ColSideColors=col.col,sepcol="white",sepwidth=c(0.001,0.001),rowsep=1:nrow(fs), scale="none",trace="none",density.info="none",labRow=NULL,labCol=NULL)
    }
    dev.off()
  }
  return(fs)
}
