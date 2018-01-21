ssGSEA<- function(x, c, plot.file, color=c("purple", "#910202", "#035005","#1273d6","#ff8207")){
  # Our up-stream matrix is already with gene symbols, so provide fake "genewprobe" and "geneset2gene" for runFAIME.
  genewprobe <- rownames(x)
  names(genewprobe) <- rownames(x)
  geneset2gene <- cbind(rownames(x), rownames(x))
  fs <- runFAIME(x, genewprobe, geneset2gene, weightRank = FALSE)
  fs.scale <- t(scale(t(fs)))
  fs.scale <- replace(fs.scale, fs.scale < -2, -2)
  fs.scale <- replace(fs.scale, fs.scale > 2, 2)
  # sn=unique(s[,1])
  # if(exists("plot.file")){
  #   fs.scale=fs.scale[,names(sort(c))]
  #   col.col=c()
  #   cluster.member=unique(sort(c))
  #   for(i in 1:length(cluster.member)){
  #     cluster.names <- names(which(c==cluster.member[i]))
  #     col.col <- c(col.col,rep(color[i],length(cluster.names)))
  #   }
  #   png(plot.file,width=1000,height=800,pointsize=18)
  #   if(require("gplots")){
  #     heatmap.2(fs.scale, Rowv=NA,Colv=NA, col=greenred(75),ColSideColors=col.col,sepcol="white",sepwidth=c(0.001,0.001),rowsep=1:nrow(fs), scale="none",trace="none",density.info="none",labRow=NULL,labCol=NULL)
  #   }
  #   dev.off()
  # }
  return(fs)
}
