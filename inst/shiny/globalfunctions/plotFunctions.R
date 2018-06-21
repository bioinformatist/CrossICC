#loading dependencies
suppressMessages(library(cluster))
suppressMessages(library(ggplot2))
#plot heat map from matrix and annotation information
plot_expression_heatmap_with_cluster<-function(df,sample.cluster, genes,cluster_row=FALSE,showRowname=FALSE){
  plot.matrix<-df
  samplename<-colnames(plot.matrix)

  annotation.list<-sample.cluster[samplename]

  annotation.list<-sort(annotation.list)
  plot.matrix<-plot.matrix[,names(annotation.list)]
  annotation.frame<-data.frame(cluster=as.factor(annotation.list))
  rownames(annotation.frame)<-names(annotation.list)
  #heatmap colors
  colorlength <- 3
  if(length(unique(annotation.frame[,1]))>3){
    colorlength <- length(unique(annotation.frame[,1]))
    color.list<-brewer.pal(colorlength, "Set2")
  }else{
    colorlength <- length(unique(annotation.frame[,1]))
    color.list<-brewer.pal(3, "Set2")[1:colorlength]
  }

  names(color.list)<-unique(annotation.frame[,1])
  #plot heatmap
  pheatmap::pheatmap(plot.matrix[genes,],
                     scale = 'none',
                     border_color = NA,
                     cluster_cols = FALSE,
                     cluster_rows = cluster_row,
                     annotation_col = annotation.frame,
                    show_rownames = showRowname,
                     show_colnames = FALSE,
                     colorRampPalette(c("blue", "white", "red"))(100))

}

# plot sihouttle ------
# sih are a sih function from crossICC object
plot_sihouttle_with_crossICCout <- function(sih){
  # max.sliw<-which.max(max(sih[,3])) + 1
  color.list<-c()
  colorlength <- 3
  if(length(unique(sih[,1]))>3){
    colorlength <- length(unique(sih[,1]))
    color.list<-brewer.pal(colorlength, "Set2")
  }else{
    colorlength <- length(unique(sih[,1]))
    color.list<-brewer.pal(3, "Set2")[1:colorlength]
  }

  plot(sih,col=color.list,border=NA)
}

plot_balanced_heatmap<-function(all.k){
  pheatmap::pheatmap(all.k,
                     border_color = NA,
                     show_rownames = FALSE,
                     show_colnames = FALSE,
                     colorRampPalette(c("#FC8D59", "#FFFFBF", "#91CF60"))(50))

}

#plot iterator
plot_iter_with_crossICC<-function(crossICC.object){
  df<-crossICC.object
  ggplot() + geom_path(data = df$iter.sig, aes(x = Iteration, y = Signatures), linetype = 2) + geom_vline(xintercept = max(df$iter.sig$Iteration) -1, size = 2.5, color = 'pink')
}

