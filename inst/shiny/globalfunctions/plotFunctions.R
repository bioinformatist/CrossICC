#loading dependencies
suppressMessages(library(cluster))
#plot heat map from matrix and annotation information
plot_expression_heatmap_with_cluster<-function(df,sample.cluster, genes){
  plot.matrix<-df
  samplename<-colnames(plot.matrix)
  if(class(sample.cluster)=="list"){
    sample.cluster<-sample.cluster[[1]]
  }
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
                     cluster_rows = FALSE,
                     annotation_col = annotation.frame,
                     annotation_colors = list(cluster=color.list),
                     show_colnames = FALSE,
                     colorRampPalette(c("blue", "white", "red"))(100))

}

# plot sihouttle ------
# sih are a sih function from crossICC object
plot_sihouttle_with_crossICCout <- function(sih){
  max.sliw<-which.max(max(sih[,3])) + 1
  color.list<-c()
  colorlength <- 3
  if(length(unique(sih[,1]))>3){
    colorlength <- length(unique(sih[,1]))
    color.list<-brewer.pal(colorlength, "Set2")
  }else{
    colorlength <- length(unique(sih[,1]))
    color.list<-brewer.pal(3, "Set2")[1:colorlength]
  }

  color.list<-brewer.pal(colorlength, "Set2")
  plot(sih,col=color.list[1:max.sliw])
}

plot_balanced_heatmap<-function(all.k){
  pheatmap::pheatmap(all.k,
                     border_color = NA,
                     show_rownames = FALSE,
                     show_colnames = FALSE,
                     colorRampPalette(c("#FC8D59", "#FFFFBF", "#91CF60"))(50))

}


