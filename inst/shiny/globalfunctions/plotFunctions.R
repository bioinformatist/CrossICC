#loading dependencies
suppressMessages(library(cluster))
suppressMessages(library(ggplot2))
suppressMessages(library(ggsci))
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
  xx<-pheatmap::pheatmap(plot.matrix[genes,],
                     scale = 'none',
                     border_color = NA,
                     cluster_cols = FALSE,
                     cluster_rows = cluster_row,
                     annotation_col = annotation.frame,
                    show_rownames = showRowname,
                     show_colnames = FALSE,
                     colorRampPalette(c("blue", "white", "red"))(100))
  return(xx)

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
 xx<- pheatmap::pheatmap(all.k,
                     border_color = NA,
                     show_rownames = FALSE,
                     show_colnames = FALSE,
                     colorRampPalette(c("#FC8D59", "#FFFFBF", "#91CF60"))(50))
 return(xx)

}

#plot iterator
plot_iter_with_crossICC<-function(crossICC.object){
  df<-crossICC.object
  ggplot() + geom_path(data = df$iter.sig, aes(x = Iteration, y = Signatures), linetype = 2) + geom_vline(xintercept = max(df$iter.sig$Iteration) -1, size = 2.5, color = 'pink')
}

# ggplot2 plottheme
#plot theme
#axis.x without rotate
plotDefaultTheme2=theme(
  panel.background = element_rect(fill = "white", color = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(angle = 0, color = "black",size = 15),
  axis.text.y = element_text(angle = 0, color = "black",size = 20),
  axis.title = element_text(face = "bold", color = "black", size = 20),
  legend.title = element_text(face = "bold", color = "black", size = 20),
  legend.text = element_text(color = "black", size = 20)
)


plotStackBarplot<-function(df,int.vect1,int.vect2,input.theme){


  plotdata<-data.frame(X=as.factor(df[,int.vect1]),Y=as.factor(df[,int.vect2]))

  g<-ggplot(plotdata,aes(X))+geom_bar(aes(fill = Y),position = "fill")+xlab(int.vect1)+ylab(int.vect2)+scale_fill_aaas()

  if(input.theme=="default"){
    g=g+theme(legend.position = "none")+ plotDefaultTheme2
  }else if(input.theme=="Tufte"){
    g=g+geom_rangeframe() + theme_tufte()
  }else if(input.theme=="Economist"){
    g=g+ theme_economist()+ scale_colour_economist()
  }else if(input.theme=="Solarized"){
    g=g+ theme_solarized()+ scale_colour_solarized("blue")
  }else if(input.theme=="Stata"){
    g=g+ theme_stata() + scale_colour_stata()
  }else if(input.theme=="Excel 2003"){
    g=g+ theme_excel() + scale_colour_excel()
  }else if(input.theme=="Inverse Gray"){
    g=g+ theme_igray()
  }else if(input.theme=="Fivethirtyeight"){
    g=g+scale_color_fivethirtyeight()+ theme_fivethirtyeight()
  }else if(input.theme=="Tableau"){
    g=g+theme_igray()+ scale_colour_tableau()
  }else if(input.theme=="Stephen"){
    g=g+theme_few()+ scale_colour_few()
  }else if(input.theme=="Wall Street"){
    g=g+theme_wsj()+ scale_colour_wsj("colors6", "")
  }else if(input.theme=="GDocs"){
    g=g+theme_gdocs()+ scale_color_gdocs()
  }else if(input.theme=="Calc"){
    g=g+theme_calc()+ scale_color_calc()
  }else if(input.theme=="Pander"){
    g=g+theme_pander()+ scale_colour_pander()
  }else if(input.theme=="Highcharts"){
    g=g+theme_hc()+ scale_colour_hc()
  }
  return(g)
}


