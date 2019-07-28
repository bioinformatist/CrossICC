#loading dependencies
suppressMessages(library(cluster))
suppressMessages(library(ggplot2))
suppressMessages(library(ggsci))
suppressMessages(library(tibble))
suppressMessages(library(ggalluvial))
suppressMessages(library(dplyr))
suppressMessages(library(limma))
suppressMessages(library(pheatmap))
# suppressMessages(library(survminer))
#plot heat map from matrix and annotation information
plot_expression_heatmap_with_cluster<-function(df,sample.cluster, genes,cluster_row=FALSE,showRowname=FALSE){
  clustercolor <- list(Cluster = c(Cluster1 = "#E41A1C", Cluster2 = "#377EB8", Cluster3 = "#4DAF4A", Cluster4 = "#984EA3",
                                   Cluster5 = "#FF7F00", Cluster6 = "#FFFF33", Cluster7 = "#A65628", Cluster8 = "#F781BF", Cluster9 = "#999999"))

  plot.matrix <- df
  samplename <- colnames(plot.matrix)
  sample.cluster <- data.frame(Cluster = sample.cluster)
  sample.cluster$Cluster <- paste("Cluster", sample.cluster$Cluster, sep = "")

  sample.cluster <- tibble::rownames_to_column(sample.cluster, var = "sampleid") %>%
    filter(sampleid %in% samplename) %>%
    arrange(Cluster) %>%
    data.frame(row.names = 1)

  plot.matrix <- plot.matrix[, rownames(sample.cluster)]
  #plot heatmap
  xx<-pheatmap(plot.matrix[genes,],
                     scale = 'none',
                     annotation_colors = clustercolor,
                     border_color = NA,
                     cluster_cols = FALSE,
                     cluster_rows = cluster_row,
                     annotation_col = sample.cluster,
                     show_rownames = showRowname,
                     show_colnames = FALSE,
                     colorRampPalette(c("blue", "white", "red"))(100))

  return(xx)

}

# plot sihouttle ------
# sih are a sih function from crossICC object
plot_sihouttle_with_crossICCout <- function(sih){
  # max.sliw<-which.max(max(sih[,3])) + 1
  color.list<-c( "#E41A1C", "#377EB8",  "#4DAF4A", "#984EA3","#FF7F00", "#FFFF33", "#A65628",  "#F781BF",  "#999999")

  plot(sih,col=color.list,border=NA)
}

plot_balanced_heatmap<-function(all.k){
 xx<- pheatmap::pheatmap(all.k,
                     border_color = NA,
                     show_rownames = TRUE,
                     show_colnames = TRUE,
                     colorRampPalette(c("#91CF60", "#FFFFBF","#FC8D59" ))(50))
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



Sankeyplot <- function(df,int.vect1,int.vect2,input.theme="default"){

  var1=df[,int.vect1]
  var2=df[,int.vect2]

  stastdata <- data.frame(table(paste(var1, var2, sep = "-")), stringsAsFactors = FALSE)
  stastdata <- dplyr::mutate(stastdata,
                             cluster1 = strsplit2(Var1, split = "-")[, 1],
                             cluster2 = strsplit2(Var1, split = "-")[, 2])

  g <- ggplot(data = stastdata,
               aes(axis1 = cluster1, axis2 = cluster2,
                   y = Freq)) +
    scale_x_discrete(limits = colnames(stastdata)[2:3], expand = c(.1, .05)) +
    geom_alluvium(aes(fill = cluster2)) + scale_fill_d3(palette = "category20c")+
    geom_stratum() + geom_text(stat = "stratum", label.strata = TRUE)

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

  return(sk)
}


plotSurvival<-function(df,int.var,Time,Status){
  
  plot.df<-df[,which(colnames(df) %in% c(int.var,Time,Status))]
  
  names(plot.df)<-c("feature","time","status")
  
  fit<-survival::survfit( survival::Surv(time, status) ~ feature,
           data = plot.df )
  
  ggsurvplot(fit, data = plot.df, pval = TRUE)
  
}


