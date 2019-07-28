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


# other function for run shiny 

#' Title Get fihser test p value for overlaps
#' @param vec1 a string vector
#' @param vec2 a string vector
#' @param universe total number of all strings that vec1 and vec2 comes from
#' @return  a P value
get_overlap_test_by_fisher <- function(vec1, vec2, universe) {
  mat <- matrix(c(universe - length(union(vec1, vec2)), length(setdiff(vec1, vec2)), length(setdiff(vec1, vec2)), length(intersect(vec1,
                                                                                                                                   vec2))), nrow = 2)
  
  fr <- fisher.test(mat, alternative = "greater")
  return(fr$p.value)
}


#' Title get jaccard index of list factors

#' @param list1 a list that contains a lot vectors
#' @param list2 a list that contains a lot vectors
#' @param universe (Optional) total number of all strings that vec1 and vec2 comes from
#' @return a data frame of Jaccard index or a list contains two dataframe (jaccard index and Fisher's test P value list )
#'
get_jarrad_index_df_fromlist <- function(list1, list2, universe = NULL) {
  # get individual jaccard index
  get_J <- function(vec1, vec2) {
    I <- length(intersect(vec1, vec2))
    S <- I/(length(vec1) + length(vec2) - I)
    return(S)
  }
  
  # get individual jaccard index to a list
  get_J_list <- function(vec1, list.temp) {
    res.list <- unlist(lapply(list.temp, get_J, vec1 = vec1))
    return(res.list)
  }
  # fisher P list
  get_P_list <- function(vec1, list.temp, universe) {
    P.list <- unlist(lapply(list.temp, get_overlap_test_by_fisher, vec1 = vec1, universe = universe))
    return(P.list)
  }
  # get matrix output
  res.m.list <- lapply(list1, get_J_list, list.temp = list2)
  S.df <- data.frame(matrix(unlist(res.m.list), nrow = length(list1), byrow = TRUE), stringsAsFactors = FALSE)
  colnames(S.df) <- names(list2)
  row.names(S.df) <- names(list1)
  if (!is.null(universe)) {
    P.m.list <- lapply(list1, get_P_list, list.temp = list2, universe = universe)
    p.df <- data.frame(matrix(unlist(P.m.list), nrow = length(list1), byrow = TRUE), stringsAsFactors = FALSE)
    colnames(p.df) <- names(list2)
    row.names(p.df) <- names(list1)
    return(list(J.index = S.df, P.fisher = p.df))
  } else {
    return(S.df)
  }
}
# get jaccard index of two annotated dataframe (with two column)
#' Title get jaccard index of list factors

#' @param df1 an annotated data frame with cluster at the seccond column
#' @param df2 an annotated data frame with cluster at the seccond column
#' @param universe (Optional) total number of all strings that vec1 and vec2 comes from
#' @return a data frame of Jaccard index or a list contains two dataframe (jaccard index and Fisher's test P value list )
#'
get_jarrad_index_df_fromDF <- function(df1, df2, universe = NULL) {
  if (is(df1, "integer")) {
    df1 <- data.frame(names(df1), df1)
  }
  if (is(df2, "integer")) {
    df2 <- data.frame(names(df2), df2)
  }
  # get individual jaccard index
  list1 <- tapply(df1[, 1], df1[, 2], list)
  list2 <- tapply(df2[, 1], df2[, 2], list)
  if (is.null(universe)) {
    return(get_jarrad_index_df_fromlist(list1, list2))
  } else {
    return(get_jarrad_index_df_fromlist(list1, list2, universe))
  }
  
}

# get Rand index
#' Title Adjust Rank Index

#' @param df input data frame
#' @param col1 name of interest variable 1 column in df
#' @param col2 name of interest variable 2 column in df
#' @return adjust ARI value
rand.index <- function(df, col1, col2) {
  group1 <- as.numeric(as.factor(df[, col1]))
  group2 <- as.numeric(as.factor(df[, col2]))
  
  x <- abs(sapply(group1, function(x) x - group1))
  x[x > 1] <- 1
  y <- abs(sapply(group2, function(x) x - group2))
  y[y > 1] <- 1
  sg <- sum(abs(x - y))/2
  bc <- choose(dim(x)[1], 2)
  ri <- 1 - sg/bc
  return(ri)
}


# get adjust Rand index
#' Title Adjust Rank Index
#' @param df input data frame
#' @param col1 name of interest variable 1 column in df
#' @param col2 name of interest variable 2 column in df
#' @return adjust ARI value
#'
Cal.ARI <- function(df, col1, col2) {
  x = df[, col1]
  y = df[, col2]
  
  if (length(x) != length(y)) {
    stop("two vectors have different lengths!\n")
  }
  
  cdsum = function(x) {
    y = x * (x - 1)/2
    return(y)
  }
  
  mkMatrix = table(x, y)
  mkMatrixSum = apply(mkMatrix, 1, cdsum)
  matrixSum = sum(mkMatrixSum)
  
  matrixColsum = colSums(mkMatrix)
  matrixRowsum = rowSums(mkMatrix)
  n = length(x)
  
  statRow = sum(cdsum(matrixRowsum))
  statCol = sum(cdsum(matrixColsum))
  
  ARI = (matrixSum - (statRow * statCol)/(cdsum(n)))/(((statRow + statCol)/2) - (statRow * statCol/cdsum(n)))
  
  return(ARI)
}



