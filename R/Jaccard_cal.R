#' Title Get fihser test p value for overlaps
#' @param vec1 a string vector
#' @param vec2 a string vector
#' @param universe total number of all strings that vec1 and vec2 comes from
#' @return  a P value
#' @export
get_overlap_test_by_fisher<-function(vec1,vec2,universe){
  mat <- matrix(
    c(
      universe - length(union(vec1, vec2)),
      length(setdiff(vec1, vec2)),
      length(setdiff(vec1, vec2)),
      length(intersect(vec1, vec2))
    ),
    nrow=2
  )

  fr <- fisher.test(mat, alternative="greater")
  return(fr$p.value)
}


#' Title get jaccard index of list factors

#' @param list1 a list that contains a lot vectors
#' @param list2 a list that contains a lot vectors
#' @param universe (Optional) total number of all strings that vec1 and vec2 comes from
#' @return a data frame of Jaccard index or a list contains two dataframe (jaccard index and Fisher's test P value list )
#' @export
#'
get_jarrad_index_df_fromlist<-function(list1,list2,universe=NULL){
  # get individual jaccard index
  get_J<-function(vec1,vec2){
    I <- length(intersect(vec1,vec2))
    S <- I/(length(vec1)+length(vec2)-I)
    return(S)
  }

  #get individual jaccard index to a list
  get_J_list<-function(vec1,list.temp){
    res.list<-unlist(lapply(list.temp, get_J,vec1=vec1))
    return(res.list)
  }
  #fisher P list
  get_P_list<-function(vec1,list.temp,universe){
    P.list<-unlist(lapply(list.temp, get_overlap_test_by_fisher,vec1=vec1,universe=universe))
    return(P.list)
  }
  #get matrxi output
  res.m.list<-lapply(list1, get_J_list,list.temp=list2)
  S.df<-data.frame(matrix(unlist(res.m.list), nrow=length(list1), byrow=T),stringsAsFactors=FALSE)
  colnames(S.df)<-names(list2)
  row.names(S.df)<-names(list1)
  if(!is.null(universe)){
    P.m.list<-lapply(list1, get_P_list,list.temp=list2,universe=universe)
    p.df<-data.frame(matrix(unlist(P.m.list), nrow=length(list1), byrow=T),stringsAsFactors=FALSE)
    colnames(p.df)<-names(list2)
    row.names(p.df)<-names(list1)
    return(list(J.index=S.df,P.fisher=p.df))
  }else{
    return(S.df)
  }
}
#get jaccard index of two annotated dataframe (with two column)
#' Title get jaccard index of list factors

#' @param df1 an annotated data frame with cluster at the seccond column
#' @param df2 an annotated data frame with cluster at the seccond column
#' @param universe (Optional) total number of all strings that vec1 and vec2 comes from
#' @return a data frame of Jaccard index or a list contains two dataframe (jaccard index and Fisher's test P value list )
#' @export
#'
get_jarrad_index_df_fromDF<-function(df1,df2,universe=NULL){
  if(class(df1)=="integer"){
    df1<-data.frame(names(df1),df1)
  }
  if(class(df2)=="integer"){
    df2<-data.frame(names(df2),df2)
  }
  # get individual jaccard index
  list1<-tapply(df1[,1],df1[,2],list)
  list2<-tapply(df2[,1],df2[,2],list)
  if(is.null(universe)){
    return(get_jarrad_index_df_fromlist(list1,list2))
  }else{
    return(get_jarrad_index_df_fromlist(list1,list2,universe))
  }

}
