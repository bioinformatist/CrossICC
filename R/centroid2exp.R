#' Title To predict Samples cluter based on their centroid
#'
#' @param centroid a centroid calculated by crossICC
#' @param vd a expression dataframe to predict clusters, with featurer in rows and samples in colomns,(no NA permitted )
#' @export
#'
#' @examples
#' \donttest{
#'crossICC.object <- CrossICC(example.matrices, max.iter = 5, use.shiny = FALSE)
#'train.centroid<-cluster.centroid(crossICC.object[[2]]$platforms[[1]],crossICC.object[[2]]$gene.signature,crossICC.object[[2]]$clusters$clusters)
#'data(test.data)
#'predict.result<-centroid2exp(train.centroid,test.data)
#' }
centroid2exp <- function(centroid,vd){

  #do nornmalization as the same as in training dat aset
  vd <- m.f.s(list(vd), fdr.cutoff = 1, filter.cutoff = 1, perform.mad = FALSE)[[2]][[1]]
  gene.sig=intersect(rownames(centroid),rownames(vd[complete.cases(vd),]))
  if(!isTRUE(all.equal(sort(rownames(centroid)), sort(gene.sig)))){
    stop("missing prediction features in your expression dataset\n p") # to valide data without missing data
  }
  vd=t(scale(t(vd[gene.sig,])))
  centroid=centroid[gene.sig,]
  vclass<- c()
  vcor <- c()
  for(i in 1:ncol(vd)){
    d=vd[,i]
    c.cor=c()
    pv=c()
    for(j in colnames(centroid)){
      centroidj=centroid[,j]
      corj=cor.test(centroidj,d,use="complete",method="pearson")
      c.cor[j]=corj$estimate
      pv[j]=corj$p.value
    }
    maxk=which.max(c.cor)
    group=names(maxk)
    vcor=rbind(vcor,c(colnames(vd)[i],group,c.cor[maxk],pv[maxk]))
    if(c.cor[maxk]>0.1){
      vclass[colnames(vd)[i]]=group
    }else{
      vclass[colnames(vd)[i]]=NA
    }
    # vclass[colnames(vd)[i]]=group
  }
  return(list(overlap.gene=gene.sig,cluster=vclass,correlation=vcor,normalized.matrix=vd))
}
