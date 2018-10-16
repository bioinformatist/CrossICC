#' @export
cluster.centroid <- function(x, gene.signature, cluster){
  x.sample.name<-colnames(x)
  # print(x.sample.name)
  cluster<-cluster[x.sample.name]
  # print(cluster)
  # gene.signature <- intersect(rownames(x),gene.signature)
  x.scale <- t(scale(t(x[gene.signature,])))
  centroids <- c()
  for(i in unique(cluster)){
    sub.mat<-x.scale[gene.signature, names(which(cluster==i))]
    if(class(sub.mat)=="numeric"){
      centroids <- cbind(centroids,sub.mat)
    }else{
      centroids <- cbind(centroids,
                         apply(sub.mat, 1, median, na.rm=FALSE)
      )
    }

  }
  colnames(centroids) <- unique(cluster)
  return(centroids)
}
