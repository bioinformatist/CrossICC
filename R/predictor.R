#' Title To calculate the correlation between the predictor centroid and the validation centroid.
#'
#' @param pre.dat a eSet object or eSet-like matrix. with features in rows and samples in columns
#' @param model a CrossICC Object
#' @return a list contains a vecter that store the predict clusters and a normalized expression matrix
#' @export
#'
#' @examples
#' \donttest{
#' predictor(example.matrices$x, crossICC.object)
#' }
#'
#'
predictor <- function(pre.dat, model) {
  predict.data<-pre.dat
  crossICC.object<-InterationResult()
  #validation.Data shoud be format features in rows and samples in columns

  crossICC.object.summary<-summary.CrossICC(crossICC.object)
  # get centroid
  train.centroid<-cluster.centroid(crossICC.object$platforms[[1]],crossICC.object$gene.signature,crossICC.object.summary$clusters)
  #prediction
  vali.predict.bycentroid<-centroid2exp(train.centroid,predict.data)
  #get prediction result
  vali.predict.bycentroid.cluter<-vali.predict.bycentroid$cluster
  vali.predict.normalized.matrix<-vali.predict.bycentroid$normalized.matrix
  return(list(cluster=vali.predict.bycentroid.cluter,matrix=vali.predict.normalized.matrix))
}

# predictor <- function(x, fs, gene.signature, geneset2gene, filter.cutoff = 1, fdr.cutoff = 1) {
#   filtered.x <- m.f.s(list(x), fdr.cutoff = fdr.cutoff, filter.cutoff = filter.cutoff, perform.mad = FALSE)[[2]][[1]]
#   genewprobe <- gene.signature
#   names(genewprobe) <- gene.signature
#   c1 <- runFAIME(filtered.x, genewprobe, geneset2gene, na.last = "keep", weightRank = FALSE)
#   c2 <- fs
#
#   FAIME.scale <- function(fs) {
#     fs.scale <- t(scale(t(fs)))
#     fs.scale <- replace(fs.scale, fs.scale < -2, -2)
#     fs.scale <- replace(fs.scale, fs.scale > 2, 2)
#     fs.scale
#   }
#
#   c1 <- FAIME.scale(c1)
#   c2 <- FAIME.scale(c2)
#   # for(j in 1:ncol(c1)){
#   #   pcentroidj <- c1[,j]
#   #   cor_row <- c()
#   #   for(i in 1:ncol(c2)){
#   #     vcentroidi <- c2[,i]
#   #     corj <- cor(pcentroidj, vcentroidi, use="complete", method="pearson")
#   #     cor_row[length(cor_row)+1] <- corj$estimate
#   #   }
#   #   cor <- rbind(cor,cor_row)
#   # }
#   # rownames(cor) <- colnames(c1)
#   # colnames(cor) <- colnames(c2)
#   # mapply(cor, c1, c2)
#   cor(c1, c2)
# }
