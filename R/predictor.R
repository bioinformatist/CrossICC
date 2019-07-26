#' Title To calculate the correlation between the predictor centroid and the validation centroid.
#'
#' @param pre.dat a eSet object or eSet-like matrix. with features in rows and samples in columns
#' @param model a CrossICC Object
#' @return a list contains a vecter that store the predict clusters and a normalized expression matrix
#' @export
#'
#' @examples
#' \donttest{
#' CrossICC.object <- CrossICC(demo.platforms, max.iter = 2)
#' predictor(demo.platforms[[1]], CrossICC.object)
#' }
#'
predictor <- function(pre.dat, model) {
    predict.data <- pre.dat
    crossICC.object <- model
    # validation.Data shoud be format features in rows and samples in columns

    crossICC.object.summary <- summaryCrossICC(crossICC.object)

    # using interset gene list for prediction
    predictFeaturelist <- row.names(pre.dat)
    overlap.feature <- intersect(predictFeaturelist, crossICC.object$gene.signature)

    differ.length <- length(crossICC.object$gene.signature) - length(overlap.feature)
    if (differ.length > 0 & differ.length < 5) {
        warning(paste("missing ", differ.length, " features in your expression data set, continute predicting any way ", sep = ""))
    } else if (differ.length >= 5) {

        stop(paste("missing too many (", differ.length, ") features in your expression data set, plz replace your predict data set",
            sep = ""))
    }
    # get centroid
    centroid.list <- lapply(crossICC.object$platforms, cluster.centroid, gene.signature = overlap.feature, cluster = crossICC.object.summary$clusters)

    # get centroid of the centroid

    train.centroid <- centroidOfcentroid(centroid.list, cluster = crossICC.object.summary$clusters)
    # prediction
    vali.predict.bycentroid <- centroid2exp(train.centroid, predict.data)
    # get prediction result
    vali.predict.bycentroid.cluter <- vali.predict.bycentroid$cluster
    vali.predict.normalized.matrix <- vali.predict.bycentroid$normalized.matrix
    return(list(cluster = vali.predict.bycentroid.cluter, matrix = vali.predict.normalized.matrix))
}

#' Title Return centroid of centroid from each platform
#'
#' @param centroid.list a list stored the centroid
#' @param cluster  a named vector with gene names as name and the cluster number as vector value
#' @return a list contains a vecter that store the predict clusters and a normalized expression matrix
centroidOfcentroid <- function(centroid.list, cluster) {
    cluster.name <- unlist(lapply(centroid.list, colnames))
    final.matrix <- do.call(cbind, centroid.list)
    res <- t(apply(final.matrix, 1, tapply, cluster.name, mean))
    return(res)
}

