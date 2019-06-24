centroid <- function(x, na.rm = FALSE) {
    apply(x, 1, median, na.rm = na.rm)
}
centroid.pairwise <- function(x, x.c, y, y.c, method = "pearson") {
    cor = c()
    x.centroid = matrix(data = NA, ncol = length(x.c), nrow = nrow(x))
    for (i in 1:length(x.c)) {
        x.classi = x.c[[i]]
        x.class = x[, x.classi]
        x.centroid[, i] = centroid(x.class)
    }
    y.centroid = matrix(data = NA, ncol = length(y.c), nrow = nrow(y))
    for (i in 1:length(y.c)) {
        y.classi = y.c[[i]]
        y.class = y[, y.classi]
        y.centroid[, i] = centroid(y.class)
    }
    for (j in 1:ncol(x.centroid)) {
        x.centroidj = x.centroid[, j]
        cor_row = c()
        for (i in 1:ncol(y.centroid)) {
            y.centroidi = y.centroid[, i]
            corj = cor.test(x.centroidj, y.centroidi, use = "complete", method = method)
            cor_row[length(cor_row) + 1] = corj$estimate
        }
        cor = rbind(cor, cor_row)
    }
    rownames(cor) = paste("x", 1:length(x.c), sep = "")
    colnames(cor) = paste("y", 1:length(y.c), sep = "")
    return(cor)
}
