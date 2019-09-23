### subset mergeExpressionSet
check <- function(x) {
    if (!is.element(class(x), c("list", "mergeExpressionSet", "ExpressionSet", "matrix")))
        stop("all data must be either a list, a mergeExpressionSet, a matrix, or an ExpressionSet")
}

mergeget <- function(x) {
    if (is.element(class(x), "mergeExpressionSet")) {
        return(exprs(x))
    }
    if (is.element(class(x), "ExpressionSet")) {
        return(x)
    }
    if (is.element(class(x), "list")) {
        if (length(x) != 4)
            stop("if you want to merge a list, a list should have at least four slots, 'expression matirx', 'phenodata', 'gene names' and 'notes'.")
        if (!is.matrix(x[[1]]))
            stop("first object of the input list must be an expression matrix.")
        if (!is.vector(x[[3]]))
            stop("third object of the input list must be a gene name vector.")
        if (is.null(x[[2]]))
            stop("second object of the input list can not be NULL.")

        rownames(x[[1]]) <- x[[3]]
        pd <- data.frame(x[[2]])
        if (is.null(colnames(x[[1]])))
            sn = as.character(c(seq_len(ncol(x[[1]])))) else sn <- colnames(x[[1]])
        rownames(pd) <- sn
        ad <- new("AnnotatedDataFrame", data = pd)
        es <- new("ExpressionSet", exprs = x[[1]], phenoData = ad)
        return(es)
    }
    if (is.element(class(x), "matrix")) {
        if (is.null(rownames(x)))
            stop("if you want to merge matrix, rownames of matrix can not be NULL.")

        pd <- data.frame(rep(NA, (ncol(x))))
        if (is.null(colnames(x)))
            sn = as.character(seq_len(ncol(x))) else sn <- colnames(x)
        row.names(pd) <- sn

        rownames(pd) <- sn
        ad <- new("AnnotatedDataFrame", data = pd)
        es <- new("ExpressionSet", exprs = x, phenoData = ad)
        return(es)
    }
    stop("If you want to merge, the input object should be 'ExpressionSet', 'list', 'matrix', or 'mergeExpressionSet'.")
}



AverageDuplicates <- function(data.exprs, data.acc) {
    data.acc <- as.character(data.acc)
    dups <- rev(duplicated(rev(data.acc))) + duplicated(data.acc)
    dups <- ifelse(dups == 0, 0, 1)

    if (sum(dups) > 0) {
        data1.exprs <- as.matrix(data.exprs[dups == 0, ])
        data1.acc <- data.acc[dups == 0]
        data2.exprs <- as.matrix(data.exprs[dups == 1, ])
        data2.acc <- data.acc[dups == 1]
        data3.acc <- unique(data2.acc)
        data3.exprs <- matrix(NA, length(data3.acc), ncol(data.exprs))

        k <- 0
        for (i in data3.acc) {
            k <- k + 1
            data3.exprs[k, ] <- apply(as.matrix(data2.exprs[data2.acc == i, ]), 2, mean, na.rm = TRUE)
        }

        data4.exprs <- rbind(data1.exprs, data3.exprs)
        data4.acc <- c(data1.acc, data3.acc)
        keep <- ifelse(as.character(data4.acc) == "", 0, 1)
        data.exprs <- as.matrix(data4.exprs[keep == 1, ])
        data.acc <- as.matrix(data4.acc[keep == 1])
    }
    return(list(data = as.data.frame(data.exprs), acc = data.acc))
}


# All statements that work with list studynames have been removed for its substitute calling
mergeExprs <- function(...) {
    arg <- list(...)

    x <- alist(... = )
    k <- 0

    # When there's only one list passed to ..., it should be flattened as matrices, such as: arg[[1]] is flattened list with certain
    # sub-lists
    if ((is(arg, "list")) & (length(arg) == 1)) {
        for (j in seq_len(length(arg[[1]]))) {
            k <- k + 1
            x[[k]] <- mergeget(arg[[1]][[j]])
        }
    } else {
        for (i in seq_len(length(arg))) {
            check(arg[[i]])
            if (is.element(class(arg[[i]]), "mergeExpressionSet")) {
                mm <- mergeget(arg[[i]])
                for (j in seq_len(length(arg[[i]]))) {
                  k <- k + 1
                  x[[k]] <- mm[[j]]
                }
            } else {
                k <- k + 1
                x[[k]] <- mergeget(arg[[i]])
            }
        }
    }



    tt <- length(x)

    nnote <- matrix(NA, tt, 2)

    for (i in seq_len(tt)) {
        if (i == 1)
            iid <- as.matrix(featureNames(x[[i]])) else iid <- rbind(iid, as.matrix(featureNames(x[[i]])))
        nnote[i, 2] <- ""
    }

    iid <- as.vector(sort(unique(iid)))

    # generate the matrices with missing value 'NA'

    for (i in seq_len(tt)) {
        y <- assayData(x[[i]])[["exprs"]]
        idy <- featureNames(x[[i]])
        y.avg <- AverageDuplicates(y, idy)
        assayData(x[[i]]) <- list(exprs = as.matrix(y.avg$data))

        featureNames(x[[i]]) <- y.avg$acc
    }

    # generate the vector with common id '1', o.w. '0'

    idmatrix <- matrix(0, length(iid), tt)
    index <- as.vector(nnote[, 2])
    for (i in seq_len(tt)) {
        idx <- featureNames(x[[i]])
        cc <- match(iid, idx)
        idmatrix[, i] <- ifelse(is.na(cc), 0, 1)
    }

    rownames(idmatrix) <- iid

    # generate the list that we want
    merged <- new("mergeExpressionSet", data = x, geneStudy = idmatrix, notes = "")
    return(merged)
}



check.length <- function(x, wh) {
    if (length(pData(x)[, wh]) != ncol(assayData(x)[["exprs"]]))
        stop("Phenodata error: phenodata should have the same length as the number of columns of expression data.")
    return()
}


.mergemodel <- function(x, outcome = NULL, event = NULL, method) {
    if (method == "linear") {
        reg <- lm(x ~ outcome, na.action = na.omit)
        sx <- summary(reg)$sigma
        xx <- x/sx
        stdbeta <- lm(outcome ~ xx, na.action = na.omit)$coeff[2]

        tmp <- lm(outcome ~ x)
        beta <- tmp$coeff[2]
        tvalue <- summary(tmp)$coeff[2, 3]

        return(list(beta = beta, stdbeta = stdbeta, tvalue = tvalue))
    }
    if (method == "cox") {
        ddtt <- list(time = outcome, status = event, expr = x)
        expr <- x
        status <- event
        time <- outcome

        cox.out <- coxph(Surv(time, status) ~ expr, data = ddtt, na.action = na.omit)

        beta <- cox.out$coeff

        tmp <- status[!is.na(status) & !is.na(outcome) & !is.na(x)]
        nobs = sum(tmp == 1)
        zscore <- cox.out$coeff/sqrt(cox.out$var)

        std <- list(time = outcome, status = event, expr = (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE))
        cox.stdout <- coxph(Surv(time, status) ~ expr, data = std, na.action = na.omit)
        stdbeta <- cox.stdout$coeff

        return(list(beta = beta, stdbeta = stdbeta, zscore = zscore))
    }
    if (method == "logistic") {
        reg <- glm(event ~ x, family = binomial, na.action = na.omit)
        beta <- summary(reg)$coefficients[2, 1]
        zscore <- summary(reg)$coefficients[2, 3]

        tmp <- event[!is.na(x) & !is.na(event)]

        a <- x[tmp == 0]
        b <- x[tmp == 1]
        na <- sum(tmp == 1)
        nb <- length(tmp) - na

        if (na == 0 || nb == 0)
            sp <- sd(x) else sp <- (na * (sd(a, na.rm = TRUE)^2) + nb * (sd(b, na.rm = TRUE)^2))/length(x)

        stdx <- x/sqrt(sp)
        regstd <- glm(event ~ stdx, family = binomial, na.action = na.omit)

        stdbeta <- summary(regstd)$coefficients[2, 1]

        return(list(beta = beta, stdbeta = stdbeta, zscore = zscore))
    }

}


.intcor.order = function(x) return((length(x) + 1) - order(x))

maxintcor <- function(A, B) {
    n1 = dim(A)[2]
    n2 = dim(B)[2]
    AB = cbind(A, B)
    ns = n1
    n = n1 + n2
    nn = ns + 1
    covs = cov(AB, use = "pairwise")
    AA = covs[seq_len(ns), seq_len(ns)]
    BB = covs[nn:n, nn:n]
    AB = covs[seq_len(ns), nn:n]
    AAinv = ginv(AA)
    BBinv = ginv(BB)
    eigens = eigen(AAinv %*% AB %*% BBinv %*% t(AB))
    return(sqrt(max(Re(eigens$values))))
}


.intcor <- function(x, method = "pearson", exact = TRUE) {
    nn <- length(x)
    if (nn <= 1)
        stop("Number of studies in the mergeExpressionSet should not less than 2.")
    if (method != "pearson" && !exact)
        stop("When exact is FALSE, you can only use the method ''pearson''.")
    if (method != "pearson" && method != "spearman" && exact)
        stop("When exact is TRUE, you can only use the methods, ''pearson'' and ''spearman''.")

    geneid <- rep(1, nrow(x@geneStudy))
    for (i in seq_len(nn)) {
        geneid <- geneid * (x@geneStudy[, i])
    }
    geneuid <- geneNames(x)[geneid == 1]
    nuid <- length(geneuid)

    pcor <- alist(... = )
    nnote <- rep(NA, nn)

    for (i in seq_len(nn)) {
        matches1 <- match(geneuid, featureNames(exprs(x)[[i]]))
        pcor[[i]] <- assayData(exprs(x)[[i]])[["exprs"]][matches1, ]
        if (method == "spearman" && exact)
            pcor[[i]] <- t(apply(pcor[[i]], 1, .intcor.order))

        nnote[i] <- paste("study", i, sep = " ")
    }

    names(pcor) <- names(x)
    np <- nn * (nn - 1)/2
    ppair <- matrix(0, np, 2)
    k <- 1
    for (i in seq_len(nn - 1)) {
        for (j in (i + 1):nn) {
            ppair[k, ] <- c(i, j)
            k <- k + 1
        }
    }

    norm = function(x) return((x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE))

    icor <- matrix(NA, nuid, np)
    canc <- rep(0, np)
    for (i in seq_len(np)) {
        CC1 <- pcor[[ppair[i, 1]]]
        CC2 <- pcor[[ppair[i, 2]]]
        if (exact)
            icor[, i] <- .integ.cal(CC1, CC2, method = method) else icor[, i] <- .icor(CC1, CC2, method = method)

        A = t(apply(CC1, 1, norm))
        B = t(apply(CC2, 1, norm))

        canc[i] <- maxintcor(A, B)
    }

    rownames(icor) <- geneuid
    labels <- names(x)
    if (is.null(labels))
        labels = ""
    result <- new("mergeCor", pairwise.cors = icor, notes = labels, maxcors = canc)
    return(result)
}

isna <- function(x) return(is.na(x))

.icor = function(A, B, method = NULL) {
    rn <- rownames(A)
    norm = function(x) return((x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE))
    A = t(apply(A, 1, norm))
    B = t(apply(B, 1, norm))

    m = dim(A)[1]
    n1 = dim(A)[2]
    n2 = dim(B)[2]
    m1 = mat1 = cov(A, use = "pairwise.complete.obs")
    m2 = mat2 = cov(B, use = "pairwise.complete.obs")
    mat12 = matrix(nrow = n1, ncol = n2)
    for (i in seq_len(n1)) {
        for (j in seq_len(n2)) {
            mat12[i, j] = cov(A[, i], B[, j], use = "pairwise.complete.obs")
        }
    }

    iicor = function(index) {
        mm = mat12
        aa = A[index, ]
        bb = B[index, ]
        whca = (seq_len(length(aa)))[is.na(aa)]
        whcb = (seq_len(length(bb)))[is.na(bb)]

        if (length(whca) != 0) {
            aa <- aa[-whca]
            mm <- mat12[-whca, ]
            m1 <- mat1[-whca, -whca]
        }

        if (length(whcb) != 0) {
            bb <- bb[-whcb]
            mm <- mm[, -whcb]
            m2 <- mat2[-whcb, -whcb]
        }
        return((aa %*% mm %*% bb)/sqrt(aa %*% m1 %*% aa)/sqrt(bb %*% m2 %*% bb))
    }
    ans <- vapply(seq_len(m), iicor, numeric(1))
    names(ans) <- rn
    return(ans)
}

.nullicornorm <- function(A, B, n = 10000) {
    rn <- rownames(A)
    norm = function(x) return((x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE))
    A = t(apply(A, 1, norm))
    B = t(apply(B, 1, norm))

    m = dim(A)[1]
    n1 = dim(A)[2]
    n2 = dim(B)[2]
    mat1 = cov(A, use = "pairwise.complete.obs")
    mat2 = cov(B, use = "pairwise.complete.obs")
    mat12 = matrix(nrow = n1, ncol = n2)
    for (i in seq_len(n1)) {
        for (j in seq_len(n2)) {
            mat12[i, j] = cov(A[, i], B[, j], use = "pairwise.complete.obs")
        }
    }
    getone = function(dum) {
        n1 = dim(mat1)[1]
        n2 = dim(mat2)[1]
        x = rnorm(n1)
        y = rnorm(n2)
        v1 = x %*% mat1 %*% x
        v2 = y %*% mat2 %*% y
        c = x %*% mat12 %*% y
        return(c/(sqrt(v1) * sqrt(v2)))
    }
    return(vapply(seq_len(n), getone, numeric(1)))
}

.integ.cal <- function(mat1, mat2, method = NULL) {
    rn <- rownames(mat1)
    integ <- rep(0, nrow(mat1))

    m1 <- t(apply(mat1, 1, .integ.norm))
    m2 <- t(apply(mat2, 1, .integ.norm))

    sziicor = function(index) {
        x <- m1[index, ]
        y <- m2[index, ]

        A <- m1[-index, ]
        B <- m2[-index, ]
        nammx = function(indx) {
            aa = A[indx, ]
            sum(aa * x, na.rm = TRUE)
        }
        rx = vapply(seq_len(nrow(m1) - 1), nammx, numeric(1))

        nammy = function(indy) {
            bb = B[indy, ]
            sum(bb * y, na.rm = TRUE)
        }
        ry = vapply(seq_len(nrow(m2) - 1), nammy, numeric(1))

        integ[index] <- cor(rx, ry, use = "pairwise.complete.obs")
    }
    ans = vapply(seq_len(nrow(mat1)), sziicor, numeric(1))
    names(ans) = rn
    return(ans)
}


.integ.norm <- function(x) {
    xm <- mean(x, na.rm = TRUE)
    ss <- sum((x - xm) * (x - xm), na.rm = TRUE)
    return((x - xm)/sqrt(ss))
}

.ic.plot <- function(A, B, mmatches, xlab = NA, ylab = NA, title = NULL, main = NULL, scale = NULL, square = FALSE, ...) {
    x <- A[mmatches, ]
    y <- B[mmatches, ]

    norm = function(x) return((x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE))
    A = t(apply(A, 1, norm))
    B = t(apply(B, 1, norm))

    m = dim(A)[1]
    n1 = dim(A)[2]
    n2 = dim(B)[2]
    mat12 = matrix(nrow = n1, ncol = n2)
    for (i in seq_len(n1)) {
        for (j in seq_len(n2)) {
            mat12[i, j] = cov(A[, i], B[, j], use = "pairwise.complete.obs")
        }
    }

    cormat <- mat12

    nx = length(x)
    ny = length(y)
    dc = dim(cormat)
    if ((dc[1] != nx && dc[1] != ny) || (dc[2] != nx && dc[2] != ny)) {
        stop("cormat is not conformable to x or y")
    }
    if (dc[1] != nx) {
        cormat = t(cormat)
        dc = dim(cormat)
    }
    rg = function(x) return(max(x) - min(x))

    x = norm(x)
    y = norm(y)
    if (is.null(scale))
        scale = 2 * max(c(rg(x), rg(y)))

    cx = as.vector(abs(cormat))/max(abs(cormat))
    if (square == TRUE)
        cx = cx^2
    sn = sign(as.vector(cormat))
    cl = sn
    cl[cl == -1] = 2

    xx = rep(x, ny)
    yy = rep(y, rep(nx, ny))
    sn2 = sign(abs(xx[sn == -1]) - abs(yy[sn == -1]))
    xx[sn == -1] = sn2 * xx[sn == -1]
    yy[sn == -1] = -sn2 * yy[sn == -1]

    plot(xx, yy, pch = 16, cex = scale * cx^2, axes = FALSE, xlab = xlab, ylab = ylab, cex.lab = 1.5, main = main, ...)

    return()
}

.dens.mergeExpressionSet <- function(x, main = NA, x.legend = NULL, y.legend = NULL, cex.legend = NULL, title = NULL, method = NULL,
    ...) {
    nn <- length(x)
    if (nn <= 1)
        stop("Number of studies in the mergeExpressionSet should not less than 2.")

    nnote <- rep(NA, nn)

    for (i in seq_len(nn)) {
        nnote[i] <- nnote[i] <- paste("study", i, sep = " ")
    }

    geneid <- rep(1, nrow(x@geneStudy))
    for (i in seq_len(nn)) {
        geneid <- geneid * (x@geneStudy[, i])
    }
    geneuid <- geneNames(x)[geneid == 1]
    nuid <- length(geneuid)
    rtn <- alist(... = )

    if (nn > 2) {
        np <- nn * (nn - 1)/2
        ppair <- matrix(0, np, 2)
        k <- 1
        for (i in seq_len(nn - 1)) {
            for (j in (i + 1):nn) {
                ppair[k, ] <- c(i, j)
                k <- k + 1
            }
        }

        ppair1 <- matrix(0, np, 2)
        k <- 1
        for (i in seq_len(nn - 1)) {
            for (j in (i + 1):nn) {
                ppair1[k, ] <- c(nn + 1 - i, nn + 1 - j)
                k <- k + 1
            }
        }


        exp1 <- matrix(NA, 10000, np)
        obs <- matrix(NA, nuid, np)

        for (i in seq_len(np)) {
            matches1 <- match(geneuid, featureNames(exprs(x)[[ppair[i, 1]]]))
            matches2 <- match(geneuid, featureNames(exprs(x)[[ppair[i, 2]]]))
            exprs1 <- assayData(exprs(x)[[ppair[i, 1]]])[["exprs"]][matches1, ]
            exprs2 <- assayData(exprs(x)[[ppair[i, 2]]])[["exprs"]][matches2, ]

            obs[, i] <- .icor(exprs1, exprs2)
            exp1[, i] <- .nullicornorm(exprs1, exprs2)
        }

        k <- m <- 1
        par(mfrow = c(nn, nn), oma = c(0, 0, 8, 0))

        for (i in seq_len(nn)) {
            for (j in seq_len(nn)) {
                if (k <= np) {
                  if (ppair[k, 1] == i & ppair[k, 2] == j) {
                    d1 <- exp1[, k]
                    d <- obs[, k]
                    xmax <- max(c(density(d1)$x, density(d)$x))
                    xmin <- min(c(density(d1)$x, density(d)$x))
                    ymax <- max(c(density(d1)$y, density(d)$y))
                    range <- xmax - xmin
                    plot(density(d), xlim = c(xmin - range/10, xmax + range/10), ylim = c(0, ymax + ymax/5), main = main, ...)
                    lines(density(d1), col = 2, ...)

                    rtn[[k]] <- d1
                    k = k + 1
                  }
                }
                if (m <= np) {
                  if (ppair1[np + 1 - m, 1] == i & ppair1[np + 1 - m, 2] == j) {
                    if (ppair1[np + 1 - m, 1] == nn & ppair1[np + 1 - m, 2] == 1) {
                      aa <- c(0.5, 0.5)
                      bb <- c(0.7, 0.8)
                      plot(aa, bb, xlim = c(0, 1), ylim = c(0, 1), xlab = NA, ylab = NA, pch = " ", xaxt = "n", yaxt = "n")
                      legend.text <- c("observed", "null distribution")
                      op <- par(bg = "white")
                      if (is.null(x.legend))
                        x.legend <- 0.1
                      if (is.null(y.legend))
                        y.legend <- 0.9
                      if (is.null(cex.legend))
                        cex.legend <- 1
                      legend(x.legend, y.legend, paste(legend.text), col = seq_len(2), lty = 1, cex = cex.legend)
                    } else {
                      aa <- c(0.5, 0.5)
                      bb <- c(0.7, 0.8)
                      plot(aa, bb, xlim = c(0, 1), ylim = c(0, 1), xlab = NA, ylab = NA, pch = " ", xaxt = "n", yaxt = "n")
                    }
                    m = m + 1
                  }
                }
                if (i == j) {
                  aa <- c(0.5, 0.5)
                  bb <- c(0.7, 0.8)
                  plot(aa, bb, xlim = c(0, 1), ylim = c(0, 1), xlab = NA, ylab = NA, pch = " ", xaxt = "n", yaxt = "n")
                  text(0.5, 0.5, cex = 1.2, paste(names(x)[i]), col = 1)
                }
            }
        }
        if (is.null(title))
            title <- "Integrated Correlation:\n\ntrue and null distributions"
        mtext(title, line = 0.5, cex = 1.2, outer = TRUE)
        par(mfrow = c(1, 1))
    }

    if (nn == 2) {
        geneid <- x@geneStudy[, 1] * x@geneStudy[, 2]

        geneuid <- geneNames(x)[geneid == 1]
        nuid <- length(geneuid)

        matches1 <- match(geneuid, featureNames(exprs(x)[[1]]))
        matches2 <- match(geneuid, featureNames(exprs(x)[[2]]))
        exprs1 <- assayData(exprs(x)[[1]])[["exprs"]][matches1, ]
        exprs2 <- assayData(exprs(x)[[2]])[["exprs"]][matches2, ]
        exprs1[is.na(exprs1)] <- 0
        exprs2[is.na(exprs2)] <- 0

        d <- .icor(exprs1, exprs2)
        d1 <- .nullicornorm(exprs1, exprs2)

        xmax <- max(c(density(d1)$x, density(d)$x))
        xmin <- min(c(density(d1)$x, density(d)$x))
        ymax <- max(c(density(d1)$y, density(d)$y))
        range <- xmax - xmin
        if (is.null(main))
            main <- "Integrated Correlation:\ntrue and null distributions"
        plot(density(d), xlim = c(xmin - range/10, xmax + range/10), ylim = c(0, ymax + ymax/5), main = main, ...)
        lines(density(d1), col = 2, ...)
        legend.text <- c("observed", "null distribution")
        op <- par(bg = "white")
        if (is.null(x.legend))
            x.legend <- xmax - range/3
        if (is.null(y.legend))
            y.legend <- ymax
        if (is.null(cex.legend))
            cex.legend <- 0.9
        legend(x.legend, y.legend, paste(legend.text), col = seq_len(2), lty = 1, cex = cex.legend)
        rtn <- d1
    }
    return(rtn)

}

.show.mergeExpressionSet <- function(object) {
    x = exprs(object)
    y = list()
    for (i in seq_len(length(x))) {
        show(x[[i]])
    }
    return()
}



.summary.mergeExpressionSet <- function(object) {
    report <- alist(... = )
    nn <- length(object)

    nnote <- rep(NA, nn)

    for (i in seq_len(nn)) {
        nnote[i] <- paste("study", i, sep = " ")
    }
    report[[1]] <- matrix(NA, 2, nn)
    report[[1]][1, ] <- names(object)
    for (i in seq_len(nn)) {
        report[[1]][2, i] <- nrow(assayData(exprs(object)[[i]])[["exprs"]])
    }

    report[[2]] <- matrix(NA, 2, nn)
    report[[2]][1, ] <- names(object)
    for (i in seq_len(nn)) {
        report[[2]][2, i] <- ncol(assayData(exprs(object)[[i]])[["exprs"]])
    }

    report[[3]] <- nnote

    names(report) <- c("Number of Genes in Each Study", "Number of Samples in Each Study", "Notes")
    return(report)
}




