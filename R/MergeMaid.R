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
            sn = as.character(c(1:ncol(x[[1]]))) else sn <- colnames(x[[1]])
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
            sn = as.character(c(1:ncol(x))) else sn <- colnames(x)
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
    if ((is.element(class(arg), "list")) & (length(arg) == 1)) {
        for (j in 1:length(arg[[1]])) {
            k <- k + 1
            x[[k]] <- mergeget(arg[[1]][[j]])
        }
    } else {
        for (i in 1:length(arg)) {
            check(arg[[i]])
            if (is.element(class(arg[[i]]), "mergeExpressionSet")) {
                mm <- mergeget(arg[[i]])
                for (j in 1:length(arg[[i]])) {
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

    for (i in 1:tt) {
        if (i == 1)
            iid <- as.matrix(featureNames(x[[i]])) else iid <- rbind(iid, as.matrix(featureNames(x[[i]])))
        nnote[i, 2] <- ""
    }

    iid <- as.vector(sort(unique(iid)))

    # generate the matrices with missing value 'NA'

    for (i in 1:tt) {
        y <- assayData(x[[i]])[["exprs"]]
        idy <- featureNames(x[[i]])
        y.avg <- AverageDuplicates(y, idy)
        assayData(x[[i]]) <- list(exprs = as.matrix(y.avg$data))

        featureNames(x[[i]]) <- y.avg$acc
    }

    # generate the vector with common id '1', o.w. '0'

    idmatrix <- matrix(0, length(iid), tt)
    index <- as.vector(nnote[, 2])
    for (i in 1:tt) {
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

.model.outcome <- function(x, method = NULL, outcome = NULL, outcome2 = NULL) {
    if (is.null(method))
        stop("Specify the method you want to use.")
    nn <- length(x)
    if (nn <= 1)
        stop("Number of studies in the mergeExpressionSet should not less than 2.")

    if (method == "linear" | method == "logistic") {
        if (!is.null(outcome)) {
            if (!is.element(class(outcome), "list") & !is.vector(outcome))
                stop("Phenodata error: the data type of phenodata is not correct, please look at the help files.")
            if (is.element(class(outcome), "list")) {
                if (length(outcome) != nn)
                  stop("Phenodata error: the length of input phenodata list should be equal to the number of studies.")
                for (i in 1:nn) {
                  if (length(outcome[[i]]) != ncol(assayData(exprs(x)[[i]])[["exprs"]]))
                    stop("Phenodata error: phenodata should have the same length as the number of columns of expression data.")
                }
                out <- outcome
            }
            if (is.vector(outcome) & !is.list(outcome)) {
                if (length(outcome) != nn)
                  stop("Phenodata error: You should specify the phenodata.") else {
                  out <- alist(... = )
                  for (i in 1:nn) {
                    check.length(exprs(x)[[i]], outcome[i])
                    out[[i]] <- pData(exprs(x)[[i]])[, outcome[i]]
                  }
                }
            }
        } else {
            stop("Phenodata error: You should specify the phenodata.")
        }
    }

    if (method == "cox") {
        if (!is.null(outcome) & !is.null(outcome2)) {
            if ((!is.element(class(outcome), "list") & !is.vector(outcome)) | (!is.element(class(outcome2), "list") & !is.vector(outcome2)))
                stop("Phenodata error: the data type of phenodata is not correct, please look at the help files.")
            if (is.element(class(outcome), "list")) {
                if (length(outcome) != nn)
                  stop("Phenodata error: the length of input phenodata list should be equal to the number of studies.")
                for (i in 1:nn) {
                  if (length(outcome[[i]]) != ncol(assayData(exprs(x)[[i]])[["exprs"]]))
                    stop("Phenodata error: phenodata should have the same length as the number of columns of expression data.")
                }
                out <- outcome
            }

            if (is.vector(outcome) & !is.list(outcome)) {
                if (length(outcome) != nn)
                  stop("Phenodata error: You should specify the phenodata.") else {
                  out <- alist(... = )
                  for (i in 1:nn) {
                    check.length(exprs(x)[[i]], outcome[i])
                    out[[i]] <- pData(exprs(x)[[i]])[, outcome[i]]
                  }
                }
            }

            if (is.element(class(outcome2), "list")) {
                if (length(outcome2) != nn)
                  stop("Phenodata error: the length of input phenodata list should be equal to the number of studies.")
                for (i in 1:nn) {
                  if (length(outcome2[[i]]) != ncol(assayData(exprs(x)[[i]])[["exprs"]]))
                    stop("Phenodata error: phenodata should have the same length as the number of columns of expression data.")
                }
                out2 <- outcome2
            }
            if (is.vector(outcome2) & !is.list(outcome)) {
                if (length(outcome2) != nn)
                  stop("Phenodata error: You should specify the phenodata.") else {
                  out2 <- alist(... = )
                  for (i in 1:nn) {
                    check.length(exprs(x)[[i]], outcome2[i])
                    out2[[i]] <- pData(exprs(x)[[i]])[, outcome2[i]]
                  }
                }
            }
        } else {
            stop("Phenodata error: You should specify the phenodata.")
        }
    }

    if (method != "linear" & method != "logistic" & method != "cox")
        stop("ARG should be one of linear, cox, logistic")

    geneid <- rep(1, nrow(x@geneStudy))
    for (i in 1:nn) {
        geneid <- geneid * (x@geneStudy[, i])
    }
    geneuid <- geneNames(x)[geneid == 1]
    nuid <- length(geneuid)

    beta <- matrix(0, nuid, nn)
    stdbeta <- matrix(0, nuid, nn)
    zscore <- matrix(0, nuid, nn)

    for (i in 1:nn) {
        matches1 <- match(geneuid, featureNames(exprs(x)[[i]]))
        exprs1 <- assayData(exprs(x)[[i]])[["exprs"]][matches1, ]
        if (method == "linear") {
            outcome1 <- out[[i]]
            result1 <- apply(exprs1, 1, .mergemodel, outcome = outcome1, method = "linear")
            for (gene in 1:nuid) {
                beta[gene, i] <- result1[[gene]][[1]]
                stdbeta[gene, i] <- result1[[gene]][[2]]
                zscore[gene, i] <- result1[[gene]][[3]]
            }
        }
        if (method == "logistic") {
            event1 <- out[[i]]
            result1 <- apply(exprs1, 1, .mergemodel, event = event1, method = "logistic")
            for (gene in 1:nuid) {
                beta[gene, i] <- result1[[gene]][[1]]
                stdbeta[gene, i] <- result1[[gene]][[2]]
                zscore[gene, i] <- result1[[gene]][[3]]
            }
        }
        if (method == "cox") {
            outcome1 <- out[[i]]
            event1 <- out2[[i]]
            result1 <- apply(exprs1, 1, .mergemodel, outcome = outcome1, event = event1, method = "cox")
            for (gene in 1:nuid) {
                beta[gene, i] <- result1[[gene]][[1]]
                stdbeta[gene, i] <- result1[[gene]][[2]]
                zscore[gene, i] <- result1[[gene]][[3]]
            }
        }
    }
    rownames(beta) <- rownames(stdbeta) <- rownames(zscore) <- geneuid

    nnote <- rep(NA, nn)

    for (i in 1:nn) {
        nnote[i] <- paste("study", i, sep = " ")
    }


    colnames(beta) <- colnames(stdbeta) <- colnames(zscore) <- names(x)

    result <- new("mergeCoeff", coeff = beta, coeff.std = stdbeta, zscore = zscore, method = method)
    return(result)
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
    AA = covs[1:ns, 1:ns]
    BB = covs[nn:n, nn:n]
    AB = covs[1:ns, nn:n]
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
    for (i in 1:nn) {
        geneid <- geneid * (x@geneStudy[, i])
    }
    geneuid <- geneNames(x)[geneid == 1]
    nuid <- length(geneuid)

    pcor <- alist(... = )
    nnote <- rep(NA, nn)

    for (i in 1:nn) {
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
    for (i in 1:(nn - 1)) {
        for (j in (i + 1):nn) {
            ppair[k, ] <- c(i, j)
            k <- k + 1
        }
    }

    norm = function(x) return((x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE))

    icor <- matrix(NA, nuid, np)
    canc <- rep(0, np)
    for (i in 1:np) {
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
    for (i in 1:n1) {
        for (j in 1:n2) {
            mat12[i, j] = cov(A[, i], B[, j], use = "pairwise.complete.obs")
        }
    }

    iicor = function(index) {
        mm = mat12
        aa = A[index, ]
        bb = B[index, ]
        whca = (1:length(aa))[is.na(aa)]
        whcb = (1:length(bb))[is.na(bb)]

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
    ans = sapply(1:m, iicor)
    names(ans) = rn
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
    for (i in 1:n1) {
        for (j in 1:n2) {
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
    return(sapply(1:n, getone))
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
        rx = sapply(1:(nrow(m1) - 1), nammx)

        nammy = function(indy) {
            bb = B[indy, ]
            sum(bb * y, na.rm = TRUE)
        }
        ry = sapply(1:(nrow(m2) - 1), nammy)

        integ[index] <- cor(rx, ry, use = "pairwise.complete.obs")
    }
    ans = sapply(1:(nrow(mat1)), sziicor)
    names(ans) = rn
    return(ans)
}


.integ.norm <- function(x) {
    xm <- mean(x, na.rm = TRUE)
    ss <- sum((x - xm) * (x - xm), na.rm = TRUE)
    return((x - xm)/sqrt(ss))
}


.plot.mergeCoeff <- function(x, y, main = NULL, oma = NULL, ...) {
    if (is.null(main)) {
        if (x[[2]] == "coeff")
            cc <- NULL
        if (x[[2]] == "coeff.std")
            cc <- "standardized coefficients"
        if (x[[2]] == "zscore")
            cc <- "zscore"
        if (x[[3]] == "linear")
            method <- "Linear Regression"
        if (x[[3]] == "cox")
            method <- "Cox Hazard Rate"
        if (x[[3]] == "logistic")
            method <- "Logistic Regression"
        main = paste(method, "Coefficients\n\n", cc, sep = " ")

    }

    if (ncol(x[[1]]) == 2)
        plot(x[[1]][, 1], x[[1]][, 2], main = main, ...)
    if (ncol(x[[1]]) > 2)
        pairs(x[[1]], main = main, ...)

    return()
}

.fdr.mergeCor <- function(x, fdr = NULL) {
    nn <- length(x@notes)
    if (is.null(labels))
        labels <- x@notes
    if (ncol(x@pairwise.cors) != nn * (nn - 1)/2)
        stop("You need to specify the names for the studies.")
    if (nn <= 1)
        stop("You need to specify more than one study.")

    np <- nn * (nn - 1)/2
    ppair <- matrix(0, np, 2)
    k <- 1
    for (i in 1:(nn - 1)) {
        for (j in (i + 1):nn) {
            ppair[k, ] <- c(i, j)
            k <- k + 1
        }
    }
    nuid <- nrow(x@pairwise.cors)
    exp1 <- matrix(NA, 10000, np)
    obs <- matrix(NA, nuid, np)

    for (i in 1:np) {
        mat1 = x@correlation.matrix[[ppair[i, 1]]]
        mat2 = x@correlation.matrix[[ppair[i, 2]]]
        mat12 = x@pairwise.covs[[i]]

        if (np == 1)
            obs[, i] <- x@pairwise.cors else obs[, i] <- x@pairwise.cors[, i]
        exp1[, i] <- .getnull(mat1, mat2, mat12)
    }

    od <- order(fdr)
    temp <- fdr[od]
    rr <- rep(0, length(fdr))
    results <- rep(0, length(fdr))

    d <- apply(obs, 1, mean)
    d1 <- apply(exp1, 1, mean)

    if (range(d)[1] > range(d1)[2] || range(d)[2] < range(d1)[1])
        stop("We can't calculate in your case, please recheck the distribution of integrative correlation using the function intcorDens().")


    fexp1 <- approxfun(density(d1)$x, density(d1)$y, method = "linear")
    fexp0 <- approxfun(density(d)$x, density(d)$y, method = "linear")

    for (tt in 1:512) {
        if (density(d)$x[512 - tt] < max(d))
            break
    }

    for (j in 1:length(fdr)) {
        aa <- 0
        cutoff1 <- 512
        for (i in tt:512) {
            if (density(d)$x[512 - i] < density(d1)$x[512]) {
                aa <- (integrate(fexp1, density(d)$x[512 - i], density(d1)$x[512])[[1]])/(integrate(fexp0, density(d)$x[512 - i], density(d)$x[512])[[1]])
            }
            if (aa > temp[j]) {
                if (i != tt)
                  cutoff1 <- 512 - i
                break
            }
        }
        if (cutoff1 == 512)
            rr[j] <- NA else rr[j] <- density(d)$x[cutoff1]
    }

    results[order(fdr)] <- rr
    return(list(fdr = fdr, cutoff = results))

}



.plot.mergeExpressionSet <- function(x, y, labels = NULL, geneid = NULL, xlab = NA, ylab = NA, title = NULL, main = NULL, scale = NULL,
    square = FALSE, plotype = 1, ...) {
    nn <- length(x)
    if (is.null(labels))
        labels <- names(x)
    if (length(labels) != nn)
        stop("You need to specify the names for the studies.")
    if (nn <= 1)
        stop("You need to specify more than one study.")

    nn <- length(x)
    if (nn <= 1)
        stop("Number of studies in the mergeExpressionSet should not less than 2.")

    geneid1 <- rep(1, nrow(x@geneStudy))
    for (i in 1:nn) {
        geneid1 <- geneid1 * (x@geneStudy[, i])
    }
    gn <- geneNames(x)[geneid1 == 1]
    nuid <- length(gn)

    if (!is.null(geneid)) {
        mmatches <- match(geneid, gn)
        if (is.na(mmatches))
            stop("No such common gene id.")
    }

    if (nn > 2) {
        pcor <- alist(... = )
        nnote <- rep(NA, nn)


        np <- nn * (nn - 1)/2
        ppair <- matrix(0, np, 2)
        k <- 1
        for (i in 1:(nn - 1)) {
            for (j in (i + 1):nn) {
                ppair[k, ] <- c(i, j)
                k <- k + 1
            }
        }

        indextmp <- ppair[, 1] * nn + ppair[, 2]

        ppair1 <- matrix(0, np, 2)
        k <- 1
        for (i in 1:(nn - 1)) {
            for (j in (i + 1):nn) {
                ppair1[k, ] <- c(nn + 1 - i, nn + 1 - j)
                k <- k + 1
            }
        }

        for (i in 1:nn) {
            matches1 <- match(gn, featureNames(exprs(x)[[i]]))
            pcor[[i]] <- assayData(exprs(x)[[i]])[["exprs"]][matches1, ]
            nnote[i] <- paste("study", i, sep = " ")
        }

        names(pcor) <- names(x)

        icor <- matrix(NA, nuid, np)
        for (i in 1:np) {
            CC1 <- pcor[[ppair[i, 1]]]
            CC2 <- pcor[[ppair[i, 2]]]
            icor[, i] <- .icor(CC1, CC2)
        }

        avg.cc <- apply(icor, 1, mean)

        if (is.null(geneid))
            mmatches <- (1:nuid)[avg.cc == max(avg.cc)]

        score <- avg.cc[mmatches]

        if (is.null(plotype) || plotype == 1) {
            k <- m <- 1
            par(mfrow = c(nn, nn), oma = c(0, 0, 8, 0))

            cp <- alist(... = )
            for (i in 1:nn) {
                m1 <- t(apply(pcor[[i]], 1, .integ.norm))
                x <- m1[mmatches, ]
                cx <- m1[-mmatches, ] %*% matrix(x, length(x), 1)

                cp[[i]] <- cx
            }

            for (i in 1:nn) {
                for (j in 1:nn) {
                  if (k <= np) {
                    if (ppair[k, 1] == i & ppair[k, 2] == j) {
                      cx <- cp[[ppair[k, 1]]]
                      cy <- cp[[ppair[k, 2]]]
                      plot(as.vector(cy), as.vector(cx), xlab = xlab, ylab = ylab, main = main, ...)
                      abline(h = 0)
                      abline(v = 0)
                      k = k + 1
                    }
                  }
                  if (m <= np) {
                    if (ppair1[np + 1 - m, 1] == i & ppair1[np + 1 - m, 2] == j) {
                      pp1 <- j * nn + i
                      pp2 <- c(1:np)[match(as.character(pp1), as.character(indextmp))]
                      score <- icor[, pp2][mmatches]
                      aa <- c(0.5, 0.5)
                      bb <- c(0.7, 0.8)
                      plot(aa, bb, xlim = c(0, 1), ylim = c(0, 1), xlab = NA, ylab = NA, pch = " ", xaxt = "n", yaxt = "n")
                      text(0.47, 0.5, cex = 1.2, paste("CORR = ", as.character(signif(score, digits = 3))), col = 4)
                      m <- m + 1
                    }
                  }
                  if (i == j) {
                    aa <- c(0.5, 0.5)
                    bb <- c(0.7, 0.8)
                    plot(aa, bb, xlim = c(0, 1), ylim = c(0, 1), xlab = NA, ylab = NA, pch = " ", xaxt = "n", yaxt = "n")
                    text(0.5, 0.5, cex = 1.2, paste(labels[i]), col = 1)
                  }
                }
            }

            if (is.null(title))
                title <- paste("Integrated Correlation\n\nGene", gn[mmatches], ": average integrated score is", as.character(signif(score,
                  digits = 3)), sep = " ")
            mtext(title, line = 0.5, cex = 1.2, outer = TRUE)
            par(mfrow = c(1, 1))
        }

        if (is.null(plotype) || plotype == 2) {
            if (is.null(plotype))
                par(ask = TRUE)
            k <- m <- 1
            par(mfrow = c(nn, nn), oma = c(0, 0, 8, 0))

            for (i in 1:nn) {
                for (j in 1:nn) {
                  if (k <= np) {
                    if (ppair[k, 1] == i & ppair[k, 2] == j) {
                      cx <- pcor[[ppair[k, 1]]]
                      cy <- pcor[[ppair[k, 2]]]
                      .ic.plot(cx, cy, mmatches, xlab = xlab, ylab = ylab, main = main, scale = scale, square = square, ...)
                      k = k + 1
                    }
                  }
                  if (m <= np) {
                    if (ppair1[np + 1 - m, 1] == i & ppair1[np + 1 - m, 2] == j) {
                      pp1 <- j * nn + i
                      pp2 <- c(1:np)[match(as.character(pp1), as.character(indextmp))]
                      score <- icor[, pp2][mmatches]
                      aa <- c(0.5, 0.5)
                      bb <- c(0.7, 0.8)
                      plot(aa, bb, xlim = c(0, 1), ylim = c(0, 1), xlab = NA, ylab = NA, pch = " ", xaxt = "n", yaxt = "n")
                      text(0.47, 0.5, cex = 1.2, paste("CORR = ", as.character(signif(score, digits = 3))), col = 4)
                      m <- m + 1
                    }
                  }
                  if (i == j) {
                    aa <- c(0.5, 0.5)
                    bb <- c(0.7, 0.8)
                    plot(aa, bb, xlim = c(0, 1), ylim = c(0, 1), xlab = NA, ylab = NA, pch = " ", xaxt = "n", yaxt = "n")
                    text(0.5, 0.5, cex = 1.2, paste(labels[i]), col = 1)
                  }
                }
            }

            if (is.null(title))
                title <- paste("Integrated Correlation\n\nGene", gn[mmatches], ": average integrated score is", as.character(signif(score,
                  digits = 3)), sep = " ")
            mtext(title, line = 0.5, cex = 1.2, outer = TRUE)
            par(mfrow = c(1, 1))
            if (is.null(plotype))
                par(ask = FALSE)
        }
    }
    if (nn == 2) {
        nn <- length(x)
        if (nn <= 1)
            stop("Number of studies in the mergeExpressionSet should not less than 2.")

        geneid1 <- rep(1, nrow(x@geneStudy))
        for (i in 1:nn) {
            geneid1 <- geneid1 * (x@geneStudy[, i])
        }
        gn <- geneNames(x)[geneid1 == 1]
        nuid <- length(gn)

        pcor <- alist(... = )
        nnote <- rep(NA, nn)

        for (i in 1:nn) {
            matches1 <- match(gn, featureNames(exprs(x)[[i]]))
            pcor[[i]] <- assayData(exprs(x)[[i]])[["exprs"]][matches1, ]

            nnote[i] <- nnote[i] <- paste("study", i, sep = " ")
        }

        CC1 <- pcor[[1]]
        CC2 <- pcor[[2]]
        icor <- .icor(CC1, CC2)

        names(icor) <- gn

        if (!is.null(geneid)) {
            mmatches <- match(geneid, gn)
            if (is.na(mmatches))
                stop("No such common gene id.")
        }

        if (is.null(geneid))
            mmatches <- (1:nuid)[icor == max(icor)]
        if (is.null(main))
            main <- paste("Integrated Correlation\n\nGene", gn[mmatches])

        score <- icor[mmatches]

        if (is.null(plotype) || plotype == 1) {
            m1 <- t(apply(pcor[[1]], 1, .integ.norm))
            m2 <- t(apply(pcor[[2]], 1, .integ.norm))

            x <- m1[mmatches, ]
            y <- m2[mmatches, ]

            cx <- m1[-mmatches, ] %*% matrix(x, length(x), 1)
            cy <- m2[-mmatches, ] %*% matrix(y, length(y), 1)

            par(oma = c(0, 0, 0, 3))
            if (is.na(xlab))
                plot(as.vector(cx), as.vector(cy), ylab = ylab, xlab = paste("CORR = ", as.character(signif(score, digits = 3))), main = main,
                  ...) else plot(as.vector(cx), as.vector(cy), ylab = ylab, xlab = paste("\n", xlab, "\nCORR = ", as.character(signif(score, digits = 3))),
                main = main, ...)
            abline(h = 0)
            abline(v = 0)
        }

        if (is.null(plotype) || plotype == 2) {
            if (is.null(plotype))
                par(ask = TRUE)
            cx = pcor[[1]]
            cy = pcor[[2]]
            par(oma = c(0, 0, 0, 3))
            if (is.na(xlab))
                .ic.plot(cx, cy, mmatches, ylab = ylab, xlab = paste("CORR = ", as.character(signif(score, digits = 3))), main = main,
                  scale = scale, square = square, ...) else .ic.plot(cx, cy, mmatches, ylab = ylab, xlab = paste("\n", xlab, "\nCORR = ", as.character(signif(score, digits = 3))),
                main = main, scale = scale, square = square, ...)
            if (is.null(plotype))
                par(ask = FALSE)
        }
    }
    return()
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
    for (i in 1:n1) {
        for (j in 1:n2) {
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


.hist.mergeCor <- function(x, labels = NULL, main = NULL, xlab = NULL, title = NULL, ...) {
    nn <- length(x@notes)
    if (is.null(labels))
        labels <- x@notes
    if (ncol(x@pairwise.cors) != nn * (nn - 1)/2)
        stop("You need to specify the names for the studies.")
    if (nn <= 1)
        stop("You need to specify more than one study.")

    if (is.null(main))
        main <- NA
    if (is.null(xlab))
        xlab <- NA
    if (nn > 2) {
        np <- nn * (nn - 1)/2
        ppair <- matrix(0, np, 2)
        k <- 1
        for (i in 1:(nn - 1)) {
            for (j in (i + 1):nn) {
                ppair[k, ] <- c(i, j)
                k <- k + 1
            }
        }
        indextmp <- ppair[, 1] * nn + ppair[, 2]

        ppair1 <- matrix(0, np, 2)
        k <- 1
        for (i in 1:(nn - 1)) {
            for (j in (i + 1):nn) {
                ppair1[k, ] <- c(nn + 1 - i, nn + 1 - j)
                k <- k + 1
            }
        }

        cp <- x@pairwise.cors
        k <- m <- 1
        par(mfrow = c(nn, nn), oma = c(0, 0, 5, 0))


        for (i in 1:nn) {
            for (j in 1:nn) {
                if (k <= np) {
                  if (ppair[k, 1] == i & ppair[k, 2] == j) {
                    hist(cp[, k], xlab = xlab, main = main, ...)
                    k = k + 1
                  }
                }
                if (m <= np) {
                  if (ppair1[np + 1 - m, 1] == i & ppair1[np + 1 - m, 2] == j) {
                    pp1 <- j * nn + i
                    pp2 <- c(1:np)[match(as.character(pp1), as.character(indextmp))]
                    hist(cp[, pp2], xlab = xlab, main = main, ...)
                    m <- m + 1
                  }
                }
                if (i == j) {
                  aa <- c(0.5, 0.5)
                  bb <- c(0.7, 0.8)
                  plot(aa, bb, xlim = c(0, 1), ylim = c(0, 1), xlab = NA, ylab = NA, pch = " ", xaxt = "n", yaxt = "n")
                  text(0.5, 0.5, cex = 1.2, paste(labels[i]), col = 1)
                }
            }
        }
        if (is.null(title))
            title <- "Integrated Correlation"
        mtext(title, line = 0.5, cex = 1.2, outer = TRUE)
        par(mfrow = c(1, 1))
    }
    if (nn == 2) {
        if (is.na(main))
            main <- "Integrated Correlation"
        cp <- x@pairwise.cors
        hist(cp, main = main, xlab = xlab, ...)
    }

    return()
}


.dens.mergeExpressionSet <- function(x, main = NA, x.legend = NULL, y.legend = NULL, cex.legend = NULL, title = NULL, method = NULL,
    ...) {
    nn <- length(x)
    if (nn <= 1)
        stop("Number of studies in the mergeExpressionSet should not less than 2.")

    nnote <- rep(NA, nn)

    for (i in 1:nn) {
        nnote[i] <- nnote[i] <- paste("study", i, sep = " ")
    }

    geneid <- rep(1, nrow(x@geneStudy))
    for (i in 1:nn) {
        geneid <- geneid * (x@geneStudy[, i])
    }
    geneuid <- geneNames(x)[geneid == 1]
    nuid <- length(geneuid)
    rtn <- alist(... = )

    if (nn > 2) {
        np <- nn * (nn - 1)/2
        ppair <- matrix(0, np, 2)
        k <- 1
        for (i in 1:(nn - 1)) {
            for (j in (i + 1):nn) {
                ppair[k, ] <- c(i, j)
                k <- k + 1
            }
        }

        ppair1 <- matrix(0, np, 2)
        k <- 1
        for (i in 1:(nn - 1)) {
            for (j in (i + 1):nn) {
                ppair1[k, ] <- c(nn + 1 - i, nn + 1 - j)
                k <- k + 1
            }
        }


        exp1 <- matrix(NA, 10000, np)
        obs <- matrix(NA, nuid, np)

        for (i in 1:np) {
            matches1 <- match(geneuid, featureNames(exprs(x)[[ppair[i, 1]]]))
            matches2 <- match(geneuid, featureNames(exprs(x)[[ppair[i, 2]]]))
            exprs1 <- assayData(exprs(x)[[ppair[i, 1]]])[["exprs"]][matches1, ]
            exprs2 <- assayData(exprs(x)[[ppair[i, 2]]])[["exprs"]][matches2, ]

            obs[, i] <- .icor(exprs1, exprs2)
            exp1[, i] <- .nullicornorm(exprs1, exprs2)
        }

        k <- m <- 1
        par(mfrow = c(nn, nn), oma = c(0, 0, 8, 0))

        for (i in 1:nn) {
            for (j in 1:nn) {
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
                      legend(x.legend, y.legend, paste(legend.text), col = 1:2, lty = 1, cex = cex.legend)
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
        legend(x.legend, y.legend, paste(legend.text), col = 1:2, lty = 1, cex = cex.legend)
        rtn <- d1
    }
    return(rtn)

}

.show.mergeExpressionSet <- function(object) {
    x = exprs(object)
    y = list()
    for (i in 1:length(x)) {
        show(x[[i]])
    }
    return()
}



.summary.mergeExpressionSet <- function(object) {
    report <- alist(... = )
    nn <- length(object)

    nnote <- rep(NA, nn)

    for (i in 1:nn) {
        nnote[i] <- paste("study", i, sep = " ")
    }
    report[[1]] <- matrix(NA, 2, nn)
    report[[1]][1, ] <- names(object)
    for (i in 1:nn) {
        report[[1]][2, i] <- nrow(assayData(exprs(object)[[i]])[["exprs"]])
    }

    report[[2]] <- matrix(NA, 2, nn)
    report[[2]][1, ] <- names(object)
    for (i in 1:nn) {
        report[[2]][2, i] <- ncol(assayData(exprs(object)[[i]])[["exprs"]])
    }

    report[[3]] <- nnote

    names(report) <- c("Number of Genes in Each Study", "Number of Samples in Each Study", "Notes")
    return(report)
}




