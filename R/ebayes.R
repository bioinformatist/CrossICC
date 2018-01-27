ebayes <- function(eSet.subset, class, cutoff = 0.1){
  k <- length(unique(class))
  design <- model.matrix(~ 0 + factor(class))
  colnames(design) <- paste("K", 1:k, sep = "")
  K <- colnames(design)
  fit <- limma::lmFit(eSet.subset, design)

  x <- c()
  for(i in 1:(k - 1)){
    for(j in (i + 1):k){
      x <- c(x, paste(K[j], K[i], sep = "-"))
    }
  }
  for(i in 1:k){
    x <- c(x, paste(K[i], paste("(", paste(K[-i], collapse="+"), ")", "/", k-1, sep=""), sep="-"))
  }

  contrast.matrix <- limma::makeContrasts(contrasts = x,levels = design)
  fit2 <- limma::contrasts.fit(fit,contrast.matrix)
  fit2 <- limma::eBayes(fit2)
  gene.signature <- limma::topTable(fit2, number = 20000, adjust.method = 'BH', p.value = cutoff)
  r <- list(gene.signature)
  for(i in (k*(k - 1) / 2 + 1):(k * (k - 1) / 2 + k)){
    geneset.i <- limma::topTable(fit2, number = 20000, coef = i, adjust.method = 'BH', p.value = cutoff)
    r[[length(r)+1]] <- geneset.i
  }
  ml <- r[-1]
  names(ml) <- K
  geneset2gene <- do.call(rbind, lapply(names(ml), function(x) data.frame(rep(x, nrow(ml[[x]])), row.names(ml[[x]]))))
  list(full.m = r[[1]], geneset2gene = geneset2gene)
}
