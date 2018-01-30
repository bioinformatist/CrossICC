#' Title
#'
#' @param x
#' @param gene.signature
#' @param geneset2gene
#' @param filter.cutoff
#' @param fdr.cutoff
#' @param fs
#'
#' @return
#' @export
#'
#' @examples
predictor <- function(x, fs, gene.signature, geneset2gene, filter.cutoff = 1, fdr.cutoff = 1) {
  filtered.x <- m.f.s(list(x), fdr.cutoff = fdr.cutoff, filter.cutoff = filter.cutoff, perform.mad = FALSE)[[2]][[1]]
  genewprobe <- gene.signature
  names(genewprobe) <- gene.signature
  c1 <- runFAIME(filtered.x, genewprobe, geneset2gene, na.last = "keep", weightRank = FALSE)
  c2 <- fs

  FAIME.scale <- function(fs) {
    fs.scale <- t(scale(t(fs)))
    fs.scale <- replace(fs.scale, fs.scale < -2, -2)
    fs.scale <- replace(fs.scale, fs.scale > 2, 2)
    fs.scale
  }

  c1 <- FAIME.scale(c1)
  c2 <- FAIME.scale(c2)
  # for(j in 1:ncol(c1)){
  #   pcentroidj <- c1[,j]
  #   cor_row <- c()
  #   for(i in 1:ncol(c2)){
  #     vcentroidi <- c2[,i]
  #     corj <- cor(pcentroidj, vcentroidi, use="complete", method="pearson")
  #     cor_row[length(cor_row)+1] <- corj$estimate
  #   }
  #   cor <- rbind(cor,cor_row)
  # }
  # rownames(cor) <- colnames(c1)
  # colnames(cor) <- colnames(c2)
  # mapply(cor, c1, c2)
  cor(c1, c2)
}
