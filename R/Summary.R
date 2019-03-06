#' Generate two Score Groups
#'
#' Creates a grouping variable which divides the sample in two groups (high and low scorers)
#' of roughly equal size, without taking into account persons with extreme scores.
#'
#' The score groups are used for tests of item homogeneity.
#' @param     dat.items A data frame with the responses to the items.
#' @param     label If TRUE the levels of the group factor are named according to the split used,
#' if FALSE (default) the group factor has levels 1 and 2.
#' @export
#' @return   Score group variable, a factor with two levels.
score_groups <- function(dat.items, label = FALSE){
  x <- apply(dat.items,1,sum,na.rm=T)
  sgrp <- rep(NA, length(x))
  m <- sum(apply(dat.items, 2, max, na.rm = TRUE))
  ss <- as.data.frame(table(x[x > 0 & x < m]))
  xx <-  abs(cumsum(ss$Freq) - sum(ss$Freq)/2)
  ind <- which(xx == min(xx))
  if (length(ind)==2) {
    if (table(x)[1] <= table(x)[length(table(x))]) {
      ind <- ind[2]
    } else {
      ind <- ind[1]
    }
  }
  ss$Var1 <- as.numeric(levels(ss$Var1))
  sgrp[x <= ss$Var1[ind]] <- 1
  sgrp[x > ss$Var1[ind]] <- 2
  sgrp <- factor(sgrp)
  if (label == TRUE) {
    levels(sgrp) <- c(paste(0,"-",ss$Var1[ind],sep=""),paste(ss$Var1[ind]+1,"-",m,sep=""))
  }
  sgrp
}




