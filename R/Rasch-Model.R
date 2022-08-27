W_thresh <- function(x){
  mi <- x
  W1 <- diag(sum(mi)-1)
  W1 <- rbind(W1,rep(0,sum(mi)-1))
  W1[sum(mi),cumsum(mi[-length(mi)])] <- -1
  W1
}

#' Computation of Item Targets for Polytomous Models
#'
#' The item target is the value of the person parameter where
#' item information is maximized.
#' @param  obj An object of class "eRm" (but not "dRm"), a fitted partial
#' credit model using  the function PCM in package eRm or of class "pcmodel"
#' (from package psychotools).
#' @return vector with item targets.
#' @author Marianne Mueller
#' @import eRm
#' @importFrom psychotools pcmodel threshpar
#' @importFrom stats optimize
#' @export
#' @examples
#' \dontrun{
#'  pc.mod <- PCM(desc2[, 5:14])
#'  item_target(pc.mod)
#' }
item_target <- function(obj){
  if(class(obj)[1]=="pcmodel") obj$model <- "pcmodel"
  if(!(obj$model%in%c("pcmodel","PCM"))) stop("Item targets are computed only for polytomous models estimated with pcmodel or PCM")
  if(obj$model=="pcmodel")
    betasum.l <- lapply(threshpar(obj), cumsum)
  else {
    thresh1 <- thresholds(obj)[[3]][[1]][, -1] - mean(thresholds(obj)[[3]][[1]][, -1], na.rm=T)
    betasum.l <- lapply(as.list(as.data.frame(t(thresh1))), function(x) cumsum(na.omit(x)))
  }
  var.X <- function(x) {
    function(theta) {
      xvec <- 0:length(x)
      pp <- (exp(theta*xvec)*exp(-c(0, x)))/rep((exp(theta*xvec)%*%exp(-c(0, x))),length(xvec))
      (xvec - rep(pp%*%xvec,length(xvec)))^2%*%pp
    }
  }
  varf <- lapply(betasum.l, var.X)
  target <- unlist(sapply(varf, optimize, interval = c(-4, 4), maximum = T)[1, ])
  cat("\n")
  cat("Item targets:","\n\n")
  print(target, 3)
  cat("\n\n")
  invisible(target)
}

#' Conditional Likelihood Ratio Tests (CLR)
#'
#' The conditional likelihood ratio tests compare item parameters in low and high score groups
#' for an overall test of homogeneity, and in groups defined by the levels of exogenous factors
#' for tests of no differential item functioning (DIF).
#' @author Marianne Mueller
#' @references Andersen, E.B. (1973). A goodness of fit test for the Rasch model. \emph{Psychometrika}, 38, 123-140.
#' @param dat.items A data frame with the responses to the items.
#' @param dat.exo  A single factor or a data frame consisting of one or more exogenous factor variables.
#' @param model If model="RM" a Rasch model will be fitted,
#' if model="PCM" a partial credit model for polytomous items is used.
#' @return matrix with test statistics, df and p values.
#' @import eRm
#' @importFrom stats symnum complete.cases pchisq rnorm runif
#' @importFrom psychotools pcmodel
#' @export
#' @examples #CLR overall test and test of  no DIF for agegrp and sex
#' clr_tests(amts[,4:13],amts[,2:3])
clr_tests <- function(dat.items, dat.exo=NULL, model = c("RM","PCM")) {
  ok <- complete.cases(dat.items)
  dat.i <- dat.items[ok, ]
  sgrp <- score_groups(dat.i)
  model <- match.arg(model)
  if (model == "RM") {
    mod0 <- RM(dat.i)
    lr0 <- LRtest(mod0,sgrp)
    if (is.null(dat.exo)){
      mm <- round(t(c(clr=lr0$LR,df=lr0$df,pvalue=lr0$pvalue)),digits=3)
      row.names(mm)[1]="overall"
    } else {
      if (!is.data.frame(dat.exo)) {
        gname <- deparse(substitute(dat.exo))
        dat.exo <- data.frame(dat.exo)
        names(dat.exo) <- gname
      }
      ok1 <-  complete.cases(cbind(dat.items, dat.exo))
      mod1 <- RM(dat.items[ok1,])
      lrexo <- function(x){
        lr <- LRtest(mod1,x)
        out <- c(lr$LR,lr$df,lr$pvalue)
      }
      mm <- round(rbind(c(clr=lr0$LR,df=lr0$df,pvalue=lr0$pvalue),t(apply(dat.exo[ok1,,drop=F],2,lrexo))),digits=3)
      row.names(mm) <- c("overall", names(dat.exo[ok1, , drop = F]))
    }
    symp <- symnum(mm[,3], cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),symbols = c(" ***", " **", " *", " .", " "))
    mm<- noquote(cbind(mm,sig=symp))
  } else {
    mod0 <- pcmodel(dat.i,hessian=F)
    clr0 <- -2*(mod0$loglik - (pcmodel(dat.i[sgrp==1, ],hessian=F)$loglik + pcmodel(dat.i[sgrp==2, ],hessian=F)$loglik))
    if (is.null(dat.exo)){
      mm <- round(t(c(clr=clr0,df=mod0$df)),digits=3)
      row.names(mm)[1]="overall"
    } else {
      if (!is.data.frame(dat.exo)) {
        gname <- deparse(substitute(dat.exo))
        dat.exo <- data.frame(dat.exo)
        names(dat.exo) <- gname
      }
      ok1 <-  complete.cases(cbind(dat.items, dat.exo))
      clrhomo <- function(exo, data){
        sum(sapply(split(data, exo, drop = TRUE),function(x){pcmodel(x,hessian=F)$loglik}))
      }
      mod1 <- pcmodel(dat.items[ok1, ],hessian=F)
      ll <- mod1$loglik
      clr <- c(clr0, -2*(ll - apply(dat.exo[ok1, ,drop = F],2,clrhomo, dat.items[ok1, ])))
      df <- c(mod0$df,mod1$df*(sapply(dat.exo[ok1, , drop = F],nlevels) - 1))
      mm <- cbind(clr, df)
      row.names(mm) <- c("overall", names(dat.exo[ok1, , drop = F]))
    }
    pvalue <- apply(mm,1,function(x){1-pchisq(x[1],x[2])})
    symp <- symnum(pvalue, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),symbols = c(" ***", " **", " *", " .", " "))
    mm <- noquote(cbind(clr = format(mm[, 1], digits = 3), df = format(mm[, 2], digits = 1), pvalue = format(pvalue, digits = 2), sig = symp))
  }
  cat("\n")
  cat("Conditional Likelihood Ratio Tests:")
  cat("\n")
  mm
}
