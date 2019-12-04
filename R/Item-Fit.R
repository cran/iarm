expvar <- function(data, coeff){
  rv <- rowSums(data, na.rm = TRUE)
  rv[rv == 0] <- NA
  k <- dim(data)[2]
  n <- dim(data)[1]
  mi <- apply(data, 2, max, na.rm = TRUE)
  m <- sum(mi)
  Er <- matrix(rep(0, m*k), ncol = k)
  Varx <- matrix(rep(0, m*k), ncol = k)
  gr <- psychotools::elementary_symmetric_functions(coeff)[[1]]
  grminus <- vector("list", k)
  for (i in 1:k) {
    grminus[[i]] <- psychotools::elementary_symmetric_functions(coeff[-i])[[1]]
  }
  for (r in 1:(m-1)) {
    for (i in 1:k) {
      pscore <- rep(0, mi[i])
      for (x in 1:mi[i]) {
        if (x <= r & m - mi[i] + x >= r) {
          pscore[x] <- exp(-coeff[[i]][x])*(grminus[[i]])[r + 1 - x]/gr[r + 1]
        }
      }
      pscore <- c(1 - sum(pscore), pscore)
      Er[r, i] <- (0:mi[i])%*%pscore
      Varx[r, i] <- ((0:mi[i])-Er[r, i])^2%*%pscore
    }
  }
  Er[m, ] <- 1
  Ehat <- Er[rv, ]
  VarX <- Varx[rv, ]
  colnames(Ehat) <- paste("Item", 1:k, sep="")
  colnames(VarX) <- paste("Item", 1:k, sep="")
  list(Ehat=Ehat, VarX=VarX)
}

#' Observed and Expected Item Mean Scores
#'
#' Homogeneity of item responses in the low and high score groups is analyzed by looking at observed and expected item mean scores
#' together with standardized residuals. If the Andersen's CLR test has shown some evidence against homogeneity,
#' this comparison can indicate which items might be responsible.
#' @param  object An object of class "Rm", a fitted Rasch model or partial
#' credit model using  the functions RM or PCM in package eRm, or an object of class "pcmodel",
#'  a fitted partial credit model using the function pcmodel in package psychotools.
#' @return list with observed and expected mean scores together with standardized residuals for the two score groups.
#' @author Marianne Mueller
#' @import eRm
#' @importFrom psychotools pcmodel threshpar
#' @importFrom stats symnum coef
#' @export
#' @examples
#' rm.mod <- RM(amts[,4:13])
#' item_obsexp(rm.mod)
#' pc.mod <- PCM(desc2[,5:14])
#' item_obsexp(pc.mod)
item_obsexp <- function(object){
  if (!any("Rm"%in%class(object),class(object) =="pcmodel")) stop("object must be of class Rm or pcmodel!")
  if(class(object)[1]=="pcmodel") object$model <- "pcmodel"
  if (object$model == "RM") {
    if (sum(is.na(object$X)) > 0){
      message("Model was refitted with only complete cases")
      dat.items <- na.omit(object$X)
      object <- RM(dat.items)
    }
    X <- object$X
    k <- dim(X)[2]
    coeff <- (-1)*coef(object)
    mi <- rep(1, k)
  } else {
    if (object$model == "PCM"){
      if (sum(is.na(object$X)) > 0){
        message("Model was refitted with only complete cases")
        dat.items <- na.omit(object$X)
        object <- PCM(dat.items)
      }
      X <- object$X
      k <- dim(X)[2]
      mi <- apply(X, 2, max, na.rm = TRUE)
      thresh1 <- thresholds(object)[[3]][[1]][, -1] - mean(thresholds(object)[[3]][[1]][, 1])
      coeff <- lapply(as.list(as.data.frame(t(thresh1))), cumsum)
    } else {
      if (sum(is.na(object$data)) > 0){
        message("Model was refitted with only complete cases")
        dat.items <- na.omit(object$data)
        object <- pcmodel(dat.items)
      }
      X <- object$data
      k <- dim(X)[2]
      mi <- apply(X, 2, max, na.rm = TRUE)
      thresh1 <- coef(threshpar(object),type="matrix")
      coeff <- lapply(as.list(as.data.frame(t(thresh1))), cumsum)
    }
  }
  m <- sum(mi)
  score <- apply(X, 1, sum, na.rm=T)
  sgrp <- score_groups(X)
  resgrp <- function(x){
    EV <- expvar(X[sgrp==x & score > 0 & score < m, ], coeff)
    resm <- cbind(apply(X[sgrp==x & score > 0 & score < m, ], 2, mean, na.rm=T), colMeans(EV$Ehat, na.rm=T),
                  length(EV$VarX[!is.na(EV$VarX[, 1]), 1]) * (apply(X[sgrp==x & score > 0 & score < m, ], 2, mean, na.rm=T)
                                                              - colMeans(EV$Ehat, na.rm=T))/sqrt(colSums(EV$VarX, na.rm = T)))
    colnames(resm)= c("mean obs", "mean exp", "std.res")
    symp <- symnum(resm[, 3], cutpoints=c(-100, -3, -2, 2, 3, 100), symbols=c("--", "-", " ", "+", "++"))
    resm <- noquote(cbind(format(resm, digits = 3), sig = symp))
    resm
  }
  result <- lapply(c(1,2),resgrp)
  cat("Score group 1:","\n")
  print(result[[1]])
  cat("\n")
  cat("Score group 2:","\n")
  print(result[[2]])
  cat("\n")
  invisible(result)
}


#' Item Outfit and Infit Statistics
#'
#' To avoid bias observed item responses are compared to expected responses under
#' the conditional distribution of responses given the total score. This leads to standardized  residuals
#' which can be summarized to outfit and infit statistics in the usual way.
#' @param  object An object of class "Rm", a fitted Rasch model or partial
#' credit model using  the functions RM or PCM in package eRm, or an object of class "pcmodel",
#'  a fitted partial credit model using the function pcmodel in package psychotools.
#' @param se If TRUE the standard errors will be included.
#' @param p.adj Correction method for multiple testing. The methods are "BH","holm", "hochberg", "hommel", "bonferroni", "BY", "none". See \code{\link{p.adjust}}.
#' @details The fit statistics and their standard errors are calculated as described in Christensen et al.
#' P values are are based on the normal distribution of the standardized fit statistics.
#' @author Marianne Mueller
#' @import eRm
#' @importFrom psychotools pcmodel threshpar
#' @importFrom stats coef pnorm na.exclude p.adjust
#' @export
#' @references Christensen, K. B. , Kreiner, S. & Mesbah, M. (Eds.)
#' \emph{Rasch Models in Health}. Iste and Wiley (2013), pp. 86 - 90.
#'
#'  Kreiner, S. & Christensen, K. B. (2011) Exact evaluation of Bias in Rasch model residuals.
#' \emph{Advances in Mathematics Research}, 12, 19-40.
#' @return an object of class outfit containing:
#' \item{outfit}{outfit statistics}
#' \item{outfit.se}{standard errors of outfit statistics}
#' \item{out.pvalue}{p values of outfit statistics}
#' \item{infit}{infit statistics}
#' \item{infit.se}{standard errors of infit statistics}
#' \item{in.pvalue}{p values of infit statistics}
#' @examples
#' rm.mod <- RM(amts[,4:13])
#' out_infit(rm.mod)
out_infit <- function(object, se=TRUE, p.adj= c("BH","holm", "hochberg", "hommel", "bonferroni", "BY", "none")){
  if (!any("Rm"%in%class(object),class(object) =="pcmodel")) stop("object must be of class Rm or pcmodel!")
  if(class(object)[1]=="pcmodel") object$model <- "pcmodel"
  padj <- match.arg(p.adj)
  if (object$model=="RM") {
    if (sum(is.na(object$X)) > 0){
      message("Model was refitted with only complete cases")
      dat.items <- na.omit(object$X)
      object <- RM(dat.items)
    }
    X <- object$X
    rv <- rowSums(X, na.rm = TRUE)
    k <- dim(X)[2]
    koeff <- (-1)*coef(object)
    mi <- rep(1,k)
  } else {
    if (object$model == "PCM") {
      if (sum(is.na(object$X)) > 0){
        message("Model was refitted with only complete cases")
        dat.items <- na.omit(object$X)
        object <- PCM(dat.items)
      }
      X <- object$X
      rv <- rowSums(X, na.rm = TRUE)
      k <- dim(X)[2]
      mi <- apply(X, 2, max, na.rm = TRUE)
      thresh1 <- thresholds(object)[[3]][[1]][, -1] - mean(thresholds(object)[[3]][[1]][, 1])
      koeff <- lapply(as.list(as.data.frame(t(thresh1))), cumsum)
    } else {
      if (sum(is.na(object$data)) > 0){
        message("Model was refitted with only complete cases")
        dat.items <- na.omit(object$data)
        object <- pcmodel(dat.items)
      }
      X <- object$data
      rv <- rowSums(X, na.rm = TRUE)
      k <- dim(X)[2]
      mi <- apply(X, 2, max, na.rm = TRUE)
      thresh1 <- coef(threshpar(object),type="matrix")
      koeff <- lapply(as.list(as.data.frame(t(thresh1))), cumsum)
    }
  }
  m <- sum(mi)
  rvstrich <- rv[rv > 0 & rv < m]
  ER <- matrix(rep(0, (m-1)*k), ncol = k)
  VarZ <- matrix(rep(0, (m-1)*k), ncol = k)
  VarX.R <- matrix(rep(0, (m-1)*k), ncol = k)
  gr <- elementary_symmetric_functions(koeff)[[1]]
  grminus <- vector("list",k)
  for (i in 1:k) {
    grminus[[i]] <- elementary_symmetric_functions(koeff[-i])[[1]]
  }
  for (r in 1:(m-1)) {
    for (i in 1:k) {
      pscore <- rep(0, mi[i])
      for (x in 1:mi[i]) {
        if (x <= r  & x >= r-m+mi[i]) {
          pscore[x] <- exp(-koeff[[i]][x]) * (grminus[[i]])[r+1-x] / gr[r+1]
        }
      }
      pscore <- c(1-sum(pscore), pscore)
      ER[r,i] <- as.numeric((0:mi[i]) %*% pscore)
      VarX.R[r,i] <- as.numeric(((0:mi[i])-ER[r,i]) ^ 2 %*% pscore)
      VarZ[r,i] <- as.numeric((((0:mi[i]) - ER[r,i])^2/VarX.R[r, i] - 1)^2 %*% pscore)
    }
  }
  Ehat <- ER[rvstrich,]
  VarX <- VarX.R[rvstrich,]
  Z.vi2 <- ((X[rv > 0 & rv < m,] - Ehat) ^ 2) / VarX
  Outfit = colMeans(Z.vi2, na.rm=T)
  X <- as.data.frame(X)
  nir <- t(sapply(split(X[rv > 0 & rv < m, ], factor(rvstrich, levels = 1:(m-1))),function(x){apply(x, 2, function(y){length(na.exclude(y))})}))
  if (se==TRUE) {
    Outfit.se <- sqrt(colSums((nir/colSums(nir)^2) * VarZ))
    fit <- cbind(Outfit,Outfit.se)
    pwert <- function(x){ifelse(x[, 1] > 1, 2*(1 - pnorm((x[, 1] - 1)/x[, 2])), 2*(pnorm((x[, 1] - 1)/x[, 2])))}
    out.pwert <- p.adjust(pwert(fit),method=padj, n= 2*k)
  }
  Wri <- VarX.R / matrix(rep(colSums(nir * VarX.R), m-1), ncol = k, byrow = T)
  Wvi <- Wri[rvstrich, ]
  Infit <-  colSums(Wvi * Z.vi2, na.rm = T)
  if (se==TRUE) {
    Infit.se <- sqrt(colSums(nir * Wri^2 * VarZ,na.rm = T))
    names(Infit)=colnames(X)
    fit <- cbind(Infit,Infit.se)
    in.pwert <- p.adjust(pwert(fit),method=padj, n=2*k)
    result <- list(Outfit=Outfit, Outfit.se=Outfit.se, out.pvalue=out.pwert, Infit=Infit, Infit.se=Infit.se, in.pvalue=in.pwert,padj)
    names(result)[3] <- paste("pvalue",padj,sep=".")
    names(result)[6] <- paste("pvalue",padj,sep=".")
  } else {
    result <- list(Outfit=Outfit, Infit=Infit)
  }
  class(result) <- "outfit"
  result
}

#' Print Method for the Output of out_infit
#' @param x object of class outfit.
#' @param ... arguments passed to other functions.
#' @export
print.outfit <- function(x, ...){
  if(length(x)==2) {
    print(data.frame(Outfit=x[[1]],Infit=x[[2]]),digits=3)
    cat("\n")
  } else {
    symp <- function(x){symnum(x,cutpoints=c(0,0.001,0.01,0.05,0.1,1),symbols=c("***   "," **  "," *   "," .   ","    "))}
    tab2 <- noquote(cbind(Outfit=round(x[[1]],digits=3),se=round(x[[2]],digits=3), pvalue=round(x[[3]],digits=3),sig=symp(x[[3]]),
                                 Infit=round(x[[4]],digits=3),se=round(x[[5]],digits=3),pvalue=round(x[[6]],digits=3) ,sig=symp(x[[6]])))
    cat("\n")
    print(tab2)
    cat("\n")
    cat("P value adjustment:", x[[7]])
  }
}

# #' Computes Outfit and Infit Statistics for a Bootstrap Sample
# #' Model fit with raschmodel or pcmodel (much faster then eRm)
# #' @importFrom psychotools elementary_symmetric_functions
# #' @importFrom psychotools raschmodel itempar pcmodel threshpar
# #' @param data A dataframe with the responses to the items
# #' @param model If model="RM" a Rasch model will be fitted,
# #' if model="PCM" or "pcmodel" a partial credit model for polytomous items is used.
# #' @return vector with Outfit and Infit
outin_boot <- function(data, model= c("RM","PCM","pcmodel")){
  rv <- rowSums(data, na.rm = TRUE)
  k <- dim(data)[2]  # Anzahl Items
  mode <- match.arg(model)
  if (mode == "RM") {
    koeff <- itempar(raschmodel(data, hessian = F))
    mi <- rep(1,k)
  }
  else {
    mi <- apply(data, 2, max, na.rm = TRUE)
    thresh1 <- threshpar(pcmodel(data, hessian = F, nullcats="ignore"))
    koeff <- lapply(thresh1, cumsum)
  }
  m <- sum(mi)
  rvstrich <- rv[rv > 0 & rv < m]
  ER <- matrix(rep(0, (m-1)*k), ncol = k)
  VarZ <- matrix(rep(0, (m-1)*k), ncol = k)
  VarX.R <- matrix(rep(0, (m-1)*k), ncol = k)
  gr <- elementary_symmetric_functions(koeff)[[1]]
  grminus <- vector("list",k)
  for (i in 1:k) {
    grminus[[i]] <- elementary_symmetric_functions(koeff[-i])[[1]]
  }
  for (r in 1:(m-1)) {
    for (i in 1:k) {
      pscore <- rep(0, mi[i])
      for (x in 1:mi[i]) {
        if (x <= r  & x >= r-m+mi[i]) {
          pscore[x] <- exp(-koeff[[i]][x]) * (grminus[[i]])[r+1-x] / gr[r+1]
        }
      }
      pscore <- c(1-sum(pscore), pscore)
      ER[r,i] <- as.numeric((0:mi[i]) %*% pscore)
      VarX.R[r,i] <- as.numeric(((0:mi[i])-ER[r,i]) ^ 2 %*% pscore)
      VarZ[r,i] <- as.numeric((((0:mi[i]) - ER[r,i])^2/VarX.R[r, i] - 1)^2 %*% pscore)
    }
  }
  Ehat <- ER[rvstrich,]
  VarX <- VarX.R[rvstrich,]
  Z.vi2 <- ((data[rv > 0 & rv < m,] - Ehat) ^ 2) / VarX
  Outfit = colMeans(Z.vi2, na.rm=T)
  nir <- apply(data[rv > 0 & rv < m, ],2, function(x){tapply(x,factor(rvstrich,levels=1:(m-1)),function(y){length(na.exclude(y))})})
  nir[is.na(nir)] <- 0
  Outfit.se <- sqrt(colSums((nir/colSums(nir)^2) * VarZ))
  Wri <- VarX.R / matrix(rep(colSums(nir * VarX.R), m-1), ncol = k, byrow = T)
  Wvi <- Wri[rvstrich, ]
  Infit <-  colSums(Wvi * Z.vi2, na.rm = T)
  Infit.se <- sqrt(colSums(nir * Wri^2 * VarZ,na.rm = T))
  c((Outfit-1)/Outfit.se, (Infit-1)/Infit.se)
}

#' Computes Bootstrapping P Values for Outfit and Infit Statistics
#'@param object  an object of class "Rm" (output of RM or PCM) or class "pcmodel"
#'@param B Number of replications.
#'@param p.adj Correction method for multiple testing. The methods are "BH","holm", "hochberg", "hommel", "bonferroni", "BY", "none". See \code{\link{p.adjust}}.
#'@export
#'@import eRm
#'@importFrom psychotools raschmodel threshpar
#'@importFrom stats p.adjust
#'@importFrom utils capture.output
#'@return object of class bootfit with outfit and infit statistics and corresponding p values.
boot_fit <- function(object,B, p.adj= c("BH","holm", "hochberg", "hommel", "bonferroni", "BY", "none")){
  if(class(object)[1]=="pcmodel") object$model <- "pcmodel"
  padj <- match.arg(p.adj)
  if (object$model=="RM") {
    if (sum(is.na(object$X)) > 0){
      message("Model was refitted with only complete cases")
      dat.items <- na.omit(object$X)
      object <- RM(dat.items)
    }
    X <- object$X
    k <- dim(X)[2]
    koeff <- (-1)*coef(object)
  } else {
    if (object$model == "PCM") {
      if (sum(is.na(object$X)) > 0){
        message("Model was refitted with only complete cases")
        dat.items <- na.omit(object$X)
        object <- PCM(dat.items)
      }
      X <- object$X
      k <- dim(X)[2]
      koeff <- thresholds(object)[[3]][[1]][, -1] - mean(thresholds(object)[[3]][[1]][, 1])
      koeff2 <- vector("list",k)
      for (i in 1:k) koeff2[[i]] <- koeff[i,]
    } else {
      if (sum(is.na(object$data)) > 0){
        message("Model was refitted with only complete cases")
        dat.items <- na.omit(object$data)
        object <- pcmodel(dat.items)
      }
      X <- object$data
      k <- dim(X)[2]
      koeff <-  coef(threshpar(object),type="matrix")
      koeff2 <- threshpar(object)
    }
  }
  Fr1 <- out_infit(object)
  Fr0 <- c((Fr1$Outfit-1)/Fr1$Outfit.se, (Fr1$Infit-1)/Fr1$Infit.se)
  invisible(capture.output(persons <- PP_gpcm(X,t(koeff),slopes=rep(1,k),type="wle")[[1]][[1]][,1]))
  outin <- matrix(rep(NA,2*k*B),ncol=2*k)
  condp <- function(x,x0){
    if (x0 <= 0) p <- sum(x <= x0)/sum(x <= 0)
    else p <- sum(x >= x0)/sum(x >= 0)
    return(p)
  }
  cat("\n Number of bootstrap samples:  ")
  b <- 1
  while (b < B+1){
    if (object$model=="RM") {
      xstar <- sim.rasch(persons,koeff)
    } else {
      xstar <- sim.poly(persons,koeff2)
      mincat <- apply(xstar,2,min,na.rm = TRUE)
      if (any(mincat > 0)) next
    }
    outin[b,] <- outin_boot(xstar, object$model)

    if (object$model=="RM") {
      if (b/50==trunc(b/50)) cat(b,", ",sep="")
    } else {
      if (b/10==trunc(b/10)) cat(b,", ",sep="")
    }
    b <- b+1
  }
  cat("\n \n")
  pvalue <- apply(rbind(outin,Fr0),2,function(x){condp(x[-(B+1)],x[B+1])})
  pvalue <- p.adjust(pvalue,method=padj, n= 2*k)
  result <- cbind(Outfit=Fr1[[1]],pvalue=pvalue[1:k], Infit=Fr1[[4]],pvalue=pvalue[(k+1):(2*k)])
  result <- list(result,padj)
  names(result)[2] <- "adjust"
  class(result) <- "bootfit"
  result
}

#' Print Method for the Output of boot_fit
#' @param x object of class bootfit.
#' @param ... arguments passed to other functions.
#' @export
print.bootfit <- function(x,...){
  symp <- function(x){symnum(x,cutpoints=c(0,0.001,0.01,0.05,0.1,1),symbols=c("***   "," **  "," *   "," .   ","    "))}
  tab2 <- noquote(cbind(Outfit=round(x[[1]][,1],digits=3), pvalue=round(x[[1]][,2],digits=3),sig=symp(x[[1]][,2]),
                        Infit=round(x[[1]][,3],digits=3),pvalue=round(x [[1]][,4],digits=3) ,sig=symp(x[[1]][,4])))
  cat("\n")
  print(tab2)
  cat("\n")
  cat("P value adjustment:", x[[2]])
}


# Original version William Revelle
sim.poly <- function (persons, thresh) {
  if (length(persons) == 1) {
    faehig <- rnorm(persons)
    n <- persons
  } else {
    faehig <- persons
    n <- length(persons)
  }
  mi <- sapply(thresh,length)
  k <- length(mi)
  mtvec <- sequence(mi)
  aa <- t(as.data.frame(lapply(split(mtvec,rep(1:k,mi)),function(ii){c(ii,rep(NA,length.out=max(mi)-length(ii)))})))
  aa <- cbind(rep(0,k),aa)

  psi.l=lapply(thresh, function(x){(-1)*cumsum(x)})
  bb=t(as.data.frame(lapply(psi.l,function(ii){c(ii,rep(NA,length.out=max(mi)-length(ii)))})))
  bb=cbind(rep(0,k),bb)

  pIRF <- function(theta, a, b){
    wts <- exp(b + a*theta)
    summe <- apply(wts, 1,sum,na.rm=TRUE)
    prob <- sweep(wts,1,summe, FUN="/")
    return(prob)
  }

  sim.data <- matrix(NA, nrow=n, ncol=k)

  for(j in (1:n)) {
    prob <- pIRF(faehig[j], aa, bb)
    cumprob <- t(apply(prob, 1, cumsum))
    jju <- runif(k)
    sim.data[j,] <-  apply(jju > cumprob,1,sum,na.rm=T)
  }
  colnames(sim.data)=paste("I",1:k,sep="")
  return(sim.data)
}


pscore_poly  <- function(i,x,r,coeff){
  pscore <- c()
  m <- sum(sapply(coeff, length))
  mi <- length(coeff[[i]])
  if (x > r | m-mi+x < r) {
    pscore <- 0
  } else {
    pscore <- ifelse(x > 0, exp(-coeff[[i]][x]), 1)*(psychotools::elementary_symmetric_functions(coeff[-i])[[1]])[r+1-x]/psychotools::elementary_symmetric_functions(coeff)[[1]][r+1]
  }
  pscore
}


#' Item Restscore Association
#'
#' The observed Gamma coefficient between the score of a single item and the total score of the remaining items
#' is compared with the corresponding expected Gamma coefficient under the Rasch model.
#' @param object An object of class "Rm", a fitted Rasch model or partial
#' credit model using  the functions RM or PCM in package eRm, or an object of class "pcmodel",
#'  a fitted partial credit model using the function pcmodel in package psychotools.
#' @param p.adj Correction method for multiple testing. The methods are "BH","holm", "hochberg", "hommel", "bonferroni", "BY", "none". See \code{\link{p.adjust}}.
#' @import eRm
#' @importFrom vcdExtra GKgamma
#' @importFrom stats p.adjust
#' @export
#' @return a matrix containing:
#' \item{observed}{observed gamma coefficients}
#' \item{expected}{expected gamma coefficients}
#' \item{se}{standard errors}
#' \item{pvalue}{p values (under normal distribution assumption)}
#' \item{sig}{significance stars: 0  " *** "  0.001  " ** "  0.01  " * "  0.05   " . "  0.1  " "  1}
#' @references Kreiner, S. (2011). A note on item-restscore association in Rasch models.
#' \emph{Applied Psychological Measurement}, 35, 557-561.
#' @author Marianne Mueller
#' @examples
#' rm.mod <- RM(amts[,4:13])
#' item_restscore(rm.mod)
item_restscore <- function(object, p.adj= c("BH","holm", "hochberg", "hommel", "bonferroni", "BY", "none")){
  if (!any("Rm"%in%class(object),class(object) =="pcmodel")) stop("object must be of class Rm or pcmodel!")
  if(class(object)[1]=="pcmodel") object$model <- "pcmodel"
  padj <- match.arg(p.adj)
  if (object$model=="RM") {
    if (sum(is.na(object$X)) > 0){
      message("Model was refitted with only complete cases")
      dat.items <- na.omit(object$X)
      object <- RM(dat.items)
    }
    X <- object$X
    k <- dim(X)[2]
    n <- dim(X)[1]
    coeff <- (-1)*coef(object)
    mi <- rep(1,k)
  } else {
    if (object$model=="PCM") {
      if (sum(is.na(object$X)) > 0){
        message("Model was refitted with only complete cases")
        dat.items <- na.omit(object$X)
        object <- PCM(dat.items)
      }
      X <- object$X
      k <- dim(X)[2]
      n <- dim(X)[1]
      mi <- apply(X, 2, max, na.rm = TRUE)
      thresh1 <- thresholds(object)[[3]][[1]][, -1] - mean(thresholds(object)[[3]][[1]][, 1])
      coeff <- lapply(as.list(as.data.frame(t(thresh1))), cumsum)
    } else {
      if (sum(is.na(object$data)) > 0){
        message("Model was refitted with only complete cases")
        dat.items <- na.omit(object$data)
        object <- pcmodel(dat.items)
      }
      X <- object$data
      k <- dim(X)[2]
      n <- dim(X)[1]
      mi <- apply(X, 2, max, na.rm = TRUE)
      thresh1 <- coef(threshpar(object),type="matrix")
      coeff <- lapply(as.list(as.data.frame(t(thresh1))), cumsum)
    }
  }
  m <- sum(mi)
  score <- apply(X,1,sum,na.rm=T)
  restscore <- matrix(rep(score,k),ncol=k)- X
  mm <- t(apply(rbind(X,restscore),2,function(x){
    gc <- GKgamma(table(x[1:n],x[(n+1):(2*n)]))
    gc <- c(gc[[1]],gc[[4]])
    gc
    }))
  #! beob gammas.se sind nicht gleich wie Svend's
  # erwartete gammas
  rvneu <- factor(score,levels = 0:m)
  npmat <- vector("list", k)
  for (i in 1:k) {
    nmat <- matrix(table(rvneu)[matrix(0:(m - mi[i]) + rep(0:mi[i], each = m - mi[i] + 1), ncol = mi[i] + 1) + 1], ncol = mi[i] + 1)
    pmat <- sapply(0:mi[i], function(x){sapply(0:(m - mi[i]),function(r){pscore_poly(i, x, r + x, coeff)})})
    npmat[[i]] <- nmat*pmat
  }
  expected <- sapply(npmat, function(x){GKgamma(x)[[1]]})
  pvalue <- ifelse((mm[, 1] - expected) > 0, 2*(1 - pnorm((mm[, 1] - expected)/mm[, 2])), 2*(pnorm((mm[, 1] - expected)/mm[, 2])))

  mm <- cbind(observed = mm[, 1], expected, se = mm[, 2], round(p.adjust(pvalue,method=padj, n=k),digits=4))
  symp <- symnum(pvalue, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols=c("*** ", "** ", "* ", "." , " "))
  mm <- noquote(cbind(format(mm, digits = 3), sig = symp))
  colnames(mm)[4] <- paste("pvalue",padj,sep=".")
  cat("\n")
  mm
}

