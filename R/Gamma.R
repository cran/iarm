#' Conditional and Partial Gamma Coefficients
#'
#' Calculates conditional and partial Gamma coefficients for x and y given z with confidence intervals.
#' @param x,y,z  Three numeric vectors or factors.
#' @param conf.level Confidence level for the returned confidence interval.
#' @return data frame with estimates, standard errors and confidence interval limits.
#' @importFrom stats xtabs qnorm
#' @seealso  \code{\link{partgam_DIF}}, \code{\link{partgam_LD}}
#' @export
#' @author Marianne Mueller
#' @references Davis, J. A. A Partial coefficient for Goodman and Kruskal's Gamma.
#'  \emph{Journal of the American Statistical Association}, 62 (317), 1967, pp. 189-193.
partgam <- function(x, y, z, conf.level = 0.95){
  xxx <- xtabs(~ x + y + z)
  n <- dim(xxx)[1]
  m <- dim(xxx)[2]
  k <- dim(xxx)[3]
  sigma2.k <- rep(NA, k)
  gamma.k  <- rep(NA, k)
  C.k <- rep(NA, k)
  D.k <- rep(NA, k)
  pi.c.k <- pi.d.k <- array(0, c(n, m, k))
  for (l in 1:k){
    xx <- xxx[, , l]
    row.x <- row(xx)
    col.x <- col(xx)
    for(i in 1:n){
       for(j in 1:m){
         pi.c.k[i, j, l]<-sum(xx[row.x < i & col.x < j]) + sum(xx[row.x > i & col.x > j])
         pi.d.k[i, j, l]<-sum(xx[row.x < i & col.x > j]) + sum(xx[row.x > i & col.x < j])
       }
    }
    C.k[l] <- sum(pi.c.k[, , l]*xx)/2
    D.k[l] <- sum(pi.d.k[, , l]*xx)/2
    psi.k <- 2*(D.k[l]*pi.c.k[, , l]-C.k[l]*pi.d.k[, , l])/(C.k[l]+D.k[l])^2
    sigma2.k[l] <- sum(xx*psi.k^2)
    gamma.k[l] <- (C.k[l] - D.k[l])/(C.k[l] + D.k[l])
  }
  C <- sum(C.k)
  D <- sum(D.k)
  psi <- 2*(D*pi.c.k[, , ]-C*pi.d.k[, , ])/(C + D)^2
  sigma2 <- c(sigma2.k, sum(xxx*psi^2))
  gamma <- c(gamma.k, (C - D)/(C + D))
  pr2 <- 1 - (1 - conf.level)/2
  CIa <- outer(sqrt(sigma2),qnorm(pr2)*c(-1, 1)) + gamma
  mm <- data.frame(gamma, se = sqrt(sigma2), CI1 = ifelse(CIa[, 1] > -1 ,CIa[, 1], -1),
  CI2 = ifelse(CIa[, 2] < 1, CIa[, 2], 1))
  row.names(mm)= c(paste("conditional", 1:k), "partial")
  mm
}

#' Partial Gamma to detect Differential Item Functioning (DIF)
#'
#' Items should function in the same way for all subgroups of persons. An item
#' shows differential item functioning (DIF) if there is a significant association between the item score and
#' an exogenous variable, controlling for the scale score. Partial Gamma coefficients are used as test statistics.
#' @param dat.items A data frame with the responses to the items.
#' @param dat.exo  A single grouping factor or a data frame consisting of several exogenous factor variables.
#' @return data frame with Gamma coefficients, standard errors, p values and confidence limits for every
#' pair of an item and an exogenous variable.
#' @importFrom stats quantile
#' @seealso  {\code{\link{partgam_LD}}}
#' @export
#' @author Marianne Mueller
#' @examples
#' partgam_DIF(amts[,4:13],amts[,2:3])
#' @references Bjorner, J., Kreiner, S., Ware, J., Damsgaard, M. and Bech, P. Differential item
#' functioning in the Danish translation of the SF-36. \emph{Journal of Clinical Epidemiology},
#' 51 (11), 1998, pp. 1189-1202.
partgam_DIF <- function(dat.items,dat.exo){
  if (!is.data.frame(dat.exo)) {
    gname <- deparse(substitute(dat.exo))
    dat.exo <- data.frame(dat.exo)
    names(dat.exo) <- gname
  }
  if (is.null(names(dat.items))) names(dat.items) <- paste("I",1:dim(dat.items)[2],sep="")
  score <- apply(dat.items,1,sum,na.rm=T)
  fz <- cut(score,unique(quantile(score,0:10/10)),include.lowest=T)
  k <- dim(dat.items)[2]
  l <- dim(dat.exo)[2]
  result <- data.frame(Item=character(), Var=character(), gamma=double(),se=double(),pvalue=double(),sig=character(),lower=double(),upper=double(),stringsAsFactors=FALSE)
  z <- 1
  for (i in 1:k){
    for (j in 1:l){
      mm <- partgam(dat.items[, i], dat.exo[, j], fz)
      pvalue <- round(ifelse(mm[dim(mm)[1],1] > 0, 2*(1 - pnorm(mm[dim(mm)[1],1]/mm[dim(mm)[1],2])), 2*(pnorm(mm[dim(mm)[1],1]/mm[dim(mm)[1],2]))), digits = 4)
      symp <- symnum(pvalue, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c(" ***", " **", " *", " .", " "))
      result[z,] <- c(names(dat.items)[i],names(dat.exo)[j],mm[dim(mm)[1],1:2],pvalue,symp,mm[dim(mm)[1],3:4])
      z <- z + 1
    }
  }
  result
}

#' Partial Gamma to detect Local Dependence (LD)
#'
#' Rasch models assume locally independent items. There should be no substantial correlation left between two items
#' once the underlying factor has been taken into account.
#' Partial Gamma coefficients between pairs of items controlled for the rest score
#' can be used to assess this requirement.
#' @param dat.items A data frame with the responses to the items.
#' @return data frame with Gamma coefficients, standard errors, p values and confidence limits for every
#' pair of items.
#' @importFrom stats quantile
#' @export
#' @author Marianne Mueller
#' @seealso  \code{\link{partgam_DIF}}
#' @examples
#' partgam_LD(amts[,4:13])
#' @references Christensen, K. B. , Kreiner, S. & Mesbah, M. (Eds.)
#' \emph{Rasch Models in Health}. Iste and Wiley (2013), pp. 133 - 135.
partgam_LD <- function(dat.items){
  score <- apply(dat.items,1,sum,na.rm=T)
  ok <- complete.cases(dat.items)
  k <- dim(dat.items)[2]
  result <- data.frame(Item1=character(),Item2=character(), gamma=double(),se=double(),pvalue=double(),sig=character(),lower=double(),upper=double(),stringsAsFactors=FALSE)
  z <- 1
  for (i in 1:(k-1)){
    for (j in (i+1):k){
      rest <- score[ok] - dat.items[ok,j]
      fac.z <- cut(rest,unique(quantile(rest,0:10/10)),include.lowest=T)
      mm <- partgam(dat.items[ok,i], dat.items[ok,j],fac.z)
      pvalue <- round(ifelse(mm[dim(mm)[1],1] > 0, 2*(1 - pnorm(mm[dim(mm)[1],1]/mm[dim(mm)[1],2])), 2*(pnorm(mm[dim(mm)[1],1]/mm[dim(mm)[1],2]))), digits = 4)
      symp <- symnum(pvalue, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c(" ***", " **", " *", " .", " "))
      result[z,] <- c(names(dat.items)[i],names(dat.items)[j],mm[dim(mm)[1],1:2],pvalue,symp,mm[dim(mm)[1],3:4])
      z <- z + 1
    }
  }
  result
}



