#' Conditional and Partial Gamma Coefficients
#'
#' Calculates conditional and partial Gamma coefficients for x and y given z with confidence intervals.
#' @param x,y,z  Three numeric vectors or factors.
#' @param conf.level Confidence level for the returned confidence interval.
#' @return matrix with estimates, standard errors and confidence interval limits.
#' @importFrom stats xtabs qnorm
#' @export
#' @author Marianne Mueller
#' @examples # Partial Gamma coefficient between an item and an exogenuous variable, given the total score
#' score <- apply(amts[, 4:13], 1, sum, na.rm = TRUE)
#' fz <- cut(score,unique(quantile(score,0:10/10)))
#' partgam(amts$firstww,amts$sex,fz)
#' @references Davis, J. A. A Partial coefficient for Goodman and Kruskal's Gamma.
#'  \emph{Journal of the American Statistical Association}, 62 (317), 1967, pp. 189-193.
partgam <- function(x, y, z, conf.level = 0.95){
  xxx <- xtabs(~ x + y + fz)
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




