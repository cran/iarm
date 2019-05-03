#' Person Estimates with MLE and WLE
#'
#' Computes Person estimates with maximum likelihood estimation (MLE) and  weighted likelihood estimation (WLE) for raw scores 0 to m.
#' @param  object An object of class "Rm", a fitted Rasch model or partial
#' credit model using  the functions RM or PCM in package eRm, or an object of class "pcmodel",
#'  a fitted partial credit model using the function pcmodel in package psychotools.
#' @param  properties If TRUE additional properties of the estimates are given (see below).
#' @param allperson If TRUE person estimates (MLE and WLE) for all persons in the data set are delivered.
#' @import eRm
#' @importFrom PP PP_gpcm
#' @importFrom psychotools personpar
#' @importFrom stats coef uniroot na.omit
#' @export
#' @return If properties = False a matrix containing:
#' \item{Raw score}{raw score}
#' \item{MLE}{MLE of person parameters}
#' \item{WLE}{WLE of person parameters}
#'
#' If properties = TRUE a list with two components, one for MLE and the other for WLE. Each component
#' contains:
#' \item{Raw score}{raw score}
#' \item{MLE or WLE}{person estimates}
#' \item{SEM}{standard error of measurement}
#' \item{Bias}{bias}
#' \item{RMSE}{root mean square error}
#' \item{Score.SEM}{score sem}
#' @references Christensen, K. B. , Kreiner, S. & Mesbah, M. (Eds.)
#' \emph{Rasch Models in Health}. Iste and Wiley (2013), pp. 63 - 70.
#' @author Marianne Mueller
#' @examples
#' rm.mod <- RM(amts[,4:13])
#' person_estimates(rm.mod)
person_estimates <- function(object, properties = F, allperson = F){
    if (!any("Rm"%in%class(object),class(object) =="pcmodel")) stop("object must be of class Rm or pcmodel!")
    if(class(object)[1]=="pcmodel") object$model <- "pcmodel"
    if (object$model == "pcmodel") {X <- object$data
    } else {X <- object$X
    }
    if (object$model == "RM") {
      k <- dim(X)[2]
      coeff <- (-1)*coef(object)
      m <- k
      respm <- rbind(rep(0, k), lower.tri(matrix(1, k, k)) + diag(k))
    } else {
      if (object$model == "PCM"){
        coeff <- thresholds(object)[[3]][[1]][, -1]- mean(thresholds(object)[[3]][[1]][, 1])
      } else {
        coeff <- coef(threshpar(object),type="matrix")
      }
      k <- dim(X)[2]
      mi <- apply(X, 2, max, na.rm = TRUE)
      m <- sum(mi)
      respm <- matrix(0, ncol = k, nrow = m + 1)
      respm[, 1] <- c(0:mi[1], rep(mi[1], nrow(respm) - mi[1] - 1))
      for (i in 2:k) respm[, i] <- c(rep(0, cumsum(mi)[i - 1] + 1), 1:mi[i], rep(mi[i], nrow(respm) - cumsum(mi)[i]  -1))
    }
    invisible(capture.output(mm <- cbind(0:m, PP_gpcm(respm, t(coeff), slopes = rep(1, k), type = "mle" )[[1]][[1]][, 1],
                PP_gpcm(respm,t(coeff),slopes=rep(1,k),type="wle")[[1]][[1]][,1])))
    rownames(mm) <- rep(" ", m + 1)
    colnames(mm) <- c("Raw Score", "MLE", "WLE")
    if (allperson){
      properties <- F
      rv <- rowSums(X, na.rm = TRUE)
      mm <- mm[rv+1,]
      mm
    } else {
      if (properties == F) {
        mm
      } else {
        koeff <- lapply(as.list(as.data.frame(t(coeff))), function(x) cumsum(na.omit(x)))
        gr <- elementary_symmetric_functions(koeff)[[1]]
        s.theta <- function(r){
          function(x){
            ((exp(x*(0:m))*gr)/as.vector(exp(x*(0:m))%*%gr))%*%(0:m) - r
          }
        }
        if (object$model == "pcmodel") mm[1, 2] <- NA else  mm[1, 2] <- person.parameter(object)$pred.list[[1]]$y[1]
        try(mm[1, 2] <- uniroot(s.theta(0.25), c(-10, 10))$root)
        mm[m + 1, 2] <- uniroot(s.theta(m - 0.25), c(-6, 6))$root
        rvec = 0:m
        pers_prop <- function(x, persons){
          pr <- (exp(x[2]*rvec)*gr)/as.vector(exp(x[2]*rvec)%*%gr)
          bias <- pr%*%persons - x[2]
          sem <- sqrt((persons - as.vector(pr%*%persons))^2%*%pr)
          rsem <- sqrt((persons - x[2])^2%*%pr)
          scoresem <- sqrt((rvec- x[1])^2%*%pr)
          c(SEM = sem, Bias = bias, RMSE = rsem, Score.SEM = scoresem)
        }
        result <- list(cbind(mm[, 1:2],t(apply(mm[, c(1, 2)], 1, pers_prop, persons = mm[, 2]))),
                     cbind(mm[, c(1,3)], t(apply(mm[, c(1, 3)], 1, pers_prop, persons = mm[, 3]))))
      }
  }
}

#' Properties of the Test
#'
#' Information summarizing measurement quality of the test and test targeting.
#' @param  object An object of class "Rm", a fitted Rasch model or partial
#' credit model using  the functions RM or PCM in package eRm, or an object of class "pcmodel",
#'  a fitted partial credit model using the function pcmodel in package psychotools.
#' @import eRm
#' @importFrom psychotools elementary_symmetric_functions personpar
#' @importFrom stats na.omit
#' @export
#' @return a list containing:
#' \item{Separation reliability}{the person separation reliability as calculated in package eRm for objects of class "Rm".}
#' \item{Test difficulty}{person value with an expected score equal to half of the maximum score.}
#' \item{Test target}{person value where test information is maximized.}
#' \item{Test information}{maximal value of the test information}
#' @references Christensen, K. B. , Kreiner, S. & Mesbah, M. (Eds.)
#' \emph{Rasch Models in Health}. Iste and Wiley (2013), pp. 63 - 70.
#' @author Marianne Mueller
#' @examples
#' rm.mod <- RM(amts[,4:13])
#' test_prop(rm.mod)
test_prop <- function(object){
  if (!any("Rm"%in%class(object),class(object) =="pcmodel")) stop("object must be of class Rm or pcmodel!")
  if(class(object)[1]=="pcmodel") object$model <- "pcmodel"
  if (object$model == "RM") {
    k <- dim(object$X)[2]
    koeff <- (-1)*coef(object)
    mi <- rep(1, k)
  } else {
    if (object$model == "PCM"){
      k <- dim(object$X)[2]
      mi <- apply(object$X, 2, max, na.rm = TRUE)
      thresh1 <- thresholds(object)[[3]][[1]][, -1] - mean(thresholds(object)[[3]][[1]][, 1])
      koeff <- lapply(as.list(as.data.frame(t(thresh1))), function(x) cumsum(na.omit(x)))
    } else {
      k <- dim(object$data)[2]
      mi <- apply(object$data, 2, max, na.rm = TRUE)
      thresh1 <- coef(threshpar(object),type="matrix")
      koeff <- lapply(as.list(as.data.frame(t(thresh1))), function(x) cumsum(na.omit(x)))
    }
  }
  m <- sum(mi)
  gr <- elementary_symmetric_functions(koeff)[[1]]
  var.R <- function(x) {
    rvec <- 0:m
    pr <- (exp(x*rvec)*gr)/as.vector(exp(x*rvec)%*%gr)
    (rvec - as.vector(pr%*%rvec))^2%*%pr
  }
  s.theta <- function(r){
    function(x){
      ((exp(x*(0:m))*gr)/as.vector(exp(x*(0:m))%*%gr))%*%(0:m) - r
    }
  }
  diffic <- round(uniroot(s.theta(m/2),c(-5,5))$root, digits = 3)
  target <- round(optimize(var.R,c(-4,4),maximum=T)$maximum, digits = 3)
  info <- max(person_estimates(object, properties=TRUE)[[2]][,6]^2)
  if (object$model == "pcmodel") {
    result <- list(diffic, target, info)
    names(result) = c("Test difficulty", "Test target", "Test information")
  } else {
    result <- list(SepRel(person.parameter(object))[[1]], diffic, target, info)
    names(result) = c("Separation Reliability", "Test difficulty", "Test target", "Test information")
  }
  print(unlist(result))
}

