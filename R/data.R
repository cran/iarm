#' Abbreviated Mental Test Score (AMTS)
#'
#' A dataset containing the responses of 197 persons to the ten questions of the Abbreviated Mental Test Score (AMTS).
#' The AMTS is used to identify patients with dementia.
#' One point is given for each correct answer,
#' a score of 6 or less suggests that the patient has some mental impairment.
#'
#' @format A data frame with 197 rows and 13 variables.
#' \describe{
#'   \item{id}{id number of the patient.}
#'   \item{agegrp}{a factor with levels 16-65, 66-75, 76-85, 86+ for the age  of the patient.}
#'   \item{sex}{a factor with levels male, female of the patient.}
#'   \item{age}{age of patient, with 1 if the respondent knows his/her own age and 0 otherwise.}
#'   \item{time}{time (nearest hour), with 1 if correct and 0 otherwise.}
#'   \item{address}{address, with 1 if correct and 0 otherwise.}
#'   \item{name}{name of hospital (or area of town if at home) , with 1 if correct and 0 otherwise.}
#'   \item{year}{current year, with 1 if correct and 0 otherwise.}
#'   \item{dob}{date of birth of patient, with 1 if correct and 0 otherwise.}
#'    \item{month}{month, with 1 if correct and 0 otherwise.}
#'    \item{firstww}{date of first world war, with 1 if correct and 0 otherwise.}
#'     \item{monarch}{name of monarch, with 1 if correct and 0 otherwise.}
#'     \item{countbac}{count backwards 20-1, with 1 if correct and 0 otherwise.}
#' }
#' @docType data
#' @keywords datasets
#' @references Slade, A., Fear, J. & Tennant, A.  (2006) Identifying patients at risk of nursing home admission: The Leeds
#' Elderly Assessment Dependency Screening tool (LEADS). \emph{BMC Health Services Research}, 6:31.
#' @name amts
#' @examples data(amts)
#' @examples str(amts)
NULL

#' Depression Screening DESC-II
#'
#' A dataset containing the responses of 799 patients (indication group psychiatry, otolaryngology, cardiology, neurology)
#' to the short form DESC-II with 10 items.
#' There are 5 response categories from 0 = never to 4 = always. A higher score is supposed to mean a higher depression.

#' @format A data frame with 799 rows and 14 variables.
#' \describe{
#'   \item{code}{id number of the patient}
#'   \item{group}{a factor with levels  psychiatry, otolaryngology, cardiology, neurology  for the indication group  of the patient.}
#'   \item{gender}{a factor with levels female, male of the patient.}
#'   \item{agegroup}{a factor with levels 18-34, 35-49, 50-59, 60-87 for the age  of the patient.}
#'   \item{DESC_2_1}{feeling not to be needed}
#'   \item{DESC_2_2}{loss of interest in other people}
#'   \item{DESC_2_3}{disheartened}
#'   \item{DESC_2_4}{no pleasure doing things}
#'   \item{DESC_2_5}{feeling to be no good}
#'    \item{DESC_2_6}{uninspired}
#'    \item{DESC_2_7}{pessimistic}
#'     \item{DESC_2_8}{discouraged}
#'     \item{DESC_2_9}{withdrawal}
#'     \item{DESC_2_10}{thinking of taking one's life}
#' }
#' @docType data
#' @keywords datasets
#' @references Forkmann et al. (2009) Development and validation of the Rasch-based Depression Screening (DESC) using
#' Rasch analysis and structural equation modelling.
#' \emph{J Behav Ther Exp Psychiatry}, 40(3): 468-78.
#' @seealso \url{http://psychometrikon.de/inhalt/suchen/test.php?id=000000000000000000000000000003}
#' @name desc2
#' @examples data(desc2)
#' @examples str(desc2)
NULL
