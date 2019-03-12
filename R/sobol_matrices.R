
# INTERNAL FUNCTION TO CREATE THE SCRAMBLED MATRIX FOR THE
# COMPUTATION OF FIRST-ORDER SOBOL' INDICES -----------------------------------

scrambled_sobol <- function(A, B) {
  X <- rbind(A, B)
  for(i in 1:ncol(A)) {
    AB <- A
    AB[, i] <- B[, i]
    X <- rbind(X, AB)
  }
  AB <- X
  return(AB)
}

# INTERNAL FUNCTION TO CREATE THE SCRAMBLED MATRIX FOR THE
# COMPUTATION OF SECOND-ORDER SOBOL' INDICES ----------------------------------

scrambled_second <- function(A, B) {
  X <- rbind(A, B)
  parms <- utils::combn(1:ncol(A), 2, simplify = FALSE)
  for(i in parms) {
    AB <- A
    AB[, c(i[[1]], i[[2]])] <- B[, c(i[[1]], i[[2]])]
    X <- rbind(X, AB)
  }
  n <- nrow(X)
  AB <- X[(2 * nrow(A) + 1):n, ]
  return(AB)
}

# INTERNAL FUNCTION TO CREATE THE SCRAMBLED MATRIX FOR THE
# COMPUTATION OF THIRD-ORDER SOBOL' INDICES ----------------------------------

scrambled_third <- function(A, B) {
  X <- rbind(A, B)
  parms <- utils::combn(1:ncol(A), 3, simplify = FALSE)
  for(i in parms) {
    AB <- A
    AB[, c(i[[1]], i[[2]], i[[3]])] <- B[, c(i[[1]], i[[2]], i[[3]])]
    X <- rbind(X, AB)
  }
  n <- nrow(X)
  AB <- X[(2 * nrow(A) + 1):n, ]
  return(AB)
}

# FUNCTION TO CREATE THE SOBOL' MATRIX FOR THE COMPUTATION OF
# FIRST AND TOTAL ORDER EFFECTS - AS WELL AS SECOND AND THIRD
# ORDER EFFECTS, IF DESIRED ---------------------------------------------------

#' Creation of the sample matrices
#'
#' It creates the sample matrices to compute Sobol' first and total-order indices.
#' If needed, it also creates the sample matrices required to compute second and
#' third-order indices. It uses Sobol' quasi-random number sequences.
#'
#' @param n Integer, sample size of the Sobol' matrix.
#' @param k Integer, number of model inputs.
#' @param second Logical. If \code{second = TRUE}, it creates the scrambled
#' matrix required to compute second-order indices. Default is \code{second = FALSE}.
#' @param third Logical. If \code{third = TRUE}, it creates the scrambled
#' matrix required to compute third-order indices. Default is \code{third = FALSE}.
#' @seealso Check the function \code{\link{sobol}} in the package \code{randtoolbox}
#' to see how the Sobol' quasi-random number sequences are constructed.
#'
#' @return A matrix.
#' @export
#' @details The function generates an \eqn{(n, 2k)} matrix using Sobol' quasi-random
#' number sequences. The first \emph{k}-matrix is the \strong{A} matrix and the
#' remaining \emph{k}-matrix, the \strong{B} matrix. It then generates \emph{k}
#' additional matrices {(\strong{A}^j_{\strong{B}})}, \eqn{j=1,2,...,k}, where
#' the \emph{k} matrix is composed of all columns of the \strong{A} matrix
#' except the \emph{j}-th column, which is the \emph{j} column of the
#' \strong{B} matrix. This approach leds to a total number of model runs of
#' \eqn{n(k + 2)} for first and total-order indices \insertCite{Saltelli2010a}{sensobol}.
#'
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{}
#' @examples
#' sobol_matrices(n = 100, k = 8, second = TRUE, third = TRUE)
sobol_matrices <- function(n, k, second = FALSE, third = FALSE) {
  # Create the Sobol quasi-random number sequence
  df <- randtoolbox::sobol(n = n,
                           dim = k * 2)
  # Create A matrix
  A <- df[, 1:k]
  # Create B matrix
  B <- df[, (k + 1) : (k * 2)]
  # Create AB matrix
  AB <- scrambled_sobol(A = A, B = B)
  if(second == TRUE & third == FALSE) {
    # Compute AB matrix for second order indices
    AB.2 <- scrambled_second(A = A, B = B)
    AB <- rbind(AB, AB.2)
  }
  if(second == FALSE & third == TRUE) {
    stop("The computation of third-order Sobol' indices requires the computation of second-order indices first")
  }
  if(second == TRUE & third == TRUE) {
    AB.2 <- scrambled_second(A = A, B = B)
    AB.3 <- scrambled_third(A = A, B = B)
    AB <- rbind(AB, AB.2, AB.3)
  }
  return(AB)
}
