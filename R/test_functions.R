
# Sobol ' G function -----------------------------------------------------------

#' Sobol' G function
#'
#' It implements the \insertCite{Sobol1998;textual}{sensobol} function. In this
#' case, the function works with 8 model inputs.
#'
#' @param X A data frame or numeric matrix.
#'
#' @return A numeric vector with the model output.
#' @export
#'
#' @references
#' \insertAllCited{}
#' @examples
#' A <- sobol_matrices(n = 100, k = 8)
#' Y <- sobol_Fun(A)
sobol_Fun <- function(X) {
  a <- c(0, 1, 4.5, 9, 99, 99, 99, 99)
  y <- 1
  for (j in 1:8) {
    y <- y * (abs(4 * X[, j] - 2) + a[j])/(1 + a[j])
  }
  return(y)
}

# Ishigami function -----------------------------------------------------------

ishigami <- function(X1, X2, X3) {
  A <- 2
  B <- 1
  sin(X1) + A * sin(X2) ^ 2 + B * X3 ^ 4 * sin(X1)
}

#' Ishigami function
#'
#' It implements the \insertCite{Ishigami1990;textual}{sensobol} function,
#' which requires 3 model inputs. The transformation of the
#' distribution of the model inputs (from U(0, 1) to U(-pi, +pi)) is conducted
#' internally.
#'
#' @param X A data frame or numeric matrix.
#'
#' @return A numeric vector with the model output.
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' A <- sobol_matrices(n = 100, k = 3)
#' Y <- ishigami_Mapply(A)
ishigami_Mapply <- function(X) {
  X <- apply(X, 2, function(x) x * (pi + pi) - pi)
  return(mapply(ishigami,
                X[, 1],
                X[, 2],
                X[, 3]))
}

# Bratley et al. function -----------------------------------------------------


#' Bratley, Fox and Niederreiter (1992) function
#'
#' It implements the \insertCite{Bratley1992;textual}{sensobol} function,
#' #' which requires \emph{n} model inputs.
#'
#' @param X A data frame or numeric matrix.
#'
#' @return A numeric vector with the model output.
#' @export
#' @references
#' \insertAllCited{}
#'
#' @examples
#' A <- sobol_matrices(n = 100, k = 4)
#' Y <- bratley1992_Fun(A)
bratley1992_Fun <- function(X) {
  # Preallocate
  mat <- tmp <- vector(mode = "list", length = nrow(X))
  Y <- vector(mode = "numeric", length = nrow(X))
  for(i in 1:nrow(X)) {
    mat[[i]] <- matrix(rep(X[i, ], times = ncol(X)),
                       nrow = ncol(X),
                       ncol = ncol(X),
                       byrow = TRUE)
    mat[[i]][upper.tri(mat[[i]])] <- 1
    tmp[[i]] <- matrixStats::rowProds(mat[[i]])
    Y[[i]] <- sum(tmp[[i]] * (-1) ^ (1:ncol(X)))
  }
  return(Y)
}

# Bratley and Fox (1988) function ---------------------------------------------

#' Bratley and Fox (1988) function
#'
#' It implements the \insertCite{Bratley1988;textual}{sensobol} function,
#' which requires \emph{n} model inputs.
#'
#' @param X A data frame or numeric matrix.
#'
#' @return A numeric vector with the model output.
#'
#' @export
#'
#' @examples
#' A <- sobol_matrices(n = 100, k = 6)
#' Y <- bratley1988_Fun(A)
bratley1988_Fun <- function(X) {
  y <- 1
  for (j in 1:ncol(X)) {
    y <- y * (abs(4 * X[, j] - 2))}
  return(y)
}

# Oakley and O'Hagan function -------------------------------------------------

#' Oakley & O'Hagan (2004) function
#'
#' It implements the \insertCite{Oakley2004;textual}{sensobol} function.
#' The function needs 15 model inputs. The transformation of the
#' distribution of the model inputs (from uniform to normally distributed
#' with mean = 0 and sd = 1) is conducted internally.
#'
#' @param X A data frame or numeric matrix.
#'
#' @return A numeric vector with the model output.
#'
#' @export
#' @references
#' \insertAllCited{}
#'
#' @examples
#' A <- sobol_matrices(n = 100, k = 15)
#' Y <- oakley_Fun(A)
oakley_Fun <- function(X) {

  a1 <- c(0.0118, 0.0456, 0.2297, 0.0393, 0.1177,
          0.3865, 0.3897, 0.6061, 0.6159, 0.4005,
          1.0741, 1.1474, 0.788, 1.1242, 1.1982)

  a2 <- c(0.4341, 0.0887, 0.0512, 0.3233, 0.1489,
          1.036, 0.9892, 0.9672, 0.8977, 0.8083,
          1.8426, 2.4712, 2.3946, 2.0045, 2.2621)

  a3 <- c(0.1044, 0.2057, 0.0774, 0.273, 0.1253,
          0.7526, 0.857, 1.0331, 0.8388, 0.797,
          2.2145, 2.0382, 2.4004, 2.0541, 1.9845)

  M <- c(-0.022482886, -0.18501666, 0.13418263, 0.36867264, 0.17172785, 0.13651143, -0.44034404,
         -0.081422854, 0.71321025, -0.44361072, 0.50383394, -0.024101458, -0.045939684, 0.21666181,
         0.055887417, 0.2565963, 0.053792287, 0.25800381, 0.23795905, -0.59125756, -0.081627077,
         -0.28749073, 0.41581639, 0.49752241, 0.083893165, -0.11056683, 0.033222351, -0.13979497,
         -0.031020556, -0.22318721, -0.055999811, 0.19542252, 0.095529005, -0.2862653, -0.14441303,
         0.22369356, 0.14527412, 0.28998481, 0.2310501, -0.31929879, -0.29039128, -0.20956898, 0.43139047,
         0.024429152, 0.044904409, 0.66448103, 0.43069872, 0.29924645, -0.16202441, -0.31479544,
         -0.39026802, 0.17679822, 0.057952663, 0.17230342, 0.13466011, -0.3527524, 0.25146896, -0.018810529,
         0.36482392, -0.32504618, -0.121278, 0.12463327, 0.10656519, 0.046562296, -0.21678617, 0.19492172,
         -0.065521126, 0.024404669, -0.09682886, 0.19366196, 0.33354757, 0.31295994, -0.083615456,
         -0.25342082, 0.37325717, -0.2837623, -0.32820154, -0.10496068, -0.22073452, -0.13708154,
         -0.14426375, -0.11503319, 0.22424151, -0.030395022, -0.51505615, 0.017254978, 0.038957118,
         0.36069184, 0.30902452, 0.050030193, -0.077875893, 0.003745656, 0.88685604, -0.26590028,
         -0.079325357, -0.042734919, -0.18653782, -0.35604718, -0.17497421, 0.088699956, 0.40025886,
         -0.055979693, 0.13724479, 0.21485613, -0.011265799, -0.09229473, 0.59209563, 0.031338285,
         -0.033080861, -0.24308858, -0.099798547, 0.034460195, 0.095119813, -0.3380162, 0.0063860024,
         -0.61207299, 0.081325416, 0.88683114, 0.14254905, 0.14776204, -0.13189434, 0.52878496, 0.12652391,
         0.045113625, 0.58373514, 0.37291503, 0.11395325, -0.29479222, -0.57014085, 0.46291592, -0.094050179,
         0.13959097, -0.38607402, -0.4489706, -0.14602419, 0.058107658, -0.32289338, 0.093139162,
         0.072427234, -0.56919401, 0.52554237, 0.23656926, -0.011782016, 0.071820601, 0.078277291,
         -0.13355752, 0.22722721, 0.14369455, -0.45198935, -0.55574794, 0.66145875, 0.34633299, 0.14098019,
         0.51882591, -0.28019898, -0.1603226, -0.068413337, -0.20428242, 0.069672173, 0.23112577,
         -0.044368579, -0.16455425, 0.21620977, 0.0042702105, -0.087399014, 0.31599556, -0.027551859,
         0.13434254, 0.13497371, 0.05400568, -0.17374789, 0.17525393, 0.060258929, -0.17914162, -0.31056619,
         -0.25358691, 0.025847535, -0.43006001, -0.62266361, -0.033996882, -0.29038151, 0.03410127,
         0.034903413, -0.12121764, 0.026030714, -0.33546274, -0.41424111, 0.05324838, -0.27099455,
         -0.026251302, 0.41024137, 0.26636349, 0.15582891, -0.18666254, 0.019895831, -0.24388652,
         -0.44098852, 0.012618825, 0.24945112, 0.071101888, 0.24623792, 0.17484502, 0.0085286769,
         0.2514707, -0.14659862, -0.08462515, 0.36931333, -0.29955293, 0.1104436, -0.75690139, 0.041494323,
         -0.25980564, 0.46402128, -0.36112127, -0.94980789, -0.16504063, 0.0030943325, 0.052792942,
         0.22523648, 0.38390366, 0.45562427, -0.18631744, 0.0082333995, 0.16670803, 0.16045688)

  M <- matrix(M, 15, 15, byrow = TRUE)

  Y <- vector()
  # transformation to normal distribution
  X <- apply(X, 2, function(x) stats::qnorm(x))
  for(i in 1:nrow(X)) {
    mat <- matrix(X[i, ])
    Y[i] <- a1 %*% mat + a2 %*% sin(mat) + a3 %*% cos(mat) + t(mat) %*% M %*% mat
  }
  return(Y)
}
