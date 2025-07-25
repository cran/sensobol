
# Sobol' indices for a dummy parameter
sobol_dummy_boot <- function(d, i, N, params, boot) {
  k <- length(params)
  if (boot == TRUE) {
    m <- d[i, ]
  } else if (boot == FALSE) {
    m <- d
  }

  # Computation of E(Y), V(Y), Si and Ti
  # ----------------------------------------------------------------

  f0 <- (1 / N) * sum(m[, 1] * m[, 2])
  VY <- 1 / (2 * N - 1) * sum(m[, 1]^2 + m[, 2]^2) - f0
  Si <- (1 / (N - 1) * sum(m[, 1] * m[, 2]) - f0) / VY
  STi <- 1 - (1 / (N - 1) * sum(m[, 2] * m[, 2]) - f0) / VY
  return(c(Si, STi))
}

#' Computation of Sobol' indices for a dummy parameter
#'
#' This function computes first and total-order Sobol' indices for a dummy
#' parameter following the formulae shown
#' in \insertCite{KhorashadiZadeh2017;textual}{sensobol}.
#'
#' @param Y A numeric vector with the model output obtained from the matrix created with
#' \code{\link{sobol_matrices}}. The numeric vector should not contain any NA or NaN values.
#' @param N Positive integer, the initial sample size of the base sample matrix created with \code{\link{sobol_matrices}}.
#' @param params A character vector with the name of the model inputs.
#' @param boot Logical. If TRUE, the function bootstraps the Sobol' indices. If FALSE, it provides point
#' estimates. Default is \code{boot = FALSE}.
#' @param R Positive integer, number of bootstrap replicas.
#' @param parallel The type of parallel operation to be used (if any).
#' If missing, the default is taken from the option "boot.parallel"
#' (and if that is not set, "no"). For more information, check the
#' \code{parallel} option in the \code{boot} function of the \code{\link[boot]{boot}} package.
#' @param ncpus Positive integer: number of processes to be used in parallel operation:
#' typically one would chose this to the number of available CPUs.
#' Check the \code{ncpus} option in the \code{boot} function of the \code{\link[boot]{boot}} package.
#' @param conf Confidence intervals, number between 0 and 1. Default is \code{conf = 0.95}.
#' @param type Method to compute the confidence intervals. Default is \code{type = "norm"}.
#' Check the \code{type} option in the \code{boot} function of the \code{\link[boot]{boot}} package.
#'
#' @importFrom Rdpack reprompt
#' @importFrom rlang :=
#' @references
#' \insertAllCited{}
#'
#' @return A \code{data.table} object.
#' @export
#'
#' @examples
#' # Define settings
#' N <- 100; params <- paste("X", 1:3, sep = ""); R <- 10
#'
#' # Create sample matrix
#' mat <- sobol_matrices(N = N, params = params)
#'
#' # Compute Ishigami function
#' Y <- ishigami_Fun(mat)
#'
#' # Compute and bootstrap Sobol' indices for dummy parameter
#' ind.dummy <- sobol_dummy(Y = Y, N = N, params = params, boot = TRUE, R = R)

sobol_dummy <- function(Y, N, params, boot = FALSE, R = NULL, parallel = "no",
                        ncpus = 1, conf = 0.95, type = "norm") {
  k <- length(params)
  d <- matrix(Y, nrow = N)

  # Compute with bootstrap
  # -----------------------------------------------------------------

  if (boot == TRUE) {
    out <- boot::boot(data = d, statistic = sobol_dummy_boot,
                      R = R, N = N, params = params,
                      parallel = parallel, ncpus = ncpus,
                      boot = boot)
    out <- bootstats(out, conf = conf, type = type)

  # Compute without bootstrap
  # -----------------------------------------------------------------

  } else if (boot == FALSE) {
    tmp <- sobol_dummy_boot(d = d, N = N, params = params, boot = FALSE)
    out <- data.table::data.table(tmp)
    data.table::setnames(out, "tmp", "original")
  }
  sensitivity <- c("Si", "Ti")
  parameters <- "dummy"
  out <- cbind(out, sensitivity, parameters)

  # Transform negative values to zeroes
  # -----------------------------------------------------------------

  colNames <- colnames(out)

  if (any(grepl("high.ci", colNames)) == TRUE) {
    cols_transform <- c("original", "low.ci", "high.ci")

  } else {
    cols_transform <- "original"
  }

  for(j in cols_transform){
    data.table::set(out, i = which(out[[j]] < 0), j = j, value = 0)
  }

  return(out)
}


sobol_boot <- function(d, i, N, params, matrices, R, first, total, order, boot) {

  # Stopping rule to check concordance between estimators and sample matrix
  # -------------------------------------------------------------------

  ms <- "Revise the correspondence between the matrices and the estimators"

  if (isTRUE(all.equal(matrices, c("A", "B", "AB")))) {
    if (!first == "saltelli" & !first == "jansen" |
        !total == "jansen" & !total == "sobol" & !total == "homma" &
        !total == "janon" & !total == "glen") {
      stop(ms)
    }

  } else if (isTRUE(all.equal(matrices, c("A", "B", "BA")))) {
    if (!first == "sobol"| !total == "saltelli") {
      stop(ms)
    }

  } else if (isTRUE(all.equal(matrices, c("A", "B", "AB", "BA")))) {

    if (!first == "azzini" | !total == "azzini" &
        !total == "jansen" & !total == "sobol" & !total == "homma" &
        !total == "janon" & !total == "glen" & !total == "saltelli") {

      if (!total == "azzini" | !first == "saltelli" & !first == "jansen" &
          !first == "azzini" & !first == "sobol") {

        stop(ms)
      }
    }
  }

  # -------------------------------------

  k <- length(params)

  if (boot == TRUE) {
    m <- d[i, ]

  } else if (boot == FALSE) {
    m <- d
  }
  if (order == "second") {
    k <- length(params) + length(utils::combn(params, 2, simplify = FALSE))

  } else if (order == "third") {
    k <- length(params) +
      length(utils::combn(params, 2, simplify = FALSE)) +
      length(utils::combn(params, 3, simplify = FALSE))

  } else if (order == "fourth") {
    k <- length(params) +
      length(utils::combn(params, 2, simplify = FALSE)) +
      length(utils::combn(params, 3, simplify = FALSE)) +
      length(utils::combn(params, 4, simplify = FALSE))
  }

  # Define vectors based on sample design
  # ------------------------------------------------------------------

  if (isTRUE(all.equal(matrices, c("A", "B", "AB")))) {
    Y_A <- m[, 1]
    Y_B <- m[, 2]
    Y_AB <- m[, -c(1, 2)]

  } else if (isTRUE(all.equal(matrices, c("A", "B", "BA")))) {
    Y_A <- m[, 1]
    Y_B <- m[, 2]
    Y_BA <- m[, -c(1, 2)]

  } else if (isTRUE(all.equal(matrices, c("A", "B", "AB", "BA")))) {
    Y_A <- m[, 1]
    Y_B <- m[, 2]
    Y_AB <- m[, 3:(k + 2)]
    Y_BA <- m[, (k + 3):ncol(m)]

  } # A warning might be needed here
  if (isTRUE(all.equal(matrices, c("A", "B", "AB"))) |
      isTRUE(all.equal(matrices, c("A", "B", "BA")))) {
    f0 <- 1 / (2 * N) * sum(Y_A + Y_B)
    VY <- 1 / (2 * N - 1) * sum((Y_A - f0)^2 + (Y_B - f0)^2)
  }

  # Define first-order estimators
  # --------------------------------------------------------------------

  # Define variance for estimators with A, B, AB; or A, B, BA matrices
  if (first == "saltelli" | first == "jansen" | first == "sobol") {
    f0 <- 1 / (2 * N) * sum(Y_A + Y_B)
    VY <- 1 / (2 * N - 1) * sum((Y_A - f0)^2 + (Y_B - f0)^2)
  }

  # ----------------------------------
  if (first == "sobol") {
    Vi <- 1 / N * Rfast::colsums(Y_A * Y_BA - f0^2)

  } else if (first == "saltelli") {
    Vi <- 1 / N * Rfast::colsums(Y_B * (Y_AB - Y_A))

  } else if (first == "jansen") {
    Vi <- VY - 1 / (2 * N) * Rfast::colsums((Y_B - Y_AB)^2)

  } else if (first == "azzini") {
    VY <- Rfast::colsums((Y_A - Y_B)^2 + (Y_BA - Y_AB)^2)
    Vi <- (2 * Rfast::colsums((Y_BA - Y_B) * (Y_A - Y_AB)))

  } else {
    stop("first should be sobol, saltelli, jansen or azzini")
  }

  if (first == "azzini") {
    Si <- Vi[1:length(params)] / VY[1:length(params)]

  } else {
    Si <- Vi[1:length(params)] / VY
  }

  # Define total-order estimators
  # --------------------------------------------------------------------

  # Define variance for estimators with A, B, AB; or A, B, BA matrices
  if (total == "azzini" | total == "jansen" | total == "sobol" |
      total == "homma" | total == "janon" | total == "glen" | total == "saltelli") {

    f0 <- 1 / (2 * N) * sum(Y_A + Y_B)
    VY <- 1 / (2 * N - 1) * sum((Y_A - f0)^2 + (Y_B - f0)^2)
  }

  # ----------------------------------
  if (total == "jansen") {
    Ti <- (1 / (2 * N) * Rfast::colsums((Y_A - Y_AB)^2)) / VY

  } else if (total == "sobol") {
    Ti <- ((1 / N) * Rfast::colsums(Y_A * (Y_A - Y_AB))) / VY

  } else if (total == "homma") {
    Ti <- (VY - (1 / N) * Rfast::colsums(Y_A * Y_AB) + f0^2) / VY

  } else if (total == "saltelli") {
    Ti <- 1 - ((1 / N * Rfast::colsums(Y_B * Y_BA - f0^2)) / VY)

  } else if (total == "janon") {
    Ti <- 1 - (1 / N * Rfast::colsums(Y_A * Y_AB) -
                 (1/ N * Rfast::colsums((Y_A + Y_AB) / 2))^2) /
      (1 / N * Rfast::colsums((Y_A ^ 2 + Y_AB^2) / 2) -
         (1/ N * Rfast::colsums((Y_A + Y_AB) / 2))^2)

  } else if (total == "glen") {
    Ti <- 1 - (1 / (N - 1) *
                 Rfast::colsums(((Y_A - mean(Y_A)) * (Y_AB - Rfast::colmeans(Y_AB))) /
                                  sqrt(stats::var(Y_A) * Rfast::colVars(Y_AB))))

  } else if (total == "azzini") {
    Ti <- Rfast::colsums((Y_B - Y_BA)^2 + (Y_A - Y_AB)^2) /
      Rfast::colsums((Y_A - Y_B)^2 + (Y_BA - Y_AB)^2)

  } else {
    stop("total should be jansen, sobol, homma saltelli, janon, glen or azzini")
  }
  Ti <- Ti[1:length(params)]

  # Define computation of second-order indices
  # ---------------------------------------------------------------------

  if (order == "second" | order == "third" | order == "fourth") {
    com2 <- utils::combn(1:length(params), 2, simplify = FALSE)
    tmp2 <- do.call(rbind, com2)
    Vi2 <- Vi[(length(params) + 1):(length(params) + length(com2))] # Second order
    Vi1 <- lapply(com2, function(x) Vi[x])
    final.pairwise <- t(mapply(c, Vi2, Vi1))
    Vij <- unname(apply(final.pairwise, 1, function(x) Reduce("-", x)))

    if (first == "azzini") {
      VY <- Rfast::colsums((Y_A - Y_B)^2 + (Y_BA - Y_AB)^2)
      Sij <- Vij / VY[(length(params) + 1):(length(params) + ncol(utils::combn(1:length(params), 2)))]

    } else {
      Sij <- Vij / VY
    }

  } else {
    Sij <- NULL
  }

  # Define computation of third-order indices
  # ---------------------------------------------------------------------

  if (order == "third" | order == "fourth") {
    com3 <- utils::combn(1:length(params), 3, simplify = FALSE)
    tmp3 <- do.call(rbind, com3)
    Vi3 <- Vi[(length(params) + length(com2) + 1):(length(params) + length(com2) + length(com3))] # Third order
    Vi1 <- lapply(com3, function(x) Vi[x])

    # Pairs
    tmp <- do.call(rbind, com2)
    Vij.vec <- as.numeric(paste(tmp[, 1], tmp[, 2], sep = ""))
    names(Vij) <- Vij.vec
    Vi2.1 <- unname(Vij[paste(tmp3[, 1], tmp3[, 2], sep  = "")])
    Vi2.2 <- unname(Vij[paste(tmp3[, 1], tmp3[, 3], sep  = "")])
    Vi2.3 <- unname(Vij[paste(tmp3[, 2], tmp3[, 3], sep  = "")])

    mat3 <- cbind(Vi3, Vi2.1, Vi2.2, Vi2.3, do.call(rbind, Vi1))
    Vijk <- unname(apply(mat3, 1, function(x) Reduce("-", x)))

    if (first == "azzini") {
      Sijl <- Vijk / VY[(length(params) + length(com2) + 1):(length(params) + length(com2) + length(com3))]

    } else {
      Sijl <- Vijk / VY
    }

  } else {
    Sijl <- NULL
  }

  # Define computation of fourth-order indices
  # ---------------------------------------------------------------------

  if(order == "fourth") {
    com4 <- utils::combn(1:length(params), 4, simplify = FALSE)
    tmp4 <- do.call(rbind, com4)
    Vi4 <- Vi[(length(params) + length(com2) + length(com3) + 1):
                (length(params) + length(com2) + length(com3) + length(com4))] # Fourth order
    Vi1 <- lapply(com4, function(x) Vi[x])

    # triplets
    tmp <- do.call(rbind, com3)
    Vijk.vec <- as.numeric(paste(tmp[, 1], tmp[, 2], tmp[, 3], sep = ""))
    names(Vijk) <- Vijk.vec
    Vi3.1 <- unname(Vijk[paste(tmp4[, 1], tmp4[, 2], tmp4[, 3], sep  = "")])
    Vi3.2 <- unname(Vijk[paste(tmp4[, 1], tmp4[, 2], tmp4[, 4], sep  = "")])
    Vi3.3 <- unname(Vijk[paste(tmp4[, 1], tmp4[, 3], tmp4[, 4], sep  = "")])
    Vi3.4 <- unname(Vijk[paste(tmp4[, 2], tmp4[, 3], tmp4[, 4], sep  = "")])

    # Pairs
    Vi2.1 <- unname(Vij[paste(tmp4[, 1], tmp4[, 2], sep  = "")])
    Vi2.2 <- unname(Vij[paste(tmp4[, 1], tmp4[, 3], sep  = "")])
    Vi2.3 <- unname(Vij[paste(tmp4[, 1], tmp4[, 4], sep  = "")])
    Vi2.4 <- unname(Vij[paste(tmp4[, 2], tmp4[, 3], sep  = "")])
    Vi2.5 <- unname(Vij[paste(tmp4[, 2], tmp4[, 4], sep  = "")])
    Vi2.6 <- unname(Vij[paste(tmp4[, 3], tmp4[, 4], sep  = "")])

    mat4 <- cbind(Vi4, Vi3.1, Vi3.2, Vi3.3,
                  Vi2.1, Vi2.2, Vi2.3, Vi2.4, Vi2.5,
                  Vi2.6, do.call(rbind, Vi1))

    Vijlm <- unname(apply(mat4, 1, function(x) Reduce("-", x)))

    if (first == "azzini") {
      Sijlm <- Vijlm / VY[(length(params) + length(com2) + length(com3) + 1):
                            (length(params) + length(com2) + length(com3) + length(com4))]

    } else {
      Sijlm <- Vijlm / VY
    }

  } else {
    Sijlm <- NULL
  }
  return(c(Si, Ti, Sij, Sijl, Sijlm))
}


bootstats <- function(b, conf = conf, type = type) {
  p <- length(b$t0)
  lab <- c("original", "bias", "std.error", "low.ci", "high.ci")
  tmp <- as.data.frame(matrix(nrow = p,
                              ncol = length(lab),
                              dimnames = list(NULL, lab)))
  for (i in 1:p) {
    # original estimation, bias, standard deviation
    tmp[i, "original"] <- b$t0[i]
    tmp[i, "bias"] <- mean(b$t[, i]) - b$t0[i]
    tmp[i, "std.error"] <- stats::sd(b$t[, i])
    # confidence interval

    if (type == "norm") {
      ci <- boot::boot.ci(b, index = i, type = "norm", conf = conf)

      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$norm[2]
        tmp[i, "high.ci"] <- ci$norm[3]
      }

    } else if (type == "basic") {
      ci <- boot::boot.ci(b, index = i, type = "basic", conf = conf)
      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$basic[4]
        tmp[i, "high.ci"] <- ci$basic[5]
      }

    } else if (type == "percent") {
      ci <- boot::boot.ci(b, index = i, type = "perc", conf = conf)

      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$percent[4]
        tmp[i, "high.ci"] <- ci$percent[5]
      }

    } else if (type == "bca") {
      ci <- boot::boot.ci(b, index = i, conf = conf)

      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$bca[4]
        tmp[i, "high.ci"] <- ci$bca[5]
      }
    }
  }
  return(tmp)
}


#' Computation of Sobol' indices
#'
#' It allows to compute Sobol' indices up to the fourth-order using state-of-the-art estimators.
#'
#'@param matrices Character vector with the required matrices. The default is \code{matrices = c("A", "B", "AB")}.
#' See \code{\link{sobol_matrices}}.
#' @param Y  Numeric vector with the model output obtained from the matrix created with
#' \code{\link{sobol_matrices}}. The numeric vector should not contain any NA or NaN values.
#' @param N Positive integer, the initial sample size of the base sample matrix created with \code{\link{sobol_matrices}}.
#' @param params Character vector with the name of the model inputs.
#' @param first Estimator to compute first-order indices. Options are:
#' * \code{first = "saltelli"} \insertCite{Saltelli2010a}{sensobol}.
#' * \code{first = "jansen"} \insertCite{Jansen1999}{sensobol}.
#' * \code{first = "sobol"}  \insertCite{Sobol1993}{sensobol}.
#' * \code{first = "azzini"} \insertCite{Azzini2020}{sensobol}.
#' @param total Estimator to compute total-order indices. Options are:
#' * \code{total = "jansen"} \insertCite{Jansen1999}{sensobol}.
#' * \code{total = "sobol"} \insertCite{Sobol2001}{sensobol}.
#' * \code{total = "homma"} \insertCite{Homma1996}{sensobol}.
#' * \code{total = "janon"} \insertCite{Janon2014}{sensobol}.
#' * \code{total = "glen"} \insertCite{Glen2012}{sensobol}.
#' * \code{total = "azzini"} \insertCite{Azzini2020}{sensobol}.
#' * \code{total = "saltelli"} \insertCite{Saltelli2008}{sensobol}.
#' @param order Whether to compute "first", "second", "third" or fourth-order Sobol' indices. Default
#' is \code{order = "first"}.
#' @param boot Logical. If TRUE, the function bootstraps the Sobol' indices. If FALSE, it provides point
#' estimates. Default is \code{boot = FALSE}.
#' @param R Positive integer, number of bootstrap replicas. Default is NULL.
#' @param parallel The type of parallel operation to be used (if any).
#' If missing, the default is taken from the option "boot.parallel"
#' (and if that is not set, "no"). For more information, check the
#' \code{parallel} option in the \code{boot} function of the \code{\link[boot]{boot}} package.
#' @param ncpus Positive integer: number of processes to be used in parallel operation:
#' typically one would chose this to the number of available CPUs.
#' Check the \code{ncpus} option in the \code{boot} function of the \code{\link[boot]{boot}} package.
#' @param conf Confidence interval if \code{boot = TRUE}. Number between 0 and 1. Default is \code{conf = 0.95}.
#' @param type Method to compute the confidence interval if \code{boot = TRUE}. Default is "norm".
#' Check the \code{type} option in the \code{boot} function of the \code{\link[boot]{boot}} package.
#' @importFrom rlang ":="
#' @importFrom Rdpack reprompt
#' @importFrom stats var
#' @references
#' \insertAllCited{}
#'
#' @return A \code{sensobol} object.
#' @seealso Check the function \code{\link[boot]{boot}} for further details on the bootstrapping
#' with regards to the methods available for the computation of confidence intervals in the \code{type} argument.
#' @export
#'
#' @details Any first and total-order estimator can be combined with the appropriate sampling design.
#' Check Table 3 of the vignette for a summary of all possible combinations, and Tables 1 and 2 for a
#' mathematical description of the estimators. If the analyst mismatches estimators and sampling designs,
#' the function will generate an error and urge to redefine the sample matrices or the estimators.
#'
#' For all estimators except \insertCite{Azzini2020;textual}{sensobol}'s and \insertCite{Janon2014;textual}{sensobol}'s,
#' \code{sobol_indices()} calculates the sample mean as \deqn{\hat{f}_0=\frac{1}{2N} \sum_{v=1}^{N}(f(\mathbf{A})_v + f(\mathbf{B})_v)\,,}
#' where \eqn{N} is the row dimension of the base sample matrix, and the unconditional sample variance as
#'
#' \deqn{\hat{V}(y) = \frac{1}{2N-1} \sum{v=1}^{N} ((f(\mathbf{A})_v - \hat{f})^2 + (f(\mathbf{B})_v - \hat{f})^2)\,,}
#' where \eqn{f(\mathbf{A})_v} (\eqn{f(\mathbf{B})_v}) indicates the model output \eqn{y} obtained after running the model \eqn{f}
#' in the \eqn{v}-th row of the \eqn{\mathbf{A}} (\eqn{\mathbf{B}}) matrix.
#'
#' For the Azzini estimator,
#' \deqn{\hat{V}(y) = \sum_{v=1}^{N} (f(\mathbf{A})_v - f(\mathbf{B})_v)^2 + (f(\mathbf{B}_A^{(i)})_v - f(\mathbf{A}_B^{(i)})_v) ^ 2}
#'
#' and for the Janon estimator,
#' \deqn{\hat{V}(y)=\frac{1}{N} \sum_{v=1}^{N} \frac{f(\mathbf{A})_v^2 + f(\mathbf{A}_B^{(i)})_v^2}{2}-f_0^2}
#'
#'where \eqn{f(\mathbf{A}_B^{(i)})_v} (\eqn{f(\mathbf{B}_A^{(i)})_v}) is the model output obtained after running the model \eqn{f} in
#'the \eqn{v}-th row of an \eqn{\mathbf{A}_B^{(i)})_v} (\eqn{\mathbf{B}_A^{(i)})_v}) matrix, where all columns come from \eqn{\mathbf{A}} (\eqn{\mathbf{B}})
#'except the \eqn{i}-th, which comes from \eqn{\mathbf{B}} (\eqn{\mathbf{A}}).
#'
#' @examples
#' # Define settings
#' N <- 1000; params <- paste("X", 1:3, sep = ""); R <- 10
#'
#' # Create sample matrix
#' mat <- sobol_matrices(N = N, params = params)
#'
#' # Compute Ishigami function
#' Y <- ishigami_Fun(mat)
#'
#' # Compute and bootstrap Sobol' indices
#' ind <- sobol_indices(Y = Y, N = N, params = params, boot = TRUE, R = R)
sobol_indices <- function(matrices = c("A", "B", "AB"), Y, N, params,
                          first = "saltelli", total = "jansen",
                          order = "first", boot = FALSE, R = NULL,
                          parallel = "no", ncpus = 1, conf = 0.95, type = "norm") {

  # Check concordance between boot and R arguments
  # ---------------------------------------------------------------------

  if (boot == FALSE & is.null(R) == FALSE | boot == TRUE & is.null(R) == TRUE) {
    stop("Bootstrapping requires boot = TRUE and an integer in R")
  }

  # Define parameters
  # ----------------------------------------------------------------------

  sensitivity <- parameters <- NULL
  k <- length(params)
  d <- matrix(Y, nrow = N)

  # Function when boot = FALSE
  # -----------------------------------------------------------------------

  if (boot == FALSE) {
    tmp <- sobol_boot(d = d, N = N, params = params, first = first, total = total,
                      order = order, boot = FALSE, matrices = matrices)
    out <- data.table::data.table(tmp)
    data.table::setnames(out, "tmp", "original")

    # Function when boot = TRUE
    # -----------------------------------------------------------------------

  } else if (boot == TRUE) {
    tmp <- boot::boot(data = d, statistic = sobol_boot, R = R, N = N, params = params,
                      first = first, total = total, order = order, matrices = matrices,
                      parallel = parallel, ncpus = ncpus, boot = TRUE)
    out <- data.table::data.table(bootstats(tmp, conf = conf, type = type))

  } else {
    stop("boot has to be TRUE or FALSE")
  }

  # Vectors of parameters and sensitivity indices when order = FIRST
  # -----------------------------------------------------------------------

  if (order == "first") {
    parameters <- c(rep(params, times = 2))
    sensitivity <- c(rep(c("Si", "Ti"), each = k))

    # Vectors of parameters and sensitivity indices when order = second
    # -----------------------------------------------------------------------

  } else if (order == "second") {
    vector.second <- unlist(lapply(utils::combn(params, 2, simplify = FALSE), function(x)
      paste0(x, collapse = ".")))
    parameters <- c(c(rep(params, times = 2)), vector.second)
    sensitivity <- c(rep(c("Si", "Ti"), each = length(params)),
                     rep("Sij", times = length(vector.second)))

    # Vectors of parameters and sensitivity indices when order = third
    # -----------------------------------------------------------------------

  } else if (order == "third") {
    vector.second <- unlist(lapply(utils::combn(params, 2, simplify = FALSE), function(x)
      paste0(x, collapse = ".")))
    parameters <- c(c(rep(params, times = 2)), vector.second)
    vector.third <- unlist(lapply(utils::combn(params, 3, simplify = FALSE), function(x)
      paste0(x, collapse = ".")))
    parameters <- c(parameters, vector.third)
    sensitivity <- c(rep(c("Si", "Ti"), each = k),
                     rep("Sij", times = length(vector.second)),
                     rep("Sijl", times = length(vector.third)))

    # Vectors of parameters and sensitivity indices when order = fourth
    # -----------------------------------------------------------------------

  } else if (order == "fourth") {
    vector.second <- unlist(lapply(utils::combn(params, 2, simplify = FALSE), function(x)
      paste0(x, collapse = ".")))
    parameters <- c(c(rep(params, times = 2)), vector.second)
    vector.third <- unlist(lapply(utils::combn(params, 3, simplify = FALSE), function(x)
      paste0(x, collapse = ".")))
    parameters <- c(parameters, vector.third)
    vector.fourth <- unlist(lapply(utils::combn(params, 4, simplify = FALSE), function(x)
      paste0(x, collapse = ".")))
    parameters <- c(parameters, vector.fourth)
    sensitivity <- c(rep(c("Si", "Ti"), each = k),
                     rep("Sij", times = length(vector.second)),
                     rep("Sijl", times = length(vector.third)),
                     rep("Sijlm", times = length(vector.fourth)))
  } else {

    stop("order has to be first, second, third or fourth")
  }

  # Create class and output
  # ----------------------------------------------------------------------

  ind <- structure(list(), class = "sensobol") # Create class
  ind$results <- cbind(out, sensitivity, parameters) # Add Sobol' indices
  original <- NULL
  ind$si.sum <- ind$results[sensitivity == "Si", sum(original)] # Sum of first-order indices
  ind$first <- first # First-order estimator
  ind$total <- total # Total-order estimator
  ind$C <- length(Y) # Total number of model runs
  return(ind)
}

#' Display the results obtained with the \code{sobol_indices} function.
#'
#' @param x A \code{sensobol} object produced by \code{sobol_indices}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The function \code{print.sensobol} informs on the first and total-order
#' estimators used in the computations, the total number of model runs and
#' the sum of first-order index. It also plots the estimated results.
#' @export
#'
print.sensobol <- function(x, ...) {
  cat("\nFirst-order estimator:", x$first, "| Total-order estimator:", x$total, "\n")
  cat("\nTotal number of model runs:", x$C, "\n")
  cat("\nSum of first order indices:", x$si.sum, "\n")
  print(data.table::data.table(x$results))
}

