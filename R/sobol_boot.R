
# FUNCTION TO COMPUTE BOOTSTRAP CONFIDENCE INTERVALS --------------------------

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

sobol_ci_temp <- function(b, params, type, conf, second = FALSE, third = FALSE) {
  V1 <- NULL
  # Extract first and total order effects
  ext1 <- b[, .SD[1:length(params)]][, V1]
  Si.STi <- lapply(ext1, function(x) bootstats(x,
                                               type = type,
                                               conf = conf)) %>%
    data.table::rbindlist()
  if(second == TRUE) {
    # Extract second order effects
     ext2 <- b[, .SD[(length(params) + 1) :(length(params) +
                                              factorial(length(params)) /
                                              (factorial(2) * factorial(length(params) - 2)))]][, V1]
     Sij <- lapply(ext2, function(x) bootstats(x,
                                               type = type,
                                               conf = conf)) %>%
       data.table::rbindlist()
  } else {
    Sij <- NULL
  }
  if(third == TRUE) {
    ext3 <- b[, utils::tail(.SD, factorial(length(params)) /
                              (factorial(3) * factorial(length(params) - 3)))][, V1]
    Sijk <- lapply(ext3, function(x) bootstats(x,
                                               type = type,
                                               conf = conf)) %>%
      data.table::rbindlist()
  } else {
    Sijk <- NULL
  }
  return(rbind(Si.STi, Sij, Sijk))
}

# CREATE TWO VECTORS WITH SENSITIVITY INDICES AND MODEL INPUTS ----------------

create_vectors <- function(params, second = FALSE, third = FALSE) {
  # Define for first and total only
  parameters <- rep(params, each = 2)
  sensitivity <- rep(c("Si", "STi"), times = length(params))
  # Define for second only
  pairs <- utils::combn(params, 2, simplify = FALSE)
  parameters.second <- lapply(pairs, function(x) paste0(x, collapse = ".")) %>%
    unlist()
  sensitivity.second <- rep("Sij", times = (length(params) * (length(params) - 1) / 2))
  # Define for third only
  if(length(params) > 2) {
    triplet <- utils::combn(params, 3, simplify = FALSE)
    parameters.third <- lapply(triplet, function(x) paste0(x, collapse = ".")) %>%
      unlist()
    sensitivity.third <- rep("Sijk", times = factorial(length(params)) /
                               (factorial(3) * factorial((length(params) - 3))))
  }
  out <- data.table::as.data.table(cbind(parameters, sensitivity))
  # Set conditions
  if(second == TRUE & third == FALSE) {
    parameters <- c(parameters, parameters.second)
    sensitivity <- c(sensitivity, sensitivity.second)
    out <- data.table::as.data.table(cbind(parameters, sensitivity))
  }
  if(second == FALSE & third == TRUE) {
    parameters <- c(parameters, parameters.third)
    sensitivity <- c(sensitivity, sensitivity.third)
    out <- data.table::as.data.table(cbind(parameters, sensitivity))
  }
  if(second == TRUE & third == TRUE) {
    parameters <- c(parameters, parameters.second, parameters.third)
    sensitivity <- c(sensitivity, sensitivity.second, sensitivity.third)
    out <- data.table::as.data.table(cbind(parameters, sensitivity))
  }
  return(out)
}

#' Bootstrap confidence intervals for Sobol' indices.
#'
#' It computes bootstrap confidence intervals for Sobol' indices.
#'
#' @param b The output of the \code{sobol_indices} function.
#' @param params A vector with the name of the model inputs.
#' @param type A vector of character strings representing
#' the type of intervals required. The value should be any subset
#' of the values \code{c("norm", "basic", "perc", "bca")}. For more information,
#' check the function \code{\link{boot.ci}}.
#' @param conf A scalar or vector containing the confidence
#' level(s) of the required interval(s).
#' @param second Logical. If \code{second = TRUE}, it computes the confidence
#' intervals for second-order indices. Default is \code{second = FALSE}.
#' @param third Logical. If \code{third = TRUE}, it computes the confidence
#' intervals for third-order indices. Default is \code{third = FALSE}.
#'
#' @return A data table.
#' @seealso \code{\link{boot}}, \code{\link{boot.ci}}.
#'
#' @export
#'
#' @examples
#' # Define settings:
#' n <- 1000; k <- 8; R <- 100
#' # Design the sample matrix:
#' A <- sobol_matrices(n = n, k = k, second = TRUE, third = TRUE)
#' # Compute the model output:
#' Y <- sobol_Fun(A)
#' # Compute the Sobol' indices:
#' \donttest{sens <- sobol_indices(Y = Y, params = colnames(data.frame(A)),
#' R = R, n = n, parallel = "no", ncpus = 1,
#' second = TRUE, third = TRUE)
#' # Compute confidence intervals:
#' sobol_ci(sens, params = colnames(data.frame(A)), type = "norm", conf = 0.95)}
sobol_ci <- function(b, params, type, conf, second = FALSE, third = FALSE) {
  out <- sobol_ci_temp(b = b, params = params, type = type,
                        conf = conf, second = second, third = third) %>%
    cbind(create_vectors(params = params, second = second, third = third))
  return(out)
}


# FUNCTION TO COMPUTE CONFIDENCE INTERVALS FOR DUMMY PARAMETER ----------------


#' Bootstrap confidence intervals for the dummy parameter
#'
#' It computes bootstrap confidence intervals for the dummy parameter.
#'
#' @param b The output of the \code{sobol_dummy} function.
#' @param type A vector of character strings representing
#' the type of intervals required. The value should be any subset
#' of the values \code{c("norm", "basic", "perc", "bca")}. For more information,
#' check the function \code{\link{boot.ci}}.
#' @param conf A scalar or vector containing the confidence
#' level(s) of the required interval(s).
#' @importFrom rlang ":="
#'
#' @return A data table.
#' @seealso \code{\link{boot}}, \code{\link{boot.ci}}.
#' @export
#'
#' @examples
#' # Define settings:
#' n <- 100; k <- 8; R <- 10
#' # Design the sample matrix:
#' A <- sobol_matrices(n = n, k = k, second = TRUE, third = TRUE)
#' # Compute the model output:
#' Y <- sobol_Fun(A)
#' # Compute the Sobol' indices for the dummy parameter:
#' s.dummy <- sobol_dummy(Y = Y, params = colnames(data.frame(A)), R = R, n = n)
#' # Compute the confidence intervals for the dummy parameter:
#' sobol_ci_dummy(s.dummy, type = "norm", conf = 0.95)
sobol_ci_dummy <- function(b, type = type, conf = conf) {
  V1 <- NULL
  sensitivity <- c("Si", "STi")
  parameters <- "dummy"
  out <- lapply(b[, V1], function(x) bootstats(x,
                                               type = type,
                                               conf = 0.95)) %>%
    data.table::rbindlist()
  return(cbind(out, parameters, sensitivity))
}
