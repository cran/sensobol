

# PLOT SOBOL' FIRST AND TOTAL-ORDER INDICES -----------------------------------

#' Plot Sobol' first and total-order indices
#'
#' @param x A data.table.
#' @param dummy The output of the \code{sobol_ci_dummy} function. If supplied and
#' \code{type = 1}, the plot includes an horizontal transparent frame showing the confidence
#' intervals of the first and total-order indices for the dummy parameter.
#' @param type An integer. If \code{type = 1}, it plots first and total effects.
#' If \code{type = 2}, it plots second-order effects. If \code{type = 3}, it plots
#' third-order effects. Default is \code{type = 1}.
#'
#' @return A ggplot object.
#' @import ggplot2
#' @export
#'
#' @examples
#' # Define settings:
#' n <- 500; k <- 8; R <- 100
#' # Design the sample matrix:
#' A <- sobol_matrices(n = n, k = k, second = TRUE, third = TRUE)
#' # Compute the model output:
#' Y <- sobol_Fun(A)
#' # Compute the Sobol' indices:
#' \donttest{sens <- sobol_indices(Y = Y, params = colnames(data.frame(A)),
#' R = R, n = n, parallel = "no", ncpus = 1, second = TRUE, third = TRUE)
#' # Compute the Sobol' indices for the dummy parameter:
#' s.dummy <- sobol_dummy(Y = Y, params = colnames(data.frame(A)), R = R, n = n)
#' # Compute confidence intervals:
#' sens.ci <- sobol_ci(sens, params = colnames(data.frame(A)), type = "norm", conf = 0.95)
#' # Compute confidence intervals for the dummy parameter:
#' s.dummy.ci <- sobol_ci_dummy(s.dummy, type = "norm", conf = 0.95)
#' # Plot Sobol' indices:
#' plot_sobol(sens.ci, dummy = s.dummy.ci, type = 1)}
plot_sobol <- function(x, dummy = NULL, type = 1) {
  sensitivity <- low.ci <- high.ci <- parameters <- original <- NULL
  if(type == 1) {
    if(is.null(dummy) == FALSE) {
      plot.dummy <- geom_rect(data = dummy,
                              aes(ymin = 0,
                                  ymax = high.ci,
                                  xmin = -Inf,
                                  xmax = Inf,
                                  fill = sensitivity),
                              alpha = 0.2,
                              inherit.aes = FALSE)
    } else {
      plot.dummy <- NULL
    }
    p <- x[sensitivity == "Si" | sensitivity == "STi"]
    gg <- ggplot2::ggplot(p, aes(parameters, original,
                                 fill = sensitivity)) +
      geom_bar(stat = "identity",
               position = position_dodge(0.6),
               color = "black") +
      plot.dummy +
      geom_errorbar(aes(ymin = low.ci,
                        ymax = high.ci),
                    position = position_dodge(0.6)) +
      scale_fill_discrete(name = "Sobol' indices",
                          labels = c(expression(S[italic(i)]),
                                     expression(S[italic(T[i])]))) +
      labs(x = "",
           y = "Variance") +
      theme_bw() +
      theme(legend.position = "top",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.background = element_rect(fill = "transparent",
                                             color = NA),
            legend.key = element_rect(fill = "transparent",
                                      color = NA))
  } else if(!type == 1) {
    if(type == 2) {
      plot.type <- "Sij"
    } else if(type == 3) {
      plot.type <- "Sijk"
    } else {
      stop("Type should be either 1, 2 or 3")
    }
    p <- x[sensitivity == plot.type]
    gg <- ggplot2::ggplot(p, aes(stats::reorder(parameters, original),
                    original)) +
      geom_point() +
      geom_errorbar(aes(ymin = low.ci,
                        ymax = high.ci)) +
      theme_bw() +
      labs(x = "",
           y = "Variance") +
      geom_hline(yintercept = 0,
                 lty = 2,
                 color = "red") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.background = element_rect(fill = "transparent",
                                             color = NA),
            legend.key = element_rect(fill = "transparent",
                                      color = NA),
            axis.text.x = element_text(angle = 45,
                                       hjust = 1))
  }
  return(gg)
}

# PLOT MODEL OUTPUT UNCERTAINTY -----------------------------------------------

#' Plot model output uncertainty
#'
#' It creates an histogram with the model output distribution.
#'
#' @param Y A numeric vector with the model output.
#' @param n An integer with the initial sample size.
#'
#' @return a ggplot2 object.
#' @import ggplot2
#' @export
#'
#' @examples
#' # Define settings:
#' n <- 100; k <- 8; R <- 10
#' # Design the sample matrix:
#' A <- sobol_matrices(n = n, k = k, second = FALSE, third = FALSE)
#' # Compute the model output:
#' Y <- sobol_Fun(A)
#' # Plot the model output distribution:
#' plot_uncertainty(Y, n = n)
plot_uncertainty <- function(Y, n = NULL) {
  warning("The argument n will be substituted by N in the next release of the package")
  if(is.vector(Y) == FALSE) {
    stop("Y should be a vector")
  }
  if(is.null(n) == TRUE) {
    stop("n should be specified")
  }
  Y <- Y[1:(2*n)]
  df <- data.frame(Y)
  gg <- ggplot2::ggplot(df, aes(Y)) +
    geom_histogram(color = "black",
                   fill = "white") +
    labs(x = "Y",
         y = "Count") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "transparent",
                                           color = NA),
          legend.key = element_rect(fill = "transparent",
                                    color = NA))
  return(gg)
}

# PLOT SCATTERPLOT OF MODEL OUTPUT VERSUS MODEL INPUTS ------------------------

#' Scatterplots of the model output against the model inputs
#'
#' @param x A data table, data frame or matrix with the
#' model inputs.
#' @param Y Numeric vector with the model output.
#' @param n Integer, sample size of the Sobol' matrix.
#' @param params Vector with the name of the model inputs.
#'
#' @return A ggplot object.
#' @import ggplot2
#' @export
#'
#' @examples
#' # Define settings:
#' n <- 100; k <- 8; R <- 10
#' # Design the sample matrix:
#' A <- sobol_matrices(n = n, k = k, second = TRUE, third = TRUE)
#' # Compute the model output:
#' Y <- sobol_Fun(A)
#' # Plot scatterplots:
#' plot_scatter(x = A, Y = Y, n = n, params = colnames(data.frame(A)))
plot_scatter <- function(x, Y, n, params) {
  warning("The argument n will be substituted by N in the next release of the package")
  value <- NULL
  dt <- data.table::data.table(cbind(x, Y))
  # Retrieve the A and B matrices only
  dt <- dt[1:(2 * n)]
  data.table::setnames(dt, 1:ncol(dt), paste(c(params, "Y")))
  dt.melted <- data.table::melt(dt, measure.vars = 1:length(params),
                                variable.name = "Parameters")
  gg <- ggplot2::ggplot(dt.melted, aes(x = value, y = Y)) +
    geom_hex() +
    labs(x = "Value",
         y = "Y") +
    facet_wrap(~Parameters,
               scales = "free_x") +
    theme_bw() +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_fill_continuous(name = "Count") +
    theme(legend.position = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "transparent",
                                           color = NA),
          legend.key = element_rect(fill = "transparent",
                                    color = NA))
  return(gg)
}
