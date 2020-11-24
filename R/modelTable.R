#' Get a nicely formatted table for including in a report
#' @details It is common to want to get a summary of a model in terms of
#'   parameter estimates, their standard errors, associated z-statistics or
#'   t-statistics, and p-values. Some functions don't return the p-values
#'   (such as \code{glm}, \code{MASS::rlm}), and p-values want formatting
#'   (like "p < 0.0001") and trailing 0s are required for any rounded numbers.
#'   This small package is intended to provide a simple way of achieving this.
#'   Methods are available for objects of class 'lm', 'glm', 'rlm', 'lmrob',
#'   'gam' (from \code{mgcv::gam}),
#'   'polr' (proportion odds, as fitted by \code{MASS::polr})
#' @param x A model object, such as an object returned by \code{lm}, \code{glm},
#'   or various others.
#' @param ci Logical, whether to return a confidence interval as well as the
#'   other stuff. Defaults to \code{ci = FALSE}.
#' @param alpha Numeric, defaulting to \code{alpha = .05} and defining the
#'   level of the confidence interval, if \code{ci = TRUE}.
#' @param fmt String to be passed to \code{sprintf}. Defaults to \code{fmt = "\%.4f"}
#'   and numbers are rounded to 4 decimal places, with trailing 0s left intact.
#' @param ... Currently unused.
#' @note If you have factors in your model, you'll likely want to write a
#'   simple function to sanitize the rownames of the returned object. For
#'   objects of class 'glm', confint returns a profile confidence interval
#'   using \code{stats?::confint.glm}. For objects of class 'censtreg', it
#'   returns quantiles of the posterior distribution. For objects of class
#'   'polr' of 'brglmFit' it returns intervals based on a Gaussian approximation.
#'   NOTE the comments in \code{brglm::confint.brglm}.
#' @return A character matrix.
#' @export modelTable
modelTable <- function(x, ci = FALSE, alpha = .05, fmt = "%.4f", ...){
  UseMethod("modelTable")
}


#' @method modelTable rlm
#' @export
modelTable.rlm <- function(x, ci = FALSE, alpha = .05, fmt = "%.4f"){
  co <- coef(summary(x))
  df <- x$df.residual

  p <- 2 * (1 - pt(abs(co[, 3]), df))
  co <- cbind(co, p)
  colnames(co)[4] <- "p-value"

  if (ci){
    ts <- qt(1 - alpha / 2, df)
    lo <- co[, 1] - ts * co[, 2]
    hi <- co[, 1] + ts * co[, 2]
    co <- cbind(co, lo, hi)
    names(co)[5:6] <- c(paste0("Lower.", 1 - alpha / 2),
                        paste0("Upper.", 1 - alpha / 2))
  }

  formatModelTable(co, fmt = fmt)
}

#' @method modelTable lmrob
#' @export
modelTable.lmrob <- function(x, ci = FALSE, alpha = .05, fmt = "%.4f"){
  co <- coef(summary(x))
  df <- x$df.residual

  if (ci){
    ts <- qt(1 - alpha / 2, df)
    lo <- co[, 1] - ts * co[, 2]
    hi <- co[, 1] + ts * co[, 2]
    co <- cbind(co, lo, hi)
    names(co)[5:6] <- c(paste0("Lower.", 1 - alpha / 2),
                        paste0("Upper.", 1 - alpha / 2))
  }

  formatModelTable(co, fmt = fmt)
}

#' @method modelTable censtreg
#' @export
modelTable.censtreg <- function(x, ci = FALSE, alpha = .05, fmt = "%.4f"){
  co <- coef(summary(x))[, c(1, 3)]
  co <- co[-nrow(co), ]
  df <- nrow(x$data) - nrow(co)

  tt <- co[, 1] / co[, 2]
  p <- 2 * (1 - pt(abs(tt), df))

  co <- cbind(co, tt, p)
  colnames(co)[3:4] <- c("t-value", "p-value")

  if (ci){
    probs <- c(alpha / 2, 1 - alpha / 2)
    ci <- rstan::summary(x$model, probs = probs)$summary[1:nrow(co), 4:5]
    co <- cbind(co, ci)
  }

  formatModelTable(co, fmt = fmt)
}

#' @method modelTable glm
#' @export
modelTable.glm <- function(x, ci = FALSE, alpha = .05, fmt = "%.4f"){
  res <- coef(summary(x))

  if (ci){
    ci <- suppressMessages(confint(x, level = 1 - alpha))
    res <- cbind(res, ci)
  }

  formatModelTable(res, fmt = fmt)
}

#' @method modelTable lm
#' @export
modelTable.lm <- function(x, ci = FALSE, alpha = .05, fmt = "%.4f"){
  co <- coef(summary(x))

  if (ci){
    ci <- confint.lm(x, level = 1 - alpha)
    co <- cbind(co, ci)
  }
  formatModelTable(co, fmt = fmt)
}

#' @method modelTable polr
#' @export
modelTable.polr <- function(x, ci = FALSE, alpha = .05, fmt = "%.4f"){
  s <- suppressMessages(summary(x))
  co <- coef(s)
  p <- 2 * (1 - pt(abs(co[, 3]), x$df.residual))
  co <- cbind(co, p)
  colnames(co)[4] <- "p-value"

  co <- co[c((s$pc + 1):nrow(co), 1:s$pc), ]

  if (ci){
    z <- qnorm(1 - alpha / 2)
    lo <- co[, 1] - z * co[, 2]
    hi <- co[, 1] + z * co[, 2]
    co <- cbind(co, cbind(lo, hi))
  }

  formatModelTable(co, fmt = fmt)
}

modelTable.brglmFit <- function(x, ci = FALSE, alpha = .05, fmt = "%.4f"){
  co <- coef(summary(x))
  colnames(co)[4] <- "p-value"
  if(ci){
    z <- qnorm(1 - alpha / 2)
    lo <- co[, 1] - z * co[, 2]
    hi <- co[, 1] + z * co[, 2]
    co <- cbind(co, cbind(lo, hi))
  }

  formatModelTable(co, fmt = fmt)
}

#' @method modelTable gam
#' @export
modelTable.gam <- function(x, ci = FALSE, alpha = .05, fmt = "%.4f"){
  s <- summary(x)
  co <- s$p.table
  colnames(co)[4] <- "p-value"

  if (ci){
    ci <- confint.lm(x)[1:nrow(co), ]
    co <- cbind(co, ci)
  }

  formatModelTable(co, fmt = fmt)
}

#' @method modelTable coxph
#' @export
modelTable.coxph <- function(x, ci = FALSE, alpha = .05, fmt = "%.4f"){
  co <- coef(summary(x))[, -2]
  colnames(co)[4] <- "p-value"

  if(ci){
    z <- qnorm(1 - alpha / 2)
    lo <- co[, 1] - z * co[, 2]
    hi <- co[, 1] + z * co[, 2]
    co <- cbind(co, cbind(lo, hi))
  }

  formatModelTable(co, fmt = fmt)
}



