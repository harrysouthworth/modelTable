#' Get a nicely formatted table for including in a report
#' @details It is common to want to get a summary of a model in terms of
#'   parameter estimates, their standard errors, associated z-statistics or
#'   t-statistics, and p-values. Some functions don't return the p-values
#'   (such as \code{glm}, \code{MASS::rlm}), and p-values want formatting
#'   (like "p < 0.0001") and trailing 0s are required for any rounded numbers.
#'   This small package is intended to provide a simple way of achieving this.
#'   Methods are available for objects of class 'lm', 'glm', 'rlm', 'lmrob',
#'   'polr' (proportion odds, as fitted by \code{MASS::polr})
#' @param x A model object, such as an object returned by \code{lm}, \code{glm},
#'   or many others.
#' @param fmt String to be passed to \code{sprintf}. Defaults to \code{fmt = "%.4f"}
#'   and numbers are rounded to 4 decimal places, with trailing 0s left intact.
#' @param ... Currently unused.
#' @note If you have factors in your model, you'll likely want to write a
#'   simple function to sanitiize the rownames of the returned object.
#' @return A character matrix.
#' @export modelTable
modelTable <- function(x, fmt = "%.4f", ...){
  UseMethod("modelTable")
}

#' @method modelTable rlm
#' @export
modelTable.rlm <- function(x, fmt = "%.4f"){
  co <- coef(summary(x))
  df <- x$df

  p <- 2 * (1 - pt(abs(co[, 3]), df))
  co <- cbind(co, p)
  colnames(co)[4] <- "p-value"

  formatModelTable(co, fmt = fmt)
}

#' @method modelTable lmrob
#' @export
modelTable.lmrob <- function(x, fmt = "%.4f"){
  co <- coef(summary(x))
  formatModelTable(co, fmt = fmt)
}

#' @method modelTable censtreg
#' @export
modelTable.censtreg <- function(x, fmt = "%.4f"){
  co <- coef(summary(x))[, c(1, 3)]
  df <- nrow(x$data) - nrow(co)

  tt <- co[, 1] / co[, 2]
  p <- 2 * (1 - pt(abs(tt), df))

  co <- cbind(co, tt, p)
  colnames(co)[3:4] <- c("t-value", "p-value")

  formatModelTable(co, fmt = fmt)
}

#' Format a summary table from a model
#' @param x An object summarizing a model, either a matrix or data frame.
#' @param fmt To be passed into \code{sprtinf}, defaults to \code{fmt = "%.4f"}.
#' @details It is assumed that the 4th column is p-values, and the function
#'   either formats them using \code{sprintf} (as it does to other columns),
#'   or turns them into "<0.0001".
#' @warning "<0.0001" is hardcoded but shouldn't be.
#' @export formatModelTable
formatModelTable <- function(x, fmt = "%.4f"){
  rn <- rownames(x)
  cn <- colnames(x)
  x <- apply(x, 2, sprintf, fmt = fmt)

  x[, 4] <- ifelse(x[, 4] == "0.0000", "<0.0001", x[, 4])

  rownames(x) <- rn
  colnames(x) <- cn

  rownames(x) <- gsub("arm", "", rownames(x))
  colnames(x)[4] <- "p-value"

  x
}

#' @method modelTable glm
#' @export
modelTable.glm <- function(x, fmt = "%.4f"){
  x <- coef(summary(x))
  formatModelTable(x, fmt = fmt)
}

#' @method modelTable lm
#' @export
modelTable.lm <- modelTable.glm

#' @method modelTable polr
#' @export
modelTable.polr <- function(x, fmt = "%.4f"){
  s <- suppressMessages(summary(x))
  co <- coef(s)
  p <- 2 * (1 - pt(abs(co[, 3]), x$df.residual))
  co <- cbind(co, p)
  colnames(co)[4] <- "p-value"

  co <- co[c((s$pc + 1):nrow(co), 1:s$pc), ]

  formatModelTable(co, fmt = fmt)
}

#' @method modelTable gam
#' @export
modelTable.gam <- function(x, fmt = "%.4f"){
  s <- summary(x)
  co <- s$p.table
  colnames(co)[4] <- "p-value"

  formatModelTable(co, fmt = fmt)
}

