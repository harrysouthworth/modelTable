#' Get a nicely formatted table for including in a report
#' @details It is common to want to get a summary of a model in terms of
#'   parameter estimates, their standard errors, associated z-statistics or
#'   t-statistics, and p-values. Some functions don't return the p-values
#'   (such as \code{glm}, \code{MASS::rlm}), and p-values want formatting
#'   (like "p < 0.0001") and trailing 0s are required for any rounded numbers.
#'   This small package is intended to provide a simple way of achieving this.
#'   Methods are available for objects of class 'lm', 'glm', 'rlm', 'lmrob'
#' @param x A model object, such as an object returned by \code{lm}, \code{glm},
#'   or many others.
#' @note If you have factors in your model, you'll likely want to write a
#'   simple function to sanitiize the rownames of the returned object.
#' @return A character matrix.
modelTable <- function(x, ...){
  UseMethod("modelTable")
}


modelTable.rlm <- function(x){
  co <- coef(summary(x))
  df <- x$df

  p <- 2 * (1 - pt(abs(co[, 3]), df))
  co <- cbind(co, p)
  colnames(co)[4] <- "p-value"

  formatCoTable(co)
}

modelTable.lmrob <- function(x){
  co <- coef(summary(x))
  formatCoTable(co)
}

modelTable.censtreg <- function(x){
  co <- coef(summary(x))[, c(1, 3)]
  df <- nrow(x$data) - nrow(co)

  tt <- co[, 1] / co[, 2]
  p <- 2 * (1 - pt(abs(tt), df))

  co <- cbind(co, tt, p)
  colnames(co)[3:4] <- c("t-value", "p-value")

  formatCoTable(co)
}

formatModelTable <- function(x){
  rn <- rownames(x)
  cn <- colnames(x)
  x <- apply(x, 2, sprintf, fmt = "%.4f")

  x[, 4] <- ifelse(x[, 4] == "0.0000", "<0.0001", x[, 4])

  rownames(x) <- rn
  colnames(x) <- cn

  rownames(x) <- gsub("arm", "", rownames(x))
  colnames(x)[4] <- "p-value"

  x
}

modelTable.glm <- function(x){
  x <- coef(summary(x))
  formatCoTable(x)
}

modelTable.lm <- getCoTable.glm

modelTable.polr <- function(x){
  s <- suppressMessages(summary(x))
  co <- coef(s)
  p <- 2 * (1 - pt(abs(co[, 3]), x$df.residual))
  co <- cbind(co, p)
  colnames(co)[4] <- "p-value"

  co <- co[c((s$pc + 1):nrow(co), 1:s$pc), ]

  formatCoTable(co)
}
