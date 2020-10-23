#' Compute a correlation matrix with bootstrap confidence limits
#' @param x A numeric matrix or data frame whose correlation matrix is required.
#' @param method Argument to the function that computes the correlation. Defaults
#'   to \code{method="spearman"}.
#' @param alpha Determines the size of the confidence interval. Defaults to
#'   \code{alpha=0.05} and a 95\% interval is computed.
#' @param B The number of bootstrap samples to run. Defaults to \code{B=1000}.
#' @param corfun The function to be used to compute the correlation. Defaults to
#'   \code{corfun=cor}. If you substitute in another (e.g. a robust correlation
#'   function) it must have an argument called 'method' so you might need to
#'   write a wrapper for it.
#' @export correlation
correlation <- function(x, method = "spearman", summary = "actual",
                        alpha = .05, B = 1000, corfun = cor, digits = 3){
  theCall <- match.call()

  if (!(summary %in% c("actual", "mean", "median"))){
    stop("summary should be 'actual', 'mean' or 'median'")
  }

  # Compute the correlation
  sco <- corfun(x, method=method, use = "pairwise")

  # Get bootstrap versions of it
  bfun <- function(x){
    i <- sample(1:nrow(x), size=nrow(x), replace=TRUE)
    x <- x[i, ]
    corfun(x, method = method, use = "pairwise")
  }

  bscor <- replicate(B, bfun(x))

  # Get confidence intervals
  lower <- apply(bscor, 1:2, quantile, prob=alpha/2)
  upper <- apply(bscor, 1:2, quantile, prob=1 - alpha/2)

  if (summary == "actual"){
    cor <- sco
  } else if (summary == "mean"){
    cor <- apply(bscor, 1:2, mean)
  } else if (summary == "median"){
    cor <- apply(bscor, 1:2, median)
  }

  ci <- paste0("(", round(lower, digits), ", ", round(upper, digits), ")")
  ch <- paste0(round(sco, digits), " (", round(lower, digits), ", ", round(upper, digits), ")")

  ci <- matrix(ch, ncol = ncol(cor))
  ch <- matrix(ch, ncol = ncol(cor))

  diag(cor) <- diag(ch) <- 1

  res <- list(correlation = cor, lower = lower, upper = upper, ci = ci, string = ch, call = theCall)
  class(res) <- "correlation"
  res
}

#' @method ggplot correlation
#' @export
ggplot.correlation <- function(data, mapping=aes(), which="real", ...,
                               environment = parent.frame()){
  d <- data$correlation
  if (which == "nullify"){
    d[data$upper > 0 & data$lower < 0] <- 0
  } else if (which != "real"){
    stop("which should be either 'nullify' or 'real'")
  }

  d <- as.data.frame(d, check.names=FALSE) %>%
    mutate(vars = rownames(.)) %>%
    gather(vars2, correlation, -vars)

  ggplot(d, aes(vars, vars2)) +
    geom_tile(aes(fill=correlation), color="white") +
    scale_fill_gradient(low="white", high="steelblue") +
    scale_x_discrete("") + scale_y_discrete("") +
    theme(axis.text.y=element_text(size=12),
          axis.text.x = element_text(size=12, angle=45, hjust=1)) +
    coord_fixed()
}

#' @method print correlation
#' @export
print.correlation <- function(x, ...){
  x$ch
}
