#' Format a summary table from a model
#' @param x An object summarizing a model, either a matrix or data frame.
#' @param fmt To be passed into \code{sprintf}, defaults to \code{fmt = "\%.4f"}.
#' @param NAstring How to deal with NAs. Defaults to the empty string.
#' @details It is assumed that the 4th column is p-values, and the function
#'   either formats them using \code{sprintf} (as it does to other columns),
#'   or turns them into "<0.0001".
#' @export formatModelTable
formatModelTable <- function(x, fmt = "%.4f", NAstring = ""){
  if (ncol(x) == 3){ ## Presumably point estimates, sd and t-values
    formatTable(x, types = c("n", "n", "n"), fmt = fmt, NAstring = NAstring)
  }
  if (ncol(x) == 4){
    formatTable(x, types = c("n", "n", "n", "p"), fmt = fmt, NAstring = NAstring)
  } else if (ncol(x) == 6){
    formatTable(x, types = c("n", "n", "n", "p", "lo", "hi"))
  } else {
    stop("I don't know what to do with that number of columns")
  }
}
# formatModelTable <- function(x, fmt = "%.4f"){
#   rn <- rownames(x)
#   cn <- colnames(x)
#   x <- apply(x, 2, sprintf, fmt = fmt)
#
#   x[, 4] <- ifelse(x[, 4] == "0.0000", "<0.0001", x[, 4])
#
#   rownames(x) <- rn
#   colnames(x) <- cn
#
#   rownames(x) <- gsub("arm", "", rownames(x))
#   colnames(x)[4] <- "p-value"
#
#   x
# }

#' Format a matrix or data frame
#' @param x A matrix or data frame.
#' @param types Character vector having the same number of elements as there are
#'   columns in \code{x}.
#' @param fmt To be passed to \code{sprintf}, defaults to \code{fmt = "\%.4f"}.
#' @param NAstring String to replace \code{NA}s. Defaults to \code{NAstring = ""}.
#' @details The allowed types are "i" (integer), "n" (numeric), "p" (p-values),
#'   and "c" (character). If the type is "i" a character representation of the
#'   values coerced to integer is returned; if type is "n" a character representation
#'   of the values after passing through \code{sprintf(fmt, x)} is retunred;
#'   if type is "p", values less than 0.0001 are described as such (i.e. "<0.0001"),
#'   otherwise the returned values from \code{sptritf(fmt, X)}; if type is "c",
#'   the values are returned unaltered.
#' @return A data frame, the columns having been formatted accordingly
#' @export
formatTable <- function(x, types = NULL, fmt = "%.4f", NAstring = ""){
  nc <- ncol(x)

  if (length(types) != nc){
    stop("types must have the same number of entries as x has columns")
  }

  cn <- colnames(x) # This isn't used, but is left here in anticipation of R deciding to remove rownames
  x <- as.data.frame(x)

  np <- as.numeric(substring(fmt, 3, nchar(fmt) - 1))
  np <- 10^(-np)

  reformat <- function(X, type, f){
    X <- if (type == "i"){
      as.character(as.integer(X))
    } else if (type == "n"){
      sprintf(f, X)
    } else if (type == "p"){
      ifelse(X < np, getPstring(f), sprintf(f, X))
    } else if (type == "c"){
      X
    } else if (type == "lo"){
      paste0("(", sprintf(f, X))
    } else if (type == "hi"){
      paste0(sprintf(f, X), ")")

    } else {
      stop("Allowed types are i, n, p and c.")
    }
  }

  for (i in 1:nc){
    x[, i] <- reformat(x[, i], type = types[i], f = fmt)
    x[, i] <- gsub("NA", NAstring, x[, i])
    x[, i][is.na(x[, i])] <- NAstring
  }

  if ("lo" %in% types){
    x[, ncol(x) - 1] <- paste0(x[, ncol(x) - 1], ", ", x[, ncol(x)])
    x <- x[, -ncol(x)]
    names(x)[ncol(x)] <- "Conf. Int"
  }

  x
}

#' Get a string to represent small p-values
#' @param f A string of the type passed into \code{sprintf}.
#' @details This function expects its input to be of the form "xxnnnx" where the
#'   xs are stripped off and teh ns are coerced to numeric.
getPstring <- function(f = "%.4f"){
  n0s <- substring(f, 3, nchar(f) - 1)

  wh <- try(n0s <- as.numeric(n0s), silent = TRUE)
  if (inherits(wh, "try-error")){
    stop("Expecting fmt input to sptintf to be of form '%.nnnX': something like '%.4f'")
  }

  n0s <- as.numeric(n0s) - 1

  zs <- paste(rep("0", n0s), collapse = "")

  paste0("<0.", zs, "1")
}

