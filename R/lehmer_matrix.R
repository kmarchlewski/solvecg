
#' Creates Lehmer matrix and its inverse
#'
#' The function calculates Lehmer matrix and its inverse according to analytical
#'     formulas.
#'
#' @param size A size of the matrix.
#'
#' @return A list of two elements: \code{$matrix} --- Lehmer matrix,
#'     \code{$inverse} --- the inversion of the matrix.
#'
#' @examples
#' lehmer <- lehmer(10)
#'
#' @export
#'
#' @useDynLib solvecg

lehmer <- function(size) {
  .Call("lehmer", size, PACKAGE = "solvecg")
}
