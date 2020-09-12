
#' Invert a symmetric matrix
#'
#' The function calculate an inverse of a symmetric matrix using the conjugate
#'     gradient method.
#'
#' @param A A symmetric matrix which should be inverted.
#' @param A_inv A symmetric matrix for a result of the invertion.
#'     The conjugate gradient method requires a starting point. The A_inv
#'     matrix should contain an initial guess about a result of the invertion.
#' @param precond A type of preconditioner to be used during the iteration
#'     process. Possible choices are: 'I' - identity --- no preconditioning,
#'     'J' --- Jacobi preconditioner, 'GS' Gauss-Seidel preconditioner.
#' @param itermax A maximum number of iterations.
#' @param resmax A maximal acceptable residuum.
#'
#' @return A list of two elements: \code{$iteration} --- the vector containing
#'     a number of iterations necessary to calculate each column of the inverse
#'     matrix, \code{residuum} --- the vector containing a minimal obtained
#'     residuum for each column.
#'
#' @examples
#' lehmer <- lehmer(10)
#' inverse <- matrix(0, 10, 10)
#' invert_cg(lehmer$matrix, inverse, "I", 100, 1e-9)
#' range(inverse - lehmer$inverse)
#'
#' @export
#'
#' @useDynLib solvecg

invert_cg <- function(A, A_inv, precond, itermax, resmax) {
  .Call("invert_cg", A, A_inv, precond, itermax, resmax, PACKAGE = "solvecg")
}
