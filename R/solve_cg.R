
#' Solves a linear and symmetric system of equations
#'
#' The function calculates a solution of the Ax = b system using the conjugate
#'     gradient method.
#'
#' @param A A symmetric matrix.
#' @param x A solution vector. The conjugate gradient method requires
#'     a starting point. The x vector should contain an initial guess about
#'     a result.
#' @param b A vector of right hand side values.
#' @param precond A type of preconditioner to be used during the iteration
#'     process. Possible choices are: 'I' - identity --- no preconditioning,
#'     'J' --- Jacobi preconditioner, 'GS' Gauss-Seidel preconditioner.
#' @param itermax A maximum number of iterations.
#' @param resmax A maximal acceptable residuum.
#'
#' @return A list of two elements: \code{$iteration} --- a value of iterations
#'     necessary to calculate a solution, \code{residuum} --- a value containing
#'     a minimal obtained residuum.
#'
#' @examples
#' lehmer <- lehmer(10)
#' x <- rep(0, 10)
#' b <- lehmer$matrix %*% seq(1, 10)
#' solve_cg(lehmer$matrix, x, b, "I", 100, 1e-9)
#' range(seq(1, 10) - x)
#'
#' @export
#'
#' @useDynLib solvecg

solve_cg <- function(A, x, b, precond, itermax, resmax) {
  .Call("solve_cg", A, x, b, precond, itermax, resmax, PACKAGE = "solvecg")
}
