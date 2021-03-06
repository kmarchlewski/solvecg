% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lehmer_matrix.R
\name{lehmer}
\alias{lehmer}
\title{Creates Lehmer matrix and its inverse}
\usage{
lehmer(size)
}
\arguments{
\item{size}{A size of the matrix.}
}
\value{
A list of two elements: \code{$matrix} --- Lehmer matrix,
\code{$inverse} --- the inversion of the matrix.
}
\description{
The function calculates Lehmer matrix and its inverse according to analytical
formulas.
}
\examples{
lehmer <- lehmer(10)

}
