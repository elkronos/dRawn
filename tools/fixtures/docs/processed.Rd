% Fixture for `Rscript tools/check-docs.R --self-test`: [foo()] appears only
% as R code, which is not an unprocessed link.
\name{processed}
\alias{processed}
\title{A topic whose brackets are code}
\usage{
processed(x = y[f()])
}
\arguments{
\item{x}{A value such as \code{y[f()]}.}
}
\description{
Links to \code{\link[=draw]{draw()}}.
}
\examples{
y <- 1:3
f <- function() 2
processed(y[f()])
}
