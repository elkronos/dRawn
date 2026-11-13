% Fixture for `Rscript tools/check-docs.R --self-test`: an Rd file written
% without markdown processing, whose only [foo()] link is in \details.
\name{unprocessed}
\alias{unprocessed}
\title{A topic written without markdown}
\usage{
unprocessed(x)
}
\arguments{
\item{x}{A value.}
}
\description{
Returns \code{x}.
}
\details{
The first paragraph is plain.

The second paragraph refers to [draw()] as markdown would.
}
\examples{
unprocessed(1)
}
