\name{summary.STNP}
\alias{summary.STNP}

\title{
Summarize a specification test \code{STNP} object
}
\description{
Prints a summary of a specification test object of class \code{STNP} with all the options used, including if options were \code{"default"}
}
\usage{
\method{summary}{STNP}(object, ...)
}

\arguments{
  \item{object}{
An object of class \code{STNP} resulting from function \code{\link{SpeTest}}
}
\item{...}{
Additional summary arguments
}
}
\value{
No return value, prints a summary of the test
}
\author{
Hippolyte Boucher <Hippolyte.Boucher@outlook.com>

Pascal Lavergne <lavergnetse@gmail.com>
}

\seealso{
\code{\link{SpeTest}} is the function which performs a specification test and records it along with all its options in an object of class \code{STNP}

\code{\link{print.STNP}} prints the specification test statistic and p-value only
}
\examples{
n <- 100
k <- 2
x <- matrix(rnorm(n*k),ncol=k)
y<-1+x\%*\%(1:k)+rnorm(n)

eq<-lm(y~x+0)


summary(SpeTest(eq=eq,type="icm",norma="np"))

}
