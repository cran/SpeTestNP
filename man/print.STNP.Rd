\name{print.STNP}
\alias{print.STNP}

\title{
Print a specification test \code{STNP} object
}
\description{
Prints the test statistic and p-value of a specification test object of class \code{STNP}
}
\usage{
\method{print}{STNP}(x, ...)
}

\arguments{
  \item{x}{
An object of class \code{STNP} resulting from function \code{\link{SpeTest}}
}
\item{...}{
Additional print arguments
}
}
\value{
No return value, prints the test statistic and p-value
}
\author{
Hippolyte Boucher <Hippolyte.Boucher@outlook.com>

Pascal Lavergne <lavergnetse@gmail.com>
}

\seealso{
\code{\link{SpeTest}} is the function which performs a specification test and records it along with all its options in an object of class \code{STNP}

\code{\link{summary.STNP}} prints a summary of the specification test with all the options used
}
\examples{
n <- 100
k <- 2
x <- matrix(rnorm(n*k),ncol=k)
y<-1+x\%*\%(1:k)+rnorm(n)

eq<-lm(y~x+0)


print(SpeTest(eq=eq,type="icm",nboot=50))

}
