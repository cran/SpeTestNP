\name{SpeTest_Dist}
\alias{SpeTest_Dist}

\keyword{test}
\keyword{nonparametric}
\keyword{parametric}

\concept{nonparametric test}
\concept{specification test}
\concept{linear model}
\concept{parametric model}


\title{
Nonparametric specification test statistic distribution
}

\description{
\code{SpeTest_Dist} generates a vector from the nonparametric specification test statistic distribution under the null hypothesis for one of five different tests using the bootstrap
}

\usage{
SpeTest_Dist(eq, type="icm", norma="no", boot="wild", nboot=50, para=FALSE, ker="normal",
knorm="sd", cch="default", hv="default", nbeta="default", direct="default", 
alphan="default")
}


\arguments{
\item{eq}{
A fitted model of class \code{lm} or \code{nls}
}
\item{type}{
Test type

If \code{type = "icm"} the vector is generated from the distribution of the test of Bierens (1982) under the null hypothesis (default)

If \code{type = "zheng"} the vector is generated from the distribution of the test of Zheng (1996) under the null hypothesis

If \code{type = "esca"} the vector is generated from the distribution of the test of Escanciano (2006) under the null hypothesis, significantly increases computing time

If \code{type = "pala"} the vector is generated from the distribution of the test of Lavergne and Patilea (2008) under the null hypothesis

If \code{type = "sicm"} the vector is generated from the distribution of the test of Lavergne and Patilea (2012) under the null hypothesis
}
\item{norma}{
Normalization of the test statistic

If \code{norma = "no"} the test statistic is not normalized (default)

If \code{norma = "naive"} the test statistic is normalized with a naive estimator of the variance of its components

If \code{norma = "np"} the test statistic is normalized with a nonparametric estimator of the variance of its components
}
\item{boot}{
Bootstrap type to generate the vector drawn from the distribution under the null hypothesis of the test statistic

If \code{boot = "wild"} the wild bootstrap of Wu (1986) is used

If \code{boot = "smooth"} the smooth conditional moments bootstrap of Gozalo (1997) is used
}
\item{nboot}{
Size of the vector drawn from the test statistic distribution under the null using the bootstrap, by default \code{nboot = 50}
}
\item{para}{
Parallel computing

If \code{para = FALSE} parallel computing is not used to generate the vector from the test statistic distribution under the null (default)

If \code{para = TRUE} parallel computing is used to generate the vector from the test statistic distribution under the null, significantly decreases computing time, makes use of all CPU cores except one
}
\item{ker}{
Kernel function used in the central matrix and for the nonparametric covariance estimator

If \code{ker = "normal"} the central matrix kernel function is the normal p.d.f (default)

If \code{ker = "triangle"} the central matrix kernel function is the triangular p.d.f

If \code{ker = "logistic"} the central matrix kernel function is the logistic p.d.f

If \code{ker = "sinc"} the central matrix kernel function is the sine cardinal function
}
\item{knorm}{
Normalization of the kernel function

If \code{knorm = "sd"} then the standard deviation using the kernel function equals 1 (default)

If \code{knorm ="sq"} then the integral of the squared kernel function equals 1
}
\item{cch}{
Central matrix kernel bandwidth 

If \code{type = "icm"} or \code{type = "esca"} then \code{cch} always equals \code{1}

If \code{type = "zheng"} the \code{"default"} bandwidth is the scaled
rule of thumb: \code{cch = 1.06*n^(-1/(4+k))} where \code{k} is
the number of regressors

If \code{type = "sicm"} or \code{type = "pala"} the \code{"default"}
bandwidth is the scaled rule of thumb: \code{cch = 1.06*n^(-1/5)}

The user may change the bandwidth when \code{type = "zheng"}, \code{type = "sicm"} or \code{type = "pala"}
}
\item{hv}{
If \code{norma = "np"} or \code{rejection = "bootstrap"} and \code{boot = "smooth"}, \code{hv} is the bandwidth of the nonparametric errors covariance estimator, by \code{"default"} the bandwidth is the scaled rule of thumb \code{hv = 1.06*n^(-1/(4+k))}
}
\item{nbeta}{
  If \code{type = "pala"} or \code{type = "sicm"}, \code{nbeta} is the number of "betas" in the unit hypersphere used to compute the statistic, computing time increases as \code{nbeta} gets larger

  By \code{"default"} it is equal to 20 times the square root of the number of exogenous control variables
}
\item{direct}{
  If \code{type = "pala"}, \code{direct} is the favored direction for beta, by \code{"default"} it is the OLS estimator if \code{class(eq) = "lm"}
  
  If \code{type = "sicm"}, \code{direct} is the initial direction for beta. This direction should be a vector of \code{0} (for no direction), \code{1} (for positive direction) and \code{-1} (for negative direction)

  For ex, \code{c(1,-1,0)} indicates that the user thinks that the first regressor has a positive effect on the dependent variable, that the second regressor has a negative effect on the dependent variable, and that he has no idea about the effect of the third regressor
  
  By \code{"default"} no direction is given to the hypersphere
}
\item{alphan}{
  If \code{type = "pala"}, \code{alphan} is the weight given to the favored direction for beta, by \code{"default"} it is equal to \code{log(n)*n^(-3/2)}
}
}


\details{
To generate a vector from the specification test statistic distribution under the null using the bootstrap the only argument needed is a model \code{eq} of class \code{lm} or of class \code{nls}.
But other options can and should be specified: the test statistic type \code{type}, the normalization of the test statistic \code{norma}, the bootstrap type \code{boot} and the size of the vector being generated which is equal to the number of bootstrap samples \code{nboot}, whether the vector is generated using parallel computing \code{para}, the central matrix kernel function \code{ker} and its standardization \code{knorm}, the bandwidths \code{cch} and \code{hv}. If the user has knowledge of the tests coined by Lavergne and Patilea he may choose a higher number of betas for the hypersphere (which may significantly increase computational time) and an initial "direction" to the hypersphere for the SICM test (none is given by \code{"default"}) or a starting beta for the PALA test (which is the OLS estimator by \code{"default"} if \code{class(eq) = "nls"}).

The statistic can be normalized with a naive estimator of the conditional covariance of its elements as in Zheng (1996), or with a nonparametric estimator of the conditional covariance of its elements as in in Yin, Geng, Li, Wang (2010). The vector is generated either from the wild bootstrap of Wu (1986) or from the smooth conditional moments bootstrap of Gozalo (1997).
}

\value{
\code{SpeTest_Dist} returns a vector of size \code{nboot} which is drawn from the distribution of the test statistic under the null hypothesis using the bootstrap
}

\seealso{
\code{\link{SpeTest}} is the function which performs a specification test and records it along with all its options in an object of class \code{STNP}

\code{\link{SpeTest_Stat}} is the function which only returns the specification test statistic
}


\note{
The data used to obtain the fitted model \code{eq} should not contain factors, factor variables should be transformed into dummy variables a priori

Requires the packages \code{stats}  (already installed and loaded by default in Rstudio), \code{foreach}, \code{parallel} and \code{doParallel} (if parallel computing is used to generate the vector) to be installed

For more information and to be able to use the package to its full potential see the references
}


\examples{
n <- 100
k <- 2
x <- matrix(rnorm(n*k),ncol=k)
y<-1+x\%*\%(1:k)+rnorm(n)

eq<-lm(y~x+0)

SpeTest_Dist(eq=eq,type="zheng",boot="wild",nboot=10)

eq<-nls(out~expla1*a+b*expla2+c,start=list(a=0,b=4,c=2),
data=data.frame(out=y,expla1=x[,1],expla2=x[,2]))

SpeTest_Dist(eq=eq,type="zheng",boot="wild",nboot=20)

}


\author{
Hippolyte Boucher <Hippolyte.Boucher@outlook.com>

Pascal Lavergne <lavergnetse@gmail.com>
}


\references{
H.J. Bierens (1982), "Consistent Model Specification Test",
\cite{Journal of Econometrics}, 20 (1), 105-134

J.C. Escanciano (2006), "A Consistent Diagnostic Test for Regression Models using Projections",
\cite{Economic Theory}, 22 (6), 1030-1051

P.L. Gozalo (1997), "Nonparametric Bootstrap Analysis with Applications to Demographic Effects in Demand Functions",
\cite{Journal of Econometrics}, 81 (2), 357-393

P. Lavergne and V. Patilea (2008), "Breaking the Curse of Dimensionality in Nonparametric Testing", \cite{Journal of Econometrics}, 143 (1), 103-122

P. Lavergne and V. Patilea (2012), "One for All and All for One: Regression Checks with Many Regressors",
\cite{Journal of Business and Economic Statistics}, 30 (1), 41-52


C.F.J. Wu (1986), "Jackknife, bootstrap and other resampling methods in regression analysis (with discussion)", \cite{Annals of Statistics}, 14 (4), 1261-1350

J. Yin, Z. Geng, R. Li, H. Wang (2010), "Nonparametric covariance model",
\cite{Statistica Sinica}, 20 (1), 469-479

J.X. Zheng (1996), "A Consistent Test of Functional Form via Nonparametric Estimation Techniques",
\cite{Journal of Econometrics}, 75 (2), 263-289


}
