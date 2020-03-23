\name{NMLmulti-package}
\alias{NMLmulti-package}
\alias{NMLmulti}
\docType{package}
\title{
  Compute Normalized Maximum Likelihood (NML) penalties for multinomial models
}
\description{
Computes NML penalties for models of joint multinomial data using Monte Carlo integration.
}
\details{

}
\author{
David Kellen and Karl Christoph Klauer

Maintainer: David  <davekellen@gmail.com>
}
\references{
Klauer, K. C., & Kellen, D. (2011). The flexibility of models of recognition memory: An analysis by the minimum-description length principle. Journal of Mathematical Psychology, 55, 430-450.

Klauer, K. C. & Kellen, D. (2015). The Flexibility of Models of Recognition Memory: The Case of Confidence Ratings. Journal of Mathematical Psychology, 67, 8-25.
 
}
\keyword{ package }
\seealso{
 
}
\examples{
  \dontrun{
# computes NML penalty for the SCR model
nml_scr <- run_nml(fun=NMLmulti::SCR, parl=10, ks=rep(3,5), 
            Ns=c(rep(100,4),200), fits = 1, cores=6, batchsize=5000,
            burn=10000, precision=0.1)

# computes NML penalty for the LT model
nml_lt  <- run_nml(fun=NMLmulti::LT, parl=10, ks=rep(3,5), 
            Ns=c(rep(100,4),200), fits = 1, cores=6, batchsize=5000,
            burn=10000, precision=0.1) 
  }
}
