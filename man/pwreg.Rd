% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{pwreg}
\alias{pwreg}
\title{Fit priority-adjusted proportional win-fractions (PW) regression model}
\usage{
pwreg(
  ID,
  time,
  status,
  Z,
  rho = 0,
  strata = NULL,
  fixedL = TRUE,
  eps = 1e-04,
  maxiter = 50
)
}
\arguments{
\item{ID}{a vector of unique subject-level identifiers.}

\item{time}{a vector of all the event times.}

\item{status}{a vector of the status for all the event. 0: censoring, 1:death
and 2: non-fatal event.}

\item{Z}{a matrix or a vector of covariates.}

\item{rho}{a non-negative number as the power of the survival function used
in the weight. Default (\code{rho=0}) is recommended. If there is a `strata` argument,
then `rho` is ignored}

\item{strata}{a vector of strata. `strata` needs to be specified when fit a stratified PW model.}

\item{fixedL}{logical variable indicating which variance estimator to be used. If `TRUE`,
the variance estimator for small strata is going to used.}

\item{eps}{precision for the convergence of Newton-Raphson algorithm.}

\item{maxiter}{maximum number of iterations allow for the Newton-Raphson
algorithm.}
}
\value{
An object of class \code{pwreg} with the following components.
\code{beta}:a vector of estimated regression coefficients. \code{Var}:estimated
covariance matrix for \code{beta}. \code{conv:} boolean variable indicating
whether the algorithm converged within the maximum number of iterations.
}
\description{
Fit priority-adjusted proportional win-fractions (PW) regression model.
}
\examples{
library(WR)
head(gbc)
id_unique <-unique(gbc$id)

# Randomly sample 200 subjects from gbc data
set.seed(2021)
id_sample <- sample(id_unique, 200)
gbc_reduce <- gbc[gbc$id \%in\% id_sample, ]

# Use the reduced gbc data for PW analysis
nr <- nrow(gbc_reduce)
p <- ncol(gbc_reduce)-3
ID <- gbc_reduce[,"id"]
time <- gbc_reduce[,"time"]
status <- gbc_reduce[,"status"]
Z <- as.matrix(gbc_reduce[,4:(3+p)],nr,p)
pwreg.obj <- pwreg(time=time,status=status,Z=Z,ID=ID)
print(pwreg.obj)

# Fit a stratified PW model
age_group <- cut(gbc_reduce$age, breaks = c(0, 35, 45, 55, 65, Inf), right = FALSE)
gbc_st <- gbc_reduce[,-5]
strata <- age_group
nr     <- nrow(gbc_st)
p      <- ncol(gbc_st)-3
ID     <- gbc_st[, "id"]
time   <- gbc_st[, "time"]
status <- gbc_st[, "status"]
Z      <- as.matrix(gbc_st[,4:(3+p)],nr,p)
st.pwreg.obj <- pwreg(time=time,status=status,Z=Z,ID=ID,strata=strata,fixedL=TRUE)
print(st.pwreg.obj)
}
\references{
Mao, L. and Wang, T. (2020). "A class of proportional win-fractions
regression models for composite outcomes". Biometrics, 10.1111/biom.13382

Wang, T. and Mao, L. (2021+). "Stratified Proportional Win-fractions
Regression Analysis".
}
\seealso{
\code{\link{score.proc}}
}
