% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Regressions.R
\name{RegressCrosssection}
\alias{RegressCrosssection}
\title{performs univariate crosssectional regressions for each metabolite at a given timepoint in DILT1D study
\deqn{M_i = X\Beta _i + \Epsilon _i}}
\usage{
RegressCrosssection(covariates, cSampleInfo, cMetDataLong,
  timepoint = "V0Pre", predictor = "factor(sex)")
}
\arguments{
\item{covariates}{The DILT1D covariate data frame}

\item{cSampleInfo}{The DILT1D dataframe for sample metadata}

\item{cMetDataLong}{The DILT1D dataframe for metabolite counts (or z-scores)}

\item{timepoint}{A visit with which to calculate the crosssectional regression}

\item{the}{predictor variable for the regression}
}
\value{
A list of coefficents, standard errors, tvalues and p values for the regressions
}
\description{
performs univariate crosssectional regressions for each metabolite at a given timepoint in DILT1D study
\deqn{M_i = X\Beta _i + \Epsilon _i}
}

