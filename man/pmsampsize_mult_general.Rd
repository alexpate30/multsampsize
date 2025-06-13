% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pmsampsize_mult.R
\name{pmsampsize_mult_general}
\alias{pmsampsize_mult_general}
\title{pmsampsize
\itemize{
\item Sample Size for Development of a Prediction Model
}}
\usage{
pmsampsize_mult_general(
  type,
  nagrsquared = NA,
  csrsquared = NA,
  rsquared = NA,
  parameters,
  shrinkage = 0.9,
  prevalence = NA,
  cstatistic = NA,
  seed = 123456,
  rate = NA,
  timepoint = NA,
  meanfup = NA,
  intercept = NA,
  sd = NA,
  mmoe = 1.1,
  K = NA,
  mult_n_events = NA,
  mult_rsquared_overall = NA,
  mult_nagrsquared_overall = NA
)
}
\arguments{
\item{type}{specifies the type of analysis for which sample size is being calculated
\itemize{
\item \code{"c"} specifies sample size calculation for a prediction model with a continuous outcome
\item \code{"b"} specifies sample size calculation for a prediction model with a binary outcome
\item \code{"s"} specifies sample size calculation for a prediction model with a survival (time-to-event) outcome
\item \code{"m"} specifies sample size calculation for a prediction model developed using multinomial logistic regression for a polytomous outcome
}}

\item{nagrsquared}{for \code{type="b" or type="s"} this specifies the expected value of
the Nagelkerke's R-squared of the new model, which is the Cox-Snell R-squared scaled to
lie in the 0,1 range. It is interpretable in the same way as the standard R-squared, i.e.
the percentage of variation in outcome values explained by the model.
For \code{type="m"} this specifies the expected value of
the Nagelkerke's R-squared of the pairwise distinct logistic regression models, and should be a vector.
Please read the description of \code{rsquared} for additional details about specifying the
expected R-squared performance}

\item{csrsquared}{for \code{type="b" or type="s"} this specifies the expected value of
the Cox-Snell R-squared of the new model. The Cox-Snell R-squared is the generalised
version of the well-known R-squared for continuous outcomes, based on the likelihood.
Please read the description of \code{rsquared} for additional details about specifying the
expected R-squared performance.
For \code{type="m"} this specifies the expected value of
the Cox-Snell R-squared of the pairwise distinct logistic regression models, and should be a vector.
The papers by Riley et al. (see references) outline how to obtain the Cox-Snell R-squared value
from published studies if they are not reported, using other information (such as the
C-statistic, see \code{cstatistic()} option below).}

\item{rsquared}{for \code{type="c"} this specifies the expected value of the R-squared of the new model,
where R-squared is the percentage of variation in outcome values explained by the model.
For example, the user may input the value of the R-squared reported for a
previous prediction model study in the same field.  If taking a value from a previous
prediction model development study, users should input the model's adjusted R-squared
value, not the apparent R-squared value, as the latter is optimistic (biased).  However,
if taking the R-squared value from an external validation of a previous model, the
apparent R-squared can be used (as the validation data was not used for development, and
so R-squared apparent is then unbiased).
Users should be conservative with their chosen R-squared value; for example, by taking the R-squared
value from a previous model, even if they hope their new model will improve performance.}

\item{parameters}{specifies the number of candidate predictor parameters for potential
inclusion in the new prediction model.  Note that this may be larger than the number of
candidate predictors, as categorical and continuous predictors often require two or more
parameters to be estimated.}

\item{shrinkage}{specifies the level of shrinkage desired at internal validation after
developing the new model. Shrinkage is a measure of overfitting, and can range from 0 to 1,
with higher values denoting less overfitting. We recommend a shrinkage = 0.9 (the
default in \code{pmsampsize}), which indicates that the predictor effect (beta coefficients) in
the model would need to be shrunk by 10\\% to adjust for overfitting. See references
below for further information.}

\item{prevalence}{(\code{type="b"} option) specifies the overall outcome proportion
(for a prognostic model) or
overall prevalence (for a diagnostic model) expected within the model development
dataset. This should be derived based on previous studies in the same population.}

\item{cstatistic}{(\code{type="b"} option) specifies the C-statistic reported in an
existing prediction model study to be used in conjunction with the expected
prevalence to approximate the Cox-Snell R-squared using the approach of Riley et al. 2020.
Ideally, this should be an optimism-adjusted C-statistic. The approximate Cox-Snell R-squared
value is used as described above for the \code{csrsquared()} option, and so is treated as a baseline
for the expected performance of the new model.
For \code{type="m"} this specifies the pairwise C-statistics reported in an
existing prediction model, and should be a vector.}

\item{seed}{(\code{type="b"} option) specifies the initial value of the random-number
seed used by the random-number functions when simulating data to approximate the
Cox-Snell R-squared based on reported C-statistic and expect prevalence as described
by Riley et al. 2020}

\item{rate}{(\code{type="s"} option) specifies the overall event rate in the population of interest,
for example as obtained from a previous study, for the survival outcome of interest. NB: rate must
be given in time units used for meanfup and timepoint options.}

\item{timepoint}{(\code{type="s"} option) specifies the timepoint of interest for prediction.
NB: time units must be the same as given for meanfup option (e.g. years, months).}

\item{meanfup}{(\code{type="s"} option) specifies the average (mean) follow-up time
anticipated for individuals in the model development dataset,
for example as taken from a previous study in the population of interest.
NB: time units must be the same as given for timepoint option.}

\item{intercept}{(\code{type="c"} options) specifies the average outcome value in the population of
interest e.g. the average blood pressure, or average pain score.
This could be based on a previous study, or on clinical knowledge.}

\item{sd}{(\code{type="c"} options) specifies the standard deviation (SD) of
outcome values in the population e.g.
the SD for blood pressure in patients with all other predictors set to the average.
This could again be based on a previous study, or on clinical knowledge.}

\item{mmoe}{(\code{type="c"} options) multiplicative margin of error (MMOE)
acceptable for calculation of the
intercept. The default is a MMOE of 10\\%. Confidence interval for the intercept will be
displayed in the output for reference. See references below for further information.}

\item{K}{for \code{type="m"} this specifies the number of outcome categories}

\item{mult_n_events}{for \code{type="m"} this specifies the number of expected events in each outcome category}

\item{mult_rsquared_overall}{for \code{type="m"} Rsquared for overall model}

\item{mult_nagrsquared_overall}{for \code{type="m"} Nagelkere Rsquared for overall model}
}
\description{
\code{pmsampsize} computes the minimum sample size required for the development of a new
multivariable prediction model using the criteria proposed by Riley \emph{et al}. 2018.
}
\examples{
## Random made up example

# Polytomous outcomes (multinomial logistic regression)

pmsampsize_mult_general(type = "m",
csrsquared = c(0.15,0.15,0.15),
parameters = 17,
shrinkage = 0.9,
K = 3,
mult_n_events = c(50,100,150),
mult_nagrsquared_overall = 0.15)

}
