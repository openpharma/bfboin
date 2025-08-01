% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bf_boin_oc.R
\name{get.oc.bf}
\alias{get.oc.bf}
\title{Get Operating Characteristics for the BF-BOIN Design}
\usage{
get.oc.bf(
  ntrial = 1000,
  seed = 3262,
  target = 0.25,
  p.true = c(0.1, 0.3, 0.5),
  ncohort = 10,
  cohortsize = 3,
  n.earlystop = 100,
  startdose = 1,
  titration = FALSE,
  p.saf = 0.6 * target,
  p.tox = 1.4 * target,
  cutoff.eli = 0.95,
  extrasafe = FALSE,
  offset = 0.05,
  boundMTD = FALSE,
  n.cap = 12,
  end.backfill = TRUE,
  n.per.month = 3,
  dlt.window = 1,
  p.response.true = c(1, 1, 1),
  three.plus.three = FALSE,
  accrual = "uniform",
  backfill.assign = "highest"
)
}
\arguments{
\item{ntrial}{the total number of trials to be simulated}

\item{seed}{the random number seed for simulation}

\item{target}{the target DLT rate}

\item{p.true}{a vector containing the true toxicity probabilities of the
investigational dose levels.}

\item{ncohort}{the total number of cohorts}

\item{cohortsize}{the cohort size}

\item{n.earlystop}{the early stopping parameter. If the number of patients
treated at the current dose reaches \code{n.earlystop},
stop the trial and select the MTD based on the observed data.
The default value \code{n.earlystop=100} essentially turns
off this type of early stopping.}

\item{startdose}{the starting dose level for the trial}

\item{titration}{set \code{titration=TRUE} to perform dose escalation with cohort size = 1 to accelerate dose escalation at the begining of the trial.}

\item{p.saf}{the highest toxicity probability that is deemed subtherapeutic
(i.e. below the MTD) such that dose escalation should be undertaken.
The default value is \code{p.saf=0.6*target}.}

\item{p.tox}{the lowest toxicity probability that is deemed overly toxic such
that deescalation is required. The default value is
\code{p.tox=1.4*target}).}

\item{cutoff.eli}{the cutoff to eliminate an overly toxic dose for safety.
We recommend the default value of (\code{cutoff.eli=0.95}) for general use.}

\item{extrasafe}{set \code{extrasafe=TRUE} to impose a more stringent stopping rule}

\item{offset}{a small positive number (between \code{0} and \code{0.5}) to control how strict the
stopping rule is when \code{extrasafe=TRUE}. A larger value leads to a more
strict stopping rule. The default value \code{offset=0.05} generally works well.}

\item{boundMTD}{set \code{boundMTD=TRUE} to impose the condition: the isotonic estimate of toxicity probability
for the selected MTD must be less than de-escalation boundary.}

\item{n.cap}{permanently close a dose for backfilling if the number of patients assigned
to the dose reaches \code{n.cap}}

\item{end.backfill}{when the dose escalation ends, the backfilling by definition also ends. Default is TRUE.}

\item{n.per.month}{patient accrual rate per month}

\item{dlt.window}{DLT assessment window (months)}

\item{p.response.true}{a vector containing the true response probabilities of the
investigational dose levels}

\item{three.plus.three}{modify the decision from de-escalation to stay when observing
1 DLT out of 3 patients}

\item{accrual}{"uniform" or "poisson", according to whether accrual distribution is uniform
(consistent with Shiny App) or a Poisson process (consistent with publication)}

\item{backfill.assign}{How to assign backfill dose given the open backfill doses. Options are
"highest" (default), "lowest", or "random".}
}
\value{
\code{get.oc.bf()} returns the operating characteristics of the BOIN design as a list,
       including:
       (1) selection percentage at each dose level (\code{$selpercent}),
       (2) the average number of patients treated at each dose level (\code{$npatients}),
       (3) the percentage of patients treated at each dose level on average (\code{$percentpatients}),
       (4) the average number of toxicities observed at each dose level (\code{$ntox}),
       (5) the average number of toxicities in total (\code{$totaltox}),
       (6) the average number of patients in total(\code{$totaln}),
       (7) the percentage of early stopping without selecting the MTD (\code{$percentstop}),
       (8) the average duration of the trial (\code{duration}).
}
\description{
Get Operating Characteristics for the BF-BOIN Design
}
\examples{

get.oc.bf(ntrial = 1000,
          seed = 9,
          target = 0.25,
          p.true = c(0.1, 0.5),
          ncohort = 10,
          cohortsize = 3,
          n.earlystop = 9,
          startdose = 1,
          titration = FALSE,
          cutoff.eli = 0.95,
          extrasafe = TRUE,
          offset = 0.1,
          boundMTD=FALSE,
          n.cap = 12,
          end.backfill = TRUE,
          n.per.month = 1,
          dlt.window = 1,
          p.response.true = c(0.001, 0.001))




}
\references{
Zhao Y, Yuan Y, Korn EL, Freidlin B. Backfilling patients in phase I dose-escalation trials using Bayesian optimal interval design (BOIN). Clinical Cancer Research. 2024 Feb 16;30(4):673-9.
}
\seealso{
Shiny app: \url{https://biostatistics.mdanderson.org/shinyapps/BF-BOIN/}
}
