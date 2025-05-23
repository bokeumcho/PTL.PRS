% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PTL_PRS.R
\name{pseudo_split}
\alias{pseudo_split}
\title{Pseudo-split target summary statistics into train and validation sets.
Optionally, if pseudo_test is TRUE, first split into test and rest,
then split the rest into train and validation.}
\usage{
pseudo_split(
  target_sumstats,
  subprop,
  ref_file_ps,
  tempfile,
  random_seed = 42,
  pseudo_test = FALSE,
  test_seed = 70,
  test_prop = NULL
)
}
\arguments{
\item{target_sumstats}{Target summary statistics to split.}

\item{subprop}{Proportion used in the second pseudo split (train vs. validation).}

\item{ref_file_ps}{Prefix of the PLINK file of the reference panel data in the target population.}

\item{tempfile}{Output path to save split summary statistics.}

\item{random_seed}{A random seed for the second split.}

\item{pseudo_test}{Logical. If TRUE, perform a two-step (three-set) split.}

\item{test_seed}{A random seed for the first split to generate test vs. rest.}

\item{test_prop}{Proportion of the full sample to assign to the test set.}
}
\value{
If pseudo_test is FALSE, returns a list of two matrices (train and validation sets).
If TRUE, returns a named list with elements "test", "train", and "validation".
}
\description{
Pseudo-split target summary statistics into train and validation sets.
Optionally, if pseudo_test is TRUE, first split into test and rest,
then split the rest into train and validation.
}
