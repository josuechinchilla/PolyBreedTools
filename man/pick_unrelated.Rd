% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pick_unrelated.R
\name{pick_unrelated}
\alias{pick_unrelated}
\title{Pick unrelated ids from a pedigree}
\usage{
pick_unrelated(ids, ped, coeff = 0.0625)
}
\arguments{
\item{ids}{vector of ids to choose from}

\item{ped}{data.frame containing pedigree information. Should contain columns "ID" "Sire" and "Dam"
in that order.}

\item{coeff}{numeric giving maximum coefficient of relatedness among chosen animals.}
}
\description{
Given a list of ids to choose from and a pedigree, returns a subset of the list that are as
distantly related as possible.
}
