% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solve_composition_poly.R
\name{solve_composition_poly}
\alias{solve_composition_poly}
\title{Compute genome-wide breed composition}
\usage{
solve_composition_poly(
  Y,
  X,
  ped = NULL,
  groups = NULL,
  mia = FALSE,
  sire = FALSE,
  dam = FALSE,
  ploidy = 2
)
}
\arguments{
\item{Y}{numeric matrix of genotypes (columns) from all animals (rows) in population
coded as dosage of allele B {0, 1, ..., ploidy}}

\item{X}{numeric matrix of allele frequencies (rows) from each reference panel (columns). Frequencies are
relative to allele B.}

\item{ped}{data.frame giving pedigree information. Must be formatted "ID", "Sire", "Dam"}

\item{groups}{list of IDs categorized by breed/population. If specified, output will be a list
of results categorized by breed/population.}

\item{mia}{logical. Only applies if ped argument is supplied. If true, returns a data.frame
containing the inferred maternally inherited allele for each locus for each animal instead
of breed composition results.}

\item{sire}{logical. Only applies if ped argument is supplied. If true, returns a data.frame
containing sire genotypes for each locus for each animal instead of breed composition results.}

\item{dam}{logical. Only applies if ped argument is supplied. If true, returns a data.frame
containing dam genotypes for each locus for each animal instead of breed composition results.}

\item{ploidy}{integer. The ploidy level of the species (e.g., 2 for diploid, 3 for triploid, etc.).}
}
\value{
A data.frame or list of data.frames (if groups is !NULL) with breed/ancestry composition
 results
}
\description{
Computes genome-wide breed/ancestry composition using quadratic programming on a
batch of animals.
}
