#' Compute genome-wide breed composition
#' 
#' Computes genome-wide breed/ancestry composition using quadratic programming on a
#' batch of animals.
#' 
#' @param Y numeric matrix of genotypes (columns) from all animals (rows) in population
#'  coded as dosage of allele B {0, 1, ..., ploidy}
#' @param X numeric matrix of allele frequencies (rows) from each reference panel (columns). Frequencies are
#'  relative to allele B.
#' @param ped data.frame giving pedigree information. Must be formatted "ID", "Sire", "Dam"
#' @param groups list of IDs categorized by breed/population. If specified, output will be a list
#'  of results categorized by breed/population.
#' @param mia logical. Only applies if ped argument is supplied. If true, returns a data.frame
#'  containing the inferred maternally inherited allele for each locus for each animal instead
#'  of breed composition results.
#' @param sire logical. Only applies if ped argument is supplied. If true, returns a data.frame
#'  containing sire genotypes for each locus for each animal instead of breed composition results.
#' @param dam logical. Only applies if ped argument is supplied. If true, returns a data.frame
#'  containing dam genotypes for each locus for each animal instead of breed composition results.
#' @param ploidy integer. The ploidy level of the species (e.g., 2 for diploid, 3 for triploid, etc.).
#' @return A data.frame or list of data.frames (if groups is !NULL) with breed/ancestry composition
#'  results
#' @import quadprog
#' @export
solve_composition_poly <- function(Y,
                              X,
                              ped = NULL,
                              groups = NULL,
                              mia = FALSE,
                              sire = FALSE,
                              dam = FALSE,
                              ploidy = 2) {

  # Functions require Y to be animals x SNPs. Transpose
  Y <- t(Y)
  
  # SNPs in Y should only be the ones present in X
  Y <- Y[rownames(Y) %in% rownames(X), ]
  
  # If ped is supplied, use QPsolve_par to compute genomic composition using
  #   only animals who have genotyped parents (by incorporating Sire genotype).
  if (!is.null(ped)) {
    mat_results <- lapply(colnames(Y),
                          QPsolve_par,
                          Y,
                          X,
                          ped,
                          mia = mia,
                          sire = sire,
                          dam = dam)
    
    mat_results_tab <- do.call(rbind, mat_results)
    return (mat_results_tab)
    
    # Else if groups supplied - perform regular genomic computation
    #   and list results by groups
  } else if (!is.null(groups)) {
    
    # When using regular genomic computation - adjust dosage based on ploidy
    Y <- Y / ploidy
    
    grouped_results <- lapply(groups, QPseparate, Y, X)
    return (grouped_results)
    
    # If neither using the ped or grouping option - just perform normal, unsegregated
    #   calculation
  } else {
    
    # Adjust dosage based on ploidy
    Y <- Y / ploidy
    
    results <- t(apply(Y, 2, QPsolve, X))
    return (results)
  }
}



