#' Harmonize LD Matrix with Summary Data
#'
#' @description This function extracts instrumental variables from the reference population, aligns them to the same effect allele as the summary dataset, and calculates the correlation between the instrumental variables.
#'
#' @param plink2_path The path to the PLINK 2 software.
#' @param bfile The binary file of the reference population.
#' @param gwas_sumdata The exposure or harmonized summary data.
#' @param thread The number of threads for PLINK. Default is 80.
#' @details The function ensures that the SNPs in the LD matrix are oriented to the same effect allele as those in the summary dataset.
#' @keywords IVs correlation matrix, LD matrix, SNP alignment
#' @importFrom snpStats read.plink ld
#' @importFrom utils read.table write.table
#' @importFrom methods as
#' @export
#' @return A matrix of estimated correlations between instrumental variants.
#'
ld_cor <- function(plink2_path, bfile, gwas_sumdata, thread = 80) {
  os_name <- Sys.info()["sysname"]
  slash <- ifelse(startsWith(os_name, "Win"), "\\", "/")

  # Read GWAS summary data
  gwas_data <- read.table(gwas_sumdata, header = TRUE)
  write.table(gwas_data[, c(1, 2)], "gwas_snps.a1", col.names = FALSE, row.names = FALSE, quote = FALSE)

  # Execute PLINK command to generate the reference genotype data
  system(paste0(plink2_path,
                " --bfile ", bfile,
                " --make-bed ",
                " --extract gwas_snps.a1",
                " --alt1-allele gwas_snps.a1",
                " --out ref_genotype",
                " --threads ", thread))

  # Load genotype data
  p_data <- read.plink("ref_genotype.bed", "ref_genotype.bim", "ref_genotype.fam")
  g_matrix <- as(p_data$genotypes, "numeric")

  # Calculate LD matrix
  d_value <- ncol(p_data$genotypes) - 1
  ld_matrix <- ld(p_data$genotypes, depth = d_value, stats = c("R"), symmetric = TRUE)
  ld_cor <- as.matrix(ld_matrix$R)
  diag(ld_cor) <- 1

  # Match and reorder LD matrix according to GWAS summary data
  gwas_snps <- gwas_data$SNP
  ld_snps <- colnames(ld_cor)
  matched_indices <- match(gwas_snps, ld_snps)
  valid_indices <- which(!is.na(matched_indices))
  gwas_data <- gwas_data[valid_indices, ]
  reordered_indices <- matched_indices[valid_indices]
  ld_cor <- ld_cor[reordered_indices, reordered_indices]

  return(ld_cor)
}
