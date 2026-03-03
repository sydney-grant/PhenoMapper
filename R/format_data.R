#' Counting frequency of genes in snv data
#'
#' @description
#' This function helps convert vcf format to readable format for
#' machine learning model
#'
#' @param gene a character variable of gene name
#' @param data a dataframe containing vcf
#'
#' @return a dataframe of reformatted snv data
#'
#' @keywords internal
#'
#' @examples
#' data("example_snv")
#' formatted_data = calc_mutations(gene = "TP53", data = snv)
#'
#' @noRd


calc_mutations <- function(gene, data){
  Hugo_Symbol <- NULL
  Tumor_Sample_Barcode <- NULL
  sub <- subset(data, Hugo_Symbol == gene)
  count <- c()
  for (sample in unique(data$Tumor_Sample_Barcode)){
    sub2 <- subset(sub, Tumor_Sample_Barcode == sample)
    if (nrow(sub2) == 0){c <- 0}
    if (nrow(sub2) > 0){c <- 1}
    count <- c(count, c)
  }
  col <- data.frame("Gene" = unlist(count))
  colnames(col) <- gene
  return(col)
}





#' Curate datasets to use with PhenoMap models
#' @description
#' This function take single nucleotide variants in vcf format with copy
#' number alterations represented with GISTIC gene-level copy number
#' threshold calls.
#'
#' @details
#' The two datasets are merged and formatted for use with
#' PhenoMap models. Sample names between two dataframes should match. If a
#' sample name appears in the cna file, but not the snv file, it will be
#' assumed this sample had no vaiant detected.
#'
#' @param snv a dataframe including single nucleotide variants in vcf format
#' containing the following columns:
#' Hugo_Symbol
#' Tumor_Sample_Barcode
#' @param cna a dataframe including copy number alterations represented
#' with GISTIC gene-level copy number threshold calls
#'
#' @return a dataframe containing merged genomic data
#'
#' @importFrom dplyr bind_cols
#' @importFrom pbapply pblapply
#' @importFrom utils data
#'
#' @examples
#' data("example_cna")
#' data("example_snv")
#' data <- format_data(snv = snv, cna = cna)
#'
#' @export


format_data <- function(snv, cna){
  Hugo_Symbol <- NULL
  data("genes")

  if (length(setdiff(snv$Tumor_Sample_Barcode, colnames(cna))) > 0){
    cat("WARNING: SampleIDs found present in SNV dataset, but not CNA.")}

  cat("Formatting CNA Data...")
  cna <- subset(cna, Hugo_Symbol %in% gene_df$Gene)

  cna_missing_genes <- setdiff(gene_df$Gene, cna$Hugo_Symbol)
  cna_missing <- data.frame(matrix(0, ncol = ncol(cna)-1,
                        nrow = length(cna_missing_genes)))
  colnames(cna_missing) <- colnames(cna)[-1]
  cna_missing$Hugo_Symbol <- cna_missing_genes
  cna <- rbind(cna, cna_missing)

  cna$Hugo_Symbol <- paste0("CNA_", cna$Hugo_Symbol)
  rownames(cna) <- cna$Hugo_Symbol
  cna <- as.data.frame(t(cna[,-1]))

  cat("Formatting SNV Data...")
  snv_count <- dplyr::bind_cols(pbapply::pblapply(gene_df$Gene, calc_mutations,
                                                  data = snv))
  colnames(snv_count) <- paste0("SNV_", colnames(snv_count))
  snv_count$PatientID <- unique(snv$Tumor_Sample_Barcode)

  cna$PatientID <- rownames(cna)

  samples_nosnv <- setdiff(cna$PatientID, snv_count$PatientID)

  if (length(samples_nosnv) > 0){
  nosnv <- as.data.frame(matrix(0, ncol = (ncol(snv_count)-1),
                                nrow = length(samples_nosnv)))
  colnames(nosnv) <- colnames(snv_count)[-ncol(snv_count)]
  nosnv$PatientID <- samples_nosnv
  snv_count <- rbind(snv_count, nosnv)
  }

  merged_data <- merge(snv_count, cna, by = "PatientID")

  rownames(merged_data) <- merged_data$PatientID
  merged_data <- merged_data[,-1]

  return(merged_data)
}





