# load both datasets for running data formatting example
data("example_cna")
data("example_snv")

# helper function for formatting snv data to gene-level format
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

# main function for formatting cna and snv data, checking for missing
# samples and merging data together for dataframe to use with models
format_data <- function(snv, cna){
  Hugo_Symbol <- NULL
  data("genes")

  if (length(setdiff(snv$Tumor_Sample_Barcode, colnames(cna))) > 0){
    print("WARNING: SampleIDs found present in SNV dataset, but not CNA.")}

  print("Formatting CNA Data...")
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

  print("Formatting SNV Data...")
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



# run testthat function to check if all features were correctly added
# to avoid errors in running randomforest models
test_that("correct number of features", {
  # check if all features for randomforest model still present in
  # formatted data after removing random data points
  snv <- snv[sample(1:nrow(snv), 0.9*nrow(snv)),]
  cna <- cna[sample(1:nrow(cna), 0.9*nrow(cna)),]

  data <- format_data(snv = snv, cna = cna)
  expect_equal(ncol(data), 44132)
})
