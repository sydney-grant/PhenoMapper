#' Calculate GSVA-based pathway scores with genomic data
#'
#' @description
#' This function uses our models pretrained on TCGA pan-cancer datasets
#' to calculate single sample pathway-level scores based on single nucleotide
#' variant and copy number alteration data.
#'
#'
#' @param models a list of models downloaded from OSF
#' @param data a dataframe of formatted snv and cna datasets
#'
#' @return a dataframe containing single sample scores for each pathway
#'
#' @import randomForest
#' @importFrom stats predict
#' @importFrom utils data
#'
#' @examples
#' data("example_cna")
#' data("example_snv")
#'
#' data <- format_data(snv = snv, cna = cna)
#'
#' my_pathways <- c("HALLMARK_DNA_REPAIR", "HALLMARK_G2M_CHECKPOINT")
#' my_models <- download_models(pathways = my_pathways, keep_files = FALSE)
#'
#' my_scores <- predict_scores(data = data, models = my_models)
#'
#' @export



predict_scores <- function(data, models){
  scores <- data.frame("SampleID" = rownames(data))
  for (i in seq(from = 1, to = length(models))){
    model <- models[i]
    cat(paste0("Predicting model ", i, "/", length(models), ": ",
                 model[[1]]$name))

    missing_inputs <- setdiff(row.names(model[[1]]$importance), colnames(data))
    missing_data <- data.frame(matrix(0, ncol = length(missing_inputs),
                                      nrow = nrow(data)))
    colnames(missing_data) <- missing_inputs
    data <- cbind(data, missing_data)

    preds <- predict(model[[1]], newdata = data)
    preds <- data.frame("Scores" = preds)
    colnames(preds) <- model[[1]]$name
    scores <- cbind(scores, preds)
  }
  return(scores)
}



