#' Download pre-trained models from OSF
#'
#' @description
#' This function downloads the models pre-trained on TCGA pan-cancer datasets
#' for calculating scores of 47 different pathways using genomic data.
#'
#' @details
#' RDS files are downloaded to current working directory. Downloaded files
#' may optionally be removed from directory after loaded into R environment.
#'
#' @param pathways a list of pathways
#' @param keyword a character key word for searching pathways
#' @param keep_files logical value to specify removal of downloaded files
#'
#' @return a list of pretrained models
#'
#' @importFrom osfr osf_retrieve_node
#' @importFrom osfr osf_ls_files
#' @importFrom osfr osf_download
#'
#' @examples
#' my_pathways <- c("HALLMARK_DNA_REPAIR", "HALLMARK_G2M_CHECKPOINT")
#' my_models <- download_models(pathways = my_pathways)
#'
#' my_models <- download_models(keyword = "MYC")
#'
#' @export

download_models <- function(pathways = NULL, keyword = NULL,
                            keep_files = FALSE){
  name <- NULL
  data("pathway_files")

  if (is.null(pathways) == FALSE){pathways <- paste0(pathways, ".rds")}
  if (is.null(pathways) == TRUE){pathways <-
    pathway_files$name[grep(keyword, pathway_files$name)]}

  pathway_files <- subset(pathway_files, name %in% pathways)

  osf_download(pathway_files, progress = TRUE, conflicts = "skip")

  models <- c()
  for (i in seq(from = 1, to = length(pathways))){
    model <- readRDS(pathways[i])
    model[[1]]$name <- substr(pathways[i], start = 1,
                              stop = (nchar(pathways[i])-4))
    models <- c(models, model)
    if (keep_files == FALSE){file.remove(pathways[i])}
  }
  if (length(models) < length(pathways)){cat("WARNING: Not all requested
                                               models found. Use view_model()
                                               to see list of all available
                                               models.")}
  return(models)
}




