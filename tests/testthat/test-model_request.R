# function for downloading requested models from osf repository
download_models <- function(pathways = NULL, keyword = NULL,
                            keep_files = FALSE){
  name <- NULL
  data("pathway_files")

  if (is.null(pathways) == FALSE){pathways <- paste0(pathways, ".rds")}
  if (is.null(pathways) == TRUE){
    pathways <- pathway_files$name[grep(keyword, pathway_files$name)]}

  pathway_files <- subset(pathway_files, name %in% pathways)

  osf_download(pathway_files, progress = TRUE)

  models <- c()
  for (i in seq(from = 1, to = length(pathways))){
    model <- readRDS(pathways[i])
    model[[1]]$name <- substr(pathways[i], start = 1,
                              stop = (nchar(pathways[i])-4))
    models <- c(models, model)
    if (keep_files == FALSE){file.remove(pathways[i])}
  }
  return(models)
}

view_models <- function(keyword = "ALL"){
  data("pathway_descriptions")
  if (keyword != "ALL"){
    pathway_descriptions <-
      pathway_descriptions[grep(keyword, pathway_descriptions$Pathway),]}
  return(pathway_descriptions)
}

test_that("requested models available", {
  # open list of available models and randomly select 10 to download from server
  available_models <- view_models()
  requested_models <- available_models$Pathway[sample(1:nrow(available_models),
                                                      10)]
  returned_models <- download_models(pathways = requested_models)
  expect_equal(length(returned_models), length(requested_models))
})
