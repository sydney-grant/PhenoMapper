#' View pre-trained models available on OSF
#'
#' @description
#' This function allows users to browse available models for download
#'
#' @param keyword a character key word for searching pathways
#'
#' @return a list of model names
#'
#' @examples
#' all_pathways <- view_models()
#' kegg_pathways <- view_models(keyword = "KEGG")
#'
#' @export

view_models <- function(keyword = "ALL"){
  data("pathway_descriptions")
  if (keyword != "ALL"){
    pathway_descriptions <-
      pathway_descriptions[grep(keyword, pathway_descriptions$Pathway),]}
  return(pathway_descriptions)
}


