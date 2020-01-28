#' synonym2common function.
#'
#' This function converts drug synonyms into common name using DrugBank db vocabulary
#' @param drug A character vector containing the drugs to be converted
#' @param manualSearch A boolean determining if you would like the function to promt manual searches for unconverted drugs
#' @return A dataframe containing the original drug name and converted names
#' @export
#' @examples
#' # non-common drug synonyms
#'   drugs <- c("E.P.O.", "TNK-tPA", "Imiglucerasa")
#'   converted <- synonym2common(drugs)

synonym2common <- function(drugs,manualSearch=T){
  drugs <- base::as.character(drugs)
  drugs <- data.frame(original=drugs,name=drugs,stringsAsFactors = F)
  drugs$name <- stringr::str_remove_all(string = drugs$name,pattern = "-")
  drugs$name <- base::toupper(stringr::str_remove_all(string = drugs$name,pattern = " "))
  drugs$dup <- duplicated(drugs$name)
  drugs_in <- drugs$original[drugs$dup==FALSE]
  `%>%` <- dplyr::`%>%`
  #load DrugBank DB vocabulary
  data("dbvoc")
  #find genes that are already in Official Gene Symbol format
  converted <- base::data.frame(
    original=base::sort(drugs_in[base::tolower(drugs_in) %in% base::tolower(dbvoc$Common.name)]),
    converted=base::sort(drugs_in[base::tolower(drugs_in) %in% base::tolower(dbvoc$Common.name)]),
    stringsAsFactors = F
    
  )
  converted <- converted %>%
    dplyr::filter(!base::duplicated(converted))
  
  #try to find aliases in the "Synonyms" column
  to_convert <- base::data.frame(
    original=drugs_in[!base::tolower(drugs_in) %in% base::tolower(dbvoc$Common.name)],
    stringsAsFactors = F
  )
  to_convert$original <- base::unique(to_convert$original)
  
  dbvoc$is_present <- ""
  dbvoc$syn <- ""
  for (i in 1:base::length(dbvoc$Synonyms)){
    synonyms <- dbvoc$Synonyms[i]
    if(synonyms == ""){
      dbvoc$is_present[i] <- FALSE
      dbvoc$syn[i] <- "-"
    }else{
      synonyms <- base::unlist(base::strsplit(x = synonyms,
                                             split = "|",
                                             fixed = T))
      dbvoc$is_present[i] <- base::any(base::tolower(to_convert$original) %in% base::tolower(synonyms))
      if(base::length(to_convert$original[base::tolower(to_convert$original) %in% base::tolower(synonyms)])==0){
        dbvoc$syn[i] <- "-"
      }else{
        dbvoc$syn[i] <- to_convert$original[base::tolower(to_convert$original) %in% base::tolower(synonyms)]
      }  
    }
  }
  
  selected_gene_info <- dbvoc %>%
    dplyr::filter(is_present==T) %>%
    dplyr::select(syn,Common.name) %>%
    dplyr::rename(original=syn, converted=Common.name)
  
  to_convert <- base::merge(to_convert,selected_gene_info,by = "original",all.x = T)
  
  #update converted dataframe
  converted <- base::as.data.frame(base::rbind(converted,to_convert[!base::is.na(to_convert$converted),]))
  
  #keep only non-converted gene aliases
  to_convert <- to_convert %>%
    dplyr::filter(is.na(converted)) %>%
    dplyr::select(original)
  
  to_convert$merge_name <- stringr::str_remove_all(string = to_convert$original,pattern = "-")
  to_convert$merge_name <- base::toupper(stringr::str_remove_all(string = to_convert$merge_name,pattern = " "))
  
  #promt user to search mannually the remaining drugs online
  if (manualSearch==T & length(to_convert$original>0)){
    for (i in 1:base::length(to_convert$original)){
      base::message(to_convert$original[i])
      to_convert$converted[i] <- base::readline(prompt = "Type Common drug name (or NA if non-existing): ")
    }
    converted <- base::as.data.frame(base::rbind(converted,to_convert[!base::is.na(to_convert$converted),]))
  }
  rownames(converted) <- 1:length(converted$original)
  #find original aliases that have conflicted Official Gene Symbols and chose one to keep
  dups <- converted[base::duplicated(converted$original),]
  if(length(dups$original)>0){
    for (i in 1:length(base::unique(dups$original))){
      original <- base::unique(dups$original)[i]
      lines <- which(converted$original==original)
      base::print(converted[lines,])
      selected_line <- base::as.numeric(base::readline(prompt = "Type line number of selected Official Gene Symbol: "))
      converted <- converted[-base::setdiff(lines,selected_line),]
      rownames(converted) <- 1:length(converted$original)
    }
  }else{
    base::return(converted)
  }  
}