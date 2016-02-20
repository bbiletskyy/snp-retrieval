#Use library RJSONIO
library(RJSONIO)
#Returns creates an url
urlFomSnpId <- function(rs){
  return(paste('http://www.ncbi.nlm.nih.gov/projects/SNP/snp_gene.cgi?rs=', rs, sep=""))  
}

#Turns empty strings to NULL, returns the original string if non-empty
exist <- function(s){
  if("" == s) {
    return(NULL)
  } else {
    return(s)
  }
}

#Turns nulls and empty strings to a different missing-value-marker ("-") that can be stored in a table
na <- function(val) {
  if(is.null(val) || "" == val) {
    #return(NA)
    return("-")
  } else {
    return(val)
  }  
}

#Tries so read SNP's json representation by given url, returns NULL in case of errors
snpFromJSON <- function(url){
  return(tryCatch({
    fromJSON(url)
  }, error = function(err) {
    return(NULL)
  }))  
}

#Returns GRCh38.p2 assemly from the json snp representation
assemblyFromSnp <- function(snp){
  return(as.list(snp$assembly$GRCh38.p2[[1]]))
}

#Extracts Gene models and returns NULL if geneModels are missing in SNP's JSON
geneModelsFromAssembly <- function(assembly){
  if("geneModel" %in% names(assembly)){
    return(assembly$geneModel)
  } else {
    return(NULL)
  }
}

#Removes "rs" prefix from snp-id and returns pure integer id 
noRsPrefix <- function(rsId) {
  prefix = "rs"
  id <- tolower(rsId)
  pos <- regexpr(prefix, id)
  if(pos != -1) {
    return(substr(id, pos + nchar(prefix), nchar(id)))
  } else {
    return(as.character(id))  
  }
}

#Returns a new interer SNP-id into which the provided rsId was merged, using the provided table of merges.
#The table of merges must contain "from", "to" and "current" columns. 
merged <- function(rsId, merges = NULL) {
  id <- noRsPrefix(rsId)
  mergedRs <- merges[merges$from == id, ]
  if(nrow(mergedRs) > 0) {
    idNew <- mergedRs[1, "current"]
    if(is.na(idNew)) {
      idNew <- mergedRs[1, "to"]
    } 
    return(idNew)
  } else {
    return(rsId)
  }
}
