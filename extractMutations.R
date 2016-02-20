#Make use of functions in other file
source("functions.R")

#Loading the list of merged rs  
merges <- read.delim2("data/RsMergeArch.bcp", stringsAsFactors=FALSE, col.names = c("from","to", "build_id", "orien", "create_time", "last_updated_time", "current", "orien2Current", "last"))

#Read diseases and SNP id's from CSV
data = read.csv("data/physical+health-related-traits.csv", stringsAsFactors=FALSE)

#Filter out the rows with missing SNP - the identifier of mutation
data <- data[which(data$SNP != ""),]

#Filter merges data, so that only the need ones remain
filteredMerges <- merges[which(as.character(merges$from) %in% lapply(data$SNP, noRsPrefix)),]

#Create a structure of the results table
result <- data.frame(Disease = character(), OSNP = character(), SNP = character(), Gene = character(), OCodon = character(), MCodon = character(), OAA = character(), MAA = character(), Clinical = character())
#Filling in the resulting table
for(i in 1:nrow(data)){
  r <- data[i,]
  disease <- na(r$Disease)
  snpId <- as.character(merged(r$SNP, filteredMerges))
  url <- urlFomSnpId(snpId)
  snp <- snpFromJSON(url)
  clinical <- na(snp$is_clinical)
  assembly <- assemblyFromSnp(snp)
  gene.models <- geneModelsFromAssembly(assembly)
  gene.model <- gene.models[[1]]
  gene.name <- na(gene.model$geneSymbol)
  allele.contig <- as.list(gene.model$contig_allele)
  allele.var <- as.list(gene.model$variation_allele[[1]])
  original.codon <- na(allele.contig$codon)
  original.aa <- na(allele.contig$aaAbbrev)
  mutated.codon <- na(allele.var$codon)
  mutated.aa <- na(allele.var$aaAbbrev)

  res <- data.frame(Disease=disease, OSNP = r$SNP, SNP=snpId, Gene=gene.name, OCodon=original.codon, MCodon=mutated.codon, OAA=original.aa, MAA=mutated.aa, Clinical=clinical)
  result <- rbind(result, res)
}

result
write.csv(result, "results/results_physical+health-related-traits.csv")



