#Read diseases and SNP id's from CSV
data = read.csv("data/autoimmune.csv")
#Make use of functions in other file
source("functions.R")

#Filter out the rows with missing SNP - the identifier of mutation
data <- data[which(data$SNP != ""),]

#Create a structure of the results table
result <- data.frame(Disease = character(), SNP = character(), Gene = character(), OCodon = character(), MCodon = character(), OAA = character(), MAA = character(), Clinical = character())

for(i in 1:nrow(data)){
  r <- data[i,]
  disease <- na(r$Disease)
  snp.id <- r$SNP
  url <- urlFomSnpId(snp.id)
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

  
  res <- data.frame(Disease=disease, SNP=snp.id, Gene=gene.name, OCodon=original.codon, MCodon=mutated.codon, OAA=original.aa, MAA=mutated.aa, Clinical=clinical)
  result <- rbind(result, res)
}

result
write.csv(result, "results.csv")



