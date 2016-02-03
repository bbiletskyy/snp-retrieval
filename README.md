# snp-retrieval
## Description
R-scripts for obtaining a SNP details given the SNP id
## Motivation
Risk of suffering from specific diseases are related with specific genetic mutations. One special class of genetic mutations is called SNP - single-nucleotide polymorphism (https://en.wikipedia.org/wiki/Single-nucleotide_polymorphism). Which is the variation of a single nucleotide.
In order to understand the mechanism of a SNP's effect one usually needs to find out which genes were affected by SNP. Therefore it is necessary to locate a chromosome, a gene, a codon and a resulting amino-acid before and after the mutation. 
## How to use
1. Import R files to some R developer environment (i.e. R-Studio)
2. Put your input data in cvs format into the _/data_ folder
3. Execute _extractMutations.R_ script
4. Check the results in _results.csv_ file

### Remarks
The input CSV file should contain at least 2 columns: _Disease_ and _SNP_, for the disease names and SNP id's respectively.
There is a nuber of example csv files in the _/data_ folder with the Disease-SNP mappings obtained from http://www.eupedia.com/genetics/medical_dna_test.shtml
