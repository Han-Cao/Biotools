# R functions for genetic analysis.

countAllele <- function(genotype, mode="AlleleTable", target=NULL){
	if (mode="AlleleTable"){
		counts <- ifelse(genotype[,1]==target,1,0) + ifelse(genotype[,2]==target,1,0)
		return(counts)
	}
	if (mode="vcf"){
		genotype[genotype=="0|0"] <- 0
		genotype[genotype=="0|1" | genotype=="1|0"] <- 1
		genotype[genotype=="1|1"] <- 2
		class(genotype) <- "numeric"
		return(genotype)
	}
}

#Count haplotypes in hapList
countHap <- function(table, hapList){
	for (type in hapList){
		table[,as.character(type)] <- CountAllele(type, table[,c("hap1", "hap2")])
	}
}

#Remove zero in sample ID\
#eg. A0010 -> A10
removeZero <- function(name){
	sub("[^\d\W]+0+","\\1",name)
}

