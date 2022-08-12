library(RNOmni)

# R functions for genetic analysis.

countAllele <- function(genotype, mode="AlleleTable", target=NULL){
	if (mode=="AlleleTable"){
		counts <- ifelse(genotype[,1]==target,1,0) + ifelse(genotype[,2]==target,1,0)
		return(counts)
	}
	if (mode=="vcf"){
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
		table[,as.character(type)] <- countAllele(table[,c("hap1", "hap2")], target=type)
	}
	table[,"others"] <- 2 - rowSums(table[,hapList])
	return(table)
}

#Remove zero in sample ID\
#eg. A0010 -> A10
remove_zero <- function(x) {
  sub("^([A-Za-z]+)0+", "\\1", x)
}

#convert p value to z score
convert_p2z <- function(beta, p){
  abs_z <- abs(qnorm(p/2))
  sign(beta) * abs_z
}

rankNorm_withNA <- function(u, k = 3/8){
  df_raw <- data.frame(ID=1:length(u),origin=u)
  df_no_NA <- df_raw[complete.cases(df_raw$origin),]
  df_no_NA$rankNorm <- rankNorm(df_no_NA$origin, k)
  df_merge <- left_join(df_raw, df_no_NA, by="ID")
  return(df_merge$rankNorm)
}


# convert plink results to METASOFT input
plink2metasoft <- function(plink_files, corhots){
	master_table <- data.frame()
	# read and combine plink results
	for (i in 1:length(plink_files)) {
		file = plink_files[i]
		cohort = cohorts[i]
		cohort_table <- read.table(file, header=TRUE)
		cohort_table$cohort <- cohort
		master_table <- rbind(master_table, cohort_table)
	}
	master_table <- filter(master_table, !is.na(BETA), !is.na(SE))

	# variant infomation
	variant_info <- master_table[,c("SNP", "CHR", "BP", "A1")]
	variant_info <- filter(variant_info, !duplicated(variant_info)) %>% arrange(BP)

	# merge beta
	master_table_beta <- master_table[,c("SNP", "cohort", "BETA")]
	master_table_beta <- dcast(master_table_beta, SNP ~ cohort, value.var="BETA")
	# merge se
	master_table_se <- master_table[,c("SNP", "cohort", "SE")]
	master_table_se <- dcast(master_table_se, SNP ~ cohort, value.var="SE")

	output_table <- data.frame(SNP=master_table_beta$SNP)
	for (cohort in cohorts) {
		output_table[,paste0(cohort, "_BETA")] <- master_table_beta[, cohort]
		output_table[,paste0(cohort, "_SE")] <- master_table_se[, cohort]
	}
	return(list(variants<-variant_info, metasoft<-output_table))
}
