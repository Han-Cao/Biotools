#convert p to z

convert_p2z <- function(beta, p){
  abs_z <- abs(qnorm(p/2))
  sign(beta) * abs_z
}

file_result <- "Input/convert_p2z.input"

result_table <- read.table(file_result, header = TRUE)
result_table$z <- convert_p2z(result_table[,1], result_table[,2])
write.table(result_table, file = "Input/convert_z.out", row.names = FALSE, quote = FALSE, sep = "\t")
