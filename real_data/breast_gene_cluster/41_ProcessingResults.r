

#################################################
#################################################
##
## Processing Simulations!
##
#################################################
#################################################
##
##

setwd("/nas/longleaf/home/jsgomez/github/CorrelationTesting/empirical/breast_gene_expression_v1_cluster")
result_sum <- matrix(0, ncol = 12, nrow = 4)



#################################################
## Processing preliminary results:
for(sim_ind in 1:100) {
  result_data <- read.csv(paste0("result_data/result_data_", sim_ind,".csv"))
  result_significant <- (result_data[, 3:14] <= 0.05)
  
  result_sum <- result_sum + result_significant
}


result_average <- result_sum / 100
result_average

colnames(result_average) <- paste0("n = ", seq(25,300, by = 25))
rownames(result_average) <- c(
  "Cai & Zhang (2016)",
  "Zheng et al. (2019) sum",
  "Zheng et al. (2019) max",
  "Permutation Test")
result_average

write.csv(result_average, "42_Results_average.csv")


#################################################
## Processing full results:
for(sim_ind in 1:400) {
  result_data <- read.csv(paste0("result_data/result_data_", sim_ind,".csv"))
  result_significant <- (result_data[, 3:14] <= 0.05)
  
  result_sum <- result_sum + result_significant
}


result_average <- result_sum / 400
result_average

colnames(result_average) <- paste0("n = ", seq(25,300, by = 25))
rownames(result_average) <- c(
  "Cai & Zhang (2016)",
  "Zheng et al. (2019) sum",
  "Zheng et al. (2019) max",
  "Permutation Test")
result_average

write.csv(result_average, "42_Results_average_400.csv")
