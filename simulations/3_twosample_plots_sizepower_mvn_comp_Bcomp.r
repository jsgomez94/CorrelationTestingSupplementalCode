
ROOT_PATH <- "/nas/longleaf/home/jsgomez/github/CorrelationTestingSupplementalCode/simulations/"
sim_folder <- "mvn_comp"
sim_type <- c("exps", "full")[2]

init <- TRUE
output <- NULL
options(max.print=1000000)
for (m in 1:5) {
    ##### FIRST STEP: GATHERING THE DATA:
    for (n in c(40,80,160,320)) {
    for (p in c(15,30,60,120,240,480,960)) {
        nsplit <- 1
        n1<-n
        n2<-n

        wd <- paste0(ROOT_PATH, sim_folder, "/txt_", sim_type, "/")
        setwd(wd)
        filelist_size = list.files(pattern = paste("*_m", m, "_n1_", n1, "_n2_", n2, "_p", p, "size", "txt", sep = ""))
        filelist_power = list.files(pattern = paste("*_m", m, "_n1_", n1, "_n2_", n2, "_p", p, "power", "txt", sep = ""))
        filelist <- c(filelist_power, filelist_size)
        print(paste("*_m", m, "_n1_", n1, "_n2_", n2, "_p", p, "txt", sep = ""))
        print(length(filelist))
        #assuming tab separated values with a header
        datalist = lapply(filelist, function(x) read.table(x, header = F))
        
        #assuming the same header/columns for all files
        datafr = do.call("rbind", datalist)
        head(datafr)
        datafr[, -(1:13)] <- (datafr[, -(1:13)] < 0.05)
        
        tem_o<-aggregate( .~V1+V2+V3+V4+V5,data=datafr,mean)
        colnames(tem_o) <- c(
            "n1_col", #1
            "n2_col", #2
            "p_col", #3
            "model_col", #4
            "sp_col", #5
            "max_sigma", #6
            "l1_sigma", #7
            "F_sigma", #8
            "l0_sigma", #9
            "max_rho", #10 
            "l1_rho", #11
            "F_rho", #12
            "l0_rho", #13
            "Li-Chen", #14
            "Cai-Liu-Xia", #15
            "Zheng-EtAl-2020", #16
            "Perm-cov1_BB", "Perm-cov2_BB","Perm-cov3_BB", "Perm-cov4_BB", "Perm-cov5_BB", "Perm-cov6_BB", # 17-22
            "Perm-cov1_B200", "Perm-cov2_B200","Perm-cov3_B200", "Perm-cov4_B200", "Perm-cov5_B200", "Perm-cov6_B200", # 23-28
            "Perm-cov1_B50", "Perm-cov2_B50","Perm-cov3_B50", "Perm-cov4_B50", "Perm-cov5_B50", "Perm-cov6_B50", # 29-34
            "Cai-Zhang", #35
            "Zheng-EtAl_sum", #36
            "Zheng-EtAl_max", #37
            "Perm-corr1_ratio_BB", "Perm-corr1_stab_BB", "Perm-corr1_dir_BB", #38,39,40
            "Perm-corr1_ratio_B200", "Perm-corr1_stab_B200", "Perm-corr1_dir_B200", #41,42,43
            "Perm-corr1_ratio_B50", "Perm-corr1_stab_B50", "Perm-corr1_dir_B50") #44,45,46
        output<-rbind(output,tem_o)    
    }
    }
}

dim(output)
output <- round(output, 2)
output_df <- as.data.frame(output)

cov_counter <- 0
cor_counter <- 0

rm(
    "datafr", "datalist", "filelist", "filelist_power", "filelist_size", "init", "m", "n",
    "n1", "n2", "nsplit", "p", "tem_o")
print(ls())


## Saving results:
## Ceating directories:
folder_new      <- paste0("results/")
if (!dir.exists(paste0(ROOT_PATH, folder_new))) {
    dir.create(paste0(ROOT_PATH, folder_new))
}
subfolder_new   <- paste0(folder_new, sim_folder, "_", sim_type, "/")
if (!dir.exists(paste0(ROOT_PATH, subfolder_new))) {
    dir.create(paste0(ROOT_PATH, subfolder_new))
}
subfolder2_new  <- paste0(subfolder_new, "Method_Comparison/")
if (!dir.exists(paste0(ROOT_PATH, subfolder2_new))) {
    dir.create(paste0(ROOT_PATH, subfolder2_new))
}

## Saving Outputs:
save.image(
    file = paste0(ROOT_PATH, subfolder_new, "results_", sim_folder, "_", sim_type,".RData"))


rm(list = ls())
print(ls())





#####################################################
#####################################################
## TRANSFORMATION COMPARISON
##
#####################################################
#####################################################

library(tidyverse)
library(magrittr)
library(gridExtra)

rm(list = ls())
cbPalette <- c("#999999", "#E69F00", "#56B4E9", 
               "#009E73", "#F0E442", "#0072B2", 
               "#D55E00", "#CC79A7")
ROOT_PATH <- "/nas/longleaf/home/jsgomez/github/CorrelationTestingSupplementalCode/simulations/"
sim_folder <- "mvn_comp"
sim_type <- c("exps", "full")[2]
load(paste0(ROOT_PATH, "results/", sim_folder,"_", sim_type, "/results_",sim_folder, "_", sim_type, ".RData"))
cov_counter <- cov_counter + 1
save.image(file = paste0(ROOT_PATH, "results/", sim_folder,"_", sim_type, "/results_", sim_folder, "_", sim_type,".RData"))
print(c(cov_counter, cor_counter))
ls()

## Saving directory
subfolder_new        <- paste0("results/", sim_folder, "_", sim_type, "/Transf_comparison/")
if (!dir.exists(paste0(ROOT_PATH, subfolder_new))) {
    dir.create(paste0(ROOT_PATH, subfolder_new))
}

dim(output_df)
colnames(output_df)

method_names <- c(
    "Perm-corr1_stab_BB",
    "Perm-corr1_dir_BB",
    "Perm-corr1_ratio_BB")
method_names_clean <-  c(
    "Perm-Cor: Variance-Stab.",
    "Perm-Cor: No Transf.",
    "Perm-Cor: Our Proposal")

output_cor <- output_df %>%
    as_tibble() %>%

    mutate(n_col = n1_col, .before = n1_col) %>%
    
    select(
        -"Cai-Zhang", -"Zheng-EtAl_sum", -"Zheng-EtAl_max",
        -"Li-Chen", -"Cai-Liu-Xia", -"Zheng-EtAl-2020",
        -starts_with("Perm-cov"),
        -"n1_col", -"n2_col",
        -max_sigma, -l1_sigma, -F_sigma, -l0_sigma, 
        -max_rho, -l1_rho, -F_rho, -l0_rho) %>%
    
    pivot_longer(
        cols      = c(starts_with("Perm-corr")),
        names_to  = "method",
        values_to = "R_Rate") %>%

    filter(
        method %in% method_names,
        p_col >= 60) %>%

    pivot_wider(
        names_from = "sp_col",
        names_prefix = "rejRate_sp",
        values_from = R_Rate) %>%

    mutate(
        Size = rejRate_sp0,
        "Power: Dense Alt." = rejRate_sp2,
        "Power: Sparse Alt." = rejRate_sp1,
        "Power Gap: Dense Alt." = rejRate_sp2 - rejRate_sp0,
        "Power Gap: Sparse Alt." = rejRate_sp1 - rejRate_sp0)


###############################
###############################
## Creating plots:

for (model_ind in 1:5) {

    file_name <- paste0(
        subfolder_new,
        "PowerCorr_model_m", model_ind, "_v", cor_counter, ".pdf")
    
    pdf(file_name, width = 9, height = 5)

    p1 <- output_cor %>%

        select(-rejRate_sp0, -rejRate_sp1, -rejRate_sp2, -"Power Gap: Dense Alt.", -"Power Gap: Sparse Alt.") %>%

        mutate(n_name = factor(paste0("n = ", n_col), levels = paste0("n = ", c(40,80,160,320)))) %>%
        mutate(method = str_replace_all(method, setNames(method_names_clean, method_names))) %>%
        mutate(method = factor(method, levels = method_names_clean)) %>%
    
        filter(model_col == model_ind) %>%

        pivot_longer(
            cols = c(Size, "Power: Dense Alt.", "Power: Sparse Alt."),
            names_to = "measure_type",
            values_to = "measure_val") %>%
        
        mutate(measure_type = factor(measure_type, levels = c("Size", "Power: Dense Alt.", "Power: Sparse Alt."))) %>%

        ggplot(aes(x = p_col, y = measure_val)) +
            geom_line(aes(col = method, linetype = method, linewidth = method), alpha = 1) +
            scale_linetype_manual(values = c(2, 3, 1)) +
            scale_discrete_manual("linewidth", values = c(0.75, 0.75, 1)) +
            geom_point(aes(col = method, shape = method), size = 2.2, alpha = 1) +
            scale_shape_manual(values = c(5, 13, 19)) +
            scale_color_manual(values = c("gray50", "gray30", "black")) +
            facet_grid(rows = vars(n_name), cols = vars(measure_type)) +
            ylab("Rejection Rate") +
            xlab("p") +
            theme(
                panel.spacing = unit(0.75, "lines"),
                legend.title = element_blank())
    print(p1)
    dev.off()

    
}


    

     
#####################################################
#####################################################
## COMPARING B-PERMUTATION SIZE RESULTS
## CORRELATION.
#####################################################
#####################################################

library(tidyverse)
library(magrittr)
library(gridExtra)

rm(list = ls())
cbPalette <- c("#999999", "#E69F00", "#56B4E9", 
               "#009E73", "#F0E442", "#0072B2", 
               "#D55E00", "#CC79A7")
ROOT_PATH <- "/nas/longleaf/home/jsgomez/github/CorrelationTestingSupplementalCode/simulations/"
sim_folder <- "mvn_comp"
sim_type <- c("exps", "full")[2]
load(paste0(ROOT_PATH, "results/", sim_folder,"_", sim_type, "/results_",sim_folder, "_", sim_type, ".RData"))
cor_counter <- cor_counter + 1
save.image(file = paste0(ROOT_PATH, "results/", sim_folder,"_", sim_type, "/results_", sim_folder, "_", sim_type,".RData"))
print(c(cov_counter, cor_counter))
ls()

## Saving directory
subfolder_new        <- paste0("results/", sim_folder, "_", sim_type, "/B_comparison/")
if (!dir.exists(paste0(ROOT_PATH, subfolder_new))) {
    dir.create(paste0(ROOT_PATH, subfolder_new))
}

method_names <- c(
    "Perm-corr1_ratio_B50",
    "Perm-corr1_ratio_B200",
    "Perm-corr1_ratio_BB"
    )  
method_names_clean <-  c(
    "Perm-Corr (B = 50)",
    "Perm-Corr (B = 200)",
    "Perm-Corr (B = 1000)") 


output_cor <- output_df %>%
    as_tibble() %>%

    mutate(n_col = n1_col, .before = n1_col) %>%
    
    select(
        -"Cai-Zhang", -"Zheng-EtAl_sum", -"Zheng-EtAl_max",
        -"Li-Chen", -"Cai-Liu-Xia", -"Zheng-EtAl-2020",
        -starts_with("Perm-cov"),
        -"n1_col", -"n2_col",
        -max_sigma, -l1_sigma, -F_sigma, -l0_sigma, 
        -max_rho, -l1_rho, -F_rho, -l0_rho) %>%
    
    pivot_longer(
        cols      = c(starts_with("Perm-corr")),
        names_to  = "method",
        values_to = "R_Rate") %>%

    filter(
        method %in% method_names,
        p_col >= 60) %>%

    pivot_wider(
        names_from = "sp_col",
        names_prefix = "rejRate_sp",
        values_from = R_Rate) %>%

    mutate(
        Size = rejRate_sp0,
        "Power: Dense Alt." = rejRate_sp2,
        "Power: Sparse Alt." = rejRate_sp1,
        "Power Gap: Dense Alt." = rejRate_sp2 - rejRate_sp0,
        "Power Gap: Sparse Alt." = rejRate_sp1 - rejRate_sp0)


###############################
###############################
## Creating plots:

for (model_ind in 1:5) {
    file_name <- paste0(
        subfolder_new,
        "PowerCorr_model_m", model_ind, "_v", cor_counter, ".pdf")
    
    pdf(file_name, width = 9, height = 5)

    p1 <- output_cor %>%

        select(-rejRate_sp0, -rejRate_sp1, -rejRate_sp2, -"Power Gap: Dense Alt.", -"Power Gap: Sparse Alt.") %>%

        mutate(n_name = factor(paste0("n = ", n_col), levels = paste0("n = ", c(40,80,160,320)))) %>%
        mutate(method = str_replace_all(method, setNames(method_names_clean, method_names))) %>%
        mutate(method = factor(method, levels = method_names_clean)) %>%
    
        filter(model_col == model_ind) %>%

        pivot_longer(
            cols = c(Size, "Power: Dense Alt.", "Power: Sparse Alt."),
            names_to = "measure_type",
            values_to = "measure_val") %>%
        
        mutate(measure_type = factor(measure_type, levels = c("Size", "Power: Dense Alt.", "Power: Sparse Alt."))) %>%

        ggplot(aes(x = p_col, y = measure_val)) +
            geom_line(aes(col = method, linetype = method, linewidth = method), alpha = 1) +
            scale_linetype_manual(values = c(2, 3, 1)) +
            scale_discrete_manual("linewidth", values = c(0.75, 0.75, 1)) +
            geom_point(aes(col = method, shape = method), size = 2.2, alpha = 1) +
            scale_shape_manual(values = c(5, 13, 19)) +
            scale_color_manual(values = c("gray50", "gray30", "black")) +
            facet_grid(rows = vars(n_name), cols = vars(measure_type)) +
            ylab("Rejection Rate") +
            xlab("p") +
            theme(
                panel.spacing = unit(0.75, "lines"),
                legend.title = element_blank())
    print(p1)
    dev.off()

    
}
     