# clear environment
rm(list = ls())

# goto working directory
setwd("~/Documents/GitHub/global_quality_assessment/")

# load user functions
source("library.R")

main_manuscript <- function(){
  # loop over analysis
  for (correct in c(FALSE, TRUE)){
    # read in pair info used for five (5) tests of the ability of chemical shifts assess the global quality of NMR structures
    pairs <- read.table("data/tests_info", col.names = c("pair", "pair_name", "ref", "threshold"))
    pair_names <- c("1R2P_2LPS","2FRL_2M22","2H2X_2M21","2KFC_2L1V","2N6Q_5KMZ")
    
    # initialize data 
    initialize_analysis(pairs, correct_shifts = correct)
    for (wt in 1:2){
      # variables for outlier analysis
      min_nslrs <- NULL
      thresholds <- c(1:10, "none")
      # Errors determined using average structure
      sm_average_mean <- summarize_tables("mean","mae", TRUE, FALSE, FALSE, names=pair_names, weight = wt )
      sm_average_larmord <- summarize_tables("larmord","mae", TRUE, FALSE, FALSE, names=pair_names, weight = wt)
      sm_average_ramsey <- summarize_tables("ramsey","mae", TRUE, FALSE, FALSE, names=pair_names, weight = wt)
      
      # Errors determined using conformatinally-averaged chemical shifts
      sm_conf_mean <- summarize_tables("mean","mae", FALSE, FALSE, TRUE, names=pair_names, weight = wt)
      sm_conf_larmord <- summarize_tables("larmord","mae", FALSE, FALSE, TRUE, names=pair_names, weight = wt)
      sm_conf_ramsey <- summarize_tables("ramsey","mae", FALSE, FALSE, TRUE, names=pair_names, weight = wt)
      
      # Errors determined for each model in so-called "composite" ensembles
      sm_all_mean <- summarize_tables("mean","mae", FALSE, FALSE, FALSE, names=pair_names, weight = wt)
      sm_all_larmord <- summarize_tables("larmord","mae", FALSE, FALSE, FALSE, names=pair_names, weight = wt)
      sm_all_ramsey <- summarize_tables("ramsey","mae", FALSE, FALSE, FALSE, names=pair_names, weight = wt)
      
      # Check impact of outliers
      for (t in seq_along(thresholds)){
        for (predictor in c("mean","larmord","ramsey")){
          ifelse(thresholds[t]=="none", threshold <- thresholds[t], threshold <- as.numeric(thresholds[t]))
          m <- summarize_tables(predictor,"mae", FALSE, FALSE, FALSE, names=pair_names, weight = wt, outliers = threshold)
          min_nslrs <- c(min_nslrs, mean(m[[1]]))
        }
      }
      
      min_nslrs <- matrix(min_nslrs, ncol=3, byrow = TRUE)
      rownames(min_nslrs) <- thresholds
      colnames(min_nslrs) <- c("mean", "larmord", "ramsey")
      min_nslrs <- cbind(min_nslrs,average=rowMeans(min_nslrs))
      

      #  **** Begin supporting information manuscript analysis ****
      # Part 1:
      # True-positive rates analysis
      # (1) Consensus
      # Error determined using chemical shifts computed from the average NMR structures
      get_flags_only(sm_average_mean)
      # Errors determined using conformatinally-averaged chemical shifts
      get_flags_only(sm_conf_mean)
      # Errors determined for each model in so-called "composite" ensembles
      get_flags_only(sm_all_mean)
      
      # (2) RAMSEY
      # Error determined using chemical shifts computed from the average NMR structures
      get_flags_only(sm_average_ramsey)
      # Errors determined using conformatinally-averaged chemical shifts
      get_flags_only(sm_conf_ramsey)
      # Errors determined for each model in so-called "composite" ensembles
      get_flags_only(sm_all_ramsey)
      
      # (3) LARMORD
      # Error determined using chemical shifts computed from the average NMR structures
      get_flags_only(sm_average_larmord)
      # Errors determined using conformatinally-averaged chemical shifts
      get_flags_only(sm_conf_larmord)
      # Errors determined for each model in so-called "composite" ensembles
      get_flags_only(sm_all_larmord)
      
      # Part 2:
      # Resolving power analysis (NSLR)
      # (1) Consensus
      m <- get_nslr_only(sm_all_mean)
      make_nslr_plots(m,  gsub("_",":",pair_names),figfile=paste("figures/nslr_mean_wt_", wt, "_cor_", as.character(correct),".pdf", sep = ""))
      # (2) RAMSEY
      m <- get_nslr_only(sm_all_ramsey)
      make_nslr_plots(m,  gsub("_",":",pair_names),figfile=paste("figures/nslr_ramsey_wt_", wt, "_cor_", as.character(correct),".pdf", sep = ""))
      # (3) LARMORD
      m <- get_nslr_only(sm_all_larmord)
      make_nslr_plots(m,  gsub("_",":",pair_names),figfile=paste("figures/nslr_larmord_wt_", wt, "_cor_", as.character(correct),".pdf", sep = ""))
      
      # Part 3:
      # Comparing NMR solution and X-ray crystal structures
      # (1) Consensus
      sm_xray_mean <- summarize_tables("mean","mae", FALSE, TRUE, FALSE, c("1SCL_430D","2LPS_1KXK","28SR_1LNT"))
      # (2) RAMSEY
      sm_xray_ramsey <- summarize_tables("ramsey","mae", FALSE, TRUE, FALSE, c("1SCL_430D","2LPS_1KXK","28SR_1LNT"))
      # (3) LARMORD
      sm_xray_larmord <- summarize_tables("larmord","mae", FALSE, TRUE, FALSE, c("1SCL_430D","2LPS_1KXK","28SR_1LNT"))
      #  **** End supporting information manuscript analysis ****
      save(list=ls(), file= paste("workspaces/workspace_wt_", wt, "_correct_", as.character(correct), ".RData", sep = ""))
    }
  }
}

main_manuscript()