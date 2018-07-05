# clear environment
rm(list = ls())

# set seed
set.seed(12345)

# goto working directory
setwd("~/Documents/GitHub/global_quality_assessment/")

# load user functions
source("library.R")


# loop over analysis
for (correct in c(TRUE)){
  # read in pair info used for five (5) tests of the ability of chemical shifts assess the global quality of NMR structures
  pairs <- read.table("data/tests_info", col.names = c("pair", "pair_name", "ref", "threshold"))
  pair_names <- c("1R2P_2LPS","2FRL_2M22","2H2X_2M21","2KFC_2L1V","2N6Q_5KMZ")
  
  # initialize data 
  initialize_analysis(pairs, correct_shifts = correct)
  for (wt in 2:2){
    # variables for outlier analysis
    min_nslrs <- NULL
    #thresholds <- c(0.5, 1:6, "none")
    thresholds <- c("none")
    # Errors determined using average structure
    sm_average_larmord <- summarize_tables("larmord","mae", TRUE, FALSE, FALSE, names=pair_names, weight = wt)

    # Errors determined using conformatinally-averaged chemical shifts
    sm_conf_larmord <- summarize_tables("larmord","mae", FALSE, FALSE, TRUE, names=pair_names, weight = wt)

    # Errors determined for each model in so-called "composite" ensembles
    sm_all_larmord <- summarize_tables("larmord","mae", FALSE, FALSE, FALSE, names=pair_names, weight = wt)

    # Check impact of outliers
    for (t in seq_along(thresholds)){
      for (predictor in c("larmord")){
        ifelse(thresholds[t]=="none", threshold <- thresholds[t], threshold <- as.numeric(thresholds[t]))
        m <- summarize_tables(predictor,"mae", FALSE, FALSE, FALSE, names=pair_names, weight = wt, outliers = threshold)
        min_nslrs <- c(min_nslrs, colMeans(m[[1]]))
      }
    }
    
    min_nslrs <- matrix(round(min_nslrs, 2), ncol=3, byrow = TRUE)

    rownames(min_nslrs) <- thresholds
    colnames(min_nslrs) <- c("proton", "carbon", "both")

    
    #  **** Begin supporting information manuscript analysis ****
    # Part 1:
    # True-positive rates analysis

    # (3) LARMORD
    # Error determined using chemical shifts computed from the average NMR structures
    get_flags_only(sm_average_larmord)
    # Errors determined using conformatinally-averaged chemical shifts
    get_flags_only(sm_conf_larmord)
    # Errors determined for each model in so-called "composite" ensembles
    get_flags_only(sm_all_larmord)
    
    # Part 2:
    # Resolving power analysis (NSLR)
    # (3) LARMORD
    m <- get_nslr_only(sm_all_larmord)
    make_nslr_plots(m,  gsub("_",":",pair_names),figfile=paste("figures/nslr_larmord_wt_", wt, "_cor_", as.character(correct),".pdf", sep = ""))
    generate_all_barplot(prediction_methods = c("larmord"), nucleus_groups = c("both"))
    
    #  **** End supporting information manuscript analysis ****
    save(list=ls(), file= paste("workspaces/workspace_wt_", wt, "_correct_", as.character(correct), ".RData", sep = ""))
  }
}


r <- nuclei_importance(pairs=c("1R2P_2LPS","2FRL_2M22","2H2X_2M21","2KFC_2L1V","2N6Q_5KMZ"), predictor="larmord", weight=1, error_type="mae", outliers="none", nuclei = c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8","H1'","H2'","H3'","H4'","H5'","H5''","H2","H5","H6","H8"))
r$importance/max(r$importance)

stop()
all <- NULL
pair_names <- c("1R2P_2LPS","2N6Q_5KMZ","2H2X_2M21","2FRL_2M22","2KFC_2L1V")
for (pair in pair_names){
  file <- sprintf("chemical_shifts_%s_average.txt", pair)
  m <- get_errors(pairfile = file, nucleus_group = "carbon")
  m$pair  <- pair
  m <- m[order(m$model), ]
  if(is.null(all)){all <- m} else {all <- rbind(all, m)}
  cat(sprintf("%s %4.2f/%4.2f\n", pair, m$error[1], m$error[2]))
}



all <- NULL
pair_names <- c("1R2P_2LPS","2N6Q_5KMZ","2H2X_2M21","2FRL_2M22","2KFC_2L1V")
for (pair in pair_names){
  file <- sprintf("chemical_shifts_%s_average.txt", pair)
  m <- get_errors(pairfile = file, nucleus_group = "both", outliers = "none", error_type = "tau")
  m$pair  <- pair
  m <- m[order(m$model), ]
  if(is.null(all)){all <- m} else {all <- rbind(all, m)}
  cat(sprintf("%s %4.3f/%4.3f\n", pair, m$error[1], m$error[2]))
}

