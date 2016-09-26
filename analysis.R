# goto analysis directory
setwd("~/GitSoftware/nar_manuscript_analysis/")

# source need library
source("analysis_scripts/make_tables.R")
source("analysis_scripts/analyze_data.R")
source("analysis_scripts/compile_data.R")
source("analysis_scripts/merge_free_bound.R")
source("analysis_scripts/summary.R")
source("analysis_scripts/correlations_scores.R")
source("analysis_scripts/make_error_barplot.R")
source("analysis_scripts/make_nslr_plots.R")

# compile the data used for analysis
pairs <- read.table("data/tests_info", col.names = c("pair", "pair_name", "ref", "threshold"))

# data for composite NMR bundle
compile_data(pairs, FALSE, FALSE)

# data for average structure
compile_data(pairs, TRUE, FALSE)

# data for nmr-xray pairs
compile_data(pairs, FALSE, TRUE)

# merge free and bound chemical shifts to ensure we only use common chemical shifts
merge_free_bound("2N7X", "2N82")
merge_free_bound("1Z2J", "2L94")

# TPR - true positive rates
pair_names <- c("1R2P_2LPS","2FRL_2M22","2H2X_2M21","2KFC_2L1V","2L94_1Z2J","2N82_2N7X","1Z2J_2L94","2N7X_2N82")

# Errors determined using chemical shifts computed average structure
summarize_tables("mean","mae", TRUE, FALSE, FALSE, names=pair_names)

# Errors determined using conformatinally-averaged chemical shifts
summarize_tables("mean","mae", FALSE, FALSE, TRUE, names=pair_names)

# Errors determined for each model in so-called "composite" ensembles
summarize_tables("mean","mae", FALSE, FALSE, FALSE, names=pair_names)

# Compare NSLR
m <- get_nslr_only(summarize_tables("mean","mae", FALSE, FALSE, FALSE, names=pair_names))
make_nslr_plots(m,  gsub("_",":",pair_names))

# Errors determined from average NMR and X-ray
summarize_tables("mean","mae", FALSE, TRUE, FALSE, c("1SCL_430D","2LPS_1KXK","28SR_1LNT"))

# Supporting Information
pair_names <- c("1R2P_2LPS","2FRL_2M22","2H2X_2M21","2KFC_2L1V","2L94_1Z2J","2N82_2N7X","1Z2J_2L94","2N7X_2N82")
# True-positive rates
# (1) RAMSEY
get_flags_only(summarize_tables("ramsey","mae", TRUE, FALSE, FALSE, names=pair_names))
get_flags_only(summarize_tables("ramsey","mae", FALSE, FALSE, TRUE, names=pair_names))
get_flags_only(summarize_tables("ramsey","mae", FALSE, FALSE, FALSE, names=pair_names))

# (2) LARMORD
get_flags_only(summarize_tables("larmord","mae", TRUE, FALSE, FALSE, names=pair_names))
get_flags_only(summarize_tables("larmord","mae", FALSE, FALSE, TRUE, names=pair_names))
get_flags_only(summarize_tables("larmord","mae", FALSE, FALSE, FALSE, names=pair_names))

# NSLR
# (1) Consensus
m <- get_nslr_only(summarize_tables("mean","mae", FALSE, FALSE, FALSE, names=pair_names))
make_nslr_plots(m,  gsub("_",":",pair_names),figfile="nslr_mean.pdf")
# (2) RAMSEY
m <- get_nslr_only(summarize_tables("ramsey","mae", FALSE, FALSE, FALSE, names=pair_names))
make_nslr_plots(m,  gsub("_",":",pair_names),figfile="nslr_ramsey.pdf")
# (3) LARMORD
m <- get_nslr_only(summarize_tables("larmord","mae", FALSE, FALSE, FALSE, names=pair_names))
make_nslr_plots(m,  gsub("_",":",pair_names),figfile="nslr_larmord.pdf")


# Errors determined from average NMR and X-ray
# (1) Consensus
summarize_tables("mean","mae", FALSE, TRUE, FALSE, c("1SCL_430D","2LPS_1KXK","28SR_1LNT"))
# (2) RAMSEY
summarize_tables("ramsey","mae", FALSE, TRUE, FALSE, c("1SCL_430D","2LPS_1KXK","28SR_1LNT"))
# (3) LARMORD
summarize_tables("larmord","mae", FALSE, TRUE, FALSE, c("1SCL_430D","2LPS_1KXK","28SR_1LNT"))


