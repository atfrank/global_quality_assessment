# load user functions
source("global_quality_assessment/library.R")

# initialize data 
initialize_analysis()

# read in pair info used for eight (8) tests of the ability of chemical shifts assess the global quality of NMR structures
pairs <- read.table(pair_info, col.names = c("pair", "pair_name", "ref", "threshold"))
pair_names <- c("1R2P_2LPS","2FRL_2M22","2H2X_2M21","2KFC_2L1V","2L94_1Z2J","2N82_2N7X","1Z2J_2L94","2N7X_2N82")

#  **** Begin main manuscript analysis ****
# Part 1:
# True-positive rates analysis (using Consensus predicted chemical shifts)
# Error determined using chemical shifts computed from the average NMR structures
summarize_tables("mean","mae", TRUE, FALSE, FALSE, names=pair_names)

# Errors determined using conformatinally-averaged chemical shifts
summarize_tables("mean","mae", FALSE, FALSE, TRUE, names=pair_names)

# Errors determined for each model in so-called "composite" ensembles
summarize_tables("mean","mae", FALSE, FALSE, FALSE, names=pair_names)

# Part 2:
# Resolving power analysis (NSLR)
m <- get_nslr_only(summarize_tables("mean","mae", FALSE, FALSE, FALSE, names=pair_names))
make_nslr_plots(m,  gsub("_",":",pair_names))

# Part 3:
# Comparing NMR solution and X-ray crystal structures
summarize_tables("mean","mae", FALSE, TRUE, FALSE, c("1SCL_430D","2LPS_1KXK","28SR_1LNT"))
#  **** End main manuscript analysis ****

#  **** Begin supporting information manuscript analysis ****
# Part 1:
# True-positive rates analysis
# (1) RAMSEY
# Error determined using chemical shifts computed from the average NMR structures
get_flags_only(summarize_tables("ramsey","mae", TRUE, FALSE, FALSE, names=pair_names))
# Errors determined using conformatinally-averaged chemical shifts
get_flags_only(summarize_tables("ramsey","mae", FALSE, FALSE, TRUE, names=pair_names))
# Errors determined for each model in so-called "composite" ensembles
get_flags_only(summarize_tables("ramsey","mae", FALSE, FALSE, FALSE, names=pair_names))

# (2) LARMORD
get_flags_only(summarize_tables("larmord","mae", TRUE, FALSE, FALSE, names=pair_names))
# Errors determined using conformatinally-averaged chemical shifts
get_flags_only(summarize_tables("larmord","mae", FALSE, FALSE, TRUE, names=pair_names))
# Errors determined for each model in so-called "composite" ensembles
get_flags_only(summarize_tables("larmord","mae", FALSE, FALSE, FALSE, names=pair_names))

# Part 2:
# Resolving power analysis (NSLR)
# (1) Consensus
m <- get_nslr_only(summarize_tables("mean","mae", FALSE, FALSE, FALSE, names=pair_names))
make_nslr_plots(m,  gsub("_",":",pair_names),figfile="nslr_mean.pdf")
# (2) RAMSEY
m <- get_nslr_only(summarize_tables("ramsey","mae", FALSE, FALSE, FALSE, names=pair_names))
make_nslr_plots(m,  gsub("_",":",pair_names),figfile="nslr_ramsey.pdf")
# (3) LARMORD
m <- get_nslr_only(summarize_tables("larmord","mae", FALSE, FALSE, FALSE, names=pair_names))
make_nslr_plots(m,  gsub("_",":",pair_names),figfile="nslr_larmord.pdf")

# Part 3:
# Comparing NMR solution and X-ray crystal structures
# (1) Consensus
summarize_tables("mean","mae", FALSE, TRUE, FALSE, c("1SCL_430D","2LPS_1KXK","28SR_1LNT"))
# (2) RAMSEY
summarize_tables("ramsey","mae", FALSE, TRUE, FALSE, c("1SCL_430D","2LPS_1KXK","28SR_1LNT"))
# (3) LARMORD
summarize_tables("larmord","mae", FALSE, TRUE, FALSE, c("1SCL_430D","2LPS_1KXK","28SR_1LNT"))
#  **** End supporting information manuscript analysis ****
