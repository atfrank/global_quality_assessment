setwd("~/documents/github/global_quality_assessment/synthetic_decoys")
source("scripts/functions.R")

# read in rmsd file
rmsd <- read.table("rmsd_files/all.txt", col.names = c("id","model","rmsd"))

# read in larmord predicted shifts
predcs <- nmR::load_cs_data(csfile = "chemical_shift_files/all.txt", names = c("id", "resid", "resname", "nucleus", "predCS", "model"))

# read in larmord corrected observed shifts
expcs <-  nmR::load_cs_data(csfile="observed_shifts_corrected/all_measured_corrected.txt", names = c("id", "resname", "resid", "nucleus", "expCS", "error"))

# merge
cs <- merge(predcs, rmsd, by = c("model","id"))
cs <- merge(expcs, cs, by = c("resname","resid","nucleus","id"))
# unweighted errors
cs$weight <- 1

# get error matrices
errors <- plyr::ddply(.dat=cs, .var=c("model", "rmsd", "id", "resname", "nucleus"), .fun=nmR::score_mae)
errors <- plyr::ddply(.dat=errors, .var=c("id"), .fun = get_error_matrices)

# save
save(cs, file = "error_files/chemical_shifts.RData")
save(errors, file = "error_files/error_matrix_unweighted.RData")

#predcs2 <- nmR::load_cs_data(csfile="chemical_shift_files/all.txt", accuracyFile = "data/larmord_accuracy_resname_nucleus.txt", atomBasedWeights = FALSE, names = c("id", "resid", "resname", "nucleus", "predCS", "model"))
#cs2 <- merge(predcs2, rmsd, by = c("model","id"))
#cs2 <- merge(expcs, cs2, by = c("resname","resid","nucleus","id"))
#errors <- plyr::ddply(.dat=cs2, .var=c("model", "rmsd", "id", "resname", "nucleus"), .fun=nmR::score_mae)
#errors <- plyr::ddply(.dat=errors, .var=c("id"), .fun = get_error_matrices)
#save(errors, file = "error_files/error_matrix_weighted.RData")