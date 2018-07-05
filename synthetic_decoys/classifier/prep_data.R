setwd("~/Documents/GitHub/global_quality_assessment/synthetic_decoys/classifier/")

# read in data
rmsd <- read.table("rmsd_files/all.txt", col.names = c("id", "model", "rmsd"))
predcs <- nmR::load_cs_data(csfile="chemical_shift_files/all.txt", accuracyFile = "data/larmord_accuracy_resname_nucleus.txt", atomBasedWeights = FALSE, names = c("id", "resid", "resname", "nucleus", "predCS", "model"))
expcs <-  nmR::load_cs_data(csfile="chemical_shift_files/all_measured_corrected.txt", names = c("id", "resname", "resid", "nucleus", "expCS", "error"))

# merge
cs <- merge(predcs, rmsd, by = c("model", "id"))
cs <- merge(expcs, cs, by = c("resname", "resid", "nucleus", "id"))
cs$flag <- 0
cs$flag[cs$rmsd < 3 ] <- 1

# get error matrices
errors <- plyr::ddply(.dat=cs, .var=c("model","flag", "rmsd", "id", "resname", "nucleus"), .fun=nmR::score_mae)
errors <- plyr::ddply(.dat=errors, .var=c("id"), .fun = get_error_matrices)

# save data
save(cs, file = "chemical_shifts.RData")
save(errors, file = "error_matrices.RData")
