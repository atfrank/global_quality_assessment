setwd("~/documents/github/global_quality_assessment/")
source("synthetic_decoys/scripts/functions.R")

# read in rmsd file
rmsd <- read.table("test/rmsd.txt", col.names = c("id","model","rmsd"))

# read in larmord predicted shifts
predcs <- nmR::load_cs_data(csfile = "test/all_computed_larmord.txt", names = c("id", "resid", "resname", "nucleus", "predCS", "model"))

# read in larmord corrected observed shifts
expcs <-  nmR::load_cs_data(csfile="test/all_measured_corrected.txt", names = c("id", "resname", "resid", "nucleus", "expCS", "error"))

# merge
cs <- merge(predcs, rmsd, by = c("model","id"))
cs <- merge(expcs, cs, by = c("resname","resid","nucleus","id"))
# unweighted errors
cs$weight <- 1

# get error matrices
errors <- plyr::ddply(.dat=cs, .var=c("model", "rmsd", "id", "resname", "nucleus"), .fun=nmR::score_mae)
errors <- plyr::ddply(.dat=errors, .var=c("id"), .fun = get_error_matrices)
names <- gsub(":","", gsub("\'","p",colnames(errors)))
colnames(errors) <- names

# save
write.table(errors, "test/test_error_unweighted.txt", col.names = T, row.names = F, quote = F)

# NORMALIZE
rnas <- unique(errors$id)
errors <- as.data.frame(errors)

# scramble data
errors <- errors[sample(1:nrow(errors)),]

# delete NA columns
errors <- errors[,colSums(is.na(errors))<nrow(errors)]

# fill in NA
errors[is.na(errors)] <- 0

# scale and center
normalize <- function(x){
  preProcValues <- preProcess(x[, !colnames(x) %in% c("id","rmsd","model")], method = c("center", "scale"), na.remove = T)
  x_tf <- predict(preProcValues, x)
  return(x_tf)
}
errors <- ddply(.dat = errors, .var = "id", .fun = normalize)
write.table(errors, "test/test_error_unweighted_normed.txt", col.names = T, row.names = F, quote = F)
