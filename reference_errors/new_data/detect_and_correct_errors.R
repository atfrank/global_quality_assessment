library(nmR)
setwd("~/Documents/GitHub/global_quality_assessment/reference_errors/new_data/")
expnames <- c("resname","resid","nucleus","expCS","error")
prednames <- c("model","resid","resname","nucleus","predCS","id")  
rnas <- unlist(strsplit("5KQE 5IEM 5LWJ 2RVO 5UZT 5UF3 5WQ1 5V17 5V16 5KH8 5LSN 6EZ0 5N5C"," "))
for (rna in rnas){
  # Detect errors
  # Read in data
  expcs <- read.table(paste("~/Documents/GitHub/global_quality_assessment/reference_errors/new_data/",rna,"/chemical_shifts.txt.gz",sep=""), col.names = expnames, stringsAsFactors = FALSE)
  expcs <- subset(expcs, expCS != 0.0)
  predcs_larmord <- read.table(paste("~/Documents/GitHub/global_quality_assessment/reference_errors/new_data/",rna,"/chemical_shifts.larmord.txt.gz",sep=""), col.names = prednames, stringsAsFactors = FALSE)
  
  # LARMORD
  reference <- detect_referencing_error(expcs, predcs_larmord,  FALSE)  
  accuracy <- accuracy_estimate(expcs, predcs_larmord, verbose = FALSE)  
  save(reference, file=paste("~/Documents/GitHub/global_quality_assessment/reference_errors/new_data/",rna,"/referencing_errors_larmord.RData",sep=""))
  save(accuracy, file=paste("~/Documents/GitHub/global_quality_assessment/reference_errors/new_data/",rna,"/accuracy_larmord.RData",sep=""))
  
  # Correct Shifts
  # LARMORD
  errors <- correct_referencing_error(expcs, predcs_larmord, 5, -2, FALSE)  
  expcs_corrected <- errors$shifts
  errors <- errors$errors
  save(expcs_corrected, file=paste("~/Documents/GitHub/global_quality_assessment/reference_errors/new_data/",rna,"/chemical_shifts_corrected_larmord.RData",sep=""))
  save(errors, file=paste("~/Documents/GitHub/global_quality_assessment/reference_errors/new_data/",rna,"/detected_errors_larmord.RData",sep=""))
  write.table(expcs_corrected, file=paste("~/Documents/GitHub/global_quality_assessment/reference_errors/new_data/",rna,"/chemical_shifts_corrected_larmord.txt",sep="") , quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(expcs_corrected, file=paste("~/Documents/GitHub/global_quality_assessment/reference_errors/new_data/corrected_shifts/observed_shifts_corrected_larmord_", rna,".txt",sep="") , quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(errors, file=paste("~/Documents/GitHub/global_quality_assessment/reference_errors/new_data/",rna,"/detected_errors_larmord.txt",sep="") , quote = FALSE, row.names = FALSE, col.names = FALSE)

}

