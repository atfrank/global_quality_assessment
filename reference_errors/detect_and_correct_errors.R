library(nmR)
setwd("/Users/atfrank/Documents/GitSoftware/chemical_shift_referencing_errors/")
expnames <- c("resname","resid","nucleus","expCS","error")
prednames <- c("model","resid","resname","nucleus","predCS","id")  
rnas <- unlist(strsplit("2LPS 28SR 1KKA 1L1W 1LC6 1LDZ 1NC0 1OW9 1PJY 1R7W 1SCL 1UUU 1XHP 1YSV 1Z2J 1ZC5 2FDT 2JWV 2K66 2KOC 2L3E 2LBJ 2LBL 2LDL 2LDT 2LHP 2LI4 2LK3 2LP9 2LPA 2LQZ 2LU0 2LUB 2LUN 2LV0 2M12 2M21 2M22 2M24 2M4W 2M5U 2M8K 2MEQ 2MFD 2MHI 2MNC 2MXL 2N2O 2N2P 2N4L 2N6S 2N6T 2N6W 2N6X 2NCI 2QH2 2QH4 2Y95 4A4S 4A4T 4A4U 5A17 5A18 5KQE 5KMZ 2L1V"," "))
for (rna in rnas){
  # Detect errors
  # Read in data
  expcs <- read.table(paste("/Users/atfrank/GitSoftware/chemical_shift_referencing_errors/",rna,"/chemical_shifts.txt.gz",sep=""), col.names = expnames, stringsAsFactors = FALSE)
  expcs <- subset(expcs, expCS != 0.0)
  predcs_larmord <- read.table(paste("/Users/atfrank/GitSoftware/chemical_shift_referencing_errors/",rna,"/chemical_shifts.larmord.txt.gz",sep=""), col.names = prednames, stringsAsFactors = FALSE)
  predcs_ramsey <- read.table(paste("/Users/atfrank/GitSoftware/chemical_shift_referencing_errors/",rna,"/chemical_shifts.ramsey.txt.gz",sep=""), col.names = prednames, stringsAsFactors = FALSE)
  
  # LARMORD
  reference <- detect_referencing_error(expcs, predcs_larmord,  FALSE)  
  accuracy <- accuracy_estimate(expcs, predcs_larmord, verbose = FALSE)  
  save(reference, file=paste("/Users/atfrank/GitSoftware/chemical_shift_referencing_errors/",rna,"/referencing_errors_larmord.RData",sep=""))
  save(accuracy, file=paste("/Users/atfrank/GitSoftware/chemical_shift_referencing_errors/",rna,"/accuracy_larmord.RData",sep=""))
  
  # RAMSEY
  reference <- detect_referencing_error(expcs, predcs_ramsey, FALSE)  
  accuracy <- accuracy_estimate(expcs, predcs_larmord, verbose = FALSE)  
  save(reference, file=paste("/Users/atfrank/GitSoftware/chemical_shift_referencing_errors/",rna,"/referencing_errors_ramsey.RData",sep=""))
  save(accuracy, file=paste("/Users/atfrank/GitSoftware/chemical_shift_referencing_errors/",rna,"/accuracy_ramsey.RData",sep=""))
  cat(sprintf("done checking referencing errors %s\n",rna))

  # Correct Shifts
  # LARMORD
  errors <- correct_referencing_error(expcs, predcs_larmord, 5, -2, FALSE)  
  expcs_corrected <- errors$shifts
  errors <- errors$errors
  save(expcs_corrected, file=paste("/Users/atfrank/GitSoftware/chemical_shift_referencing_errors/",rna,"/chemical_shifts_corrected_larmord.RData",sep=""))
  save(errors, file=paste("/Users/atfrank/GitSoftware/chemical_shift_referencing_errors/",rna,"/detected_errors_larmord.RData",sep=""))
  write.table(expcs_corrected, file=paste("/Users/atfrank/GitSoftware/chemical_shift_referencing_errors/",rna,"/chemical_shifts_corrected_larmord.txt",sep="") , quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(expcs_corrected, file=paste("/Users/atfrank/GitSoftware/chemical_shift_referencing_errors/corrected_shifts/observed_shifts_corrected_larmord_", rna,".txt",sep="") , quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(errors, file=paste("/Users/atfrank/GitSoftware/chemical_shift_referencing_errors/",rna,"/detected_errors_larmord.txt",sep="") , quote = FALSE, row.names = FALSE, col.names = FALSE)

  # RAMSEY
  errors <- correct_referencing_error(expcs, predcs_ramsey, 5, -2, verbose = FALSE)  
  expcs_corrected <- errors$shifts
  errors <- errors$errors
  save(expcs_corrected, file=paste("/Users/atfrank/GitSoftware/chemical_shift_referencing_errors/",rna,"/chemical_shifts_corrected_ramsey.RData",sep=""))
  save(errors, file=paste("/Users/atfrank/GitSoftware/chemical_shift_referencing_errors/",rna,"/detected_errors_ramsey.RData",sep=""))
  write.table(expcs_corrected, file=paste("/Users/atfrank/GitSoftware/chemical_shift_referencing_errors/",rna,"/chemical_shifts_corrected_ramsey.txt",sep="") , quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(expcs_corrected, file=paste("/Users/atfrank/GitSoftware/chemical_shift_referencing_errors/corrected_shifts/observed_shifts_corrected_ramsey_", rna,".txt",sep="") , quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(errors, file=paste("/Users/atfrank/GitSoftware/chemical_shift_referencing_errors/",rna,"/detected_errors_ramsey.txt",sep="") , quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat(sprintf("done correcting referencing errors %s\n",rna))
  
}

