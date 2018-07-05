setwd("~/Documents/GitHub/global_quality_assessment/")
train_set <- unlist(strsplit("1KKA 1L1W 1LC6 1LDZ 1NC0 1OW9 1PJY 1R7W 1R7Z 1UUU 1YSV 2FDT 2GM0 2JXQ 2JXS 2K3Z 2K41 2KOC 2KYD 2LBJ 2LBL 2LDL 2LDT 2LHP 2LI4 2LK3 2LP9 2LPA 2LU0 2LUB 2LV0 2RN1 2Y95 4A4S 4A4T 4A4U", " "))
rmsd <- read.table("synthetic_decoys/rmsd_files/all.txt", col.names = c("id", "model", "rmsd"))
predcs <- nmR::load_cs_data(csfile="synthetic_decoys/chemical_shift_files/all.txt", accuracyFile = "data/larmord_accuracy_resname_nucleus.txt", atomBasedWeights = FALSE, names = c("id", "resid", "resname", "nucleus", "predCS", "model"))
expcs <-  nmR::load_cs_data(csfile="synthetic_decoys/chemical_shift_files/all_measured_corrected.txt", names = c("id", "resname", "resid", "nucleus", "expCS", "error"))
cs <- merge(predcs, rmsd, by = c("model", "id"))
cs <- merge(expcs, cs, by = c("resname", "resid", "nucleus", "id"))

cs <- cs[ cs$nucleus %in% c("H2", "C8", "C6") , ]
cs_cleaner <- clean_data(cs)
cs_cleaner$flag <- 0
cs_cleaner$flag[cs_cleaner$rmsd < 3 ] <- 1
errors <- ddply(.dat=cs_cleaner, .var=c("model","flag", "id"), .fun=nmR::score_mae)
colnames(errors) <- c("model","flag", "id", "error")
errors <- merge(errors, rmsd, by = c("model", "id"))
errors <- errors[order(errors$error, decreasing = FALSE), ]
nslrs <- ddply(.dat=errors, .var=c("id"), .fun= function(x){data.frame(nslr=nmR::nslr(x$flag), rmsd=x$rmsd[1])})
nslrs <- nslrs[!(nslrs$id %in% train_set), ]
rownames(nslrs) <- 1:nrow(nslrs)
stop()

r <- nuclei_importance(pairs=c("1R2P_2LPS","2FRL_2M22","2H2X_2M21","2KFC_2L1V","2N6Q_5KMZ"), predictor="larmord", weight=1, error_type="mae", outliers="none", nuclei = c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8","H1'","H2'","H3'","H4'","H5'","H5''","H2","H5","H6","H8"))
r$importance/max(r$importance)

rmsd_frame <- NULL
nslrs_frame <- NULL
rnas <- unlist(strsplit("1XHP 1Z2J 1ZC5 28SR 2JWV 2K66 2L1V 2L3E 2LPS 2LQZ 2LUN 2M12 2M21 2M22 2M24 2M4W 2M5U 2M8K 2MEQ 2MHI 2MNC 2MXL 2N2O 2N2P 2N4L 2N6S 2N6T 2N6W 2N6X 2NCI 2QH2 2QH4 2RVO 5A17 5A18 5IEM 5KH8 5KMZ 5KQE 5LSN 5LWJ 5N5C 5UF3 5UZT 5V16 5V17 5WQ1 6EZ0", " "))
for (rna in rnas){
  file <- sprintf("~/Documents/GitHub/global_quality_assessment/synthetic_decoys/electio_selections/%s_real.txt", rna)
  w <- read.table(file, col.names = c("id", "model", "rmsd", "w"))
  w$status <- 0
  w$status[w$rmsd < 3] <- 1
  w <- w[order(w$w, decreasing = TRUE), ]
  rmsd_frame <- c(rmsd_frame, w$rmsd[1])
  nslrs_frame <- c(nslrs_frame, nmR::nslr(w$status))
  cat(sprintf("%s %4.2f %4.2f\n", rna, w$rmsd[1], nmR::nslr(w$status)))
}
result <- data.frame(rnas, rmsd_frame, nslrs_frame)

