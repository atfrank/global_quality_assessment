# cutoff selection
errors <- read.table("synthetic_decoys/error_files/errors_unweighted_normed.txt", header = T)
rnas <- as.character(unique(errors$id))
cutoffs = c(1.5, 2, 2.5, 3, 3.5, 4)

M = matrix(NA, length(rnas), length(cutoffs))

for(i in 1:length(cutoffs)){
  for(j in 1:length(rnas)){
    rna = rnas[j]
    cutoff = cutoffs[i]
    x = read.table(paste0("synthetic_decoys/data/cutoff_test/et_clf_result_",rna,"_",cutoff,".txt"),
                   col.names = c("id","rmsd","flag","pred","prob"))
    x = x[order(x$prob),]
    M[j, i] = x[1, "rmsd"]
  }
}
colnames(M) = cutoffs
rownames(M) = rnas
write.table(M, "synthetic_decoys/data/cutoff_selection_based_on_best_decoy.R", 
            col.names = T, row.names = T, quote = T)
#
#1.5               2               2.5               3               3.5        
#Min.   :0.0141   Min.   :0.0140   Min.   :0.0140   Min.   :0.0140   Min.   :0.0140  
#1st Qu.:0.1391   1st Qu.:0.1395   1st Qu.:0.1389   1st Qu.:0.1371   1st Qu.:0.2266  
#Median :0.3920   Median :0.3933   Median :0.3892   Median :0.2303   Median :0.4150  
#Mean   :0.4629   Mean   :0.4818   Mean   :0.4202   Mean   :0.4395   Mean   :0.4692  
#3rd Qu.:0.5323   3rd Qu.:0.5940   3rd Qu.:0.5323   3rd Qu.:0.5365   3rd Qu.:0.6071  
#Max.   :2.5577   Max.   :3.1093   Max.   :2.2448   Max.   :3.1093   Max.   :1.6192  
#4         
#Min.   :0.0137  
#1st Qu.:0.1391  
#Median :0.3912  
#Mean   :0.4627  
#3rd Qu.:0.5365  
#Max.   :2.5102  

# apply on testing systems
cutoff = 2.5
errors = read.table("test/test_error_unweighted_normed.txt",header = T)








