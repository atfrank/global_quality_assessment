rnas = unique(errors$id)
selected_models = matrix(NA,67,4)
errors$mae = rowMeans(errors[,5:62],na.rm=T)
for(i in 1:67){
  rna = rnas[i]
  cs_errors = subset(errors, id==rna)
  selected_models[i,1] = rna
  selected_models[i,2] = cs_errors[which.min(cs_errors$mae), "rmsd"] # select using mae
}
boxplot(selected_models[,3], horizontal = T, axes = F, col = "grey", staplewex = 1,main="select using cs sum of error (normalized)")
text(x=fivenum(selected_models[,4]), labels =fivenum(selected_models[,2]), y=1.25)

selected_models = matrix(NA,67,2)
for(i in 1:67){
  rna = rnas[i]
  cs_errors = subset(errors, id==rna)
  cs_errors = cs_errors[,c("id","model","rmsd")]
  pred = read.table(paste0("RF_classifier/pred/rf_test_predscore_",rna,"_3.0.txt"),header = T)
  #selected_models[i,1] = rna
  merged = merge(cs_errors, pred, by=c("id","model"))
  selected_models[i,3] = merged[which.max(merged$score), "rmsd"]
}
for(i in 1:67){
  rna = rnas[i]
  cs_errors = subset(errors, id==rna)
  cs_errors = cs_errors[,c("id","model","rmsd")]
  pred = read.table(paste0("ERT_classifier/pred/ert_test_predscore_",rna,"_3.0.txt"),header = T)
  #selected_models[i,1] = rna
  merged = merge(cs_errors, pred, by=c("id","model"))
  selected_models[i,4] = merged[which.max(merged$score), "rmsd"]
}
