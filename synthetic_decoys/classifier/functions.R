get_nucleus_errors <- function(errors, nuclei = c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8","H1'","H2'","H3'","H4'","H5'","H5''","H2","H5","H6","H8"), outliers = "none", unweighted = TRUE){
  errors <- errors[errors$nucleus %in% nuclei, ]
  tmp <- expand.grid(unique(errors$model), nuclei, "GUA")
  tmp <- rbind(tmp, expand.grid(unique(errors$model), nuclei, "ADE"))
  tmp <- rbind(tmp, expand.grid(unique(errors$model), nuclei, "CYT"))
  tmp <- rbind(tmp, expand.grid(unique(errors$model), nuclei, "URA"))  
  colnames(tmp) <- c("model", "nucleus","resname")
  cnames <- unique(paste(tmp$nucleus, tmp$resname, sep=":"))
  errors <- merge(errors, tmp, all = TRUE)
  errors <- errors[order(errors$model, errors$nucleus, errors$resname), ]  
  mat <- matrix(errors$V1, byrow = TRUE, ncol = length(cnames))
  colnames(mat) <- cnames
  model_info <- unique(errors[,c("model", "flag", "rmsd")])
  model_info <- model_info[complete.cases(model_info),]
  mat <- cbind(model_info, mat)
  return(mat)
}
