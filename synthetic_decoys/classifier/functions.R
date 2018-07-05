# user function
norm_errors <- function(errors){
  # normalizes the errors based on the mean
  errors$total <- rowSums(errors[,!(colnames(errors) %in% c("id", "model", "reference_flag" ,"rmsd"))], na.rm = TRUE)
  cols <- colnames(errors)[!(colnames(errors) %in% c("id", "model", "reference_flag","rmsd"))]
  for (col in cols){
    #errors[, col] <- errors[, col]/median(errors[, col])
    errors[, col] <- errors[, col]/median(errors[, col])
  }
  return(errors)
}

get_balanced_set_1 <- function(train, N=10, imbalance=1.0, replace=FALSE){
  train_true <- subset(train, reference_flag==1)
  train_false <- subset(train, reference_flag==0)
  n <- min(nrow(train_false),nrow(train_true))
  ifelse(n<N, Nlocal <- n, Nlocal <- N)
  return(rbind(train_false[sample(1:nrow(train_false), Nlocal, replace = replace),], train_true[sample(1:nrow(train_true), floor(Nlocal*imbalance), replace = replace),]))
}

get_balanced_set_2 <- function(train, ngroups=10, n=20){
  # get data that even spans rmsd range
  train <- train[order(train$rmsd),]
  train$group <- cut(train$rmsd, breaks = ngroups, labels = c(1:ngroups))
  train <- train[sample(1:nrow(train)),]
  train <- ddply(.data=train, .variables=c("group"), .fun=function(x)head(x,n))
  return(train[,!(colnames(train) %in% "group")])
}

scores <- function(predictions, N=10){
  # function to quanitfy the ability to correctly classifying models as native or non-native
  # returns the NSLR and the fraction of the high probability models that are considered native
  require(nmR)
  #print(head(predictions))
  return(data.frame(NSLR=nslr(predictions$reference_flag), RMSD=round(mean(predictions$rmsd[1:N]),3), TOPN=mean(as.numeric(predictions$reference_flag[1:N])-1)))
}

nslr <- function(X){
  #' Sum of Logarithmic Ranks Function
  #'
  #' This function allows you to compute the normalized sum of logarithmic ranks
  #' @param X vector of 0 (inactives) and 1 (actives) that was sorted based on some scores (e.g., agreement between measured and predicted shifts)
  #' @export
  #' @examples
  #' random_nslr(sample(c(rep(0,100),rep(1,10))))
  ri <- which(X==1)
  N <- length(X)
  i <-  1:length(ri)
  SLRmax <- -sum(log(i/N)) # logo rank
  return(-sum(log(ri/N))/SLRmax)
} 


get_error_matrices <- function(errors, nuclei = c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8","H1'","H2'","H3'","H4'","H5'","H5''","H2","H5","H6","H8"), remove_na = FALSE){
  # function to write out matrix of chemical shift errors
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
  if(remove_na){mat[is.na(mat)] <- 0}
  return(mat)
}
