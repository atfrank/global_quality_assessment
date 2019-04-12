# user function
norm_errors <- function(errors){
  # normalizes the errors based on the mean
  #errors$total <- rowSums(errors[,!(colnames(errors) %in% c("id", "model", "reference_flag" ,"rmsd"))], na.rm = TRUE)
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
  # colnames order should match errors
  tmp$resname <- as.character(tmp$resname)
  tmp <- tmp[order(tmp$model, tmp$nucleus, tmp$resname),]
  cnames <- unique(paste(tmp$nucleus, tmp$resname, sep=":"))
  
  errors <- merge(errors, tmp, all = TRUE)
  #errors <- ddply(.data = errors, .variables = "id", .fun = function(x){merge(x, tmp, all = T)})
  errors <- errors[order(errors$model, errors$nucleus, errors$resname), ]  
  mat <- matrix(errors$V1, byrow = TRUE, ncol = length(cnames))
  colnames(mat) <- cnames
  model_info <- unique(errors[,c("model", "rmsd")])
  #model_info <- ddply(.data = errors, .variables = "id", .fun = function(x){unique(x[,c("model","rmsd")])})
  model_info <- model_info[complete.cases(model_info),]
  mat <- cbind(model_info, mat)
  if(remove_na){mat[is.na(mat)] <- 0}
  return(mat)
}

get_error_matrices_test <- function(errors, nuclei = c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8","H1'","H2'","H3'","H4'","H5'","H5''","H2","H5","H6","H8"), remove_na = FALSE){
  # function to write out matrix of chemical shift errors
  errors <- errors[errors$nucleus %in% nuclei, ]
  tmp <- expand.grid(unique(errors$id), nuclei, "GUA")
  tmp <- rbind(tmp, expand.grid(unique(errors$id), nuclei, "ADE"))
  tmp <- rbind(tmp, expand.grid(unique(errors$id), nuclei, "CYT"))
  tmp <- rbind(tmp, expand.grid(unique(errors$id), nuclei, "URA"))  
  colnames(tmp) <- c("id", "nucleus","resname")
  cnames <- gsub("\'","p",unique(paste(tmp$nucleus, tmp$resname, sep="")))
  errors <- merge(errors, tmp, all = TRUE)
  errors <- errors[order(errors$id, errors$nucleus, errors$resname), ]  
  mat <- matrix(errors$V1, byrow = TRUE, ncol = length(cnames))
  colnames(mat) <- cnames
  if(remove_na){mat[is.na(mat)] <- 0}
  return(mat)
}

prep_data_test <- function(expcs_file, predcs_file, rna = "2LPS", errors_file){
  # read in data
  rmsd <- read.table("../../data/rmsd_info", col.names = c("id", "model", "rmsd"), stringsAsFactors = F)
  rmsd <- rmsd[grepl(rna, rmsd$model),]
  predcs <- nmR::load_cs_data(csfile = predcs_file, accuracyFile = "data/larmord_accuracy_resname_nucleus.txt", atomBasedWeights = FALSE, names = c("id", "resid", "resname", "nucleus", "predCS", "model"))
  expcs <-  nmR::load_cs_data(csfile = expcs_file, names = c("resname", "resid", "nucleus", "expCS", "error"))
  
  # merge
  cs <- merge(predcs, rmsd, by = c("model","id"))
  cs <- merge(expcs, cs, by = c("resname", "resid", "nucleus"))
  cs$flag <- 0
  cs$flag[cs$rmsd < 3] <- 1

  # get error matrix
  errors <- plyr::ddply(.dat=cs, .var=c("model","flag", "rmsd", "id", "resname", "nucleus"), .fun=nmR::score_mae)
  errors <- plyr::ddply(.dat=errors, .var=c("id"), .fun = get_error_matrices_test)

  # save data
  write.table(errors, errors_file, col.names = T, row.names = F, quote = F)
}

score_sequence_similarity <- function(similarity_file = "sequence_similarity.txt", rnas, threshold, output){
  ###########################################################################################
  # prepare input sequence similarity file:                                                 #
  # vsearch --allpairs_global fasta.txt --iddef 0 --id 0 --alnout output.aln                #
  # grep 'PDB' output.aln > tmp.txt                                                         #
  # sed '/nt/d' ./tmp.txt > output.txt                                                      #
  ###########################################################################################
  library(stringr)
  # read in vsearch output
  data <- read.table(similarity_file, fill = T, header = F)
  data$V3 <- gsub("[[:space:]]", "", str_replace_all(gsub("*:A|PDBID|CHAIN|SEQUENCE","",as.character(data$V3)), "[^[:alnum:]]", " "))
  data$V1 <- as.character(data$V1)
  
  # construct half similarity matrix
  M <- matrix(NA, length(rnas), length(rnas))
  query_rnas <- NULL
  colnames(M) <- rnas
  rownames(M) <- rnas
  for(i in 1:nrow(data)){
    if(data$V1[i] == "Query"){
      query_rna <- gsub("[^[:alnum:]]", "", gsub("*:A|PDBID|CHAIN|SEQUENCE","",data$V2[i]))
      if(query_rna %in% rnas){
        query_rnas <- c(query_rnas, query_rna)
      }
    } else{
      if(data$V3[i] %in% rnas){M[rownames(M) == query_rna, data$V3[i]] <- data$V1[i]}
    }
  }
  for(col in 1:ncol(M)){
    M[,col] <- as.numeric(gsub("%", "", M[,col]))/100
  }
  
  # complete M to a full matrix
  M[lower.tri(M)] <- t(M)[lower.tri(M)]
  diag(M) <- 1
  
  # create a list to store which rnas should be removed
  x <- vector(mode = "list", length = length(rnas))
  names(x) <- rnas
  for(rna in rnas){
    a = M[rownames(M)==rna, ]
    remove = as.vector(names(a[a>0.8]))
    x[[rna]] <- remove
  }
  save(x, file = output)
}