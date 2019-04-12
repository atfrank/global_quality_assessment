# load errors matrices
library(plyr)
library(caret)
source("scripts/functions.R")
errors <- get(load("error_files/error_matrix_unweighted.RData"))
rnas <- unique(errors$id)

# normalize the error for each RNA (specificed by id) separately
#errors <- ddply(.data = errors, .variables = c("id"), .fun = norm_errors)
# rf, ert and xgboost should be insensitive to feature scaling
errors <- as.data.frame(errors)

# scramble data just to be sure
errors <- errors[sample(1:nrow(errors)),]

# make sure the reference_flag is a factor so that when used in randomForest will recognize this as a classification problem
#errors <- errors[!(errors$rmsd <= 3.5 & errors$rmsd >= 1.5), ]
#errors$flag <- ifelse(errors$rmsd < 2, 0, 1)
#errors$flag <- ifelse(errors$rmsd < 3, 0, 1)
#errors$flag <- as.factor(errors$flag)



# choosing training and testing RNAs
test_rnas <- c("1R2P","2LPS","2N6Q","5KMZ","2H2X","2M21","2FRL","2M22","2KFC","2L1V")
train_rnas <- unique(errors$id)[!unique(errors$id) %in% test_rnas]
rnames <- c("model")
errors <- errors[(errors$id %in% train_rnas), !(colnames(errors) %in% rnames)]

# delete empty columns
errors <- errors[,colSums(is.na(errors))<nrow(errors)]

# fill in NA
#for(i in 4:ncol(errors)){
#  errors[is.na(errors[,i]), i] <- mean(errors[,i], na.rm = TRUE)
#}
#errors[is.na(errors)] <- 0

# set column names for training and testing dataframes
names <- gsub(":","", gsub("\'","p",colnames(errors)))
colnames(errors) <- names

normalize <- function(x){
  preProcValues <- preProcess(x[, !colnames(x) %in% c("id","rmsd","flag")], method = c("center", "scale"), na.remove = T)
  x_tf <- predict(preProcValues, x)
  return(x_tf)
}

# normalize rna-wise
#errors <- ddply(.dat = errors, .var = "id", .fun = normalize)
errors <- ddply(.dat = errors, .var = "id", .fun = norm_errors)

# fill in NA
for(i in 4:ncol(errors)){
  errors[is.na(errors[,i]), i] <- mean(errors[,i], na.rm = TRUE)
}

# save
#save(errors, "error_files/errors_unweighted_norm.RData")
write.table(errors, "error_files/errors_unweighted_median_scaled.txt", col.names = T, row.names = F, quote = F)