setwd("~/projects/global_quality_assessment/NEW/")
options(java.parameters = "- Xmx1024m")
library(caret)
library(plyr)
library(dplyr)
library(glmnet)
args <- commandArgs(TRUE)

# load in normalized error matrix
errors <- read.table("error_files/errors_unweighted_normed.txt", header = T)
errors$flag <- as.factor(errors$flag)
rnas <- as.character(unique(errors$id))

# load sequence similarity list
M <- get(load("sequence_files/remove_list_0.8.RData"))

# leave one out cross validation
rna <- rnas[as.integer(args[1])+1]
print(rna)
# split train and test
train <- errors[!errors$id %in% M[[rna]], ]
#train$id <- NULL
levels(train$flag) <- c("new","old") # Please use factor levels that can be used as valid R variable names.

# balance training set
train <- downSample(train[,colnames(train)!="flag"], train$flag)

# create cv index s.t. decoys of one RNA will not be both in train and cv
train$cv_id = as.numeric(train$id)
folds <- groupKFold(train$cv_id, k = 5)
fitControl <- trainControl(method = "cv", number = 5, verboseIter = TRUE, index = folds,
                           summaryFunction = twoClassSummary, classProbs = TRUE)
train = train[, !colnames(train) %in% c("id","rmsd","cv_id")]

# logistic regression
logit_grid <- expand.grid(alpha = 1, lambda = seq(0.001,0.1,by = 0.001))
set.seed(42)
system.time({logit_fit = train(Class ~ ., data = train,
               method = "glmnet",
               trControl = fitControl,
               tuneGrid = logit_grid,
               metric = "ROC")})
print(logit_fit)

# use GOOD decoys to center and scale all decoys
#errors_raw <- read.table("error_files/error_matrix_unweighted.txt", header = T)
#test <- subset(errors_raw, id == rna)
#test <- test[, colnames(errors)[!colnames(errors) %in% c("flag")]]
#test[is.na(test)] <- 0
#test_in <- subset(test, rmsd < 1.5 | rmsd > 3.5)
#test_in[, c("id","rmsd")] <- NULL
#preProcValues <- preProcess(test_in, method = c("center", "scale"), na.remove = T)
#test <- predict(preProcValues, test)


# NOT FILTERING
test <- subset(errors, id == rna)

# predict 
test$pred <- predict(logit_fit, test)
test$pred <- ifelse(test$pred == "new", 0, 1)
test$pred <- as.factor(test$pred)
proba <- predict(logit_fit, test, type = "prob")
test$prob <- proba[,2]

test_score = as.vector(confusionMatrix(test$pred, test$flag)$byClass)
cnames = names(confusionMatrix(test$pred, test$flag)$byClass)

# write output 
write.table(test[, c("id","rmsd","flag","pred","prob")], paste0("data/logit_unfiltered_clf_result_", rna, ".txt"), col.names = F, row.names = F, quote = F)
write.table(data.frame(cnames, test_score), paste0("data/logit_unfiltered_clf_score_", rna, ".txt"), col.names = F, row.names = F, quote = F)



