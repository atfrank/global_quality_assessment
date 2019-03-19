setwd("~/documents/GitHub/global_quality_assessment/synthetic_decoys/")
options(java.parameters = "- Xmx1024m")
library(caret)
library(plyr)
library(dplyr)
library(rJava)
library(extraTrees)

args <- commandArgs(TRUE)


# load in normalized error matrix
errors <- get(load("error_files/errors_unweighted_norm.RData"))
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
test <- errors[errors$id == rna, ]

#inTraining <- createDataPartition(errors$flag, p = .75, list = FALSE)
#training <- errors[inTraining,]
#testing <- errors[-inTraining,]

# balance training set
train <- downSample(train[,colnames(train)!="flag"], train$flag)

# create cv index s.t. decoys of one RNA will not be both in train and cv
train$cv_id = as.numeric(train$id)
folds <- groupKFold(train$cv_id, k = 5)
fitControl <- trainControl(method = "cv", number = 5, verboseIter = TRUE, index = folds,
                           summaryFunction = twoClassSummary, classProbs = TRUE)
train = train[, !colnames(train) %in% c("id","rmsd","cv_id")]

# extraTrees 
et_grid <- expand.grid(mtry = 1:10, numRandomCuts = 1:10)
set.seed(42)
system.time({et_fit = train(Class ~ ., data = train,
               method = "extraTrees",
               trControl = fitControl,
               tuneGrid = et_grid,
               numThreads = 4,
               metric = "ROC")})
print(et_fit)

# predict 
test$pred <- predict(et_fit, test)
test$pred <- ifelse(test$pred == "new", 0, 1)
proba <- predict(et_fit, test, type = "prob")
test$prob <- proba[,2]

test$pred <- as.factor(test$pred)
test_score = as.vector(confusionMatrix(test$pred, test$flag)$byClass)
cnames = names(confusionMatrix(test$pred, test$flag)$byClass)

# write output 
write.table(test[, c("id","rmsd","flag","pred","prob")], paste0("data/et_norm_clf_result_", rna, ".txt"), col.names = F, row.names = F, quote = F)
write.table(data.frame(cnames, test_score), paste0("data/et_norm_clf_score_", rna, ".txt"), col.names = F, row.names = F, quote = F)



