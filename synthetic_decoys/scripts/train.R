setwd("~/projects/global_quality_assessment/NEW/")
library(caret)
library(plyr)
library(dplyr)
library(rJava)
library(extraTrees)
library(ranger)
library(xgboost)

# load errors matrices
errors = get(load("error_files/error_matrix_unweighted.RData"))

# normalize the error for each RNA (specificed by id) separately
#errors <- ddply(.data = errors, .variables = c("id"), .fun = norm_errors)
# rf, ert and xgboost should be insensitive to feature scaling
errors <- as.data.frame(errors)

# scramble data just to be sure
errors <- errors[sample(1:nrow(errors)),]

# make sure the reference_flag is a factor so that when used in randomForest will recognize this as a classification problem
errors$flag <- ifelse(errors$rmsd < 3.0, 1, 0)
errors$flag <- as.factor(errors$flag)

# choosing training and testing RNAs
test_rnas <- c("1R2P","2LPS","2N6Q","5KMZ","2H2X","2M21","2FRL","2M22","2KFC","2L1V")
train_rnas <- unique(errors$id)[!unique(errors$id) %in% test_rnas]
rnames <- c("id", "model", "rmsd")
train <- errors[(errors$id %in% train_rnas), !(colnames(errors) %in% rnames)]

# delete empty columns
train <- train[,colSums(is.na(train))<nrow(train)]

# fill in NA
train[is.na(train)] <- 1

# set column names for training and testing dataframes
names <- gsub(":","", gsub("\'","p",colnames(train)))
colnames(train) <- names

# split training data into train and cv
set.seed(998)
levels(train$flag) <- c("new","old") # Please use factor levels that can be used as valid R variable names.
inTraining <- createDataPartition(train$flag, p = .75, list = FALSE)
training <- train[inTraining,]
testing <- train[-inTraining,]

# grid search tuning in caret
fitControl <- trainControl(method = "cv", number = 5, verboseIter = TRUE,
                           summaryFunction = twoClassSummary, classProbs = TRUE)

# ranger
rf_grid =  expand.grid(mtry = 1:10, splitrule = 'gini', min.node.size = 1)
set.seed(42)
system.time({rf_fit = train(flag ~ ., data = training, 
               method = 'ranger', 
               trControl = fitControl,
               tuneGrid = rf_grid, 
               num.threads = 4,
               metric = "ROC")})
tiff("rf.tiff",width=4,height = 4,units = 'in',res = 300)
plot(rf_fit, main = "random forest")
dev.off()
sink("log.txt")
cat("Best params for RF:\n")
cat(paste(names(rf_fit$bestTune), rf_fit$bestTune, sep = ":", collapse = ","))
sink()

# extraTrees 
et_grid = expand.grid(mtry = 1:10, numRandomCuts = 1:10)
set.seed(42)
system.time({et_fit = train(flag ~ ., data = training,
               method = "extraTrees",
               trControl = fitControl,
               tuneGrid = et_grid,
               numThreads = 4,
               metric = "ROC")})
tiff("et.tiff",width=4,height = 4,units = 'in',res = 300)
plot(et_fit, main = "extra trees")
dev.off()
sink("log.txt")
cat("Best params for ERT:\n")
cat(paste(names(et_fit$bestTune), et_fit$bestTune, sep = ":", collapse = ","))
sink()

# xgboost
set.seed(42)
system.time({gbm_fit = train(flag ~ ., data = training,
               method = 'gbm',
               trControl = fitControl,
               verbose = FALSE,
               tuneLength = 10,
               metric = "ROC")})
tiff("gbm.tiff",width=4,height = 4,units = 'in',res = 300)
plot(gbm_fit, main = "xgboost")
dev.off()
sink("log.txt")
cat("Best params for XGBoost:\n")
cat(paste(names(gbm_fit$bestTune), gbm_fit$bestTune, sep = ":", collapse = ","))
sink()

# make predictions using three models
testing$rf_binary = predict(rf_fit, newdata = testing)
testing$rf_prob = predict(rf_fit, newdata = testing, type = "prob") 
testing$et_binary = predict(et_fit, newdata = testing)
testing$et_prob = predict(et_fit, newdata = testing, type = "prob") 
testing$gbm_binary = predict(gbm_fit, newdata = testing)
testing$gbm_prob = predict(gbm_fit, newdata = testing, type = "prob") 

# save prediction for further analysis
save(testing, file = "cv_preds.RData")