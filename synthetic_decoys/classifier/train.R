# load errors matrices
load("error_matrices.RData")

# scramble data just to be sure
errors <- errors[sample(1:nrow(errors)),]

# normalize the error for each RNA (specificed by id) separately
library(plyr)
library(dplyr)
errors <- ddply(.data = errors, .variables = c("id"), .fun = norm_errors)
errors <- as.data.frame(errors)

# make sure the reference_flag is a factor so that when used in randomForest will recognize this as a classification problem
errors$flag <- ifelse(errors$rmsd < 3.0, 1, 0)
errors$flag <- as.factor(errors$flag)


# SPLIT INTO TRAINING AND TESTING
# choosing training and testing RNAs
train_rnas <- "1KKA 1L1W 1LC6 1LDZ 1NC0 1OW9 1PJY 1R7W 1UUU 1XHP 1YSV 1Z2J 1ZC5 2FDT 2JWV 2K66 2KOC 2L3E 2LBJ 2LBL 2LDL 2LDT 2LPA 2LQZ 2LU0 2LUB"
train_rnas <- unlist(strsplit(train_rnas, " "))
rnames <- c("id", "model", "rmsd")

train <- errors[(errors$id %in% train_rnas), !(colnames(errors) %in% rnames)]
train <- ddply(.data = errors[(errors$id %in% train_rnas), !(colnames(errors) %in% c("model"))], .variables = c("id"), .fun = get_balanced_set_1, N=600, imbalance=0.25, replace=FALSE)
train <- train[(train$id %in% train_rnas), !(colnames(train) %in% rnames)]

test <- errors[!(errors$id %in% train_rnas), !(colnames(errors) %in% rnames)]
test_info <- errors[!(errors$id %in% train_rnas), (colnames(errors) %in% rnames)]

train[is.na(train)] <- 1.0
test[is.na(test)] <- 1.0

# set column names for training and testing dataframes
names <- unlist(strsplit("reference_flag C1pGUA C2pGUA C3pGUA C4pGUA C5pGUA C2GUA C5GUA C6GUA C8GUA H1pGUA H2pGUA H3pGUA H4pGUA H5pGUA H5ppGUA H2GUA H5GUA H6GUA H8GUA C1pADE C2pADE C3pADE C4pADE C5pADE C2ADE C5ADE C6ADE C8ADE H1pADE H2pADE H3pADE H4pADE H5pADE H5ppADE H2ADE H5ADE H6ADE H8ADE C1pCYT C2pCYT C3pCYT C4pCYT C5pCYT C2CYT C5CYT C6CYT C8CYT H1pCYT H2pCYT H3pCYT H4pCYT H5pCYT H5ppCYT H2CYT H5CYT H6CYT H8CYT C1pURA C2pURA C3pURA C4pURA C5pURA C2URA C5URA C6URA C8URA H1pURA H2pURA H3pURA H4pURA H5pURA H5ppURA H2URA H5URA H6URA H8URA total", " "))
colnames(train) <- names
colnames(test) <- names

write.table(train, file="train.txt", col.names = F, row.names = F, quote = F)
write.table(test, file="test.txt", col.names = F, row.names = F, quote = F)
write.table(test_info, file="test_info.txt", col.names = F, row.names = F, quote = F)

# BUILD CLASSIFIER
library(randomForest)
library(caret)
library(RRF)
rf <- randomForest(formula = reference_flag~., data = train, na.action = na.exclude, replace = FALSE, nodesize=5, do.trace = TRUE, ntree = 500)

# Make predictions using classifier
test$reference_flag_binary <- predict(rf, test, type = "response")
test$reference_flag_prob <- predict(rf, test, type = "prob")
t <- confusionMatrix(test$reference_flag, predict(rf, test))
test <- cbind(test_info, test)

# get performance after sorting based on total error
test <- test[order(test[,c("total")], decreasing = FALSE), ]
test_nslr_total_error <- ddply(.data = test, .variables = c("id"), .fun = scores, N=10)

# get nmR package ready to use
#install.packages("nmR")
library("nmR")

# get performance after sorting based on predicted probability of being native
test <- test[order(test[,c("reference_flag_prob")][,1]), ]
test_nslr_class <- ddply(.data = test, .variables = c("id"), .fun = scores, N=10)

# compare results obtained when using the total error and classifier
comp <- cbind(test_nslr_total_error, test_nslr_class)
comp$diff <- round(test_nslr_total_error$NSLR - test_nslr_class$NSLR,4)
print(comp)