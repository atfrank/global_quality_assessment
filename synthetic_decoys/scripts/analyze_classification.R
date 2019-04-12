setwd("~/documents/GitHub/global_quality_assessment/synthetic_decoys/")
library(caret)
rnas <- unique(read.table("error_files/errors_unweighted_norm.txt", header = T)$id)
M_logit <- NULL
for(rna in rnas){
  print("one rna done")
  pred <- read.table(paste0("data/logit_unfiltered_clf_result_", rna, ".txt"),
                           col.names = c("id","rmsd","class","pred","prob"))
  
  x_in <- pred[pred$rmsd > 3.5 | pred$rmsd < 1.5, ]
  x_in$ref <- as.factor(ifelse(x_in$rmsd > 3.5, 1, 0))
  x_in$pred <- as.factor(x_in$pred)
  tmp <- as.vector(confusionMatrix(x_in$pred, x_in$ref)$byClass)
  if(is.null(M_logit)){
    M_logit <- tmp
  } else {
    M_logit <- rbind(M_logit, tmp)
  }
}

# draw density plot of two classes
x = NULL
for(rna in rnas){
  print("one rna done")
  pred <- read.table(paste0("data/logit_unfiltered_clf_result_", rna, ".txt"),
                     col.names = c("id","rmsd","class","pred","prob"))
  
  x_out <- pred[pred$rmsd <= 3.5 & pred$rmsd >= 1.5, ]
  if(is.null(x)){
    x <- x_out
  } else {
    x <- rbind(x, x_out)
  }
}

x$class = as.factor(x$class)
means = ddply(x, "class", summarise, rmsd.mean = mean(rmsd))
ggplot(x, aes(x = rmsd, fill = class)) +
       geom_density(alpha = .3) + #alpha used for filling the density
       geom_vline(data = means, aes(xintercept = rmsd.mean, colour = class),
                  linetype = "longdash", size=1) + 
       ggtitle("decoys with RMSD between 1.5 and 3.5 (logit unfiltered, use RMSD = 3 as threshold)")

#