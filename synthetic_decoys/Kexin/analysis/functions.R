load_model_accuracy <- function(file = "loo_acc_3.0.txt", colnames = c("id", "sen", "spec", "overall"), header = TRUE){
  data_test <- read.table(file = file, header = T)
  group1 = get(load(file = "larmord_train_rnas.RData"))
  group2 = get(load(file = "not_in_larmord.RData"))
  data_test$type <- "not_in_larmord"
  data_test$type[data_test$rna %in% group1] <- "in_larmord"
  data_test$rna <- as.character(data_test$rna)
  return(data_test)
}

plot_color_bar <- function(data, metric = "spec"){
  require(tidyverse)
  data$group <- data$type
  data$group <- as.factor(data$group)
  data$value <- as.numeric(as.character(data[, metric]))
  data$individual <- data$rna
  data$rna <- 1:nrow(data)
  data$colour <- "red"
  data$colour[data$group == "not_in_larmord"] <- "blue"
  # Set a number of 'empty bar' to add at the end of each group
  empty_bar=10
  to_add = data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
  colnames(to_add) = colnames(data)
  to_add$group=rep(levels(data$group), each=empty_bar)
  data=rbind(data, to_add)
  data = data %>% arrange(group, value)
  data$rna=seq(1, nrow(data))
  # Get the name and the y position of each label
  label_data=data
  number_of_bar=nrow(label_data)
  angle= 90 - 360 * (label_data$rna-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust <- ifelse( angle < -90, 1, 0)
  label_data$angle<-ifelse(angle < -90, angle+180, angle)
  
  # prepare a data frame for base lines
  base_data=data %>% 
    group_by(group) %>% 
    summarize(start=min(rna), end=max(rna) - empty_bar) %>% 
    rowwise() %>% 
    mutate(title=mean(c(start, end)))
  base_data$mean <- mean(data$value[data$group == "not_in_larmord"], na.rm = TRUE)
  base_data$mean[base_data$group == "in_larmord"] <- mean(data$value[data$group == "in_larmord"], na.rm = TRUE)
  
  print(base_data)
  # prepare a data frame for grid (scales)
  grid_data = base_data
  grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
  grid_data$start = grid_data$start - 1
  grid_data=grid_data[-1,]
  
  # Make the plot
  p = ggplot(data, aes(x=as.factor(rna), y=value, fill=data$colour)) +       # Note that rna is a factor. If x is numeric, there is some space between the first bar
    
    geom_bar(aes(x=as.factor(rna), y=value, fill=data$colour), stat="identity", alpha=0.5) +
    
    # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
    
    geom_segment(data=grid_data, aes(x = end, y = 0.80, xend = start, yend = 0.80), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 0.60, xend = start, yend = 0.60), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 0.40, xend = start, yend = 0.40), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 0.20, xend = start, yend = 0.20), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 0.00, xend = start, yend = 0.00), colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    
    # Add text showing the value of each 100/75/50/25 lines
    annotate("text", x = rep(max(data$rna), 6), y = c(0.00, 0.20, 0.40, 0.60, 0.80, 1.00), label = c("0.00", "0.20", "0.40", "0.60", "0.80", "1.00") , color="black", size=3.5 , angle=0, fontface="bold", hjust=1) +
    
    geom_bar(aes(x=as.factor(rna), y=value, fill=group, colour = "black"), stat="identity", colour = "black") +
    ylim(-0.5,1.5) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(0,4), "cm") 
    ) +
    coord_polar() +
    geom_text(data=label_data, aes(x=rna, y=value+0.1, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
    
    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = 0.5, xend = end, yend = 0.5 ), colour = "red", alpha=0.8, size=0.5 , inherit.aes = FALSE, linetype= 2) +
    geom_segment(data=base_data, aes(x = start, y = mean, xend = end, yend = mean ), colour = "white", alpha=0.8, size=0.5 , inherit.aes = FALSE, linetype= 2) +
    geom_segment(data=base_data, aes(x = start, y = 0, xend = end, yend = 0 ), colour = "black", alpha=0.8, size=0.7 , inherit.aes = FALSE )
  p
}

sort_inf_files <- function(rnas, inf_path = "INF/", rewrite = FALSE){
  for(rna in rnas){
    inf = read.table(paste0(inf_path,rna,"_inf.csv"),header=T,sep=",")
    inf$model = gsub(paste0(rna,"_*"),"",inf$fn)
    inf$model = as.numeric(gsub("*.pdb.outCR","",inf$model))
    inf = inf[order(inf$model),]
    inf$id = rna
    if(rewrite){
      write.table(inf[,c("id","model","inf_all")], file = paste0(inf_path, rna, "_inf_sorted.csv"),
                  quote = F, col.names = T, row.names = F)
    }
  }
}

calc_DI <- function(rnas, inf_path){
  errors = read.table("errors_normalized.txt", header = T)
  errors = errors[order(errors$id, errors$model),]
  inf_all = NULL
  for(rna in rnas){
    inf = read.table(paste0(inf_path, rna, "_inf_sorted.csv"), header = T)
    if(is.null(inf_all)){
      inf_all = inf
    } else {
      inf_all = rbind(inf_all, inf)
    }
  }
  errors_merged = merge(errors, inf_all, by = c("id","model"))
  errors_merged$DI = errors_merged$rmsd/errors_merged$inf_all
  return(errors_merged)
}

calc_rmsd_of_selected_structures <- function(rnas, file_path, errors){
  for(rna in rnas){
    cs_errors = subset(errors, id==rna)
    cs_errors = cs_errors[,c("id","model","rmsd")]
    pred = read.table(paste0(file_path, rna,"_3.0.txt"), header = T)
    merged = merge(cs_errors, pred, by=c("id","model"))
    selected = merged[which.max(merged$score), "rmsd"]
  }
  return(selected)
}