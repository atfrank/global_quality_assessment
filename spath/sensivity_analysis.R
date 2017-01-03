combine_features_shifts <- function(rna_1="1R2P", rna_2="2LPS",include_these=NULL){
  # column names for feature file
  names <- unlist(strsplit("type name nucleus model id resid dexpCS RG RA RC RU RT H1p.0 H2p.1 H3p.2 H4p.3 H5p.4 H5pp.5 H2.6 H5.7 H6.8 H8.9 C1p.10 C2p.11 C3p.12 C4p.13 C5p.14 C2.15 C5.16 C6.17 C8.18 randomCoil hbond_b_b hbond_b_bk hbond_bk_bk ma po cc alpha beta gamma delta epsilon zeta chi nu0 nu1 nu2 nu3 nu4 stacking rc rc_po rc_hm", " "))
  # read in features
  features <- read.table(paste("features_", rna_1, "_", rna_2, ".txt", sep = ""), col.names = names)
  names <- unlist(strsplit("type name nucleus model id resid dexpCS cc rc ma po hbond_b_b hbond_b_bk hbond_bk_bk stacking alpha beta gamma delta epsilon zeta chi nu0 nu1 nu2 nu3 nu4", " "))
  features <- features[,names]
  # adjust model number in features to match that in the predicted chemical shift files
  features$model <- features$model + 1
  # make hbond binary
  features$hbond_b_b <- as.numeric(features$hbond_b_b > 0)
  features$hbond_b_bk <- as.numeric(features$hbond_b_bk > 0)
  features$hbond_bk_bk <- as.numeric(features$hbond_bk_bk > 0)
  # read in larmord and ramsey
  predicted_larmord <- read.table(paste("scratch/larmord_", rna_1, "_", rna_2, ".txt", sep = ""), col.names = c("model","resid", "resname", "nucleus", "predCS_larmord", "id") )
  predicted_ramsey <- read.table(paste("scratch/ramsey_", rna_1, "_", rna_2, ".txt", sep = ""), col.names = c("model","resid", "resname", "nucleus", "predCS_ramsey", "id") )
  # merge with with predicted shifts
  merge_names <- c("model", "resid", "nucleus", "id")
  structure_features <- merge(features, predicted_larmord, by = merge_names)
  merge_names <- c("model", "resid", "nucleus", "id", "resname")
  structure_features <- merge(structure_features, predicted_ramsey, by = merge_names)
  structure_features$predCS_consensus <- (structure_features$predCS_larmord+structure_features$predCS_ramsey)/2
  # remove unnecessary columns
  # do analysis 
  ignore_me <- c("type name model id resid dexpCS randomCoil RG RA RC RU RT H1p.0 H2p.1 H3p.2 H4p.3 H5p.4 H5pp.5 H2.6 H5.7 H6.8 H8.9 C1p.10 C2p.11 C3p.12 C4p.13 C5p.14 C2.15 C5.16 C6.17 C8.18")
  ignore_me <- unlist(strsplit(ignore_me, " "))
  ignore_me <- ignore_me[!(ignore_me %in% include_these)]
  return(list(consensus=structure_features[, !(colnames(structure_features) %in% c(ignore_me, "predCS_ramsey","predCS_larmord"))], larmord=structure_features[, !(colnames(structure_features) %in% c(ignore_me, "predCS_ramsey","predCS_consensus"))], ramsey=structure_features[, !(colnames(structure_features) %in% c(ignore_me, "predCS_larmord", "predCS_consensus"))]))
}

combine_features_shifts_all <- function(include_these=NULL){
  rnas <- data.frame(rna_1=unlist(strsplit("1R2P 2FRL 2H2X 2KFC 2N6Q", " ")), rna_2=unlist(strsplit("2LPS 2M22 2M21 2L1V 5KMZ ", " ")))
  for (i in seq_along(rnas$rna_1)){
    l <- combine_features_shifts(rnas$rna_1[i], rnas$rna_2[i], include_these)
    if (i==1){
      larmord <- l$larmord
      ramsey <- l$ramsey
      consensus <- l$consensus
    } else {
      larmord <- rbind(larmord, l$larmord)
      ramsey <- rbind(ramsey, l$ramsey)
      consensus <- rbind(consensus, l$consensus)
    }
  }
  names(consensus)[names(consensus) == 'predCS_consensus'] <- 'predCS'
  names(ramsey)[names(ramsey) == 'predCS_ramsey'] <- 'predCS'
  names(larmord)[names(larmord) == 'predCS_larmord'] <- 'predCS'
  return(list(larmord=larmord, ramsey=ramsey, consensus=consensus))
}

get_importance <- function(f, nuc="C1'", method="larmord", verbose = FALSE){
  if (method=="larmord"){tmp <- subset(f$larmord, nucleus==nuc)}
  if (method=="ramsey"){tmp <- subset(f$ramsey, nucleus==nuc)}
  if (method=="consensus"){tmp <- subset(f$consensus, nucleus==nuc)}
  imp <- compute_importance(tmp)
  imp$names <- rownames(imp)
  return(imp)
}

get_importance_by_group <- function(f, group="proton", method="larmord", verbose = FALSE){
  require(randomForest)
  require(plyr)
  if (method=="larmord"){tmp <- f$larmord}
  if (method=="ramsey"){tmp <- f$ramsey}
  if (method=="consensus"){tmp <- f$consensus}
  
  if (group == "carbon"){nuclei <- c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8")}
  if (group == "proton"){nuclei <- c("H1'","H2'","H3'","H4'","H5'","H5''","H2","H5","H6","H8")}
  if (group == "carbon_sugar"){nuclei <- c("C1'","C2'","C3'","C4'","C5'")}
  if (group == "proton_sugar"){nuclei <- c("H1'","H2'","H3'","H4'","H5'","H5''")}
  if (group == "carbon_base"){nuclei <- c("C2","C5","C6","C8")}
  if (group == "proton_base"){nuclei <- c("H2","H5","H6","H8")}
  tmp <- tmp[(tmp$nucleus %in% nuclei),]
  print(unique(tmp$nucleus))
  imp <- compute_importance(tmp, verbose)
  return(imp)
}

get_importance_by_group_all <- function(f, methods = c("consensus","larmord","ramsey"), groups = c("carbon", "proton", "carbon_sugar", "proton_sugar", "carbon_base", "proton_base"), verbose = TRUE){
  for (method in methods){
    for (group in groups){
      rfname <- paste("predictive_model_", group, "_", method, ".RData", sep = "")
      rf <- get_importance_by_group(f, group, method, verbose)
      save(rf, file = rfname)
    }
  }
}

compute_importance <- function(tmp, ntree = 100, verbose=FALSE){
  require(randomForest)
  tmp <- tmp[,!(colnames(tmp) %in% c("nucleus","type"))]
  rf <- randomForest(predCS~., data=tmp, importance=TRUE, do.trace=verbose, ntree = ntree)
  return(rf)
}

get_importance_by_nuclei <- function(f, nuclei="C1'", method="larmord", ntree = 100, verbose = FALSE){
  require(randomForest)
  require(plyr)
  if (method=="larmord"){tmp <- f$larmord}
  if (method=="ramsey"){tmp <- f$ramsey}
  if (method=="consensus"){tmp <- f$consensus}
  tmp <- tmp[(tmp$nucleus %in% nuclei),]
  print(unique(tmp$nucleus))
  tmp <- tmp[,!(colnames(tmp) %in% c("nucleus", "resname", "type"))]
  imp <- compute_importance(tmp, ntree, verbose)
  return(imp)
}

get_importance_by_nuclei_all <- function(f, methods = c("consensus","larmord","ramsey"), nuclei = c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8","H1'","H2'","H3'","H4'","H5'","H5''","H2","H5","H6","H8"), ntree = 100, verbose = TRUE){
  for (method in methods){
    for (nucleus in nuclei){
      rfname <- paste("predictive_model_", nucleus, "_", method, ".RData", sep = "")
      rf <- get_importance_by_nuclei(f, nucleus, method, ntree, verbose)
      save(rf, file = rfname)
    }
  }
}

summarize_importance_by_nuclei_all <- function(methods = c("consensus","larmord","ramsey"), nuclei = c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8","H1'","H2'","H3'","H4'","H5'","H5''","H2","H5","H6","H8"), ntree = 100, verbose = TRUE){
  for (method in methods){
    for (nucleus in nuclei){
      rfname <- paste("predictive_model_", nucleus, "_", method, ".RData", sep = "")
      load(rfname)
      tmp <- as.matrix(rf$importance[,1])
      tmp <- tmp/max(tmp)
      ifelse(exists("mat"), mat <- cbind(mat, tmp), mat <- tmp) 
      names <- rownames(mat)
    }
    #print(round(mat, 3))
    mat <- matrix(as.numeric(mat), ncol = ncol(mat), byrow = FALSE)
    #print(round(mat, 3))
    means <- rowMeans(mat)
    
    ifelse(exists("all_mat"), all_mat <- cbind(all_mat, means), all_mat <- means) 
    rm("mat")
  }
  rownames(all_mat) <- names
  colnames(all_mat) <- methods
  return(round(all_mat, 3))
}

summarize_importance_by_group <- function(group="proton"){
  if (group == "carbon"){nuclei <- c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8")}
  if (group == "proton"){nuclei <- c("H1'","H2'","H3'","H4'","H5'","H5''","H2","H5","H6","H8")}
  if (group == "carbon_sugar"){nuclei <- c("C1'","C2'","C3'","C4'","C5'")}
  if (group == "proton_sugar"){nuclei <- c("H1'","H2'","H3'","H4'","H5'","H5''")}
  if (group == "carbon_base"){nuclei <- c("C2","C5","C6","C8")}
  if (group == "proton_base"){nuclei <- c("H2","H5","H6","H8")}
  imp <- summarize_importance_by_nuclei_all(nuclei = nuclei)
  for(i in 1:ncol(imp)){imp[,i] <- imp[,i]/max(imp[,i])}
  return(round(imp, 2))
}





get_mutual_info <- function(f, nuc="C1'", method="larmord"){
  require(infotheo)
  if (method=="larmord"){tmp <- subset(f$larmord, nucleus==nuc)}
  if (method=="ramsey"){tmp <- subset(f$ramsey, nucleus==nuc)}
  if (method=="consensus"){tmp <- subset(f$consensus, nucleus==nuc)}
  names <- colnames(tmp)
  names <- names[!(names %in% c("predCS","nucleus"))]
  inf <- NULL
  for (name in names){
    inf <- c(inf, mutinformation(X=as.factor(tmp[,name]), Y=as.factor(tmp[,"predCS"])))
  }
  return(data.frame(variables=names, mutual_information=inf, nucleus=nuc))
}

make_sensivity_plots <- function(imp, figfile="tmp.pdf", width = 0.75*9.430556, height = 0.75*6.055556){
  pdf(file = paste("c_",figfile,sep = ""), width = width, height = height)
  par(lwd=0.2,mfrow=c(1,1),mgp=c(1.6, 0.4, 0),tcl=-0.3,oma=c(0.1,0.2,0.2,0.2))
  mat <- subset(imp, nucleus=="C1'")$`%IncMSE`/max(subset(imp, nucleus=="C1'")$`%IncMSE`)
  mat <- cbind(mat, subset(imp, nucleus=="C2'")$`%IncMSE`/max(subset(imp, nucleus=="C2'")$`%IncMSE`))
  mat <- cbind(mat, subset(imp, nucleus=="C3'")$`%IncMSE`/max(subset(imp, nucleus=="C3'")$`%IncMSE`))
  mat <- cbind(mat, subset(imp, nucleus=="C4'")$`%IncMSE`/max(subset(imp, nucleus=="C4'")$`%IncMSE`))
  mat <- cbind(mat, subset(imp, nucleus=="C5'")$`%IncMSE`/max(subset(imp, nucleus=="C5'")$`%IncMSE`))
  mat <- cbind(mat, subset(imp, nucleus=="C2")$`%IncMSE`/max(subset(imp, nucleus=="C2")$`%IncMSE`))
  mat <- cbind(mat, subset(imp, nucleus=="C5")$`%IncMSE`/max(subset(imp, nucleus=="C5")$`%IncMSE`))
  mat <- cbind(mat, subset(imp, nucleus=="C6")$`%IncMSE`/max(subset(imp, nucleus=="C6")$`%IncMSE`))
  mat <- cbind(mat, subset(imp, nucleus=="C8")$`%IncMSE`/max(subset(imp, nucleus=="C8")$`%IncMSE`))
  mat <- mat/max(mat)
  names <- names_carbons <- c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8")
  l_names <- c("hydrogen bonding","α","β","γ","δ", "ε","ζ","ν0","ν1","ν2","ν3","ν4","stacking","ring current")
  barplot(sqrt(mat), beside = T, names.arg = names, col=rainbow(n=nrow(mat)), ylab="normalized sensivity", xlab="nucleus", cex.axis = 1.0, cex.names = 1.0, ylim=c(0, 1.0))
  # add segment
  n <- ncol(mat)
  # width of each boxplot is 0.8
  m <- nrow(mat)+1
  x0s <- seq(1,n*m,m)
  x1s <- seq(m,n*m,m)
  # these are the y-coordinates for the horizontal lines
  # that you need to set to the desired values.
  y0s <- y_carbons <- colMeans(mat)
  # add segments
  segments(x0 = x0s, x1 = x1s, y0 = y0s, col = "black", lwd="1")
  dev.off()
  
  pdf(file = paste("p_",figfile,sep = ""), width = width, height = height)
  par(lwd=0.2,mfrow=c(1,1),mgp=c(1.6, 0.4, 0),tcl=-0.3,oma=c(0.1,0.2,0.2,0.2))
  mat <- subset(imp, nucleus=="H1'")$`%IncMSE`/max(subset(imp, nucleus=="H1'")$`%IncMSE`)
  mat <- cbind(mat, subset(imp, nucleus=="H2'")$`%IncMSE`/max(subset(imp, nucleus=="H2'")$`%IncMSE`))
  mat <- cbind(mat, subset(imp, nucleus=="H3'")$`%IncMSE`/max(subset(imp, nucleus=="H3'")$`%IncMSE`))
  mat <- cbind(mat, subset(imp, nucleus=="H4'")$`%IncMSE`/max(subset(imp, nucleus=="H4'")$`%IncMSE`))
  mat <- cbind(mat, subset(imp, nucleus=="H5'")$`%IncMSE`/max(subset(imp, nucleus=="H5'")$`%IncMSE`))
  mat <- cbind(mat, subset(imp, nucleus=="H5''")$`%IncMSE`/max(subset(imp, nucleus=="H5'")$`%IncMSE`))
  mat <- cbind(mat, subset(imp, nucleus=="H2")$`%IncMSE`/max(subset(imp, nucleus=="H2")$`%IncMSE`))
  mat <- cbind(mat, subset(imp, nucleus=="H5")$`%IncMSE`/max(subset(imp, nucleus=="H5")$`%IncMSE`))
  mat <- cbind(mat, subset(imp, nucleus=="H6")$`%IncMSE`/max(subset(imp, nucleus=="H6")$`%IncMSE`))
  mat <- cbind(mat, subset(imp, nucleus=="H8")$`%IncMSE`/max(subset(imp, nucleus=="H8")$`%IncMSE`))
  mat <- mat/max(mat)
  names <- names_protons <- c("H1'","H2'","H3'","H4'","H5'","H5''","H2","H5","H6","H8")
  l_names <- c("hydrogen bonding","magnetic anisotropy","bond polarization","α","β","γ","δ", "ε","ζ","ν0","ν1","ν2","ν3","ν4","stacking","ring current")
  barplot(sqrt(mat), beside = T, names.arg = names, col=rainbow(n=nrow(mat)), ylab="normalized sensivity", xlab="nucleus", cex.axis = 1.0, cex.names = 1.0, ylim=c(0, 1.0))
  # add segment
  n <- ncol(mat)
  # width of each boxplot is 0.8
  m <- nrow(mat)+1
  x0s <- seq(1,n*m,m)
  x1s <- seq(m,n*m,m)
  # these are the y-coordinates for the horizontal lines
  # that you need to set to the desired values.
  y0s <- y_protons <- colMeans(mat)
  
  # add segments
  segments(x0 = x0s, x1 = x1s, y0 = y0s, col = "black", lwd="1")
  dev.off()
  return(data.frame(names=c(names_carbons, names_protons), sensivities=c(y_carbons, y_protons)))
}

make_sensivity_plots_all <- function(sensivity="sensitivity.RData"){
  load(sensivity)
  s_l <- make_sensivity_plots(imp_l, figfile = "figure_larmord.pdf")
  s_r <- make_sensivity_plots(imp_r, figfile = "figure_ramsey.pdf")
  s_c <- make_sensivity_plots(imp_c, figfile = "figure_consensus.pdf")
  s <- cbind(cbind(s_c, s_r[,2]),s_l[,2])
  colnames(s) <- c("nucleus", "consensus", "ramsey", "larmord")
  s[,-1] <- round(s[,-1],3)
  save(s, file = "mean_sensivities.RData")
  return(s)
}

get_fluctuations <- function(x, feature="rc", circular_data = FALSE){
  require(circular)
  if (circular_data){
    return(data.frame(sd=circular::sd.circular(x[,feature])))
  } else {
    return(data.frame(sd=sd(x[,feature])))
  }
}


get_difference <- function(x, feature="rc", circular_data = FALSE){
  endpoints <- c(min(x$model), max(x$model))
  x <- x[x$model %in% endpoints, ]
  require(circular)
  if (circular_data){
    diff <- round(1 - cos(rad(x[1,feature]) - rad(x[2,feature])),3)
  } else {
    diff <- abs(x[1,feature] - x[2,feature])
  }
  return(data.frame(df=diff))
}
