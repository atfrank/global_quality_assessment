initialize_analysis <- function(workdir="~/GitSoftware/global_quality_assessment/"){
  # initialize data and variables need for analysis
  # goto analysis directory
  setwd(workdir)
  
  # data for composite NMR bundle
  compile_data(pairs, FALSE, FALSE)
  
  # data for average structure
  compile_data(pairs, TRUE, FALSE)
  
  # data for nmr-xray pairs
  compile_data(pairs, FALSE, TRUE)
  
  # merge free and bound chemical shifts to ensure we only use common chemical shifts
  merge_free_bound("2N7X", "2N82")
  merge_free_bound("1Z2J", "2L94")

}

compile_data <- function(pairs, average_data=FALSE, nmr_xray=FALSE, corrected_shifts = FALSE){
  # function to compile data used in this analysis
  # get weight files
  weight_larmord_1 <- read.table("data/larmord_accuracy_nucleus.txt",col.names = c("nucleus","weight_larmord_1"))
  weight_larmord_2 <- read.table("data/larmord_accuracy_resname_nucleus.txt",col.names = c("resname","nucleus","weight_larmord_2"))
  weight_ramsey_1 <- read.table("data/ramsey_accuracy_nucleus.txt",col.names = c("nucleus","weight_ramsey_1"))
  weight_ramsey_2 <- read.table("data/ramsey_accuracy_resname_nucleus.txt",col.names = c("resname","nucleus","weight_ramsey_2"))
  
  for (i in seq_along(pairs$pair)){
    pair <- pairs$pair[i]
    ref <- pairs$ref[i]
    pair_name <- pairs$pair_name[i]
    threshold <- pairs$threshold[i]
    
    # file names
    if (average_data){
      larmord_file <- paste("data/larmord_",pair,"_average.txt",sep = "")
      ramsey_file <- paste("data/ramsey_",pair,"_average.txt",sep = "")
      outfile <- paste("data/chemical_shifts_",pair_name,"_average.txt",sep = "")
    } 
    if (nmr_xray){
      larmord_file <- paste("data/larmord_",pair,"_nmr_xray.txt",sep = "")
      ramsey_file <- paste("data/ramsey_",pair,"_nmr_xray.txt",sep = "")
      outfile <- paste("data/chemical_shifts_",pair_name,"_nmr_xray.txt",sep = "")
    } 
    
    if (!average_data && !nmr_xray){
      larmord_file <- paste("data/larmord_",pair,".txt",sep = "")
      ramsey_file <- paste("data/ramsey_",pair,".txt",sep = "")
      outfile <- paste("data/chemical_shifts_",pair_name,".txt",sep = "")
    }
    
    if (ref %in% c("1Z2J","2L94","2N7X","2N82")){
      measured_file <- paste("data/measured_shifts_merged_",ref,".dat",sep = "")
    } else {
      if (corrected_shifts){
        measured_file <- paste("data/measured_shifts_corrected_clean_",ref,".dat",sep = "")
      } else {
        measured_file <- paste("data/measured_shifts_",ref,".dat",sep = "")
      }
    }
    
    if(file.exists(larmord_file) && file.exists(larmord_file)){
      # column names
      larmord_names = c("model", "resid", "resname", "nucleus", "larmord_predCS", "id" )
      ramsey_names = c("model", "resid", "resname", "nucleus", "ramsey_predCS", "id" )
      measured_names <- c("resname", "resid", "nucleus", "expCS", "expError")
      
      # load larmord
      cs_l <- read.table(larmord_file, col.names = larmord_names)
      cs_l$resname <- as.character(cs_l$resname)
      # load ramsey
      cs_r <- read.table(ramsey_file, col.names = ramsey_names)
      cs_r$resname <- as.character(cs_r$resname)
      cs_r$resname[cs_r$resname=="URI"] <- "URA"
      
      # load measure shifts
      expcs <- read.table(measured_file, col.names = measured_names)
      expcs$resname <- as.character(expcs$resname)
      
      # merge files
      cs <- merge(cs_l, cs_r)
      cs <- merge(cs, expcs)
      
      names <- c("resid", "resname", "nucleus", "model", "larmord_predCS",  "ramsey_predCS", "expCS")
      cs <- cs[,names]
      
      cs <- cs[order(c(cs$model, cs$resid, cs$nucleus)),]
      cs <- cs[complete.cases(cs),]
      
      # add reference flag
      cs$reference_flag <- 0
      if (threshold<0){
        cs$reference_flag[cs$model<=abs(threshold)] <- 1
      } else {
        cs$reference_flag[cs$model>threshold] <- 1
      }
      
      # specical handling if working with averaged data
      if(average_data){
        cs$reference_flag <- cs$model-1
      }
      
      # add type info
      cs$type <- "carbon"
      cs$type[grep("H",cs$nucleus)] <- "proton"
      
      # add weights
      cs <- merge(cs, weight_larmord_1)
      cs <- merge(cs, weight_larmord_2)
      cs <- merge(cs, weight_ramsey_1)
      cs <- merge(cs, weight_ramsey_2)
      
      # write out file
      write.table(cs, file=outfile, quote = F, col.names = T, row.names = F)
      cat(sprintf("done with %s\n", pair))
    }
  }
}

merge_free_bound <- function(free, bound, corrected_shifts=FALSE){
  # function to merge chemical shifts for the free vs. bound and bound vs. free analysis
  # ensures that the comparison is carried out using a common set of chemical shifts
  # such that the results are not sensitive to differences in the level of peak assignments in the free and bound states
  # merge and write out common chemical shifts
  
  # input chemical shifts file name
  if (corrected_shifts){
    freefile <- paste("data/measured_shifts_corrected_clean_",free,".dat",sep = "")
    boundfile <- paste("data/measured_shifts_corrected_clean_",bound,".dat",sep = "")
  } else {
    freefile <- paste("data/measured_shifts_",free,".dat",sep = "")
    boundfile <- paste("data/measured_shifts_",bound,".dat",sep = "")
  }
  
  # output chemical shifts file name
  freefile_out <- paste("data/measured_shifts_merged_",free,".dat",sep = "")
  boundfile_out <- paste("data/measured_shifts_merged_",bound,".dat",sep = "")
  
  # read in chemical shifts
  free <- read.table(freefile, col.names=c("resname", "resid", "nucleus", "free", "expError_f"))
  bound <- read.table(boundfile, col.names=c("resname", "resid", "nucleus", "bound", "expError_b"))  
  
  # merge
  cs <- merge(free, bound, by=c("resname", "resid", "nucleus"))
  cs$error <- 0
  
  # select out need frames
  free <- cs[,c("resname", "resid", "nucleus", "free", "error")]
  bound <- cs[,c("resname", "resid", "nucleus", "bound", "error")]
  
  # write out chemical shift file
  write.table(free, file = freefile_out, quote = F, col.names = F, row.names = F)
  write.table(bound, file = boundfile_out, quote = F, col.names = F, row.names = F)
}


get_cs <- function(pairfile, nucleus_group = "both", prediction_method = "larmord", weight = 2, error_type = "rmse", conformational_averaging = FALSE){
  # function to get chemical shift data for a pair of NMR RNA structures
  library(nmR)
  library(plyr)
  # read in shifts
  cs <- read.table(paste(pairfile, sep=""), header = T)
  outfile <- paste(paste("errors/",nucleus_group,sep=""), prediction_method, error_type, pairfile, sep="_")

  # select subset
  if (nucleus_group=="proton" || nucleus_group=="carbon"){
    cs <- subset(cs, type == nucleus_group)
  } else {
    nucleus_group <- "both"
  }

  # set weight
  if (prediction_method=="larmord"){
    cs$predCS <- cs$larmord_predCS
    if (weight==1){
      cs$weight <- 1/cs$weight_larmord_1
    } else{
      cs$weight <- 1/cs$weight_larmord_2
    }
  }

  if (prediction_method=="ramsey"){
    cs$predCS <- cs$ramsey_predCS
    if (weight==1){
      cs$weight <- 1/cs$weight_ramsey_1
    } else{
      cs$weight <- 1/cs$weight_ramsey_2
    }
  }

  if (prediction_method=="mean"){
    cs$predCS <- (cs$ramsey_predCS + cs$larmord_predCS)/2
    if (weight==1){
      cs$weight <- 1/((cs$weight_ramsey_1 + cs$weight_larmord_1)/2)
    } else{
      cs$weight <- 1/((cs$weight_ramsey_2 + cs$weight_larmord_2)/2)
    }
  }

  if (conformational_averaging){
		cs <- ddply(.dat=cs, .var=c("resid","resname","nucleus","expCS","weight","reference_flag","type"), .fun=function(x){mean(x$predCS)})
		cs$model <- cs$reference_flag
		cs$predCS <- cs$V1
  }
  return(cs)
}


get_slrs <- function(pairfile, nucleus_group = "both", prediction_method = "larmord", weight = 1, error_type = "rmse", conformational_averaging = FALSE){
  # function to compute the resolving scores (i.e., the NSLR) for one of the test (i.e., chemical shifts in th pairfile) present in the manuscript
  # 
  library(nmR)
  library(plyr)
  # read in shifts
  cs <- read.table(paste("data/",pairfile, sep=""), header = T)
  outfile <- paste(paste("errors/",nucleus_group,sep=""), prediction_method, error_type, pairfile, sep="_")

  # select subset
  if (nucleus_group=="proton" || nucleus_group=="carbon"){
    cs <- subset(cs, type == nucleus_group)
  } else {
    nucleus_group <- "both"
  }

  # set weight
  if (prediction_method=="larmord"){
    cs$predCS <- cs$larmord_predCS
    if (weight==1){
      cs$weight <- 1/cs$weight_larmord_1
    } else{
      cs$weight <- 1/cs$weight_larmord_2
    }
  }

  if (prediction_method=="ramsey"){
    cs$predCS <- cs$ramsey_predCS
    if (weight==1){
      cs$weight <- 1/cs$weight_ramsey_1
    } else{
      cs$weight <- 1/cs$weight_ramsey_2
    }
  }

  if (prediction_method=="mean"){
    cs$predCS <- (cs$ramsey_predCS + cs$larmord_predCS)/2
    if (weight==1){
      cs$weight <- 1/((cs$weight_ramsey_1 + cs$weight_larmord_1)/2)
    } else{
      cs$weight <- 1/((cs$weight_ramsey_2 + cs$weight_larmord_2)/2)
    }
  }

  if (conformational_averaging){
    cs <- ddply(.dat=cs, .var=c("resid","resname","nucleus","expCS","weight","reference_flag","type"), .fun=function(x){mean(x$predCS)})
    cs$model <- cs$reference_flag
    cs$predCS <- cs$V1
  }

  if (error_type=="rmse"){
    errors <- ddply(.dat=cs, .var=c("model","reference_flag"), .fun=score_rmse)
  }
  if (error_type=="mae") {
    errors <- ddply(.dat=cs, .var=c("model","reference_flag"), .fun=score_mae)
  }
  if (error_type=="geo_mae") {
    errors <- ddply(.dat=cs, .var=c("model","reference_flag"), .fun=score_geo_mae)
  }
  if (error_type=="tau") {
    tmp <- ddply(.dat=cs, .var=c("model","reference_flag", "type"), .fun=correlation_kendall)
    tmp$weight <- tmp$N/sum(tmp$N)
    errors <- ddply(.dat=tmp, .var=c("model","reference_flag"), .fun=function(x){sum(tmp$weight*x$cor)})
  }
  if (error_type=="r") {
    tmp <- ddply(.dat=cs, .var=c("model","reference_flag", "type"), .fun=correlation_pearson)
    tmp$weight <- tmp$N/sum(tmp$N)
    errors <- ddply(.dat=tmp, .var=c("model","reference_flag"), .fun=function(x){sum(tmp$weight*x$cor)})
  }
  if (error_type=="rho") {
    tmp <- ddply(.dat=cs, .var=c("model","reference_flag", "type"), .fun=correlation_spearman)
    tmp$weight <- tmp$N/sum(tmp$N)
    errors <- ddply(.dat=tmp, .var=c("model","reference_flag"), .fun=function(x){sum(tmp$weight*x$cor)})
  }

  colnames(errors) <- c("model","flag","error")
  write.table(errors,outfile,col.names=TRUE, quote=FALSE, row.names=FALSE)

  if (error_type %in% c("mae","rmse","geo_mae")){
    errors <- errors[order(errors$error, decreasing = FALSE),]
  } else {
    errors <- errors[order(errors$error, decreasing = TRUE),]
  }


  return(data.frame(nslr=nslr(errors$flag),flag=errors$flag[1], error=errors$error[1]))
}


make_error_barplot <- function(pairs=c("1R2P_2LPS","2L94_1Z2J", "2FRL_2M22","1Z2J_2L94", "2H2X_2M21","2N82_2N7X", "2KFC_2L1V", "2N7X_2N82"), nucleus_group = "both", prediction_method = "larmord", error_type = "rmse"){
  # make barplot shown in figure 2
  figfile <- paste(paste("figures/",nucleus_group, sep=""), prediction_method, error_type, "fig.pdf", sep="_")
  pdf(file = figfile, width = 20, height = 20)
  par(lwd=3,mfrow=c(4,2),mgp=c(1.6, 0.4, 0),tcl=-0.3,oma=c(0.1,0.2,0.2,0.2))
  for (pair in pairs){
    pairfile <- paste("chemical_shifts_",pair,".txt",sep="")
    infile <- paste(paste("errors/",nucleus_group,sep=""), prediction_method, error_type, pairfile, sep="_")
    errors <- read.table(infile, header = T)
    
    cols <- rep("white",nrow(errors))
    cols[errors$flag==0] <- "blue"
    
    border_cols <- rep("black",nrow(errors))
    border_cols[which.min(errors$error)] <- "red"
    barplot(errors$error, col = cols, names.arg=errors$model, xlim=c(1,50),border = border_cols, lwd="2", las=2, cex.names = 2.1, cex.axis = 2.1)
  }
  dev.off()
}

generate_all_barplot <- function(){
  # make barplots like that figure 2
  for (nucleus_group in c("proton","carbon","both")){
    for (prediction_method in c("ramsey","mean","larmord")){
      for (error_type in c("rmse","mae","tau")){
        make_error_barplot(nucleus_group = nucleus_group, prediction_method = prediction_method, error_type = error_type)
      }
    }
  }
}

make_error_scatterplot <- function(pairs=c("2KFC_2L1V"), nucleus_group = "both", prediction_method = "mean", weight = 2, error_type = "rmse", conformational_averaging = TRUE){
  # used to make correlation plots shown in the figure that demonstrates the use of chemical shift errors to assessed the quality of NMR structures
  figfile <- paste(paste("figures/correlation",nucleus_group, sep=""), prediction_method, error_type, "fig.pdf", sep="_")
  pdf(file = figfile, width = 5, height = 5)
  par(lwd=1,mfrow=c(2,2),mgp=c(1.6, 0.4, 0),tcl=-0.3,oma=c(0.1,0.2,0.2,0.2), pty="s",cex=0.6)
  for (pair in pairs){
    pairfile <- paste("data/chemical_shifts_",pair,".txt",sep="")
    cs <- get_cs(pairfile, nucleus_group, prediction_method, weight, error_type, conformational_averaging)
    print(head(cs$type))
    protons <- cs$type=="proton"&cs$reference_flag=="0"
    carbons <- cs$type=="carbon"&cs$reference_flag=="0"
    plot(cs$expCS[protons], cs$predCS[protons],xlab = "obs. chemical shifts (ppm)", ylab = "comp. chemical shifts(ppm)", pch=21, bg = "orange")
    abline(0,1, lty="dashed")
    plot(cs$expCS[carbons], cs$predCS[carbons],xlab = "obs. chemical shifts (ppm)", ylab = "comp. chemical shifts(ppm)", pch=21, bg = "blue")
    abline(0,1, lty="dashed")
    protons <- cs$type=="proton"&cs$reference_flag=="1"
    carbons <- cs$type=="carbon"&cs$reference_flag=="1"
    plot(cs$expCS[protons], cs$predCS[protons],xlab = "obs. chemical shifts (ppm)", ylab = "comp. chemical shifts(ppm)", pch=21, bg = "orange")
    abline(0,1, lty="dashed")
    plot(cs$expCS[carbons], cs$predCS[carbons],xlab = "obs. chemical shifts (ppm)", ylab = "comp. chemical shifts(ppm)", pch=21, bg = "blue")
    abline(0,1, lty="dashed")
  }
  dev.off()
}

make_nslr_plots <- function(m, labels=NULL, figfile="test.pdf"){
  # using to make the box-plots shown in figure 2
  pdf(file = figfile, width = 10, height = 10)
  m <- m[1:8,]
  cols <- c("orange","blue","red")
  data <- data.frame(nslr=c(m[1:8,1],m[1:8,2],m[1:8,3]))
  data$label <- as.factor(c(rep("A",8),rep("B",8),rep("C",8)))
  par(lwd=3,mfrow=c(2,2),mgp=c(1.6, 0.4, 0),tcl=-0.3,oma=c(0.1,0.2,0.2,0.2),cex=1.5)
  print(as.matrix(m))
  barplot(t(as.matrix(m)), beside = TRUE, col = cols, las=2, ylim=c(0,1),  names.arg=labels)
  abline(h=random_nslr(20,40),lwd="2",lty="dashed")
  abline(h=mean(as.matrix(m)[,1]),lwd="2",lty="dashed",col=cols[1])
  abline(h=mean(as.matrix(m)[,2]),lwd="2",lty="dashed",col=cols[2])
  abline(h=mean(as.matrix(m)[,3]),lwd="2",lty="dashed",col=cols[3])
  p <- boxplot(nslr~., data, ylim = c(0, 1),col=cols,names=c("proton","carbon","both"))
  abline(h=random_nslr(20,40),lwd="2",lty="dashed")
  dev.off()
}

make_table <- function(pair, predictor="larmord", average_data=FALSE, nmr_xray=FALSE, conformational_averaging=FALSE, weight=1, error_types=c("mae","rmse","r","tau","rho")){
  # make summary tables 
  # list contains: 1) NSLR -- resolving score, 2) Flags -- specifying whether the lowest error model was a reference model, 3) Lowest Errors -- error for the lowest error model in a given test
  # the three columns correspond to the results obtained when using 1H, 13C, and both 1H and 13C
  nslrs <- NULL
  flags <- NULL
  errors <- NULL
  nucleus_groups <- c("proton","carbon","both")
  for (error_type in error_types){
    for (nucleus_group in nucleus_groups){
      if (average_data){
        file <- paste("chemical_shifts_",pair,"_average.txt",sep="")
      }
      if (nmr_xray){
        file <- paste("chemical_shifts_",pair,"_nmr_xray.txt",sep="")
      }
      if(!average_data && !nmr_xray) {
        file <- paste("chemical_shifts_",pair,".txt",sep="")
      }
      out <- get_slrs(file, error_type = error_type, weight = weight, nucleus_group = nucleus_group, prediction_method = predictor, conformational_averaging=conformational_averaging)
      nslrs <- c(nslrs, out$nslr)
      flags <- c(flags, out$flag)
      errors <- c(errors, out$error)
    }
  }
  nslrs <- matrix(nslrs,ncol=length(nucleus_groups),byrow = T)
  flags <- matrix(flags,ncol=length(nucleus_groups),byrow = T)
  errors <- matrix(errors,ncol=length(nucleus_groups),byrow = T)
  return(list(nslrs, flags,errors))
}


summarize_tables <- function(predictor="larmord", error_type = "mae", averaged_data=FALSE, nmr_xray=FALSE, conformational_averaging=FALSE, names=c("2L94_1Z2J","1Z2J_2L94","2N82_2N7X","2N7X_2N82","1R2P_2LPS","2FRL_2M22","2H2X_2M21","2KFC_2L1V"), weight=1){
  # function summarizes the results presented in the paper
  # this is the main driver function
  nslrs <- flags <- errors <- NULL
  error_types  <- c("mae","rmse","tau","r","rho","geo_mae")
  error_type_index <- which(error_types==error_type)

  for (i in seq_along(names)){
    m <- make_table(names[i],predictor,averaged_data, nmr_xray, conformational_averaging, weight, error_types)
    nslrs <- c(nslrs, m[[1]][error_type_index,])
    flags <- c(flags,m[[2]][error_type_index,])
    errors <- c(errors,m[[3]][error_type_index,])
  }

  nslrs <- matrix(nslrs,ncol=ncol(m[[1]]),byrow = T)
  flags <- matrix(flags,ncol=ncol(m[[1]]),byrow = T)
  errors <- matrix(errors,ncol=ncol(m[[1]]),byrow = T)

  rownames(nslrs) <- names
  rownames(flags) <- names
  rownames(errors) <- names
  return(list(nslrs, flags,errors))
}

get_nslr_only <- function(summary){
  m <- summary[[1]]
  names <- c(rownames(m),"Mean")
  m <- rbind(m,colMeans(m))
  rownames(m) <- names
  return(round(m,2))
}

get_flags_only <- function(summary){
  m <- summary[[2]]
  names <- c(rownames(m),"Total", "Success-Rate")
  total <- colSums(m)
  success <- colMeans(m)
  m <- rbind(m,total)
  m <- round(rbind(m,success),2)
  rownames(m) <- names
  return(m)
}

get_errors_only <- function(summary){
  m <- summary[[3]]
  names <- c(rownames(m),"Mean")
  m <- rbind(m,colMeans(m))
  rownames(m) <- names
  return(round(m,2))
}

create_tpr_tables <- function(predictors, averaged_data=FALSE, nmr_xray=FALSE, conformational_averaging=FALSE, names=c("2L94_1Z2J","1Z2J_2L94","2N82_2N7X","2N7X_2N82","1R2P_2LPS","2FRL_2M22","2H2X_2M21","2KFC_2L1V")){
  m <- NULL
  for (predictor in predictors){
    table <- summarize_tables(predictor,"mae", averaged_data, nmr_xray, conformational_averaging, names)
    success_mae <- get_flags_only(table)['Success-Rate',]
    table <- summarize_tables(predictor,"rmse", averaged_data, nmr_xray, conformational_averaging, names)
    success_rmse <- get_flags_only(table)['Success-Rate',]
   m <- c(m, paste(success_mae, " (", success_rmse,")", sep=""))
  }
  return(matrix(m, ncol = length(predictors), byrow=F))
}

create_nslr_tables <- function(predictors, averaged_data=FALSE, nmr_xray=FALSE, conformational_averaging=FALSE, names=c("2L94_1Z2J","1Z2J_2L94","2N82_2N7X","2N7X_2N82","1R2P_2LPS","2FRL_2M22","2H2X_2M21","2KFC_2L1V")){
  m <- NULL
  for (predictor in predictors){
    table <- summarize_tables(predictor,"mae", averaged_data, nmr_xray, conformational_averaging, names)
    success_mae <- get_nslr_only(table)['Mean',]
    table <- summarize_tables(predictor,"rmse", averaged_data, nmr_xray, conformational_averaging, names)
    success_rmse <- get_nslr_only(table)['Mean',]
    m <- c(m, paste(success_mae, " (", success_rmse,")", sep=""))
  }
  return(matrix(m, ncol = length(predictors), byrow=F))
}

correlation_kendall <- function(x){
  #' Kendall Correlation Scoring Function
  #'
  #' This function computes the 1 - tau
  #' @param x input dataframe. Should contain field: predCS and expCS.
  #' @export
  #' @examples
  #' score_kendall(x)
  return(data.frame(cor=cor(x$expCS,x$predCS,method="kendall"),N=nrow(x)))
}

correlation_pearson <- function(x){
  #' Pearson Correlation Scoring Function
  #'
  #' This function computes the 1 - R
  #' @param x input dataframe. Should contain field: predCS and expCS.
  #' @export
  #' @examples
  #' score_pearson(x)
  return(data.frame(cor=cor(x$expCS,x$predCS,method="pearson"),N=nrow(x)))
}

correlation_spearman <- function(x){
  #' Spearman Correlation Scoring Function
  #'
  #' This function computes the 1 - rho
  #' @param x input dataframe. Should contain field: predCS and expCS.
  #' @export
  #' @examples
  #' score_spearman(x)
  return(data.frame(cor=cor(x$expCS,x$predCS,method="spearman"),N=nrow(x)))
}



