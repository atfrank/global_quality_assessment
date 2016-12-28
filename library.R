initialize_analysis <- function(pairs, correct_shifts=FALSE, workdir="~/GitSoftware/global_quality_assessment/"){
  # initialize data and variables need for analysis
  # goto analysis directory
  setwd(workdir)
  
  # data for composite NMR bundle
  compile_data(pairs, FALSE, FALSE, correct_shifts)
  
  # data for average structure
  compile_data(pairs, TRUE, FALSE, correct_shifts)
  
  # data for nmr-xray pairs
  compile_data(pairs, FALSE, TRUE, correct_shifts)
}

compile_data <- function(pairs, average_data=FALSE, nmr_xray=FALSE, corrected_shifts = FALSE, sensivities_flag=TRUE){
  # function to compile data used in this analysis
  # get weight files
  weight_larmord_1 <- read.table("data/larmord_accuracy_nucleus.txt",col.names = c("nucleus","weight_larmord_1"))
  weight_larmord_2 <- read.table("data/larmord_accuracy_resname_nucleus.txt",col.names = c("resname","nucleus","weight_larmord_2"))
  weight_ramsey_1 <- read.table("data/ramsey_accuracy_nucleus.txt",col.names = c("nucleus","weight_ramsey_1"))
  weight_ramsey_2 <- read.table("data/ramsey_accuracy_resname_nucleus.txt",col.names = c("resname","nucleus","weight_ramsey_2"))
  
  if (sensivities_flag){
    assign("sensivities",get(load("data/mean_sensivities.RData")))
    sensivities[,2] <- sensivities[,2]/max(sensivities[,2])
    sensivities[,3] <- sensivities[,3]/max(sensivities[,3])
    sensivities[,4] <- sensivities[,4]/max(sensivities[,4])
  }
  
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
        measured_file <- paste("reference_errors/",ref,"/chemical_shifts_corrected_larmord.txt",sep = "")
        
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
      
      # add sensivities
      if (sensivities_flag){
        colnames(sensivities) <- c("nucleus", "sensi_mean", "sensi_ramsey", "sensi_larmord")
        cs <- merge(cs, sensivities, by = c("nucleus"))
      }
      
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

get_slrs <- function(pairfile, nucleus_group = "both", prediction_method = "larmord", weight = 1, error_type = "rmse", conformational_averaging = FALSE, outliers =FALSE){
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
    } 
    if (weight==2){
      cs$weight <- 1/cs$weight_larmord_1
    } 
    if (c("sensi_larmord") %in% colnames(cs)){
      if (weight==3){
        cs$weight <- 1/sqrt(cs$sensi_larmord)*cs$weight_larmord_1
      } 
      if (weight==4){
        cs$weight <- 1/sqrt(cs$sensi_larmord)*cs$weight_larmord_2
      } 
      if (weight==5){
        cs$weight <- sqrt(cs$sensi_larmord)/cs$weight_larmord_1
      } 
      if (weight==6){
        cs$weight <- sqrt(cs$sensi_larmord)/cs$weight_larmord_2
      } 
      
    }
  }
  if (prediction_method=="ramsey"){
    cs$predCS <- cs$ramsey_predCS
    if (weight==1){
      cs$weight <- 1/cs$weight_ramsey_1
    } 
    if (weight==2){
      cs$weight <- 1/cs$weight_ramsey_1
    } 
    if (c("sensi_ramsey") %in% colnames(cs)){
      if (weight==3){
        cs$weight <- 1/sqrt(cs$sensi_ramsey)*cs$weight_ramsey_1
      } 
      if (weight==4){
        cs$weight <- 1/sqrt(cs$sensi_ramsey)*cs$weight_ramsey_2
      } 
      if (weight==5){
        cs$weight <- sqrt(cs$sensi_ramsey)/cs$weight_ramsey_1
      } 
      if (weight==6){
        cs$weight <- sqrt(cs$sensi_ramsey)/cs$weight_ramsey_2
      } 
      
    }
  }
  if (prediction_method=="mean"){
    cs$predCS <- (cs$ramsey_predCS + cs$larmord_predCS)/2
    if (weight==1){
      cs$weight <- 1/((cs$weight_ramsey_1 + cs$weight_larmord_1)/2)
    } 
    if(weight==2){
      cs$weight <- 1/((cs$weight_ramsey_2 + cs$weight_larmord_2)/2)
    }
    if (c("sensi_mean") %in% colnames(cs)){
      if (weight==3){
        cs$weight <- 1/sqrt(cs$sensi_mean)*((cs$weight_ramsey_1 + cs$weight_larmord_1)/2)
      } 
      if(weight==4){
        cs$weight <- 1/sqrt(cs$sensi_mean)*((cs$weight_ramsey_2 + cs$weight_larmord_2)/2)
      }
      if (weight==5){
        cs$weight <- sqrt(cs$sensi_mean)/((cs$weight_ramsey_1 + cs$weight_larmord_1)/2)
      } 
      if(weight==6){
        cs$weight <- sqrt(cs$sensi_mean)/((cs$weight_ramsey_2 + cs$weight_larmord_2)/2)
      }
    }
  }

  if (conformational_averaging){
    cs <- ddply(.dat=cs, .var=c("resid","resname","nucleus","expCS","weight","reference_flag","type"), .fun=function(x){mean(x$predCS)})
    cs$model <- cs$reference_flag
    cs$predCS <- cs$V1
  }
  
  # remove outliers
  if(outliers){cs <- remove_outliers(cs)}

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

remove_outliers <- function(cs, threshold=5.0){
  diff <- ddply(.data = cs, .variables = c("resid", "nucleus"), .fun = function(x){data.frame(diff=mean(x$weight*abs(x$expCS-x$predCS)))})
  outliers <- subset(diff, diff>threshold)
  outliers$tag <- paste(outliers$resid, outliers$nucleus, sep = ":")
  cs$tag <- paste(cs$resid, cs$nucleus, sep = ":")
  if (nrow(outliers)!=0){
    #print(outliers)
    cs <- cs[ !(cs$tag %in% outliers$tag), ]
  }
  return(cs)
}

make_error_barplot <- function(pairs=c("1R2P_2LPS", "2FRL_2M22", "2H2X_2M21", "2KFC_2L1V", "2N6Q_5KMZ"), nucleus_group = "both", prediction_method = "larmord", error_type = "rmse"){
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
    # add segment
    n <- length(errors$error)
    a <- sum(errors$flag==0)
    b <- a + 1
    x0s <- c(0, b+2)
    x1s <- c(a+2, n+2)
    # these are the y-coordinates for the horizontal lines
    # that you need to set to the desired values.
    y0s <- c(mean(errors$error[errors$flag==0]), mean(errors$error[errors$flag==1]))
    # add segments
    segments(x0 = x0s, x1 = x1s, y0 = y0s, col = "cyan", lwd="3")
  }
  dev.off()
}

generate_all_barplot <- function(){
  # make barplots like that figure 2
  for (nucleus_group in c("proton","carbon","both")){
    for (prediction_method in c("ramsey","mean","larmord")){
      for (error_type in c("mae")){
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
  m <- m[1:5,]
  cols <- c("orange","blue","red")
  data <- data.frame(nslr=c(m[1:5,1],m[1:5,2],m[1:5,3]))
  data$label <- as.factor(c(rep("A",5),rep("B",5),rep("C",5)))
  par(lwd=3,mfrow=c(2,2),mgp=c(1.6, 0.4, 0),tcl=-0.3,oma=c(0.1,0.2,0.2,0.2),cex=1.5)
  barplot(t(as.matrix(m)), beside = TRUE, col = cols, las=2, ylim=c(0,1),  names.arg=labels)
  abline(h=random_nslr(20,40),lwd="2",lty="dashed")
  abline(h=mean(as.matrix(m)[,1]),lwd="2",lty="dashed",col=cols[1])
  abline(h=mean(as.matrix(m)[,2]),lwd="2",lty="dashed",col=cols[2])
  abline(h=mean(as.matrix(m)[,3]),lwd="2",lty="dashed",col=cols[3])
  p <- boxplot(nslr~., data, ylim = c(0, 1),col=cols,names=c("proton","carbon","both"))
  abline(h=random_nslr(20,40),lwd="2",lty="dashed")
  dev.off()
}

make_table <- function(pair, predictor="larmord", average_data=FALSE, nmr_xray=FALSE, conformational_averaging=FALSE, weight=1, error_types=c("mae","rmse","r","tau","rho"), outliers = FALSE){
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
      out <- get_slrs(file, error_type = error_type, weight = weight, nucleus_group = nucleus_group, prediction_method = predictor, conformational_averaging=conformational_averaging, outliers = outliers)
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

summarize_tables <- function(predictor="larmord", error_type = "mae", averaged_data=FALSE, nmr_xray=FALSE, conformational_averaging=FALSE, names=c("2L94_1Z2J","1Z2J_2L94","2N82_2N7X","2N7X_2N82","1R2P_2LPS","2FRL_2M22","2H2X_2M21","2KFC_2L1V"), weight=1, outliers = FALSE){
  # function summarizes the results presented in the paper
  # this is the main driver function
  nslrs <- flags <- errors <- NULL
  error_types  <- c("mae","rmse","tau","r","rho","geo_mae")
  error_type_index <- which(error_types==error_type)

  for (i in seq_along(names)){
    m <- make_table(names[i],predictor,averaged_data, nmr_xray, conformational_averaging, weight, error_types, outliers)
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

summarize_observed_cs <- function(rna, nuclei = c("H1'","H2'","H3'","H4'","H5'","H5''","H2","H5","H6","H8","C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8")){
  require(plyr)
  tmp <- data.frame(nucleus = nuclei)
  
  # load measure shifts
  measured_names <- c("resname", "resid", "nucleus", "expCS", "expError")
  measured_file <- paste("data/measured_shifts_",rna,".dat",sep = "")
  expcs <- read.table(measured_file, col.names = measured_names, stringsAsFactors = FALSE)
  
  expcs$type <- "carbon"
  expcs$type[grepl("H", expcs$nucleus)] <- "proton"
  expcs$type[grepl("N", expcs$nucleus)] <- "nitrogen"
  
  counts <- ddply(.data = expcs, .variables = c("type","nucleus"), .fun = function(x){return(data.frame(N=nrow(x)))})
  counts <- counts[counts$nucleus %in% nuclei, ]
  counts$C <- sum(counts$N[counts$type == "carbon"])
  counts$H <- sum(counts$N[counts$type == "proton"])
  counts <- merge(tmp, counts, by = c("nucleus"), all = TRUE)
  rownames(counts) <- counts$nucleus
  counts <- counts[nuclei, ]
  return(counts)
}

get_observed_shifts_matrix <- function(rnas=c("2LPS", "2M22", "2M21", "2L1V", "5KMZ")){
  
  for (i in seq_along(rnas)){
    rna <- rnas[i]
    tmp <- summarize_observed_cs(rna)
    if (i==1){
      mat <- tmp$N
    } else {
      mat <- rbind(mat, tmp$N)
    }
  }
  colnames(mat) <- rownames(tmp)
  rownames(mat) <- rnas
  return(mat)
}

get_nucleus_errors <- function(pairfile, nucleus_group = "both", prediction_method = "larmord", weight = 1, error_type = "rmse", nuclei = c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8","H1'","H2'","H3'","H4'","H5'","H5''","H2","H5","H6","H8")){
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
    } 
    if (weight==2){
      cs$weight <- 1/cs$weight_larmord_1
    } 
    if (c("sensi_larmord") %in% colnames(cs)){
      if (weight==3){
        cs$weight <- 1/sqrt(cs$sensi_larmord)*cs$weight_larmord_1
      } 
      if (weight==4){
        cs$weight <- 1/sqrt(cs$sensi_larmord)*cs$weight_larmord_2
      } 
      if (weight==5){
        cs$weight <- sqrt(cs$sensi_larmord)/cs$weight_larmord_1
      } 
      if (weight==6){
        cs$weight <- sqrt(cs$sensi_larmord)/cs$weight_larmord_2
      } 
      
    }
  }
  if (prediction_method=="ramsey"){
    cs$predCS <- cs$ramsey_predCS
    if (weight==1){
      cs$weight <- 1/cs$weight_ramsey_1
    } 
    if (weight==2){
      cs$weight <- 1/cs$weight_ramsey_1
    } 
    if (c("sensi_ramsey") %in% colnames(cs)){
      if (weight==3){
        cs$weight <- 1/sqrt(cs$sensi_ramsey)*cs$weight_ramsey_1
      } 
      if (weight==4){
        cs$weight <- 1/sqrt(cs$sensi_ramsey)*cs$weight_ramsey_2
      } 
      if (weight==5){
        cs$weight <- sqrt(cs$sensi_ramsey)/cs$weight_ramsey_1
      } 
      if (weight==6){
        cs$weight <- sqrt(cs$sensi_ramsey)/cs$weight_ramsey_2
      } 
      
    }
  }
  if (prediction_method=="mean"){
    cs$predCS <- (cs$ramsey_predCS + cs$larmord_predCS)/2
    if (weight==1){
      cs$weight <- 1/((cs$weight_ramsey_1 + cs$weight_larmord_1)/2)
    } 
    if(weight==2){
      cs$weight <- 1/((cs$weight_ramsey_2 + cs$weight_larmord_2)/2)
    }
    if (c("sensi_mean") %in% colnames(cs)){
      if (weight==3){
        cs$weight <- 1/sqrt(cs$sensi_mean)*((cs$weight_ramsey_1 + cs$weight_larmord_1)/2)
      } 
      if(weight==4){
        cs$weight <- 1/sqrt(cs$sensi_mean)*((cs$weight_ramsey_2 + cs$weight_larmord_2)/2)
      }
      if (weight==5){
        cs$weight <- sqrt(cs$sensi_mean)/((cs$weight_ramsey_1 + cs$weight_larmord_1)/2)
      } 
      if(weight==6){
        cs$weight <- sqrt(cs$sensi_mean)/((cs$weight_ramsey_2 + cs$weight_larmord_2)/2)
      }
    }
  }
  
  if (error_type=="rmse"){
    errors <- ddply(.dat=cs, .var=c("model","reference_flag","nucleus"), .fun=score_rmse)
  }
  if (error_type=="mae") {
    errors <- ddply(.dat=cs, .var=c("nucleus", "model","reference_flag"), .fun=score_mae)
  }
  if (error_type=="geo_mae") {
    errors <- ddply(.dat=cs, .var=c("model","reference_flag","nucleus"), .fun=score_geo_mae)
  }
  if (error_type=="tau") {
    errors <- ddply(.dat=cs, .var=c("model","reference_flag","nucleus"), .fun=correlation_kendall)
  }
  if (error_type=="r") {
    errors <- ddply(.dat=cs, .var=c("model","reference_flag", "nucleus"), .fun=correlation_spearman)
  }
  if (error_type=="rho") {
    errors <- ddply(.dat=cs, .var=c("model","reference_flag","nucleus"), .fun=correlation_pearson)
  }
  
  tmp <- expand.grid(unique(errors$model), nuclei)
  colnames(tmp) <- c("model", "nucleus")
  
  errors <- merge(errors, tmp, all = TRUE)
  errors <- errors[order(errors$model,errors$nucleus), ]
  cnames <- subset(errors, model == 1)$nucleus
  mat <- matrix(errors$V1, byrow = TRUE, nrow=length(unique(errors$model)), ncol = length(nuclei))
  colnames(mat) <- cnames
  mat <- mat[ ,nuclei]
  
  model_info <- unique(errors[,c("model", "reference_flag")])
  model_info <- model_info[complete.cases(model_info),]
  mat <- cbind(model_info, mat)
  return(mat)
}

make_table_errors_nucleus <- function(pairs=c("1R2P_2LPS","2FRL_2M22","2H2X_2M21","2KFC_2L1V","2N6Q_5KMZ"), predictor="larmord", weight=1, error_type="mae", outliers = FALSE){
  # make summary tables of errors 
  for (i in seq_along(pairs)){
    pair <- pairs[i]
    file <- paste("chemical_shifts_",pair,".txt",sep="")
    out <- get_nucleus_errors(file, error_type = error_type, weight = weight, nucleus_group = "both", prediction_method = predictor)
    out$id <- pair
    if (i==1){
      error_matrix <- out
    } else {
      error_matrix <- rbind(error_matrix, out)
    }
  }
  return(error_matrix)
}

nuclei_importance <- function(pairs=c("1R2P_2LPS","2FRL_2M22","2H2X_2M21","2KFC_2L1V","2N6Q_5KMZ"), predictor="larmord", weight=1, error_type="mae"){
  # see: https://www.r-bloggers.com/variable-importance-plot-and-variable-selection/
  # see: http://freakonometrics.hypotheses.org/19458
  # see: 
  require(randomForest)
  t <- make_table_errors_nucleus(pairs, predictor, weight, error_type)
  data <- t[, !(colnames(t) %in% c("model", "id"))]
  data$reference_flag <- as.factor(data$reference_flag)
  names <- unlist(strsplit("reference_flag C1p C2p C3p C4p C5p C2 C5 C6 C8 H1p H2p H3p H4p H5p H5pp H2 H5 H6 H8", " "))
  colnames(data) <- names
  rf <- randomForest(formula = reference_flag~., data = data, na.action = na.exclude)
  varImpPlot(rf)
  return(rf)
}

plot_importance <- function(rf, figfile = "test.pdf"){
  tmp <- as.data.frame(rf$importance)
  tmp$nucleus <- c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8","H1'","H2'","H3'","H4'","H5'","H5''","H2","H5","H6","H8") 
  names_tmp <- tmp[order(tmp$MeanDecreaseGini, decreasing = F),]
  pdf(file = figfile, width = 4/2, height = 6/2)
  varImpPlot(rf, scale = TRUE, lcolor = "blue",  bg = "red", pt.cex = 0.8, labels = names_tmp$nucleus, main = "", cex = 0.5)  
  abline(h=15, lty="dashed", lwd = "2.5")
  dev.off()
}




