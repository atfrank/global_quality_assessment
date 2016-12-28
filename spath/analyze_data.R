get_cs <- function(rna_1="1R2P", rna_2="2LPS"){
  
  cs_file_ramsey <- paste("scratch/ramsey_",rna_1,"_",rna_2,".txt", sep="")
  cs_file_larmord <- paste("scratch/larmord_",rna_1,"_",rna_2,".txt", sep="")
  structure_file <- paste("scratch/structure_info_",rna_1,"_",rna_2,".txt", sep="")
  struct_names <- c("model","rmsd","TM","GDT")
  
  cs_ramsey <- load_cs_data(cs_file_ramsey, names = c("model","resid","resname","nucleus","predCS_ramsey","id"), accuracyFile = "scratch/ramsey_accuracy_resname_nucleus.txt", atomBasedWeights = TRUE)
  cs_ramsey$weight_ramsey <- cs_ramsey$weight
  cs_larmord <- load_cs_data(cs_file_larmord, names = c("model","resid","resname","nucleus","predCS_larmord","id"), accuracyFile = "scratch/larmord_accuracy_resname_nucleus.txt", atomBasedWeights = TRUE)
  cs_larmord$weight_larmord <- cs_larmord$weight
  expcs <- read.table(paste("scratch/measured_shifts_",rna_2,".dat", sep=""), col.names = c("resname", "resid", "nucleus", "expCS","error"))
  cs <- merge(expcs, cs_larmord)
  cs <- merge(cs, cs_ramsey, by=c("model","resid","nucleus","resname","type","id"))

  ## add reference flag from struct file
  struct <- read.table(structure_file, col.names = struct_names)
  struct$model <- struct$model
  cs <- merge(cs, struct, by=c("model"))
  cs$reference_flag <- 0
  cs$reference_flag[cs$GDT>0.70] <- 1
  
  cs <- cs[complete.cases(cs),]
  return(cs)
}

get_errors <- function(cs, nucleus_group = "both", prediction_method = "larmord", error_type = "mae", save_data = FALSE){
  require(nmR)
  require(plyr)

  # select subset
  if (nucleus_group=="proton" || nucleus_group=="carbon"){
    cs <- subset(cs, type == nucleus_group)
  } else {
    nucleus_group <- "both"
  }
  
  # LARMORD
  if (prediction_method=="larmord"){
    cs$predCS <- cs$predCS_larmord
    cs$weight <- 1/(as.numeric(cs$weight_larmord))
  }
  
  # RAMSEY
  if (prediction_method=="ramsey"){
    cs$predCS <- cs$predCS_ramsey
    cs$weight <- 1/(as.numeric(cs$weight_ramsey))
  }
  
  # CONSENSUS
  if (prediction_method=="consensus"){
    cs$predCS <- (cs$predCS_larmord+cs$predCS_ramsey)/2
    cs$weight <- 1/((as.numeric(cs$weight_larmord)+as.numeric(cs$weight_ramsey))/2)
  }

  if (error_type=="rmse"){
    errors <- ddply(.dat=cs, .var=c("model","reference_flag","rmsd","GDT"), .fun=score_rmse)
  }
  if (error_type=="mae") {
    errors <- ddply(.dat=cs, .var=c("model","reference_flag","rmsd","GDT"), .fun=score_mae)
  }
  if (error_type=="geo_mae") {
    errors <- ddply(.dat=cs, .var=c("model","reference_flag","rmsd","GDT"), .fun=score_geo_mae)
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
  
  colnames(errors) <- c("model","flag","rmsd","GDT","error")
  
  if (save_data){
    #write.table(errors,outfile,col.names=TRUE, quote=FALSE, row.names=FALSE)
    #save(errors, file=outfile)
  }
  
  #if (error_type %in% c("mae","rmse","geo_mae")){
  #  errors <- errors[order(errors$error, decreasing = FALSE),]
  #} else {
  #  errors <- errors[order(errors$error, decreasing = TRUE),]
  #}
  return(errors)
}

plot_trace_bins <- function(cs, width = 0.5){
  require(plyr)
  require(Hmisc)
  errors <- get_errors(cs, prediction_method = "larmord",error_type = "rmse")
  data <- errors[,c("GDT","error")]
  errors <- get_errors(cs, prediction_method = "ramsey",error_type = "rmse")
  data <- cbind(data, errors[,c("error")])
  errors <- get_errors(cs, prediction_method = "consensus",error_type = "rmse")
  data <- cbind(data, errors[,c("error")])
  
  min_y <- min(as.vector(data[,-1]))
  max_y <- max(as.vector(data[,-1]))
  
  min_x <- min(as.vector(data[,1]))
  max_x <- max(as.vector(data[,1]))
  nbreaks <- abs(min_x - max_x)/width
  data <- cbind(data,cut(data[,1], nbreaks, labels = FALSE))
  colnames(data) <- c("rmsd","larmord","ramsey","consensus","label")
  r_data <- ddply(.data = data, .variables = c("label"), .fun = function(x){data.frame(rmsd_mean=mean(x$rmsd), rmsd_sd=sd(x$rmsd), larmord_mean=mean(x$larmord), larmord_sd=sd(x$larmord), ramsey_mean=mean(x$ramsey), ramsey_sd=sd(x$ramsey), consensus_mean=mean(x$consensus), consensus_sd=sd(x$consensus))})

  d <- r_data[,c("rmsd_mean","larmord_mean","larmord_sd")]
  colnames(d) <- c("x","y","sd")
  u_lwd <- 2.5
  plot(d$x, d$y, type="n", ylim=c(min_y,max_y), lwd=u_lwd, cex.axis = 1.5, xlab = "", ylab = "")
  with (data = d, expr = errbar(x, y, y+sd, y-sd, add=T, pch=21, cap=.015, errbar.col="black", bg = "black", lwd=u_lwd))

  d <- r_data[,c("rmsd_mean","ramsey_mean","ramsey_sd")]
  colnames(d) <- c("x","y","sd")
  points(d$x, d$y, type="n")
  with (data = d, expr = errbar(x, y, y+sd, y-sd, add=T, pch=21, cap=.015, errbar.col="blue", bg = "blue", lwd=u_lwd))

  d <- r_data[,c("rmsd_mean","consensus_mean","consensus_sd")]
  colnames(d) <- c("x","y","sd")
  points(d$x, d$y, type="n")
  with (data = d, expr = errbar(x, y, y+sd, y-sd, add=T, pch=21, cap=.015, errbar.col="red", bg = "red", lwd=u_lwd))
}


plot_trace <- function(cs){
  errors <- get_errors(cs, prediction_method = "larmord")
  data <- errors[,c("rmsd","error")]
  errors <- get_errors(cs, prediction_method = "ramsey")
  data <- cbind(data, errors[,c("error")])
  errors <- get_errors(cs, prediction_method = "consensus")
  data <- cbind(data, errors[,c("error")])
  
  min_y <- min(as.vector(data[,-1]))
  max_y <- max(as.vector(data[,-1]))

  min_x <- min(as.vector(data[,1]))
  max_x <- max(as.vector(data[,1]))
  nbreaks <- abs(min_x - max_x)/0.25
  
  plot(data[,1], data[,2], type="l", col="black", ylim=c(min_y,max_y),lwd="2")
  errors <- get_errors(cs, prediction_method = "ramsey")
  lines(data[,1], data[,3], type="l", col="blue",lwd="2")
  errors <- get_errors(cs, prediction_method = "consensus")
  lines(data[,1], data[,4], type="l", col="red",lwd="2")
  data <- cbind(data,cut(data[,1], nbreaks, labels = FALSE))
  colnames(data) <- c("rmsd","larmord","ramsey","consensus","label")
  r_data <- ddply(.data = data, .variables = c("label"), .fun = function(x){data.frame(rmsd_mean=mean(x$rmsd), rmsd_sd=sd(x$rmsd), larmord_mean=mean(x$larmord), larmord_sd=sd(x$larmord), ramsey_mean=mean(x$ramsey), ramsey_sd=sd(x$ramsey), consensus_mean=mean(x$consensus), consensus_sd=sd(x$consensus))})
  return(list(data,r_data))
}

get_errors_residuewise <- function(cs, nucleus_group = "both", prediction_method = "larmord", error_type = "mae", save_data = FALSE){
  require(nmR)
  require(plyr)
  
  # select subset
  if (nucleus_group=="proton" || nucleus_group=="carbon"){
    cs <- subset(cs, type == nucleus_group)
  } else {
    nucleus_group <- "both"
  }
  
  # LARMORD
  if (prediction_method=="larmord"){
    cs$predCS <- cs$predCS_larmord
    cs$weight <- 1/(as.numeric(cs$weight_larmord))
  }
  
  # RAMSEY
  if (prediction_method=="ramsey"){
    cs$predCS <- cs$predCS_ramsey
    cs$weight <- 1/(as.numeric(cs$weight_ramsey))
  }
  
  # CONSENSUS
  if (prediction_method=="consensus"){
    cs$predCS <- (cs$predCS_larmord+cs$predCS_ramsey)/2
    cs$weight <- 1/((as.numeric(cs$weight_larmord)+as.numeric(cs$weight_ramsey))/2)
  }
  
  if (error_type=="mae") {
    errors <- ddply(.dat=cs, .var=c("model","resid","reference_flag","rmsd","GDT"), .fun=score_mae)
  }

  colnames(errors) <- c("model","resid","flag","rmsd","GDT","error")
  
  if (save_data){
    #write.table(errors,outfile,col.names=TRUE, quote=FALSE, row.names=FALSE)
    #save(errors, file=outfile)
  }
  
  #if (error_type %in% c("mae","rmse","geo_mae")){
  #  errors <- errors[order(errors$error, decreasing = FALSE),]
  #} else {
  #  errors <- errors[order(errors$error, decreasing = TRUE),]
  #}
  return(errors)
}

sensitivity <- function(rna_1 = "2FRL",  rna_2 = "2M22"){
  cors_s <- cors_k <- cors_p <- NULL
  cs <- get_cs(rna_1 = rna_1, rna_2 =  rna_2 )
  nuclei <- c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8","H1'","H2'","H3'","H4'","H5'","H5''","H2","H5","H6","H8","proton_baseCarbon" ,"sugar_baseCarbon","baseProton_sugarCarbon","sugarProton_baseCarbon" ,"base" ,"sugar" ,"carbon" ,"proton" ,"baseCarbon" ,"baseProton" ,"sugarCarbon" ,"sugarProton")
  for (nuc in nuclei){
    compare <- read.table(paste("scratch/", tolower(rna_1), "_", tolower(rna_2), "_compare.txt", sep = ""), col.names = c("i","id","resid","local_rmsd"))
    tmp_cs <- subset(cs, model==1)
    tmp_cs <- get_cs_subset(tmp_cs, nuc)
    tmp_errors <- get_errors_residuewise(tmp_cs, nucleus_group = "both", prediction_method = "larmord")
    data <- merge(tmp_errors, compare, by = c("resid"))
    plot(data$resid, data$local_rmsd, type = "h", lwd = "2.5", ylim = c(-4,4))
    lines(data$resid, -5*data$error, type = "h", lwd = "2.5", col = "red")
    # plot_trace(cs)
    cors_s <- c(cors_s, cor(data$local_rmsd,data$error, method = "spearman"))
    cors_k <- c(cors_k, cor(data$local_rmsd,data$error, method = "kendall"))
    cors_p <- c(cors_p, cor(data$local_rmsd,data$error, method = "pearson"))
  }
  cors <- data.frame(nuclei,cors_s, cors_k, cors_p)
  print(tail(cors[order(abs(cors$cors_p)),],5))
}

plot_sensitivity <- function(rna_1 = "2FRL",  rna_2 = "2M22", nuc){
  cors_s <- cors_k <- cors_p <- NULL
  cs <- get_cs(rna_1 = rna_1, rna_2 =  rna_2 )
  compare <- read.table(paste("scratch/", tolower(rna_1), "_", tolower(rna_2), "_compare.txt", sep = ""), col.names = c("i","id","resid","local_rmsd"))
  tmp_cs <- subset(cs, model==1)
  tmp_cs <- get_cs_subset(tmp_cs, nuc)
  tmp_errors <- get_errors_residuewise(tmp_cs, nucleus_group = "both", prediction_method = "larmord")
  data <- merge(tmp_errors, compare, by = c("resid"))
  plot(data$resid, data$local_rmsd, type = "h", lwd = "2.5", ylim = c(-4,4))
  lines(data$resid, -5*data$error, type = "h", lwd = "2.5", col = "red")
}

get_trace <- function(rna_1 = "2FRL",  rna_2 = "2M22", width = 0.25){
  cors_s <- cors_k <- cors_p <- NULL
  cs <- get_cs(rna_1 = rna_1, rna_2 =  rna_2 )
  data <- plot_trace_bins(cs, width)
  return(data)
}


