load_basepair_interaction <- function(file){
  data <- read.table(file = file, header = F, fill = T, stringsAsFactors = F)
  data <- data[,c(1,3,4,7)]
  colnames(data) <- c("resid","resname","edge","orientation")
}

load_stack_interaction <- function(file){
  data <- read.table(file = file, header = F, fill = T, stringsAsFactors = F)
  data <- data[,c(1,4)]
  colnames(data) <- c("resid","type")
}

