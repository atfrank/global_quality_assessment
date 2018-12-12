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
}