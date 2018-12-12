sort_inf_files <- function(rnas, inf_path = "INF/", rewrite = FALSE){
  for(rna in rnas){
    inf = read.table(paste0(inf_path,rna,"_inf.csv"),header=T,sep=",")
    inf$id = gsub(paste0(rna,"_*"),"",inf$fn)
    inf$id = as.numeric(gsub("*.pdb.outCR","",inf$id))
    inf = inf[order(inf$id),]
    if(rewrite){
      write.table(inf[,c("id","inf_all")], file = paste0(inf_path, rna, "_inf_sorted.csv"),
                  quote = F, col.names = T, row.names = F)
    }
  }
}
