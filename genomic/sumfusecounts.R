chemin<-"C:/Users/maxence.lechemia/Documents/Génomique/Genomiccounts formated2/"
list_files <- list.files(path =chemin,full.names = T)
print(list_files)
for (i in list_files)
{
  print(i)
  list_files_into <- list.files(path =i,full.names = T)
  table1<-read.table(list_files_into[1])
  table2<-read.table(list_files_into[2])
  tablefinal<-table1
  tablefinal$V2<-table1$V2+table2$V2
  write.table(tablefinal,paste0(i,".txt"), sep = "  ", dec = ".",
            row.names = F, col.names = F)
}