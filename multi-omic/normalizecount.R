
chemin<-"C:/Users/maxence.lechemia/Documents/Multi-omics/formated fused genomics and new and alod transcripto/"
output_scale<-"C:/Users/maxence.lechemia/Documents/Multi-omics/formated fused genomics and new and old transcriptomic normalized/"
  
  


if (!dir.exists(output_scale)){
  dir.create(output_scale)
} else {
  print("Dir already exists!")
}  

ToNormalizeFiles<-list.files(chemin)

for (i in ToNormalizeFiles)
{
  print(i)
  a<-read.table(paste0(chemin,i))
  a$V2<-scale(a$V2)
  write.table(a, paste0(output_scale,i), append = T, sep = "  ", dec = ".",
              row.names = F, col.names = F)
}


print(paste0(output_scale,i))
