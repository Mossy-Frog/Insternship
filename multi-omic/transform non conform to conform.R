chemin<-"C:/Users/maxence.lechemia/Documents/Multi-omics/proteo transformation v2/"
output<-"C:/Users/maxence.lechemia/Documents/Multi-omics/ordered formated fused genomic and all transcripto/"
file<-list.files(chemin,full.names = T)
file2<-list.files(chemin)
DfOfDf<- list()
count<-1
for (i in file)
{
  
  DfOfDf[[count]]<-as.data.frame(read.table(i))
  count<-count+1
}

ref<-read.table("C:/Users/maxence.lechemia/Documents/Multi-omics/formated fused genomics and new and alod transcripto/A2-1.txt")

for (i in 1:length(DfOfDf))#open each files et parcourt
{
  DfOfDf[[i]]<- merge(DfOfDf[[i]],ref, by="V1",all=T)#Merge la ref avec la table des counts à normaliser
  for (j in 1:length(DfOfDf[[i]][,1]))#Parcourt chaque cases
       {
         if(is.na(DfOfDf[[i]]$V2.x[j]))#Si une ligne de la ref n'est pas dans l'autre, rajoute la ligne avec un zéro
         {
           DfOfDf[[i]]$V2.x[j]<-0
         }
         if(is.na(DfOfDf[[i]]$V2.y[j]))#si la ligne n'existe pas dans la ref, on retire
          {
          DfOfDf[[i]]<-DfOfDf[[i]][-j,]
    }
  }
  DfOfDf[[i]]<-DfOfDf[[i]][,-3]#on retire la ref
}
# Vieille façon, non efficace
# print(i)
# for (j in 1:length(ref$V1))
# {
#   count<-0
#   for (k in 1:length(DfOfDf[[i]]$V1))
#     if (ref$V1[j] ==DfOfDf[[i]]$V1[k])
#     {
#       count<-1
#       # print(k)
#     }
#   if (count==0)
#   {
#     DfOfDf[[i]]<-rbind(DfOfDf[[i]],c(ref$V1[j],0)) 
#     # print(0)
#   }
# }




for (i in 1:length(DfOfDf))#On sauvegarde ça
{
  DfOfDf[[i]] <- DfOfDf[[i]][order(DfOfDf[[i]]),]
  write.table(DfOfDf[[i]][which(is.na(DfOfDf[[i]]$V1)==F),], paste0(output,file2[i]), append = T, sep = "  ", dec = ".",
              row.names = F, col.names = F)
}