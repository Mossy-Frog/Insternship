
library(data.table)
library("dplyr")
library(stringr)
ProtDf<-fread("C:/Users/maxence.lechemia/Documents/Proteomique/MS19-014_Table S2_Burkho-K96243.csv")


ProtDf<-as.data.frame(ProtDf)

NotreDbAnnotee<- read.csv("C:/Users/maxence.lechemia/Documents/Proteomique/ProteinIDtoLocusTag.csv")

ProtDf<-ProtDf[which(is.na(ProtDf$accession)==F),]

for (i in 1:length(ProtDf[,1]))
{
  print(i)
  for (j in 1:length(NotreDbAnnotee[,"Name.x"]))
  {

    if (ProtDf[i,1]==NotreDbAnnotee[j,"Name.x"] & (is.na(NotreDbAnnotee[j,"Name.x"])==F))
    {
      ProtDf[i,1]<-NotreDbAnnotee[j,"locus_tag"]
      ProtDf<-ProtDf[which(is.na(ProtDf$accession)==F),]
    }
  }
}


ProtDf<-ProtDf[-288,]


for (i in 16:length(ProtDf))
{
  tempDF<-ProtDf[,c(1,i)]
  write.table(tempDF, paste0("C:/Users/maxence.lechemia/Documents/Proteomique/1-raw/",colnames(tempDF)[2],".txt"), append = T, sep = "  ", dec = ".",
                row.names = F, col.names = F)
}


  