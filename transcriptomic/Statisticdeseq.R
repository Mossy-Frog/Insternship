library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(glmpca)
library(stringr)
library(vsn)

statistic_deseq<-function(df,chemin,formule,UserInputPvalueThreshold,UserInputMostImportantAttribute,output)
  #Prend en valeur la table des attributs avec les facteurs, le chemin des fichiers, la formule, et la valeur de pvalue minimale
{
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = df,
                                         directory = chemin,
                                         design = as.formula(formule))#Create the DESeq2 table
  keep <- rowSums(counts(ddsHTSeq)>=15)>=5#filter it 
  dds_filtered<-ddsHTSeq[keep,]
  DdsDES<- DESeq(ddsHTSeq)#change to suitable format
  
  
  ResName<-resultsNames(DdsDES) #avec shrinkage
  
  nameLevel<-""
  for (i in 1:length(df))
  {
    print(levels(df[,i])[1])
    if (is.null(levels(df[,i])[1])==F)
    {
      nameLevel<-paste0(nameLevel,toString(levels(df[,i])[1]),"_") 
    }
  }
  MostImportantAttributeName<-""
  for (i in unique(df[UserInputMostImportantAttribute]))
  {
    MostImportantAttributeName<-paste(MostImportantAttributeName,toString(i))
  }
  for (j in 2:length(ResName))#pour toutes les comparaisons possibles dans la formule.
  {
    res <- results(DdsDES,
                   name = ResName[j],
                   lfcThreshold = 0.5,
                   parallel = T)
    res <- subset(res,padj<UserInputPvalueThreshold)
    # PlotDes<-plotMA(res,ylim=c(-5,5)) 
    # print(PlotDes)
    x <- "a1~!@#$%^&*(){}_+:\"<>?,./;'[]-="
    
    
    write.csv(res,file=paste(output,nameLevel,MostImportantAttributeName,str_replace_all(formule, "[[:punct:]]", "_"),str_replace_all(ResName[j], "[[:punct:]]", " "),".csv", sep= " "))
  }
  return(ResName)
}