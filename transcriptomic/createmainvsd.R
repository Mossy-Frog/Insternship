library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(glmpca)
library(stringr)
library(vsn)

create_main_vsd<-function(df,chemin,UserInputMostImportantAttribute)
{
  FormulaDES<-paste0("~", UserInputMostImportantAttribute)#Create the general formula for DESeq2
  
  
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = df,
                                         directory = chemin,
                                         design = as.formula(FormulaDES))#Create the DESeq2 table
  
  keep <- rowSums(counts(ddsHTSeq)>=15)>=5#keep only the genes with high enough counts
  dds_filtered<-ddsHTSeq[keep,]
  
  if (length(dds_filtered<1000))
  {
    vsd<- varianceStabilizingTransformation(dds_filtered,blind = F)#Create the vsd
    
  }else
  {
    vsd<-vst(dds_filtered,blind = F)#Create the vsd
  }
  return(list(vsd,dds_filtered))
}

