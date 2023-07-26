library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(glmpca)
library(stringr)
library(vsn)

releveling<-function(df,ShowLevel,Iter)
{
  df <- as.data.frame(lapply(df, unlist))
  print(Iter)
  ListSubset<-c()
  if (length(ShowLevel)<Iter)
  {
    return(list(df))
  }
  for (i in 1:(length(ShowLevel[[Iter]])-1)) 
  {
    print(c(i,"boucle for"))
    Iter<-strtoi(Iter)
    df[,Iter+2] <- relevel(factor(df[,Iter+2]),ShowLevel[[Iter]][i])
    ListSubset<-c(ListSubset,releveling(df,ShowLevel = ShowLevel,Iter = Iter+1))
  }
  return(ListSubset)
}