### Packages necessaires au script
library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(glmpca)
library(stringr)
library(vsn)


create_matrix_level<-function(UserInput,UserInputColNames,UserInputRelevel,UserInputIfEmpty,UserInputMostImportantAttribute)
{
  MatrixLevel<-UserInput#copy the User Input to have the structure
  for (i in 1:length(UserInput))
  {
    if (UserInputColNames[i]!=UserInputMostImportantAttribute)
    {
      MatrixLevel[[i]]<-list(unlist(MatrixLevel[[i]]),UserInputIfEmpty[i])#Add the if empty to the structure except if it is in the most important col
    }
    MatrixLevel[[i]]<-unlist(MatrixLevel[[i]])
    MatrixLevel[[i]]<-c(MatrixLevel[[i]][which(MatrixLevel[[i]]==UserInputRelevel[i])],MatrixLevel[[i]][-which(MatrixLevel[[i]]==UserInputRelevel[i])])#Put the relevel at first pos
  }
  return(MatrixLevel)
}