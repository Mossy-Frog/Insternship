library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(glmpca)
library(stringr)
library(vsn)


create_possibility_table<-function(UserInputNames)#Create the possibility table, prend en input les noms des colonnes des attributs et sort la table des possibilités POUR LES FORMULEs
{
  DfPossibility<-data.frame("Pos1"=c(""),
                            "Pos2"=c(""),
                            "Pos3"=c(""),
                            "Formula"=c(""))#Create the structure
  
  for ( e in 1:length(UserInputNames))#Parcourt une première fois le nom des cols
  {
    for ( f in 1:length(UserInputNames))#Parcourt une deuxième fois en asynchrone 
    {
      if (e==f)
      {
        TempRow<-c(UserInputNames[e],"","",paste0("~",UserInputNames[e]))#Formule simple
        DfPossibility<-rbind(DfPossibility,TempRow)
      }
      else
      {
        TempRow<-c(UserInputNames[e],UserInputNames[f],"",paste0("~",UserInputNames[e],"+",UserInputNames[f]))
        DfPossibility<-rbind(DfPossibility,TempRow)#formules complexes
        TempRow<-c(UserInputNames[e],UserInputNames[f],paste0(UserInputNames[e],":",UserInputNames[f]),paste0("~",UserInputNames[e],"+",UserInputNames[f],"+",UserInputNames[e],":",UserInputNames[f]))
        DfPossibility<-rbind(DfPossibility,TempRow)
      }
    }
  }
  
  DfPossibility<-DfPossibility[-1,]#on enlève la première ligne créée lors de la création de la structure 
  rownames(DfPossibility) = seq(length=nrow(DfPossibility))#on remet les lignes bien
  return(DfPossibility)
}
