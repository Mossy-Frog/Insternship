library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(glmpca)
library(stringr)
library(vsn)

#Création du dataframe utilisable pour DESEQ2
create_df_deseq<- function(chemin,UserInput,UserInputColNames,UserInputIfEmpty, UserInputRelevel)
{
  sampleFiles <- list.files(chemin)#[which(sampleFiles==1)]
  sampleFiles <-str_remove(sampleFiles,".txt") #On retire le txt pour ""parser"" la liste de fichiers txt et les mettre dans la colonne
  ColuNames<-c()
  ColuNames<-append(ColuNames,"Filename")
  for (i in 1:length(UserInputColNames))#créer autant de colonnes attributs qu'il n'y a des différents Attributs possibles donnés par l'utilisateurs
  {
    ColuNames<-append(ColuNames,UserInputColNames[i])
  }
  
  ColuNames<-append(ColuNames,"Number")#traçabilité des numéro de txt
  # print(ColuNames)
  df<- data.frame(c(sampleFiles))
  colnames(df)[1]<-"Samplename"
  for (i in ColuNames)
  {
    df[i]<-c("")  
  }
  
  
  
  for (j in 1:length(sampleFiles))#on parcourt la liste de txt
  {
    tempvar<-toString(sampleFiles[j])#on transforme le nom des txt en string, pour être sur pour le strsplit
    templist <- unlist(strsplit(tempvar,"-",fixed = TRUE))#strsplit pour extraire chaque élément entre les -
    df[[length(colnames(df))]][[j]]<-templist[length(templist)]#le dernier élément -le numéro- est automatiquement mis dans les lignes et colonnes qui correspondent
    templist <- templist[1:length(templist)-1]#puisqu'on a déjà traité le numéro, on peut le virer
    for (k in 1:length(templist))#on parcourt la liste qui composait le nom du fichier
    {
      for (l in 1:length(UserInput))#on parcourt la liste de liste donnée par l'utilisateur
      {
        for (m in 1:length(UserInput[[l]]))#on parcourt la liste de liste donnée par l'utilisateur
        {
          # print(UserInput[[l]][[m]])
          # print(sampleFiles[j])
          if (UserInput[[l]][[m]]==templist[k])#Si l'élément du nom du fichier correspond à ce que l'utilisateur a donné comme attribut, on met cet élément dans la bonne colonne et ligne de notre matrice
          {
            df[[l+2]][[j]]<-templist[k]#plus É parce qu'on est obligé de mettre le samplename et filename
          }
        }
      }
    } 
  }
  for (i in 1:length(df$Samplename))#on parcourt le nom des fichiers
  {
    df$Samplename[i]<-list(sampleFiles[i])#on met le nom du fichier sans le txt dans le sample file
    df$Filename[i]<-paste0(list(sampleFiles[i]),".txt")#on met le nom du fichier avec le txt dans le sample file
  }
  
  ###Check les valeurs vides dans le tableau et les remplace par les valeurs données par le user dans ce cas
  countemptycol<-1#Compteur du userinputifempty
  for (i in 3:(length(df)-1))
  {
    for (j in 1:length(df[,i]))
    {
      if (is.na(df[j,i])|df[j,i]=="")
        df[j,i]<-UserInputIfEmpty[countemptycol]
    }
    countemptycol<-countemptycol+1
  }
  df<-df[,-length(df)]
  df<-as.data.frame(df)
  for (i in 1:length(colnames(df)))
  {
    # if (grepl("Attribute",colnames(df)[i],fixed=TRUE))
    # {
    temp<-colnames(df)[i]
    # print("condition remplie,factorisation de la colonne")
    colnames(df)[i]<-factor(colnames(df)[i])#change de format la colonne
    colnames(df)[i]<-temp
    # }
  }
  # View(df)
  #Relevel les colonnes en fonction de ce que le USER a mis
  countlevelcol<-1
  for (i in 3:(length(df)))#Part de la première colonne des attributs à la dernière des attributs
  {
    # print(df[,i])
    # print(UserInputRelevel[countlevelcol])
    df[,i] <- relevel(factor(df[,i]),UserInputRelevel[[countlevelcol]])#Permet de mettre le facteur sur l'élément important dans les colonnes
    countlevelcol<-countlevelcol+1
  }
  # print(df)
  
  return(df)
}


# df<-create_df_deseq(chemin = "C:/Users/maxence.lechemia/Documents/Proteomique/1-raw", UserInput = list(list("A","B","K","pDUC")),UserInputColNames = c("strain"),UserInputIfEmpty = c("Souche Inconnue"),UserInputRelevel= c("K"))



