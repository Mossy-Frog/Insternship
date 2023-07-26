library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(glmpca)
library(stringr)
library(vsn)



row_to_subset<-function(df,row,UserInputMostImportantAttribute,UserInputRelevel,UserInputColNames,UserInputIfEmpty,UserInput)
  #Prend en valeur la table des attributs, la ligne de la formule intéressante dans la table des possibilités,
  #L''attribut le plus important selon l'utilisateur, les level les plus importants, les noms des colonnes d'attributs,
  #la liste des éléments si vides et l'userinput.
  
  #Retourne une liste contenant d'une part une liste contenant tous les subset optimaux, et d'autre part une liste avec le nom des subset optimaux pour la formule
  #Rappel: une formule peut utiliser plusieurs subsets et donner des résultats différents.
{
  
  TempDf<-data.frame(df)#création des subset
  ListDf<-list(TempDf)#Création de la liste des subsets, avec le subset non traité
  NameDf<-c("Untreated")#création de la liste des noms avec le nom du premier df
  if (toString(row[1,"Pos1"])==UserInputMostImportantAttribute&row["Pos2"]=="")#Premier cas: si il s'agit juste d'une formule avec 1 seul élément et que l'élément est le plus important
  {
    colnamestemp<-UserInputColNames[UserInputColNames != UserInputMostImportantAttribute]#On créé une liste des colonnes à trier, sans la colonne des attributs qui contient l'attribut le plus important
    for (i in colnamestemp)#On parcourt la liste créé
    {
      TempDf<- TempDf[which(TempDf[,i]==UserInputRelevel[which(UserInputColNames==i)]),]#On retire les lignes dont les colonnes ne contiennent pas les controles
    }
    ListDf<-append(ListDf,list(TempDf))#On ajoute le subset créé
    NameDf<-append(NameDf,"")#ON nomme le subset
  }
  else if(toString(row[1,"Pos1"])!=UserInputMostImportantAttribute&row["Pos2"]=="")#Deuxième cas, formule simple à 1 élément, mais ce n'est pas l'attribut le plus important
  {
    TypesOfMainAttribute<-unlist(UserInput[which(UserInputColNames==UserInputMostImportantAttribute)])#On créé une liste qui contient les différents éléments de l'attribut le plus important
    TypesOfMainAttribute<-append(TypesOfMainAttribute,UserInputIfEmpty[which(UserInputColNames==UserInputMostImportantAttribute)])#on rajoute l'élément si dans le nom du fichier l'emplacement est vide
    TypesOfMainAttribute<-unique(TypesOfMainAttribute)#On vire les duplicats
    colnamestemp<-UserInputColNames[UserInputColNames != toString(UserInputMostImportantAttribute)]#On créé une liste des colonnes à trier, sans la colonne des attributs qui contient l'attribut le plus important
    colnamestemp<-colnamestemp[colnamestemp != toString(row[1,"Pos1"])]#On retire aussi l'élement dans la formule que nous avons donc pas envie de trier
    
    for (i in TypesOfMainAttribute)#Parcourt liste des différents éléments de l'attribut
    {
      TempDf<-data.frame(df)#on créé le subset
      TempDf<- TempDf[which(TempDf[,UserInputMostImportantAttribute]==i),]#On ne met dans le subset qu'un seul élément de l'attribut important et on vire les autres
      for (j in colnamestemp) #On parcourt la liste des éléments à trier
      {
        TempDf<- TempDf[which(TempDf[,j]==UserInputRelevel[which(UserInputColNames==j)]),]#Dans les éléments à trier on ne garde que les controles
      }
      # print(TempDf)
      ListDf<-append(ListDf,list(TempDf))#On ajoute le subset créé
      NameDf<-append(NameDf,i)#On ajute le nom du subset, qui est donc le nom de l'élément de l'attribut choisi
    }
  }
  else if((toString(row[1,"Pos1"])==UserInputMostImportantAttribute|toString(row[2,"Pos1"])==UserInputMostImportantAttribute)&row["Pos2"]!="")#3ème cas:S'il y a 2 attributs dans la formule, dont l'attribut important
  {
    colnamestemp<-UserInputColNames[UserInputColNames != toString(row[1,"Pos1"])]#On créé la liste des attributs à trier en enlevant les attributs de la formule
    colnamestemp<-colnamestemp[colnamestemp != toString(row[1,"Pos2"])]#on enlève le deuxième attribut contenu dans la formule
    for (i in colnamestemp)#On parcourt la liste des attributs à trier
    {
      TempDf<- TempDf[which(TempDf[,i]==UserInputRelevel[which(UserInputColNames==i)]),]#On retire tous les éléments qui ne sont pas des controles dans la liste à trier
    }
    ListDf<-append(ListDf,list(TempDf))#On ajoute le subset
    NameDf<-append(NameDf,"")#On ajoute le nom du subset
  }
  else if (toString(row[1,"Pos1"])!= UserInputMostImportantAttribute & row["Pos2"] != "" & toString(row[1,"Pos2"]) != UserInputMostImportantAttribute & toString(row[1,"Pos1"]) != "")#4ème cas:
    #2 attributs dans la formule mais aucun ne fait partie de l'attribut important
  {
    TypesOfMainAttribute<-unlist(UserInput[which(UserInputColNames==UserInputMostImportantAttribute)])#On créé une liste qui contient les différents éléments de l'attribut le plus important
    TypesOfMainAttribute<-append(TypesOfMainAttribute,UserInputIfEmpty[which(UserInputColNames==UserInputMostImportantAttribute)])
    TypesOfMainAttribute<-unique(TypesOfMainAttribute)
    colnamestemp<-UserInputColNames[UserInputColNames != toString(row[1,"Pos1"])]#On créé la liste des attributs à trier en enlevant les attributs de la formule
    colnamestemp<-colnamestemp[colnamestemp != toString(row[1,"Pos2"])]#on enlève le deuxième attribut contenu dans la formule
    for (i in TypesOfMainAttribute)#Parcourt liste des différents éléments de l'attribut
    {
      TempDf<-data.frame(df)#Création du subset
      TempDf<- TempDf[which(TempDf[,UserInputMostImportantAttribute]==i),]#On vire les éléments de l'attribut autre que celui qu'on a dans notre itération
      for (j in colnamestemp)#On parcourt la liste des attributs à trier
      {
        if (j!=UserInputMostImportantAttribute)
        {
          TempDf<- TempDf[which(TempDf[,j]==UserInputRelevel[which(UserInputColNames==j)]),]#On retire tous les éléments qui ne sont pas des controles dans la liste à trier
        }
      }
      # print(TempDf)
      ListDf<-append(ListDf,list(TempDf))#On ajoute le subset
      NameDf<-append(NameDf,i)#On ajoute le nom du subset
    }
    
  }
  
  return(c(list(ListDf),list(NameDf)))#On retourne une liste qui contient les 2 listes
}