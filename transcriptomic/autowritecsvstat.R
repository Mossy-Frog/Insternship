#Create statistical tables with DESEQ2 
#Auteur: Maxence LECHEMIA


#Initialisation des libraires nécessaires
library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(glmpca)
library(stringr)
library(vsn)


#Initialisation des fonctions nécessaires, il faudra peut être retrouver les chemins corrects.

source("C:/Users/maxence.lechemia/Documents/scripts formalisés/createdfdeseq.R")
source("C:/Users/maxence.lechemia/Documents/scripts formalisés/createpossibilitytable.R")
source("C:/Users/maxence.lechemia/Documents/scripts formalisés/releveling.R")
source("C:/Users/maxence.lechemia/Documents/scripts formalisés/rowtosubset.R")
source("C:/Users/maxence.lechemia/Documents/scripts formalisés/Statisticdeseq.R")
source("C:/Users/maxence.lechemia/Documents/scripts formalisés/creatematrixlevel.R")


#Initialisation de la fonction, Prend en compte dans l'ordre :
# - Une liste de liste contenant les différents attributs de l'échantillon
# - Une liste c() contenant le nom des types d'attributs différents
# - une liste c() contenant le nom de l'attribut si l'algorithme ne trouve pas l'attribut dans le nom
# - une liste c() contenant l'attribut le plus important pour la factorisation pour chaque type d'attribut
# - une liste c() contenant l'index des samples pouvant être problématiques lors de l'analyse, ex : contenant des valeurs nulles ou des samples n'ayant pas la même dimesionnalité
# - une chaine de caractère donnant le type d'attribut auquel on donne le plus d'importance
# - Une chaine de caractère dirigeant vers le dossier dans lequel les counts sont
# - Une chaine de caractère pour le fichier de sorti en PDF


#Sort dans un fichier toutes les possibilités de comparaison
auto_write_csv_stat<-function(UserInput,
                               UserInputColNames,
                               UserInputIfEmpty,
                               UserInputRelevel,
                               UserInputProblmaticRowsIndex,
                               UserInputPvalueThreshold,
                               UserInputMostImportantAttribute,
                               chemin,
                               output)
{
  DfPossibility<-create_possibility_table(UserInputNames = UserInputColNames)
  df<-create_df_deseq(chemin = chemin,
                    UserInput = UserInput,
                    UserInputColNames = UserInputColNames,
                    UserInputIfEmpty =UserInputIfEmpty,
                    UserInputRelevel= UserInputRelevel)
  MatrixLevel<-create_matrix_level(UserInput = UserInput,
                                   UserInputColNames=UserInputColNames,
                                   UserInputRelevel = UserInputRelevel,
                                   UserInputIfEmpty = UserInputIfEmpty,
                                   UserInputMostImportantAttribute = UserInputMostImportantAttribute)
  ListSubset<-releveling(df = df,ShowLevel = MatrixLevel,Iter = 1)
  
  Dfcomprehension<-data.frame("Formule"=c(""),#Comprendre comment le resultname marche !
                              "Comparaison"=c(""))
  for (i in ListSubset)
  {
    for (j in 1:length(DfPossibility$Formula))
    {
      SubDf<-row_to_subset(df = i,
                           row = DfPossibility[j,],
                           UserInputMostImportantAttribute = UserInputMostImportantAttribute,
                           UserInputRelevel = UserInputRelevel,
                           UserInputColNames = UserInputColNames,
                           UserInputIfEmpty = UserInputIfEmpty,
                           UserInput = UserInput)
      for ( k in SubDf[[1]])
      {
        attempt<- try(statistic_deseq(df = k,
                                      chemin = chemin,
                                      formule =DfPossibility$Formula[j],
                                      UserInputPvalueThreshold = UserInputPvalueThreshold,
                                      UserInputMostImportantAttribute = UserInputMostImportantAttribute,
                                      output=output))#On test chaque formule avec un dataset
        print(length(attempt))
        Dfcomprehension<-rbind(Dfcomprehension,c(DfPossibility$Formula[j],toString(attempt)))
        if (length(attempt)==1)#Et si ça marche pas on fait un autre dataset en croisant les doigts ça marche
        {
          attempt<- try(statistic_deseq(df = k[-UserInputProblmaticRowsIndex,],
                                        chemin = chemin,
                                        formule =DfPossibility$Formula[j],
                                        UserInputPvalueThreshold = UserInputPvalueThreshold,
                                        UserInputMostImportantAttribute = UserInputMostImportantAttribute,
                                        output=output ))
          print("ARREUUUHHH")#ou pas de temps en temps
          Dfcomprehension<-rbind(Dfcomprehension,c(DfPossibility$Formula[j],toString(attempt)))
        }
      }
    }
  }
}

auto_write_csv_stat(UserInput = list(list("A","B","K","pDUC")),
                     UserInputColNames = c("strain"),
                     UserInputIfEmpty = c("Souche Inconnue"),
                     UserInputRelevel= c("K"),
                     UserInputProblmaticRowsIndex =c(9:12) ,
                     UserInputPvalueThreshold = 0.1,
                     UserInputMostImportantAttribute ="strain" ,
                     chemin = "C:/Users/maxence.lechemia/Documents/Proteomique/1-raw",
                     output = "C:/Users/maxence.lechemia/Documents/Proteomique/ComparaisonAutomatiqueTest/")
  