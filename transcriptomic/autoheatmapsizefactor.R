#Creation de la heatmap avec les différents samples
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
library(magrittr)
library(apeglm)
library(genefilter)
library(data.table)
library(readr)
library(gridExtra)

#Initialisation des fonctions nécessaires, il faudra peut être retrouver les chemins corrects.

source("C:/Users/maxence.lechemia/Documents/scripts formalisés/createdfdeseq.R")
source("C:/Users/maxence.lechemia/Documents/scripts formalisés/createmainvsd.R")
source("C:/Users/maxence.lechemia/Documents/scripts formalisés/printheatmapsizefactor.R")


#Initialisation de la fonction, Prend en compte dans l'ordre :
# - Une liste de liste contenant les différents attributs de l'échantillon
# - Une liste c() contenant le nom des types d'attributs différents
# - une liste c() contenant le nom de l'attribut si l'algorithme ne trouve pas l'attribut dans le nom
# - une liste c() contenant l'attribut le plus important pour la factorisation pour chaque type d'attribut
# - une chaine de caractère donnant le type d'attribut auquel on donne le plus d'importance
# - Une valeur qui edicte comment les lignes des graphiques sont séparées
# - Une valeur qui edicte comment les colonnes des graphiques sont séparées
# - Une chaine de caractère dirigeant vers le dossier dans lequel les counts sont
# - Une chaine de caractère pour le fichier de sorti en PDF
# - Une chaine de caractère donnant le chemin du fichier contenant les annotations
# - Une liste c() contenant les mots-clefs d'annotation
# - 

auto_heatmap_sizefactor<-function(UserInput,
                                  UserInputColNames,
                                  UserInputIfEmpty,
                                  UserInputRelevel,
                                  UserInputMostImportantAttribute,
                                  SepRowHeatmap,
                                  SepColHeatmap,
                                  chemin,
                                  output,
                                  PathAnnotated,
                                  listcol,
                                  SelectGene=1,
                                  TableSelection=c(),
                                  NumberSelect=100,
                                  ColTarget=c(),
                                  Flatten=TRUE,
                                  ColDf,
                                  RowDf=c(1),
                                  BoolClusterCol=T,
                                  BoolClusterRow=T)
{
  df<-create_df_deseq(chemin = chemin,
                      UserInput = UserInput,
                      UserInputColNames = UserInputColNames,
                      UserInputIfEmpty =UserInputIfEmpty,
                      UserInputRelevel= UserInputRelevel)#creation of the df, summing up the infos of the user for DESeq2
  
  vsd_dds<-create_main_vsd(df = df,
                           chemin = chemin,
                           UserInputMostImportantAttribute = UserInputMostImportantAttribute)
  vsd<-vsd_dds[[1]]#Creation of the VSD, table used for the heatmap
  dds_filtered<-vsd_dds[[2]] 
  
  a<-print_heatmap_sizefactor(df = df,
                           vsd = vsd,
                           PathAnnotated = PathAnnotated,
                           listcol= listcol,
                           SelectGene = SelectGene,
                           TableSelection = TableSelection,
                           NumberSelect = NumberSelect,
                           ColTarget = ColTarget,
                           Flatten = Flatten,
                           UserInputMostImportantAttribute = UserInputMostImportantAttribute,
                           ColDf = ColDf,
                           RowDf = RowDf,
                           BoolClusterCol = BoolClusterCol,
                           BoolClusterRow = BoolClusterRow)#Print the heatmap
  return(a)
}



a<-auto_heatmap_sizefactor(UserInput = list(list("A2","B1","C"),list("ori","R5","R10")),
                        UserInputColNames = c("strain","subculture"),
                        UserInputIfEmpty = c("A2","noSC"),
                        UserInputRelevel= c("B1","ori"),
                        UserInputMostImportantAttribute ="strain",
                        SepRowHeatmap = 4,
                        SepColHeatmap = 4,
                        chemin = "C:/Users/maxence.lechemia/Documents/Génomique/Genomics count formated and fused/",
                        output = "C:/Users/maxence.lechemia/Documents/Génomique/",
                        PathAnnotated="C:/Users/maxence.lechemia/Documents/Transcriptomique/AnnotedDBaseWithGoRealFinal.csv",
                        listcol=c("efflux","RND","multidrug","resistance","pump"),
                        SelectGene=1,
                        TableSelection=c(),
                        NumberSelect=15,
                        ColTarget=c(),
                        Flatten=TRUE,
                        ColDf=c(),
                        RowDf=c(),
                        BoolClusterCol=F,
                        BoolClusterRow=F)


jpeg(filename = "C:/Users/maxence.lechemia/Documents/Multi-omics/plotautomated.jpeg")
print(a)
dev.off()

