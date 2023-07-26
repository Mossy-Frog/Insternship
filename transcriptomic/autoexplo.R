#Automatisation de l'analyse exploratoire de counts via DESEQ2
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
source("C:/Users/maxence.lechemia/Documents/scripts formalisés/createmainvsd.R")
source("C:/Users/maxence.lechemia/Documents/scripts formalisés/printexploplot.R")


#Initialisation de la fonction, Prend en compte dans l'ordre :
# - Une liste de liste contenant les différents attributs de l'échantillon
# - Une liste c() contenant le nom des types d'attributs différents
# - une liste c() contenant le nom de l'attribut si l'algorithme ne trouve pas l'attribut dans le nom
# - une liste c() contenant l'attribut le plus important pour la factorisation pour chaque type d'attribut
# - Une valeur qui edicte comment les lignes des graphiques sont séparées
# - Une valeur qui edicte comment les colonnes des graphiques sont séparées
# - Une chaine de caractère dirigeant vers le dossier dans lequel les counts sont
# - Une chaine de caractère pour le fichier de sorti en PDF

#Sort un PDF avec les différents graphiques possibles avec DESEQ2 pour l'exploratoire
auto_explo<-function(UserInput,
                               UserInputColNames,
                               UserInputIfEmpty,
                               UserInputRelevel,
                               UserInputMostImportantAttribute,
                               SepRowHeatmap,
                               SepColHeatmap,
                               chemin,
                               output)
{
  df<-create_df_deseq(chemin = chemin,
                      UserInput = UserInput,
                      UserInputColNames = UserInputColNames,
                      UserInputIfEmpty =UserInputIfEmpty,
                      UserInputRelevel= UserInputRelevel)
  vsd_dds<-create_main_vsd(df = df,
                       chemin = chemin,
                       UserInputMostImportantAttribute = UserInputMostImportantAttribute)
  vsd<-vsd_dds[[1]]
  dds_filtered<-vsd_dds[[2]]
  print_explo_plot(vsd = vsd,
                   dds_filtered=dds_filtered,
                   UserInputColNames=UserInputColNames,
                   UserInputMostImportantAttribute=UserInputMostImportantAttribute,
                   SepRowHeatmap = 4,
                   SepColHeatmap = 4,
                   output)
}

auto_explo(UserInput = list(list("A2","B1","C"),list("ori","R5","R10")),
           UserInputColNames = c("strain","subculture"),
           UserInputIfEmpty = c("A2","ori"),
           UserInputRelevel= c("B1","ori"),
           UserInputMostImportantAttribute ="strain",
           SepRowHeatmap = 4,
           SepColHeatmap = 4,
           chemin = "C:/Users/maxence.lechemia/Documents/Génomique/Genomics count formated and fused/",
           output = "C:/Users/maxence.lechemia/Documents/Génomique/")