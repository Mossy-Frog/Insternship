#Automatisation de la création de heatmap avec les analyses statistiques de DESEQ2
#Auteur: Maxence LECHEMIA



#Initialisation des libraires nécessaires
library(DESeq2)
library(magrittr)
library(apeglm)
library(genefilter)
library(pheatmap)
library(data.table)
library(readr)
library(ggplot2)
library(gridExtra)
library(readr)
library(devtools)
library(ComplexHeatmap)


#Initialisation de la fonction nécessaire, il faudra peut être retrouver le chemins correct.
source("C:/Users/maxence.lechemia/Documents/scripts formalisés/rowlocustogene.R")




#Initialisation de la fonction, Prend en compte dans l'ordre :
# - Une chaine de caractère dirigeant vers le dossier dans lequel les fichiers de comparaisons stats sont
# - Une chaine de caratère donnant le chemin vers le dataset des annotations 
# - Une liste c() donnant les mots clefs contenus dans les annotations choisies pour la heatmap
# - Une variable pour le nombre de  gènes présents dans la heatmap
# - Booléen pour la clusterisation des lignes
# - Booléen pour la clusterisation des colonnes
# - Une chaine de caractère pour le fichier de sorti en JPEG

#Sort un JPEG avec la heatmap et les annotations selon les paramètres donnés

auto_heatmap_logfoldchange<-function(chemin,
                                     PathAnnotated,
                                     ListCol,
                                     NbGenes,
                                     ClusterRow=F,
                                     ClusterCol=F,
                                     output)
{
  list_csv_files <- list.files(path =chemin,full.names = T)
  DfOfDf<- list()
  count<-1
  for (i in list_csv_files)
  {
    
    DfOfDf[[count]]<-as.data.frame(read.csv(i))
    count<-count+1
  }
  
  for (i in 1:length(list_csv_files))
  {
    DfOfDf[[i]]<-DfOfDf[[i]][,c("X","log2FoldChange")]
    print(gsub(chemin,"",toString(list_csv_files[i])))
    colnames(DfOfDf[[i]])[2]<- gsub(".csv","",gsub(chemin,"",toString(list_csv_files[i])))
  }
  dflogfold<-data.frame(c(""))
  colnames(dflogfold)[1]<-"X"
  for (i in 1:length(DfOfDf))
  {
    dflogfold <- merge(DfOfDf[[i]],dflogfold,c("X"),all = TRUE)
  }
  
  dflogfold<-dflogfold[-c(1),]
  row.names(dflogfold)<-dflogfold$X

  print(dflogfold)
  dflogfold <- dflogfold[,!names(dflogfold) %in% c("X")]
  
  

  dflogfold["Max"]<-c(0)
  for (i in 1:length(dflogfold[,1]))
  {
    print(i)
    TempMax<-max(dflogfold[i,],na.rm = T)
    print(TempMax)
    if(is.numeric(TempMax))
    {
      dflogfold[i,"Max"]<-TempMax
    }
  }

  dflogfold <-dflogfold[order(-dflogfold$Max),]
  print(dflogfold)
  
  
  anno_row <- as.data.frame(NotreDbAnnotee[which(NotreDbAnnotee$go_locus_tag%in% row.names(dflogfold)),ListCol])
  for (i in 1:length(anno_row))
  {
    anno_row[,i]<-as.numeric(anno_row[,i])
    print(i)
  }
  rownames(anno_row) = row.names(dflogfold)[which(row.names(dflogfold)%in% NotreDbAnnotee$go_locus_tag)]
  dflogfold<-row_locus_to_gene(dflogfold,NotreDbAnnotee = NotreDbAnnotee)
  anno_row<-row_locus_to_gene(anno_row,NotreDbAnnotee = NotreDbAnnotee)
  print(length(anno_row[,1]))
  print(length(dflogfold[,1]))
  dflogfold<-dflogfold[,-c(length(dflogfold))]
  dflogfold<-as.matrix(dflogfold)
  jpeg(file=output,width = 1000,height = 1000)
  assign("dflogfold",dflogfold,envir = globalenv())
  ComplexAnnoRow<-rowAnnotation(foo1=unlist(anno_row[1]),annotation_label = colnames(anno_row)[1],show_legend=F)
  if (length(colnames(anno_row))>1)
  for ( i in 2:length(colnames(anno_row)))
  {
    ComplexAnnoRow<-c(ComplexAnnoRow,rowAnnotation(foo1=unlist(anno_row[i]),annotation_label = colnames(anno_row)[i],show_legend=F))
  }
  View(ComplexAnnoRow)
  ht<-Heatmap(dflogfold[c(1:NbGenes),],col = c("blue","pink","red"),name = "logfoldchange",na_col = "white",cluster_rows = F,cluster_columns = F, heatmap_width = unit(30, "cm"), heatmap_height = unit(30, "cm"),right_annotation = ComplexAnnoRow[1:NbGenes])
  draw(ht)
  dev.off()
  
  
  #pheatmap(CompareDf[1:NbGenes,-c("Max")],annotation_row = anno_row,cluster_rows = ClusterRow,cluster_cols = ClusterCol)
  
  
}





auto_heatmap_logfoldchange("C:/Users/maxence.lechemia/Documents/Proteomique/ComparaisonAutomatiqueTest/",
                           "C:/Users/maxence.lechemia/Documents/Transcriptomique/AnnotedDBaseWithGoRealFinal.csv",
                           ListCol = c("efflux","membrane"),
                           NbGenes=50,
                           output = "C:/Users/maxence.lechemia/Documents/Transcriptomique/heatmap.jpeg")
