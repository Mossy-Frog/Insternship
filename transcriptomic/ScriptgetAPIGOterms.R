library(httr)
library(jsonlite)
library(xml2)
library(stringr)
                                                                                                                                  # str_glue('https://www.ebi.ac.uk/QuickGO/services/ontology/go/search?query="{final_gff_gene$gene[1]}""\&"limit=25&page=1')

# https://www.ebi.ac.uk/QuickGO/services/ontology/go/search?query=def&limit=25&page=1

#Chercher les GO terms à partir des gene names

gotcha <- list()
for (i in 1:length(sub_final_gff_joined$gene.x))
{
  if (sub_final_gff_joined$gene.x[i] != "" & FALSE) #Retirer le &False pour pouvoir lancer le prog
  {
    requestURL <- paste0('https://www.ebi.ac.uk/QuickGO/services/ontology/go/search?query=',sub_final_gff_joined$gene.x[i],"&limit=25&page=1")
    r <- GET(requestURL, accept("application/json"))
    stop_for_status(r)
    json <- toJSON(content(r))
    json_returned <- head(fromJSON(json))
    if (json_returned[1]!=0)
    {
      gotcha<-append(gotcha,json_returned)
    } 
  }
}



############################Search go terms with Locus-tag-> Go term   ###########################


#https://rest.uniprot.org/docs/?urls.primaryName=idmapping#/uniprotkb/searchCursor
#https://www.uniprot.org/id-mapping

for (i in 1:length(final_gff_joined$Ontology_term))#Parcourt la colonne des locus
{
  gotcha_accession <- list()#créé une liste ou mettre l'output de la requête
  if (is.na(final_gff_joined$old_locus_tag[i])==FALSE & FALSE)###Retirer le & FALSE pour lancer le code #Si la case de la colonne n'est pas vide
  {
    if (final_gff_joined$old_locus_tag[i] != "")#si la case de la colonne ne contient pas rien
  {
      print(i)#checker l'avancement de la requête
      requestURL <- paste0('https://rest.uniprot.org/uniprotkb/search?query=',final_gff_joined$old_locus_tag[i])#créer la requête avec l'URL
      r <- GET(requestURL, accept("application/json"))#Mettre la requête dans un get pour la lancer
      stop_for_status(r)
      json <- toJSON(content(r))#transformer le json en list potable
      json_returned <- head(fromJSON(json))#pareil ?
      json_returned <- unlist(json_returned)#l'unlist pour la parcourir et mettre des si
      if (json_returned[1]!=0)#si le premier élément de json n'est pas 0 (si la requête n'est pas nulle)
      {
        empty_go<-c()#créer une liste pour mettre les go terms de chaque ligne
        for (j in 1:length(json_returned))#parcourir la liste json
        {
          if (grepl("GO:",json_returned[j]))#on regarde si il y a des go terms
          {
            empty_go<- append(empty_go,json_returned[j])#s'il y a des go terms, on les met dans la liste
            
          }
        }
        for (k in 1:length(empty_go))#on parcourt la liste des go terms retrouvés
        {
          final_gff_joined$Ontology_term[i]<-paste0(final_gff_joined$Ontology_term[i],",",empty_go[k])#on ajoute les go terms à la liste de la database
        }
        if (is.na(final_gff_joined$old_locus_tag[i]))#si la case était vide, on remplace la case vide par la liste de go term directement
        {
          final_gff_joined$Ontology_term[i]<- empty_go
        }
      }
    }
  }
}

# for (i in 1:length(sub_final_gff_joined$gene.y))
# {
#   if (sub_final_gff_joined$gene.y[i] != "")
#   {
#     requestURL <- paste0('https://www.ebi.ac.uk/QuickGO/services/ontology/go/search?query=',sub_final_gff_joined$gene.y[i],"&limit=25&page=1")
#     r <- GET(requestURL, accept("application/json"))
#     stop_for_status(r)
#     json <- toJSON(content(r))
#     json_returned <- head(fromJSON(json))
#     if (json_returned[1]!=0)
#     {
#       gotcha<-append(gotcha,json_returned)
#     } 
#   }
# }
# 
# ################                    NE PAS LANCER DE SUITE !!!! 
# 
# requestURL <- paste0('https://www.ebi.ac.uk/QuickGO/services/ontology/go/search?query=',final_gff_gene$gene[1],"&limit=25&page=1")
# print(requestURL)
# r <- GET(requestURL, accept("application/json"))
# 
# stop_for_status(r)
# 
# json <- toJSON(content(r))
# json_returned <- head(fromJSON(json))
# print(json_returned)
# 
# ################                    NE PAS LANCER DE SUITE !!!! 


# BiocManager::install("UniProt.ws")
# install.packages("UniprotR")
# BiocManager::install("GenomicAlignments")
# library(UniprotR)
# 
# Accessions <-GetAccessionList("https://s3.amazonaws.com/csvpastebin/uploads/9571fa356c67a0c7c95e8431799a051a/Accessions.csv") 
# GeneOntologyObj <- GetProteinGOInfo(Accessions) 
# print(Accessions)
# acc <- list("Q63Z36")
# acc <- unlist(acc)
# GeneOntologyObj <- GetProteinGOInfo(acc) 
# 
# closeAllConnections()
# 






# acc <- list("RS20545")
# acc <- unlist(acc)
# GeneOntologyObj <- GetProteinGOInfo(acc) 




# #Tentative ratée de trouvé le primary accession de uniprot
# suppressPackageStartupMessages({
#   library(UniProt.ws)
# })
# 
# 
# # BiocManager::install("UniProt.ws")
# # up <- UniProt.ws(272560)
# # 
# # up <- UniProt.ws(taxId=9606)
# # 
# # closeAllConnections()









###Count go-terms####


#Sert à compter le nombre de go terms dans une colonne, là c'est sur le produit final
Go_count_final<-list()
for (i in 1:length(final_gff_joined$Ontology_term))
{
  temp<- unlist(strsplit(final_gff_joined$Ontology_term[i] , split=',', fixed=TRUE))
  temp<-temp[!duplicated(temp)]
  for (j in 1:length(temp))
  {
    if (length(temp)!=0)
    {
      if (  is.na(temp[j])==FALSE)
      {
        if (grepl("GO:",temp[j]))
        {
          Go_count_final<-append(Go_count_final,temp[j]) 
        }
      } 
    }
  }
}



#Sert à compter le nombre de go terms dans une colonne, là sur uniquement le refseq, pour comparer entre le résultat final et avant
Go_count_refseq<-list()
for (i in 1:length(Refseq_CDS_gff$Ontology_term))
{
  temp<- unlist(strsplit(Refseq_CDS_gff$Ontology_term[i] , split=',', fixed=TRUE))
  temp<-temp[!duplicated(temp)]
  for (j in 1:length(temp))
  {
    if (length(temp)!=0)
    {
      if (  is.na(temp[j])==FALSE)
      {
        if (grepl("GO:",temp[j]))
        {
          Go_count_refseq<-append(Go_count_refseq,temp[j]) 
        }
      } 
    }
  }
}




# write.csv(final_gff_joined, "C:/Users/maxence.lechemia/Documents/Transcriptomique/parsed-gff-with-GOterms.csv", row.names=FALSE)


test_csv<- read.csv("C:/Users/maxence.lechemia/Documents/Transcriptomique/parsed-gff-with-GOterms.csv")



ppt_df<- test_csv[,c("gene.y","Ontology_term")]
ppt_df2<- Refseq_CDS_gff[,c("gene","Ontology_term")]
