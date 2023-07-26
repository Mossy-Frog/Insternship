#création de la base de donnée reliant le gff et donc le locus tag et les  go terms associés avec des mots clés spécifiques 
#en la comparant avec une base de données qui met en corresppondance les go terms à la fonction


source("C:/Users/maxence.lechemia/Documents/Transcriptomique/Script R parsing Gff.R")


library(httr)
library(jsonlite)
library(xml2)
library(stringr)
library(reshape2)
library(plyr)
library(rlist)
library(readr)

OurDf<- read.csv("C:/Users/maxence.lechemia/Documents/Transcriptomique/GFF-before-GO-trad.csv")


OboDF<- read.csv("C:/Users/maxence.lechemia/Documents/Transcriptomique/GOoboEquals.csv")


# OboDF <- OboDF[-1]
# 
# for (i in 1:length(OboDF))
# {
#   OboDF[i]<-toString(OboDF[i])
#   OboDF[i]<-paste0(";",OboDF[i])
#   charstosupr<-paste0("/[][!#$%()*.:;<>@^`|~.{}]",'"',"'\'","'")
#   OboDF[i]<- gsub(charstosupr, "", OboDF[i])
#   OboDF[i]<- gsub(charstosupr, "", OboDF[i])
#   OboDF[i]<- gsub(charstosupr, "", OboDF[i])
#   OboDF[i]<- gsub(",", ";", OboDF[i])
#   OboDF[i]<-noquote(OboDF[i])
#   
# }
# 
# OboDF <- data.frame(OboDF)
# 
# 
# gff_col<- gff_create_col(gff =  OboDF)
gff_ordered<- gff_ordered_parsing(gff_file = OboDF)#parsing du obo transformé en csv auparavant
final_gff <- equal_remover(good_gff = gff_ordered)#fin du parsing

write.csv(final_gff, "C:/Users/maxence.lechemia/Documents/Transcriptomique/Parsed-good-obo.csv", row.names=FALSE)

# NewColDf<- c("")
# 
# OurDf$is_efflux<- NewColDf
# 
# 
# OurDf$is_influx<- NewColDf
# 
# OurDf$is_membrane<- NewColDf



# for (i in 1:length(OurDf$Ontology_term))
# {
#   print(i)
#   temp<- unlist(strsplit(OurDf$Ontology_term[i] , split=',', fixed=TRUE))
#   for (k in 1:length(temp))
#   {
#     for (l in 1:length(OboDF$id))
#     {
#       if (identical(temp[k],OboDF$id[l]))
#       {
#         OurDf$go_process[i]<- paste0(OurDf$go_process[i],",",OboDF$ide[l])
#         
#         if(grepl("efflux",OboDF$name[l]))
#         {
#           OurDf$is_efflux[i]<-TRUE
#         }
#         if(grepl("influx",OboDF$name[l]))
#         {
#           OurDf$is_influx[i]<-TRUE
#         }
#         if(grepl("membrane",OboDF$name[l]))
#         {
#           OurDf$is_membrane[i]<-TRUE
#         }
#       }
#     }
#   }
# }

#fonction qui permet de rajouter des colonnes à notre DB avec les GO terms qui vont correspondre aux mots clés cherchés, 
#avec True si présence et False si absence de ce mot clé, retourne aussi une table de correspondance si besoin
#Input:df contenant le locus tag et go term, gocol: colonne contenant les go term, dfcontainer: le fichier parsé obo en csv, listCharsearch: mot clés recherchés
#Output: une liste contenant une la dataframe avec les nouvelles colonnes et une dataframe de correspondance, avec tous les go terms qui correpondent aux mots clé recherchés
SearchAndDestroy<-function(df,GOcol,dfcontainer,ListCharSearch)
{
  TrueLenDf<-length(df)
  CorrespondanceDf = data.frame(id=c(""))
  print(CorrespondanceDf)
  for (i in ListCharSearch)
  {
    # print(1)
    df[i]<-c(FALSE)
    # print(2)
    CorrespondanceDf[i]<-c(FALSE)
    # print(CorrespondanceDf)
  }
  print(CorrespondanceDf)
  for (i in 1:length(dfcontainer[,1]))
  {
    temprow<-dfcontainer[i,c("name","def","synonym")]
    # print(sapply(temprow, grep, pattern = "membrane"))
    # print(temprow)
    templist<-c(dfcontainer[i,"id"])
    countSearch<-0
    for (j in ListCharSearch)
    {
      countrow<-0
      for (k in 1:length(temprow))
      {
        if (grepl(j,temprow[1,k]))
        {
          countSearch<-1
          countrow<-1
        }
      }
      if(countrow==1)
      {
        templist<-append(templist,TRUE)
      }
      else
      {
        templist<-append(templist,FALSE)
      }
    }
    
    # print(templist)
    if (countSearch==1)
    {
      CorrespondanceDf[nrow(CorrespondanceDf) + 1,] <- templist
    }
  }
  for (i in 1:length(df[,GOcol]))
  {
    temp<- unlist(strsplit(df[i,GOcol] , split=',', fixed=TRUE))
    for(j in temp)
    {
      print(temp)
      TempSub<- subset(CorrespondanceDf, id == j,select=c(2:length(CorrespondanceDf)))
      if (is.na(TempSub[1,1])==F)
      for (k in 1:length(ListCharSearch))
      {
        # print(TempSub)
        if (TempSub[1,k])
        {
          df[i,k+TrueLenDf]<-TRUE
        }
      }
    }
  }
  
  for (i in 1:TrueLenDf)
  {
    print(i)
    for (j in 1:length(df[,i]))
    {
      for (k in 1:length(ListCharSearch))
      {
        if (grepl(ListCharSearch[k],df[j,i])) 
        {
          df[j,TrueLenDf+k]<-TRUE
          print("ping!")
        }
      }
    }
  }
  out<-list(df,CorrespondanceDf)
  return(out)
}

Output <- SearchAndDestroy(OurDf,"Ontology_term",final_gff,c("efflux","influx","membrane","porin","channel","pump","resistance","MFS","RND","ABC","TolC","multidrug","quorum","drug","metabolism","transport","motility","transmembrane","transporter","chromatin","nucleus","DNA","RNA","repair","export","mitosis","cytoplasm","trna","cycle","maltose","ribosome","fatty acid","docking","signaling","autophagy","amino","nucleotide","chromosome","transcription","assembly","histone","toxin","filament"))
NotreDbAnnotee<-Output[[1]]


#permet de rajouter des données via d'autres data (ici mets toutes les efflux en true qui sont contenues dans un txt. )
txtenrichment<-read_table("C:/Users/maxence.lechemia/Documents/Transcriptomique/regulator.txt")



for (i in 1:length(NotreDbAnnotee$go_locus_tag))
{
  for (j in 1:length(txtenrichment$map_location))
  {
    if (grepl(NotreDbAnnotee$go_locus_tag[i],txtenrichment$map_location[j]))
        {
          NotreDbAnnotee$efflux[i]<-TRUE
          print("ping")
    }
  }
}

txtenrichment<-read_table("C:/Users/maxence.lechemia/Documents/Transcriptomique/efflux_gene_result_efflux.txt")

for (i in 1:length(NotreDbAnnotee$go_locus_tag))
{
  for (j in 1:length(txtenrichment$map_location))
  {
    if (grepl(NotreDbAnnotee$go_locus_tag[i],txtenrichment$map_location[j]))
    {
      NotreDbAnnotee$efflux[i]<-TRUE
      print("ping")
    }
  }
}


# on peut écrire notre base de donnée :)
write.csv(NotreDbAnnotee, "C:/Users/maxence.lechemia/Documents/Transcriptomique/AnnotatedDbWithNewKeyWords.csv", row.names=FALSE)
