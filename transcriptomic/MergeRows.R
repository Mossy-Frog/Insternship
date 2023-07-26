Go_df<- read.csv("C:/Users/maxence.lechemia/Documents/Transcriptomique/parsed-joindfed-apied-joincoled-cleancoled-gff.csv")

#Besoin de merge les lignes dupliquées lors du merge des colonnes

sub_go_df <- subset(Go_df,Go_df$go_locus_tag  != "")#création sub intéressant en prenant toutes les lignes non vides à l'endroit des locus tag
sub_go_df_trash <- subset(Go_df,Go_df$go_locus_tag  == "")#création de sub de lignes qui iron à la poubelle
Merged_rows_df <- aggregate(.~go_locus_tag, sub_go_df, paste0)#merge des lignes 
Merged_rows_df <- clean_cols(Merged_rows_df,unlist(colnames(Merged_rows_df)))#clean

for (i in colnames(Merged_rows_df))#mise en place d'un format plus lisible avec item ",' item2 dans les cases
{
  for (j in 1:length(Merged_rows_df[[i]]))
  {
    temp <- unlist(Merged_rows_df[j,i])
    
    if (length(temp)!=1)
    {
      Merged_rows_df[j,i]<-""
     for (k in 1:length(temp))
     {
       print(k)
      if (temp[k]!=""& temp[k]!=",")
      {
        if (Merged_rows_df[j,i]=="")
        {
          Merged_rows_df[j,i]<-temp[k]
        }
        else if (Merged_rows_df[j,i]!=temp[k])
        {
          Merged_rows_df[j,i]<- paste0(Merged_rows_df[j,i],",",temp[k])
        }
      }
     }
    }
  }
}


col_num_position_converter <-function(df,ColNamesToClean,superior)# permet de mettre la place des gènes en supérieur, par ex choisir entre start à 1000 ou à 1005, on prendrait 1000 car start mais si c'était la fin on prendrait 1005
{
  for (i in ColNamesToClean)
  {
    for (j in 1:length(df[[i]]))
    {
      temp<-unlist(df[j,i])
      temp<- unlist(strsplit(temp , split=',', fixed=TRUE))
      TempSupOrInf<-as.integer(temp[1])
      for (k in 1:length(temp))
      {
        if ((TempSupOrInf< as.integer(temp[k]))&superior)
        {
          TempSupOrInf<-as.integer(temp[k]) 
        }
        else if ((TempSupOrInf> as.integer(temp[k]))&(superior==FALSE))
        {
          TempSupOrInf<-as.integer(temp[k]) 
        }
      }
      df[j,i]<-TempSupOrInf
    }
  }
  return(df)
}


V4Parsed<-col_num_position_converter(Merged_rows_df,"V4",FALSE)
Df_withGO_cleaned<-col_num_position_converter(V4Parsed,"V5",T)


Df_withGO_cleaned2 = data.frame(lapply(Df_withGO_cleaned, as.character), stringsAsFactors=FALSE)

write.csv(Df_withGO_cleaned2, "C:/Users/maxence.lechemia/Documents/Transcriptomique/GFF-before-GO-trad.csv", row.names=FALSE)

