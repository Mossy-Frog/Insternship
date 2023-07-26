Go_df<- read.csv("C:/Users/maxence.lechemia/Documents/Transcriptomique/parsed-gff-with-GOterms-joined-cols.csv")


clean_cols <- function(df,ColsToClean)# sert à enlever les duplicats, les "," en trop les "na" etc
{
  for (i in ColsToClean)#Parcourt colonnes à nettoyer
  {
    for (j in 1:length(df[[i]]))#Parcourt les cases parmi les colonnes
    {
      if (typeof(df[j,i])=="character")
      {
        temp<-df[j,i]
        temp<- unlist(strsplit(temp , split=',', fixed=TRUE))#transfrome case en liste
        temp<-temp[!duplicated(temp)]
        
        for (k in 1:length(temp))#parcourt la case en liste 
        {
          if (length(temp)!=0)#empêche bug
          {
            if ((is.na(temp[k])==FALSE)&(identical(temp[k],"")==FALSE)&(temp[k]!="NA")&(temp[k]!=","))
            {
              # print(temp[k])
              if (k==1)#Prend en compte le premier élément de la case 
              {
                df[j,i]<-temp[k]
              }
              else
              {
                df[j,i]<-paste0(df[j,i],",",temp[k]) #Si l'élément est pas un "" ou na ou une virgule, rajoute l'élément
              }
            }
            else
            {
              df[j,i]<-""
            }
          }
          else
          {
            df[j,i]<-""
          }
        } 
      }
    }
  }
  return(df)
}




NewDF <- clean_cols(Go_df,unlist(colnames(Go_df)))


write.csv(NewDF, "C:/Users/maxence.lechemia/Documents/Transcriptomique/parsed-joindfed-apied-joincoled-cleancoled-gff.csv", row.names=FALSE)
