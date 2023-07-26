library(reshape2)
library(stringr)
library(plyr)
library(rlist)


#Chemin du gff
chemin = "C:/Users/maxence.lechemia/Documents/Transcriptomique/BP_corr_nomChr.gff"

#Premier ""parsing"" et convertit le gff en dataframe 

# read.delim(chemin, header=F, comment.char="#") -> gff
# library(dplyr)
# 
# gff %>% 
#   group_by(V4, V5) %>% 
#   summarize(across(V8:V9, ~ifelse(all(is.na(.x)), NA, paste0(na.omit(.x), collapse = ","))), .groups = "drop")
# 
# 
# gff <- gff %>% 
#   group_by(V5) %>% 
#   summarise(V9 = paste(V9, collapse = ","))
# 


# 
# for (i in gff[-1]){
#   i <- str_split(i,";")
#   i <- str_split(i, "=")
#   print(i)
# }

##Remove the uninteresting rows



# print(gff)
# print(gff$V9[2])

##Cr?er une petite dataset contenant juste les genes
# gff_gene <- subset(gff,gff$V3  == 'gene')
#print(gff_gene$V9[2])

##cr?er les colonnes importantes au parsing, en cr?ant le maximum des colonnes

# gff_gene[c("ID", "Dbxref", "Name","gbkey","gene", "gene_biotype","gene_synonym","locus_tag")] <- str_split_fixed(gff_gene$V9, ';', 8)

##gff_gene$locus_tag[grepl( "locus", gff_gene$gene_synonym[i], fixed = TRUE)] <- gff_gene$gene_synonym






#print(length(gff_gene$gene_synonym))

##Regarde dans une colonne (ici gene synonym) pour savoir si les locustag qui auraient pu ?tre d?plac? lors de la cr?ation des colonnes
## et les replace l? ou il faut.
# for (i in 1:length(gff_gene$gene_synonym))
# {
#   #print(grepl( "locus", gff_gene$gene_synonym[i], fixed = TRUE))
#   if (grepl( "locus", gff_gene$gene_synonym[i], fixed = TRUE)==TRUE)
#   {
#     gff_gene$locus_tag[i]<-gff_gene$gene_synonym[i]
#     gff_gene$gene_synonym[i] <- ""
#     #print("AAAAAAAAAAAAAA")
# 
#   }
# }


##Regarde dans une colonne (ici gene biotype) pour savoir si les locustag qui auraient pu ?tre d?plac? lors de la cr?ation des colonnes
## et les replace l? ou il faut.


# for (i in 1:length(gff_gene$gene_biotype))
# {
#   #print(grepl( "locus", gff_gene$gene_synonym[i], fixed = TRUE))
#   if (grepl( "locus", gff_gene$gene_biotype[i], fixed = TRUE)==TRUE)
#   {
#     gff_gene$locus_tag[i]<-gff_gene$gene_biotype[i]
#     gff_gene$gene_biotype[i] <- ""
#     #print("AAAAAAAAAAAAAA")
#     
#   }
# }

##Regarde dans une colonne (ici gene) pour savoir si les biotype qui auraient pu ?tre d?plac? lors de la cr?ation des colonnes
## et les replace l? ou il faut.

# for (i in 1:length(gff_gene$gene))
# {
#   #print(grepl( "locus", gff_gene$gene_synonym[i], fixed = TRUE))
#   if (grepl( "bio", gff_gene$gene[i], fixed = TRUE)==TRUE)
#   {
#     gff_gene$gene_biotype[i]<-gff_gene$gene[i]
#     gff_gene$gene[i] <- ""
#     #print("AAAAAAAAAAAAAA")
#     
#   }
# }

# gff <- gff[gff$V3 != "sequence_feature", ]
# gff <- gff[gff$V3 != "region", ]



############################################EQUAL REMOVER#########################################
equal_remover <-function(good_gff)
{
  ##Supprime tout ce qui est derri?re un "=" dans le dataset
  
  for (i in 1:length(good_gff))#On parcourt le dataset en colonne
  {
    for (j in 1:length(good_gff[,i]))#en ligne maintenant
    {
      if (is.character(good_gff[j,i]))#si la valeur de la case est un str (obligatoire pour le strsplit)
      {
        tempvar = good_gff[j,i]#tempvar prend la valeur de la case
        # print(tempvar)
        s1 = unlist(strsplit(tempvar , split='=', fixed=TRUE))[2] #on split la valeur de la case en 2, et on garde la partie apr?s le "="
        if (is.na(s1)==FALSE)#S'il n'y a pas d'?gal, il n'y aura pas de partie apr?s l'?gal, donc ?a ferait vide, donc on emp?che ?a l?
        {
          good_gff[j,i]<-s1  #La deuxi?me partie (sans le titre ni le =, juste la valeur "v?ritable" de la case) remplace la valeur de la case
        }
      }
    }
  }
return(good_gff)  
}


# gff_gene_test <- gff_gene[,!names(gff_gene) %in% c("V9","locus_tag")]
# 
# for (j in gff_gene_test)
# {
#   if (grepl( "locus", j, fixed = TRUE)==TRUE)
#   {
#     print(j)
#   }
#     }
# 



# print(gff_gene$V8)
# print(gff_gene$V6)

##Retire les colonnes vides de donn?es et la V9 qui est maintenant bien pars?e.

# gff_gene <- gff_gene[,!names(gff_gene) %in% c("V9","V6","V8")]

#######################################################Create columns GFF############################################

gff_create_col <- function(gff)
{
  # txt= "MA0051;IRF2=xml aaaa;abab=ccccc;acac=aaaaa" 
  list_lign <- list()
  # b <-regmatches(txt, gregexpr('(?<=;).*?(?==)', txt, perl=T))[[1]]   ###It's alive it's alive !
  # for (i in 1:length(b))
  # {
  #   print(b[i])
  #   list_lign[i]<-b[i]
  # }
  
  #test premi?re partie
  
  gff_col_name<- list()#on cr?? une liste qui contiendra le nom des colonnes
  for (i in 1:length(gff[,length(gff)]))#on parcourt la derni?re colonne de la dataset, qui contient tout ce qui doit ?tre encore pars?
  {
    if (i%%100 == 0)
    {
      print(gff_col_name,quote=F)
    }
    # print(i)
    tempvar <-gff[i,length(gff)[1]]#Pour chaque case de la colonne tempvar prend la valeur de la case
    tempvar<-paste0(";", tempvar)#on ajoute au tout d?but un ; pour pouvoir faire un parsing plus simple
    list_lign <- list() #on cr?? une liste de ligne qui sera utile plus tard
    tempvar <-regmatches(tempvar, gregexpr('(?<=;).*?(?==)', tempvar, perl=T))[[1]]# ?a extrait toutes les donn?es comprises entre un ; et un =   ###It's alive it's alive !
    for (j in 1:length(tempvar))#on parcourt toutes les donn?es extraites de la case
    {
      list_lign[j]<-tempvar[j]#on mets ces donn?es dans une liste
      # print(list_lign)
    }
    for (k in 1:length(list_lign))#on parcourt cette liste
    {
      if (length(gff_col_name)==0)#Si la liste de colonne finale n'a pas encore d'?l?ment, on rajoute automatiquement un nom de colonne
      {
        gff_col_name<-append(gff_col_name,list_lign[k])
      }
      else
      {a <-0
        # if (length(grep(list_lign[k], gff_col_name)) == 0)
        for (b in 1:length(gff_col_name))#on parcourt la liste finale de colonne
        {
          if (toString(list_lign[k])== gff_col_name[b])#Si un ?l?ment de la liste contenue dans la case n'est pas dans la liste finale, on la rajoute
          {
            a<-1
          }
        
        }
      if (a ==0)
      {
        gff_col_name<-append(gff_col_name,list_lign[k])
      }
     }
    }
  }
  
  # col_tofix<-gff_CDS[,ncol(gff_CDS)]
  # print(col_tofix)
  
  final_col <-unlist(c(gff_col_name))#notre liste finale de colonne n'est pas dans un bon format, on la remet bien
  
  gff[final_col] <- str_split_fixed(gff[,ncol(gff)], ';', length(gff_col_name))#on peut split la serni?re en colonne gr?ce ? la liste de colonne pr?c?demment faite
  
  return(gff)
}









#print(gff_CDS$V9)

# gff_CDS[c("ID","Parent", "Dbxref", "Name","Note","gbkey", "gene","product","protein_id","transl_table")] <- str_split_fixed(gff_CDS$V9, ';', 10)
# gff_CDS <- gff_CDS[,!names(gff_CDS) %in% c("V9","V6","V8")]

#print(length(gff_CDS))




#####################FONCTION GFF PARSING OREDRING###########################
gff_ordered_parsing <- function(gff_file){
  
  for (i in 1:length(gff_file[,1]))#toutes les lignes de 1 ? 5k
  {
    # print(length(gff_file[,1]))
  
      temp_list <- gff_file[i,] #?a marche !
      if (i%%100 == 0)
      {
        print(i)
      }
      #print(j)
      #print(temp_list)
    for (k in 1:length(temp_list))#parcourir la liste cr??e
    {
      for (l in 1:length(colnames(gff_file))) #parcourt toutes les colonnes pour chaque variable
      {
        tempvar = toString(temp_list[k]) #Variable de la liste parcourue
        if (is.character(tempvar))#check si c'est un charact?re pour le strsplit, normalement ?a l'est avec le toString, mais ?a marche pas tout le temps 
        {
          
          #print(tempvar)
          s1 = unlist(strsplit(tempvar , split='=', fixed=TRUE))[1]#strsplit la variable, pour ne garder que la partie avant le "="
          if (is.na(s1)==FALSE) #pas vraiment utile
          {
            if (colnames(gff_file)[l]== s1) #si la premi?re partie de la variable (donc le nom) correspond au nom de la colonne
            {
              gff_file[i,l]<-temp_list[k] #change dans le dataframe la valeur de la colonne ? la bonne valeur qui correspond
            }
          }
        }
      }
    }
  }
  #retirer des donn?es d?doubl?es, dans les cases sens?es ?tre vide, il faut les vider, donc on check dans la matrice si le nom avant est ?gale au nom de la case. 
  for (i in 1:length(gff_file))#parcourir les colonnes du gff
  {
    for (j in 1:length(gff_file[,i])) #parcourir les lignes du gff
    {
      
      tempvar = toString(gff_file[j,i]) # mettre la case parcourue dans une variable
      # print(colnames(gff_file)[i])
      # print(tempvar)
      if (is.character(tempvar)) #savoir si c'est bien un string (affiche une erreur si on met un string dans le strsplit)
      {
        
        #print(tempvar)
        s1 = unlist(strsplit(tempvar , split='=', fixed=TRUE))[1]#split pour savoir le nom de la case
        if (is.na(s1)==FALSE & s1!=tempvar)#puisque s1 est la premi?re partie de la case, s'il n'y a pas de "=" s1 retourne toute la case, donc ?a virerait des cases non nomm?es qu'on a envie de garder 
            {
  
          if (grepl( colnames(gff_file)[i], s1, fixed = TRUE)==FALSE) #si le nom de la colonne pas pareil que le nom de la case, on vire
          {
            gff_file[j,i]<-"" 
          }
        }
      }
    }
  }
  return(gff_file)#on retourne le gff !!
}
# 
# 
# 
# 
# gff <- gff[gff$V3 != "sequence_feature", ]
# gff <- gff[gff$V3 != "region", ]
# 
# 
# gff_gene <- subset(gff,gff$V3  == 'gene')
# gff_CDS <- subset(gff,gff$V3  == 'CDS')
# 
# 
# gff_col_gene<- gff_create_col(gff =  gff_gene)
# gff_ordered_gene<- gff_ordered_parsing(gff_file = gff_col_gene)
# final_gff_gene <- equal_remover(good_gff = gff_ordered_gene)
# 
# gff_col_CDS<- gff_create_col(gff =  gff_CDS)
# gff_ordered_CDS<- gff_ordered_parsing(gff_file = gff_col_CDS)
# final_gff_CDS <- equal_remover(good_gff = gff_ordered_CDS)
# 

# write.csv(final_gff_CDS, "C:/Users/maxence.lechemia/Documents/Transcriptomique/GFF parsed CSV/parsed_gff_CDS_K96243_2015_refseq.csv", row.names=FALSE)
# 
# write.csv(final_gff_gene, "C:/Users/maxence.lechemia/Documents/Transcriptomique/GFF parsed CSV/parsed_gff_gene_K96243_2015_refseq.csv", row.names=FALSE)
# 
# ##Cr?er une dataset avec tous les r?sidus
# 
# gff_filter <- subset(gff,gff$V3 !="CDS" & gff$V3 !="gene" )
# 
# ##Voir si avec les r?sidus, le CDS et le gene on a tout
# a = length(gff_CDS$V1)+length(gff_filter$V1)+length(gff_gene$V1)
# print(length(gff_CDS))
# print(a)

# num_col_note <- lengths(regmatches(gff$V9[2], gregexpr(";", gff$V9[2])))
# print(num_col_note)
# for (i in num_col_note){
#   gff[c(paste0("Col",i), "V9")] <- str_split_fixed(gff$V9, ';', 2)
# }


####################################################Give parsed gff##########################################
give_parsed_gffs <- function(path, option) #On appelle les fonctions pr?c?demment faites pour cr?er rapidement les dataframes
{
  read.delim(path, header=F, comment.char="#") -> gff

  
  
  sub_gff <- subset(gff,gff$V3  == option)
  
  
  gff_col<- gff_create_col(gff =  sub_gff)
  gff_ordered<- gff_ordered_parsing(gff_file = gff_col)
  final_gff <- equal_remover(good_gff = gff_ordered)
  
  return(final_gff)
}




