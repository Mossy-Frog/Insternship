source("C:/Users/maxence.lechemia/Documents/Transcriptomique/Script R parsing Gff.R")

join_df <- function(df1,df2,columns,remove_col)
{
  if (remove_col)#Remove specific columns que si on veut
  { 
    df1 <- df1[,!names(df1) %in% c("V2", "CDS", "V6","gbkey","transl_table","start_range","end_range","V9","V8","Parent","partial","V3","V7","ID")]
    df2 <- df2[,!names(df2) %in% c("V2", "gene", "V6","gbkey","transl_table","start_range","end_range","V9","gene_biotype","V8","V3","V7","ID")]
  }
  df3<- merge(df1,df2,by= columns,all = TRUE) #merge les 2 dataframes en fonction des colonnes donn?es, le all = true pour faire un merge inclusif et non exclusif
  return(df3)
}






# final_gff_CDS <- final_gff_CDS[,!names(final_gff_CDS) %in% c("V2", "CDS", "V6","gbkey","transl_table","start_range","end_range","V9","V8","Parent","partial","V3","V7","ID","Dbxref","Name")]
# final_gff_gene <- final_gff_gene[,!names(final_gff_gene) %in% c("V2", "gene", "V6","gbkey","transl_table","start_range","end_range","V9","gene_biotype","V8","V3","V7","ID","Dbxref","Name")]
# 
# 
# 
# 
# joined_df <- join_df(final_gff_CDS,final_gff_gene,c("V4", "V5","V1"))


chemin_EMBL <-"C:/Users/maxence.lechemia/Documents/Transcriptomique/BP_corr_nomChr.gff" #Donne le chemin du gff
# chemin_Refseq <-"C:/Users/maxence.lechemia/Documents/Transcriptomique/References/K96243_2004/GCF_000011545.1_ASM1154v1_genomic.gff" #donne le chemin du deuxi?me gff



EMBL_CDS_gff <-give_parsed_gffs(chemin_EMBL,"CDS")#On parse les 2 gff le plus possible
EMBL_gene_gff <- give_parsed_gffs(chemin_EMBL,"gene")


#merge des 2 subsets (CDS et gene) en fonction des colonnes donn?es, donc l? c'est le d?but, la fin de chaque genes et en fonction des chromosomes
joined_df1 <- join_df(EMBL_CDS_gff,EMBL_gene_gff,c("V4", "V5","V1"),remove_col =  TRUE)




# Refseq_CDS_gff <- give_parsed_gffs(chemin_Refseq,"CDS") #on fait la m?me chose pour l'autre fichier gff
# Refseq_gene_gff <- give_parsed_gffs(chemin_Refseq,"gene")
# 
# 
# joined_df2 <- join_df(Refseq_CDS_gff,Refseq_gene_gff,c("V4", "V5","V1"),remove_col =  TRUE)#on fait la m?me chose pour l'autre fichier gff



write.csv(joined_df1,"./ProteinIDtoLocusTag.csv")


#################################################REMPLACER NOM DE CHROMOSOME PAR CHROMOSOME N#######################################
#puisque les fichiers ont des noms de chromosomes diff?rents, on doit rendre le nom des chromosomes identiques


#On doit reprendre les databases originelles dans les fichiers gff pour d?terminer les r?gions

df1 <-   read.delim(chemin_EMBL, header=F, comment.char="#")
df2 <-   read.delim(chemin_Refseq, header=F, comment.char="#")


#on prend un tout petit subset, qui contient juste la description de chaque r?gion, donc des chromosomes
reg_df1 <- subset(df1,df1$V3  == "region")
reg_df2 <- subset(df2,df2$V3  == "region")

#V5 contient la longueur des r?gions, on calcule le nombre de r?gions pour chaque fichier
len_reg1 <- length(reg_df1$V5)
len_reg2 <- length(reg_df2$V5)
# print(len_reg1)


#On regarde qu'il y a bien le m?me nombre de chromosomes
if (len_reg2 == len_reg1)
 {
  col_df1 <- reg_df1$V5#On prend la colonne contenant les longueurs de r?gions
  # print("ahah")
  col_df1 <- lapply(col_df1,sort,decreasing=TRUE)#on trie de mani?re d?croissante la liste(le plus grand en premier)
  for (i in 1:length(col_df1))#on parcourt la liste des chromosomes
  {
    tempindex <- match(col_df1[i],reg_df1$V5)#on prend l'index(donc la position) de chaque r?gion
    col_df1[i]<-reg_df1$V1[tempindex]#et avec l'index ?a nous permet de remplacer la longueur de la r?gion par le nom de la r?gion
  }
  for (i in 1:length(joined_df1$V1)) #on parcourt la dataset joined, dont on a envie de remplacer le nom
    
  {
    for (j in 1:length(col_df1))#on parcourt notre liste de noms de r?gions
    {
      # print(joined_df1$V1[i] ==col_df1[j])
      
      if (joined_df1$V1[i] ==col_df1[j])#quand le nom de la r?gion dans la liste est ?gal au nom de la r?gion sur la ligne de notre dataset
      {
        joined_df1$V1[i] <- paste0("Chromosome",j)#on remplace le nom de la r?gion par "chromosome" et l'index dans la liste, ce qui fait que le chromosome 1 sera le plus grand, le 2 le  deuxi?me plus grand et ainsi de suite
      }
    }
  }
  ###Ensuite c'est la m?me chose mais pour l'autre fichier
  col_df2 <- reg_df2$V5
  # print("ahah")
  col_df2 <- lapply(col_df2,sort,decreasing=TRUE)
  for (i in 1:length(col_df2))
  {
    tempindex <- match(col_df2[i],reg_df2$V5)
    col_df2[i]<-reg_df2$V1[tempindex]
  }
  for (i in 1:length(joined_df2$V1))
    
  {
    for (j in 1:length(col_df2))
    {
      if (joined_df2$V1[i] ==col_df2[j])
      {
        joined_df2$V1[i] <- paste0("Chromosome",j)
      }
    }
  }
}


#On peut join les 2 datasets !


####################################C'EST A TRIER !!!!!###########################################
final_gff_joined <- join_df(joined_df1,joined_df2,c("V4", "V5","V1"), remove_col =  FALSE) ##join les 2 datasets (EMBL et Refseq)


final_gff_joined <- final_gff_joined[,!names(final_gff_joined) %in% c("pseudo.x", "partial", "is_ordered","pseudo.y","transl.except")]


sub_final_gff_joined <- subset(final_gff_joined,(((final_gff_joined$gene.x  != ""|final_gff_joined$gene.y  != "") & (final_gff_joined$gene.y  != final_gff_joined$gene.x))| (final_gff_joined$Ontology_term!="")))



print(final_gff_joined$old_locus_tag)




