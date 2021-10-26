#Filtering 1: Filtering out housekeeping genes etc.
######################
```{r setup, include=FALSE}
knitr::opts_chunk$set(eval = FALSE)

#get the directory where the conserved markers are located,  there should not be other csv files in there
path <- "~/trans_aligned_alevin_15000_output_transindex/DX_H/consv_2515"
files <- list.files(path= path, pattern="*.csv", full.names=TRUE, recursive=FALSE) #contains all the files found in the directory

#library(sjmisc)
#a list of housekeeping genes
HK_genes <- list("RL","RS", "RRN", "ACTB", "GAPDH", "PGK1", "PPIA", "RPL", "ARBP", "B2M", "YWHAZ",
                 "SDHA", "TBP", "RS", "POL","COX", "GFP", "SRSF", "F13A", "MYH", "COX")

#To create a dataframe later, create two empty lists (for the 2 separate functions inside the 
#file iteration loop)
data_function1 <- list()
data_function2 <- list()
# define wd
folderWD <- path
# input output files
for (f in files){
  print(f)
fileInput <- f
fileOutput <- paste0(path, tools::file_path_sans_ext(basename(fileInput)), "_filtered.csv")


#then filter
#take the csv file, get stats: number of "informative" genes, difference between
# pct values --- how many have greater than e.g. 0.4 diff, 0.6 and 0.8

InFile <- read.csv(fileInput)

#Checks the number of genes originally, how many removed, how many left after removing
#housekeeping genes
GenesFiltering <- function(f, outfile){
  removed_genes <- list()
  #original number of genes
  genes_no_orig <- sum(grepl("TRI",f$X))
  writeLines(c("FILTERING"))
  filtered_f <- data.frame()
  for (gene in HK_genes){
    # Remove row based on condition
    f<-f[!grepl(gene,f$X),]
    r<-f[grepl(gene, f$X), ]
    removed_genes <- append(removed_genes, r)
  }
  filtered_f<-f
  write.csv(filtered_f, outfile)
  #number of filtered 
  genes_no_filt <- sum(grepl("TRI",filtered_f$X))
  #the number of genes removed
  gene_diff <- genes_no_orig - genes_no_filt
  #print out the information
  writeLines(c("Original number of genes: ", genes_no_orig, "Number of genes removed: ", gene_diff, "Genes left: ", genes_no_filt, " "))
  infos_1<-c(genes_no_orig, gene_diff, genes_no_filt, NA)
  out <- list(filtered_f, infos_1)
  return(out)
}
#call the the function GenesFiltering
GF1 <- GenesFiltering(InFile, fileOutput)
inFile_filtered <- read.csv(fileOutput)

infos_1 <- GF1[2]


#Filecompare gives information about the pct values of the genes
Filecompare <- function(f,outfile){
  #print(f)
  #files <- list(f1,f2)
  writeLines(c("FILE STATS", " "))
  #for (f in files){
    genes_no <- sum(grepl("TRI",f$X))
    f$Diff <- array(f$X0_pct.1 - f$X0_pct.2)
    #print(f$Diff)
    diff_list <- list()
    a<-0; b<-0; c<-0

    infocolname = list("DX_H", "DI")
    #info <- rbind(pct.score.diff = rownames(info), info)
    
    i=1
    #go over the genes, keep count of how many of them fulfill the criteria below
    for (v in 1:nrow(f)){
      a <- sum(f$Diff >= 0.8)
      b <- sum(f$Diff >= 0.6 & f$Diff <= 0.8)
      c <- sum(f$Diff >= 0.4 & f$Diff <= 0.6)
      d <- sum(f$Diff <= 0.4)
    }
    r<-c(a,b,c)
   #print out the info on how many of the genes fall under the categories below
    writeLines(c("pct difference >= 0.8: ", a, "0.8 => pct difference >= 0.6: ",b,"0.6 => pct difference >= 0.4: ", c, "less than 0.4: ", d, " ")) 
    infos_2<- c(a,b,c,d)
    data_function2 <- append(data_function2, a)
    data_function2 <- append(data_function2, b)
    data_function2 <- append(data_function2, c)
    data_function2 <- append(data_function2, d)

    
    print(" ")
    
    #print(infocols)
    #name <- infocolname[i]
    #info[,i] <- infocols
    i=i+1
    #infocols<list()
    #diff_list<-append(diff_list, init_list)
    #inforows <- c("pct-difference >= 0.8", "0.8 => pct-difference >= 0.6","0.6 => pct-difference >= 0.4") 
    # return(out_df2)
    return(infos_2)
  }
#}
out <- Filecompare(inFile_filtered, fileOutput)
infos_2 <- out
write.csv(data_function2, fileOutput)
all_infos <- c(infos_1,infos_2)
data_function1 <- append(data_function1, all_infos) 

}

```
######################

