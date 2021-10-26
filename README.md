# Transcriptomics-project
A spring project where I analysed transcriptomic data extracted from salamander spleen. 
title: "Project Readme"
author: "Atte Rasanen"
date: "5/23/2021"
---

#all relevant packages were installed via conda
#R version: 4.0.5
#Package versions used
#salmon                    v.1.4.0
#r-seurat                  4.0.1
#r-seuratobject            4.0.0
#r-seuratwrappers          v20210208
#r-monocle3                1.0.0
#bioconductor-tximport     1.18.0
#r-harmony                 0.1
#r-cowplot                 1.1.1
#r-sjmisc                  2.8.6
#bioconductor-celda        1.6.1
#bioconductor-slingshot    1.8.0

#Performed the initial Salmon step on Lunarc HPC cluster via job submission scripts
#created the genes2genes file when running salmon's alevin
```

grep ">" trinity_pwal_v2.supertrans.w.fluors.fa > genes
sed -i 's/>//g' genes
paste genes genes > gene2gene.tsv
```

#indexing the supertranscript file on salmon before the transcript level abundancies can be quantified. 
#the index helps salmon to to map the reads to the reference transcriptome used via quasi-mapping which should reduce
#the occurence mapping the same read multiple times
#The job submission script for Lunarc
#######
```
#!/bin/bash
#SBATCH --tasks-per-node=1                               # Request 1 cores
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -J indexing_salmon                       #name job
#SBATCH -t 1-00:00                         # Runtime in D-HH:MM format
#SBATCH -A {lunarc code here}
#SBATCH -p lu
#SBATCH --mem-per-cpu=31000                  # Memory per cpu 200GB total
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=email_address  # Email to which notifications will be sent

#enter wd
/home/as0400/Documents/salmon-latest_linux_x86_64/bin


#make the index
salmon index -i index -k 31 -p 4 -t trinity_pwal_v2.supertrans.w.fluors.fa

```
######################
#The Job submission script on Lunarc for running Alevin to quantify the scRNA data of the salamander spleen
```
#!/bin/bash
#SBATCH --tasks-per-node=10                # Request 1 cores
#SBATCH -N 1                # Request one node (if you request more than one core with -c, also using
                      # -N 1 means all cores will be on the same node)
#SBATCH -J Alevin_as0400                       #name job
#SBATCH -t 3-00:00                         # Runtime in D-HH:MM format
#SBATCH -A LU2021-2-7
#SBATCH -p lu
#SBATCH --mem-per-cpu=9000                  # Memory per cpu 200GB total
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=at0400ra-s@student.lu.se   # Email to which notifications will be sent

#enter wd
cd /home/as0400/Documents/salmon-latest_linux_x86_64/bin

salmon alevin -lIU -1 {CB+UMI files} -2 {read sequence files} --chromiumV3 -i /lunarc/nobackup/projects/regen_immuno/tools/Alevin/salmon-latest_linux_x86_64/bin/index_trans -p 20 -o trans_aligned_alevin_15000_output --tgMap /lunarc/nobackup/projects/a_own/txp2gene.tsv --expectCells 15000


```
######################

#-l = library -1=CB+UMI file(s), alevin requires the path to the FASTQ file containing CB+UMI raw sequences to be given under this #command line flag   -2=read-sequence #file(s)  -i = index(created before via salmon) -p=threads -o=outputfile --chromiumV3 =the type #of single cell protocol of the seq library of the input file   tgMap = #transcript to gene map file (tsv), has 2 columns mapping each #transcript in the reference to the corresponding gene 
#expectCells = "All barcodes whose total UMI counts exceed m/10 are called as cells”, where m is the frequency of the top 1% #cells #which is specified by the parameter #of this command line flag.

#need to annotate the columns(contigs) before moving onto seurat:
```
######################
#Download the script for annotation
wget https://raw.githubusercontent.com/trinityrnaseq/trinityrnaseq/master/Analysis/DifferentialExpression/rename_matrix_feature_identifiers.pl

sed -i 's/_/-/g' quants_mat_cols.txt 
perl rename_matrix_feature_identifiers.pl quants_mat_cols.txt iNewt-gene-annot-mapping > new.quants.txt

#iNewt-gene-annot-mapping contains the annoteted contigs to check if the contigs we have in quants_mat_cols.txt match
#new.quants.txt was transferred to local when creating a directory structure for the seurat step and renamed
# new.quants.txt as quants_mat_cols.txt
```
######################

#in order to create the seurat object in R, the following directory structure was used (so the subfolders are alevin, aux_info, #libParams, logs:

#Working directory:
#  alevin:
#    -modifiedinputclusters.csv (data from souporcell)
#    -featureDump.txt
#    -Alevin.log
#    -quants_mat.gz
#    -quants_mat_cols.txt
#    -quants_mat_rows.txt
#    -quants_tier_mat.gz
#    -whitelist.txt
#  aux_info:
#    -alevin_meta_info.json
#    -ambig_info.tsv
#    -expected_bias.gz
#    -fld.gz
#    -meta_info.json
#    -observed_bias_3p.gz 
#    -observed_bias.gz
#  cmd_info.json
#  lib_format_counts.json
#  
#  libParams:
#    -flenDist.txt
#  
#  logs:
#    -salmon_quant.log




