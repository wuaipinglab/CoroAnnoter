# CoroAnnoter
A semi-automatic coronavirus genome annotation tool

#pipline

####Data Preparation####

#1. Coronavirus genome

#2. blastdb

####soft Preparation####

#R, MEME, blast, ORFfinder

####annotation####

####ORF finder####

ORFfinder -in ${samplename}.fasta -g 1 -s 0 -ml 60 -out ${samplename}.ORF -outfmt 0

####blast####

blastp -query ${samplename}.ORF -db ${blastdb} -evalue 1e-2 -max_target_seqs 3 -out ${samplename}.blast_out.xls -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen stitle qcov"

####Manually remove redundant blast results####

source("./R/merge_blast.R")

####The blastdir contain genome and BLAST results：SARS_CoV-2.fasta, SARS_CoV-2.blast.*.xls

source("./R/merge_blast.R")

merge_blast("./inst/extdata")

####The SARS_CoV-2.xls file were created. 

####Manually process blast results to remove redundant ID.

####get TRS####

source("./R/getTRS.R")

GetTRS(blastfile = "./inst/extdata/Manual-SARS_CoV-2.xls", 
       genomefile = "./inst/extdata/SARS_CoV-2.fasta")
       
####MEME prediction motif####

system("sh ./R/meme.sh ./inst/extdata/TRS")

####draw protein####

source("./R/pairwise_alignment.R")

pairwise_alignment(core_seq = "ACGAAC",
                   genomefile = "./inst/extdata/TRS/Manual-SARS_CoV-2.TRS.fasta")

source("./R/get_motif_location.R")

get_motif_location(genomefile = "./inst/extdata/SARS_CoV-2.fasta","ACGAAC")

source("./R/plot_protein_single.R")

plot_protein_single(anno_R = "./inst/extdata/Manual-SARS_CoV-2_anno.csv",
                    TRSlocation = "./inst/extdata/TRS/SARS_CoV-2.TRS")
