get_motif_location <-
function(genomefile,motif){
  options(stringsAsFactors = FALSE)
  library(Biostrings)
  genome_name <- gsub("\\.fasta","",basename(genomefile))
  pathfile <- dirname(genomefile)
  genome <- readDNAStringSet(genomefile,format = "fasta")
  CS_m <- vmatchPattern(pattern = motif,subject = genome,
                       max.mismatch = 0,fixed =FALSE )
  m <- as.data.frame(CS_m)

  CS <- subseq(genome[1], start = CS_m@ends[[1]][1]-CS_m@width0-2, end = CS_m@ends[[1]][1] +3)
  names(CS) <- "CS" 
  m$group_name <- as.data.frame(CS)$x
  
  TRS_score <- c("subject","score")
  for (i in 1:nrow(m)) {
    temp_body <- subseq(genome[1], start = CS_m@ends[[1]][i]-CS_m@width0-2, end = CS_m@ends[[1]][i] +3)
    temp_PA <- pairwiseAlignment(pattern = CS,subject = temp_body,type="global-local",gapOpening = 12,gapExtension = 2,
                                 substitutionMatrix = "BLOSUM50")  
    temp_score <- c(as.character(temp_PA@subject),as.character(temp_PA@score))
    TRS_score <- rbind(TRS_score,temp_score)
  }
  colnames(TRS_score) <- TRS_score[1,]
  TRS_score <- TRS_score[-1,]
  
  TRS_final <- cbind(m,TRS_score)
  write.csv(TRS_final,paste0(pathfile,"/TRS/",genome_name,".TRS"),
            quote = FALSE,row.names = FALSE)
}
