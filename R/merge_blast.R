

merge_blast <-function(blastdir) {

  fastafils <- list.files(blastdir,pattern = "fasta")
  fastafils <- gsub("\\.fasta","",fastafils)
  for (i in 1:length(fastafils)) {
    xlsfiles <- list.files(blastdir,pattern = fastafils[i])
    xlsfiles <- grep(xlsfiles,pattern = "xls",value = TRUE)

    if (file.info(paste0(blastdir,"/",xlsfiles[1]))[1] > 0) {
    finalfile <- read.delim(paste0(blastdir,"/",xlsfiles[1]),sep = "\t",header = FALSE)
    }
    for (j in 2: length(xlsfiles)) {
      if (file.info(paste0(blastdir,"/",xlsfiles[j]))[1] > 0) {
      tempfile <-  read.delim(paste0(blastdir,"/",xlsfiles[j]),sep = "\t",header = FALSE)
      finalfile <- rbind(finalfile,tempfile)
      }
    }
    colnames(finalfile) <- c("qseqid","sseqid", "pident","length", "mismatch","gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "slen", "stitle")
    finalfile <- cbind(species=rep(fastafils[i],nrow(finalfile)),finalfile)
    write.table(finalfile,paste0(blastdir,"/",fastafils[i],".xls"),sep = "\t",quote = TRUE,row.names = FALSE,col.names = FALSE)
  }
}
