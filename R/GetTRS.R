

GetTRS <-function(blastfile,genomefile) {

  library(Biostrings)
  library(gdata)
  options(stringsAsFactors = FALSE)

  savename <- gsub("\\.xls","",basename(blastfile))
  pathfile <- dirname(blastfile)
  if(!dir.exists(paste0(pathfile,"/TRS"))) dir.create(paste0(pathfile,"/TRS"))

  ####get location
  blastout <- read.xls(blastfile,header = TRUE)
  start_end <- t(as.data.frame(strsplit(blastout$qseqid,":")))[,-1]
  start_end <- apply(start_end, 2, as.numeric)
  colnames(start_end) <- c("start","end")
  orf <- gsub("_.*","",blastout$qseqid)
  orf <- gsub("lcl\\|*","",orf)
  orf_l <- (start_end[,2]-start_end[,1])+1
  blastout <- cbind(blastout,orf,start_end,orf_l)
  write.csv(blastout,paste0(pathfile,"/",savename,"_anno.csv"),quote = TRUE,row.names = FALSE)

  ####get TRS
  blastout0 <- blastout[!duplicated(blastout$orf),]
  blastout0 <- blastout0[order(blastout0$start,decreasing = FALSE),]
  genome <- readDNAStringSet(genomefile,format = "fasta")

  TRS <-   leader <- subseq(genome[1], start = 1, end = 250)
  for (i in 2:nrow(blastout0)) {
    temp_body <- subseq(genome[1], start = blastout0$start[i]-99, end = blastout0$start[i])
    TRS <- c(TRS,temp_body)

  }
  names(TRS) <- blastout0$stitle
  names(TRS)[1] <- "leader"
  writeXStringSet(TRS,paste0(pathfile,"/TRS/",savename,".TRS.fasta"))

}
