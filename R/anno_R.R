anno_R <-
function(blastfile) {
  library(gdata)
  options(stringsAsFactors = FALSE)
  
  savename <- gsub("\\.xls","",basename(blastfile))
  pathfile <- dirname(blastfile)
  
  blastout <- read.csv(blastfile)
  start_end <- t(as.data.frame(strsplit(blastout$location,":")))[,-1]
  start_end <- apply(start_end, 2, as.numeric)
  colnames(start_end) <- c("start","end")
  orf <- gsub("_.*","",blastout$location)
  orf <- gsub("lcl\\|*","",orf)
  orf_l <- (start_end[,2]-start_end[,1])+1
  blastout <- cbind(genus=rep(locate,nrow(blastout)),blastout,orf,start_end,orf_l)
  
  write.csv(blastout,paste0(pathfile,"/",savename,"_anno_R.csv"),quote = TRUE,row.names = FALSE)
}
