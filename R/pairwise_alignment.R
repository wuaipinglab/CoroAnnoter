
pairwise_alignment <-function(genomefile,core_seq) {
  options(stringsAsFactors = FALSE)
  library(Biostrings)

    TRS_name <- gsub("\\.fasta","",genomefile)
    TRS <- readDNAStringSet(genomefile)
    TRS.m <- vmatchPattern(pattern = core_seq,subject = TRS[1],max.mismatch = 0,fixed =FALSE )

    CS <- subseq(TRS[1], start = TRS.m@ends[[1]][1]-TRS.m@width0-2, end = TRS.m@ends[[1]][1] +3)
    names(CS) <- "CS"

    TRS_score <- c("pattern","subject","score")
    for (i in 1:length(TRS)) {
    temp_PA <- pairwiseAlignment(pattern = CS,subject = TRS[i],type="global-local",gapOpening = 12,gapExtension = 2,
                      substitutionMatrix = "BLOSUM50")
    temp_score <- c(as.character(temp_PA@pattern),as.character(temp_PA@subject),as.character(temp_PA@score))
    TRS_score <- rbind(TRS_score,temp_score)
    }
    colnames(TRS_score) <- TRS_score[1,]
    TRS_score <- TRS_score[-1,]
    rownames(TRS_score) <- names(TRS)
    write.csv(TRS_score,paste0(TRS_name,".PAscore.csv"),quote = TRUE,row.names = TRUE)
  }
