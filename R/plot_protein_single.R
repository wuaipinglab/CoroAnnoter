

plot_protein_single <-function(anno_R,TRSlocation) {
  library(ggplot2)

  blastout <-read.csv(anno_R,header = TRUE)
  pathfile <- dirname(anno_R)
  savename <- gsub("\\.csv","",basename(anno_R))

  pdf(paste0(pathfile,"/",savename,".protein.pdf"),
      width = 35,height = 4)

    strain <-blastout
    strain <- strain[order(as.numeric(strain$start)),]
    lseq <- max(strain[,c(17,18)])
    strain$lseq <- lseq
    p <- ggplot()+ylim(0.5,1.5)+ labs(x = "",y="")
    ####add genome ontology
    p <- p + ggplot2::geom_segment(mapping=ggplot2::aes(x = 1, y = 1, xend = lseq+1000, yend = 1))
    p <- p + annotate("text", x = -50, y = 1,label = "5'UTR", hjust = 1, size = 8)
    p <- p + annotate("text", x = lseq+1000, y = 1,label = "3'UTR", hjust = 1, size = 8)
    p <-p+ theme_bw(base_size = 20) + # white background
      theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
      theme(axis.ticks = element_blank(), axis.text.y = element_blank()) +
      theme(panel.border = element_blank())
    ####add protein
    p <- p + geom_rect(aes(xmin = start, xmax = end,ymin = 0.9, ymax = 1.1,fill = stitle),
                       data = strain,show.legend = TRUE)+
      geom_label(aes(x = start + (end -start)/2, y = 1.22, label = stitle),
                 data = strain,size = 8)
    ####add TRS
    trs <- read.csv(TRSlocation)
    p <- p + geom_point(data = trs,aes(x = start, y = 0.85), shape = 21,
                        colour = "black", fill = "black", size = 6, show.legend = FALSE)
    p <- p+geom_text(data = trs,aes(x=start + (end -start)/2,y=0.85,label=score),colour="white")
    p <- p+labs(title = savename)
    plot(p)
  dev.off()
}
