

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
    ####Layered
    strain$order <- 1
    for(i in 2:nrow(strain)) {
      temp.location <- c(strain$end[1:(i-1)],strain$start[i])
      temp.location <- sort(temp.location)
      dif <- i- which(temp.location==strain$start[i])
      strain$order[i] <- 1+dif
    }
    ####
    order.l <- strain$order[!duplicated(strain$order)]
    if (any(order.l %in% 1)) {
      p <- p + geom_rect(aes(xmin = start, xmax = end,ymin = 0.9, ymax = 1.1,fill = orf),
                         data = strain[strain$order==1,],show.legend = TRUE)+
        geom_text(aes(x = start + (end -start)/2, y = 1, label = orf),
                   data = strain[strain$order==1,],size = 8)

    }

    if (any(order.l %in% 2)) {
      p <- p + geom_rect(aes(xmin = start, xmax = end,ymin = 0.65, ymax = 0.85,fill = orf),
                         data = strain[strain$order==2,],show.legend = TRUE)+
        geom_text(aes(x = start + (end -start)/2, y = 0.75, label = orf),
                   data = strain[strain$order==2,],size = 8)

    }
    if (any(order.l %in% 3)) {
      p <- p + geom_rect(aes(xmin = start, xmax = end,ymin = 0.4, ymax = 0.6,fill = orf),
                         data = strain[strain$order==3,],show.legend = TRUE)+
        geom_text(aes(x = start + (end -start)/2, y = 0.5, label = orf),
                  data = strain[strain$order==3,],size = 8)

    }

    ####add TRS
    trs <- read.csv(TRSlocation)
    max.order <- max(order.l)

    p <- p + geom_point(data = trs,aes(x = start, y = 1.1-0.25*max.order), shape = 21,
                        colour = "NA", fill = "black", size = 6, show.legend = FALSE)+
            geom_text(data = trs,aes(x=start + (end -start)/2,y=1.1-0.25*max.order,
                                     label=score),colour="white")
    p <- p+labs(title = savename)
    plot(p)
  dev.off()
}
