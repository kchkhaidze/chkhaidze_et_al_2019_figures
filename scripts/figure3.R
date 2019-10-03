GeneratePlotsFig3 <- function(sim.id,
                              bulk.files,
                              bulk.ind,
                              axis.name,
                              bulk.cellcounts,
                              nclonal.muts,
                              depth,
                              path.to.data,
                              path.to.output) {
  
  scatter.plots <- list()
  nsamples <- length(bulk.files)
  for (i in 1:nsamples) {
    for (j in 1:nsamples) {
      if (j+i <= nsamples) {
        
        sample.x <- fread(paste0(path.to.data, bulk.files[i]))
        sample.y <- fread(paste0(path.to.data, bulk.files[i + j]))
        
        sample.x <- GetVAF(sample.x,
                           nclonal.muts,
                           depth,
                           bulk.cellcounts[i])
        sample.y <- GetVAF(sample.y,
                           nclonal.muts,
                           depth,
                           bulk.cellcounts[i + j])
        
        fHrs <- sum(sample.x$vaf > .2)/sum(sample.x$vaf > .08) + 
          sum(sample.y$vaf > .2)/sum(sample.y$vaf > .08)
        
        ksd <- max(abs(sapply(
          seq(0,1,.01), function(x) sum(sample.x$vaf <= x)) - 
            sapply(seq(0,1,.01), function(x) sum(sample.y$vaf <= x))))
        
        sample.xy <- merge(sample.x, sample.y, by = 'id', all = T)
        for (s in 1:nrow(sample.xy)) {
          
          if (is.na(sample.xy$clone.x[s])) {
            
            sample.xy$vaf.x[s] <- 0
            sample.xy$clone.x[s] <- sample.xy$clone.y[s]
            
          } else if (is.na(sample.xy$clone.y[s])) {
            
            sample.xy$vaf.y[s] <- 0
            sample.xy$clone.y[s] <- sample.xy$clone.x[s]
          }
        }
        bs.sfs <- data.table(x = sample.xy$vaf.x,
                             y = sample.xy$vaf.y,
                             clone = sample.xy$clone.x)
        #bs.sfs <- bs.sfs[x > 0 & y > 0]
        
        scatter.plots[[length(scatter.plots) + 1]] <-
          
          ggplot(bs.sfs, aes(x, y, fill = clone, color = clone)) +
          geom_point(alpha = .5, size = 3) +
          scale_colour_manual(values = GetCloneColours()) +
          scale_fill_manual(values = GetCloneColours()) +
          xlab(paste(axis.name, bulk.ind[i])) +
          ylab(paste(axis.name, bulk.ind[j + i])) +
          theme_bw() + 
          theme(legend.position = 'none') +
          theme(axis.text = element_text(size = 15),
                axis.title = element_text(size = 15))
      }
    }
  }
  scatter.plots.glist <- lapply(scatter.plots, ggplotGrob)
  ggsave(paste0(path.to.output, sim.id,'.pdf'),
         marrangeGrob(scatter.plots.glist, nrow = 3, ncol = 3),
         width = 25.4, height = 20.28, units = 'cm')
}