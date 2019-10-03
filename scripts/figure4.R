GeneratePlotsFig4 <- function(sim.id,
                              path.to.data) {
  
  par(mfrow = c(1, 3))
  for (nsc in c(10, 400)) {
    
    tree <- read.tree(paste0(
      path.to.data,'singlecells_sample_',nsc,'sc_',sim.id,'.tre'))
    
    tree <- ape::collapse.singles(tree)
    
    clone.ids <- as.numeric(
      gsub('type', '', sapply(strsplit(tree$tip.label, '_'), '[', 2)))
    
    clone.ids <- ifelse(clone.ids %in% 1, 'n', 's')
    
    tree$tip.label <- clone.ids
    
    plot(tree, show.tip.label = F, edge.width = 2)
    tiplabels(pch = 21, col = 'black',  cex = 2,
              bg = ifelse(tree$tip.label == 's', 'red', 'blue'))
  }
  
  genot.tree <- read.nexus(paste0(
    path.to.data, 'mutations_genot_needle8_400sc_',sim.id,'.txt.nex.tre'))
  
  tips <- genot.tree$tip.label[-1]
  
  clone.ids <- as.numeric(substr(tips, nchar(tips), nchar(tips)))
  clone.ids <- ifelse(clone.ids %in% 1, 'n', 's')
  
  genot.tree$tip.label <- c('n', clone.ids)
  
  plot(genot.tree, show.tip.label = F, edge.width = 2, rotate.tree = T)
  tiplabels(pch = 21, col = 'black',  cex = 2,
            bg = ifelse(genot.tree$tip.label == 's', 'red', 'blue'))
  
}