GeneratePlotsFig5 <- function(sim.id,
                              grid,
                              nbulks,
                              radius,
                              mut.rates,
                              birth.rates,
                              death.rates,
                              start.times,
                              kill.times,
                              push.powers,
                              seed,
                              path.to.data,
                              path.to.output) {
  
  space.clone <- fread(paste0(
    path.to.data,'/space_',sim.id,'.csv'))
  
  x <- y <- grid
  z <- 1
  xc <- c(x/4, x/2, 3*x/4)
  yc <- c(.7*grid, .3*grid)
  
  bulk.ids <- 1:6
  singlecells.strategy <- data.frame()
  bulk.singlecells <- list()
  bulk.count <- 0
  
  for (i in 1:length(xc)) {
    for (j in 1:length(yc)) {
      
      bulk.count <- bulk.count + 1
      if (bulk.count %in% bulk.ids) {
        
        curr.strategy <- SingleCellSamplesFromBulk(
          x1 = xc[i] - radius,
          x2 = xc[i] + radius,
          y1 = yc[j] - radius, 
          y2 = yc[j] + radius, 
          x = x, 
          y = y,
          z = z,
          nsamples = 4,
          space = space.clone)
        singlecells.strategy <- rbind(singlecells.strategy, curr.strategy)
        bulk.singlecells[[length(bulk.singlecells) + 1]] <- curr.strategy
      }
    }
  }

  bulks.tree.sample <- SimulateTumourTree(
    x = grid,
    y = grid,
    z = 1,
    mutation_rates = mut.rates,
    birthrates = birth.rates,
    deathrates = death.rates,
    clone_start_times = start.times,
    kill_regrow_times = kill.times,
    push_power = push.powers,
    aggressions = c(1, 1), 
    father = c(0, 0), 
    seed = seed, 
    strategy = singlecells.strategy)
  
  tree <- ape::collapse.singles(bulks.tree.sample$tree[[1]])
  
  clone.ids <- sapply(strsplit(tree$tip.label, '_'), '[', 2)
  coords.x <- sapply(strsplit(tree$tip.label, '_'), '[', 3)
  coords.y <- sapply(strsplit(tree$tip.label, '_'), '[', 4)
  coords.x <- as.numeric(substr(coords.x, 2, nchar(coords.x))) + 1
  coords.y <- as.numeric(substr(coords.y, 2, nchar(coords.y))) + 1
  
  tip.labels <- list()
  for (i in 1:length(coords.x)) {
    for (b in 1:length(bulk.singlecells)) {
      
      curr.bulk.coords <- bulk.singlecells[[b]]
      
      if (coords.x[i] %in% curr.bulk.coords$x & 
          coords.y[i] %in% curr.bulk.coords$y) {
        tip.labels[[i]] <- paste0('B', bulk.ids[b])
      }
    }
  }
  
  tip.labels.dt <- data.frame(tree$tip.label)
  tip.labels.dt$bulk <- unlist(tip.labels)
  tip.labels.dt$clone <- ifelse(clone.ids %in% c('type1'), 'WT', 'Mutant')
  
  p <- ggtree(tree, size = 1.5)
  p <- p %<+% 
    tip.labels.dt + 
    geom_tippoint(aes(shape = clone,
                      color = bulk), 
                  alpha = 1,
                  size = 10) + 
    theme(legend.position = 'right', 
          legend.title = element_text(size = 20,
                                      face = 'bold'), 
          legend.text = element_text(size = 15))
  
  ggsave(paste0(path.to.output, 'bulk_tree_', sim.id, '.pdf'),
         p, width = 10, height = 7)
}





