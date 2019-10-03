library(data.table)
library(CHESS)
library(ape)

source('scripts/utilities.R')

sim.params <- fread('data/sim_params.csv')

grid <- 400

punch.coords <- PunchBiopsyCoords(radius = 27, grid = grid)
needle.coords <- NeedleBiopsyCoords(radius = 4, grid = grid)
sampling.strategy <- GetSamplingStrategy(grid, punch.coords, needle.coords)

for (sim.id in 1:4) {
  
  sim.tum <- sim.params[sim_id == sim.id]
  print(sim.tum)
  grid <- sim.tum$grid[1]
  seed <- sim.tum$seed[1]
  birth.rates <- sim.tum$birth_rate
  death.rates <- sim.tum$death_rate
  clone.start.times <- sim.tum$start_time
  push.powers <-sim.tum$push_power
  kill.regrow.times <- sim.tum$kill_time
  
  path.to.data <- 'data/supp/neutral_death/'
  
  # Bulk samples ------------------------------------------------------------
  sim.data <- SimulateTumourSample(
    x = grid,
    y = grid,
    z = 1,
    birthrates = birth.rates,
    deathrates = death.rates,
    aggressions = c(1, 1),
    push_power = push.powers,
    mutation_rates = c(10, 10), 
    clone_start_times = clone.start.times,
    kill_regrow_times = kill.regrow.times,
    father = c(0, 0), 
    seed = seed,
    depth = 100, 
    min_reads = 5, 
    depth_model = 1,
    sample_strategy = sampling.strategy)
  
  fwrite(sim.data[[1]]$mutation_data,
         paste0(path.to.data, 'vaf_wholetumour_', sim.id, '.csv'))
  
  for (p in 2:7) {
    fwrite(sim.data[[p]]$mutation_data,
           paste0(path.to.data, 'punch_sample_p', p - 1, '_', sim.id, '.csv'))
  }
  
  for (n in 8:9) {
    fwrite(sim.data[[n]]$mutation_data,
           paste0(path.to.data, 'needle_sample_n', n - 1, '_', sim.id, '.csv'))
  }
  
  # Single cell samples -----------------------------------------------------
  space.clone <- as.data.frame(fread(
    paste0(path.to.data,'space_',sim.id,'.csv')))
  space.clone <- t(apply(space.clone, 2, rev))
  space.clone[is.na(space.clone)] <- 0
  
  for (nsc in c(10, 400)) {
    
    sampling.strategy <- SingleCellSampleCoords(
      x = grid, y = grid, z = 1, nsamples = nsc, space = space.clone)
    
    tree.sample <- SimulateTumourTree(
      x = grid, 
      y = grid,
      z = 1,
      birthrates = birth.rates,
      deathrates = death.rates,
      aggressions = c(1, 1),
      push_power = push.powers,
      mutation_rates = c(10, 10),
      clone_start_times = clone.start.times,
      kill_regrow_times = kill.regrow.times,
      father = c(0, 0),
      seed = seed,
      strategy = sampling.strategy)
    
    sim.tree <- ape::collapse.singles(tree.sample$tree[[1]])
    ape::write.tree(phy = sim.tree, file = paste0(
      path.to.data, 'singlecells_sample_',nsc, 'sc_',sim.id,'.tre'))
    
    ncells.to.genot <- nsc
    sc.coords.matrix <- sampling.strategy
    
    sampling.strategy <- needle.coords[2]
    for(i in 1:ncells.to.genot){
      sampling.strategy[[length(sampling.strategy) + 1]] <- 
        as.numeric(c(sc.coords.matrix[i, ], sc.coords.matrix[i,]))
    }
    
    genot.cells.mutdata <- SimulateTumourSample(
      x = grid, 
      y = grid,
      z = 1,
      birthrates = birth.rates,
      deathrates = death.rates,
      aggressions = c(1, 1),
      push_power = push.powers,
      mutation_rates = c(10, 10),
      clone_start_times = clone.start.times,
      kill_regrow_times = kill.regrow.times,
      father = c(0, 0),
      seed = seed,
      depth = 100, 
      min_reads = 5, 
      depth_model = 1,
      sample_strategy = sampling.strategy
    )
    
    genot.ref.muts <- genot.cells.mutdata[[1]]$mutation_data
    genot.matrix <- data.frame(muts = as.character(genot.ref.muts$id))
    genot.sample.names <- c()
    
    clone.ids <- rep(0, ncells.to.genot - 1)
    
    for (sc in 2:length(genot.cells.mutdata)) {
      
      currcell.mutdata <- genot.cells.mutdata[[sc]]$mutation_data
      
      genot.sample <- ifelse(
        as.character(genot.matrix$muts) %in% currcell.mutdata$id, 1, 0)
      
      cur.sc.clone.id <- unique(currcell.mutdata$clone)
      
      clone.ids[sc - 1] <- cur.sc.clone.id[length(cur.sc.clone.id)]
      
      genot.matrix <- cbind(genot.matrix, genot.sample)
      genot.sample.names <- c(genot.sample.names,
                              paste0('s', sc - 1, '_c', clone.ids[sc - 1]))
    }
    
    colnames(genot.matrix) <- c('muts', genot.sample.names)
    genot.matrix$muts <- 1:nrow(genot.matrix)
    genot.matrix <- data.frame(apply(genot.matrix, 2, 
                                     function(x) as.integer(as.character(x))))
    
    genot.tree.filename <- paste0(
      'mutations_genot_needle8_',ncells.to.genot,'sc_',sim.id,'.txt')
    
    fwrite(genot.matrix, paste0(path.to.data, genot.tree.filename))
    
    ConvertMutsBinaryToNexus(paste0(path.to.data, genot.tree.filename))
  }
  
}


# Stochasticity - Sampling bias -------------------------------------------

for (tum in 1:4) {
  
  sim.tum <- sim.params[sim_id == tum]
  
  grid <- sim.tum$grid[1]
  birth.rates <- sim.tum$birth_rate
  death.rates <- sim.tum$death_rate
  clone.start.times <- sim.tum$start_time
  push.powers <-sim.tum$push_power
  kill.regrow.times <- sim.tum$kill_time
  
  path.to.data <- paste0('data/supp/stocheffect/T',tum,'/')
  
  for (sim.id in 1:100) {
    
    print(sim.id)
    
    seed <- sim.id #runif(1)*1e10
    
    sim.data <- SimulateTumourSample(
      x = grid,
      y = grid,
      z = 1,
      birthrates = birth.rates,
      deathrates = death.rates,
      aggressions = c(1, 1),
      push_power = push.powers,
      mutation_rates = c(10, 10), 
      clone_start_times = clone.start.times,
      kill_regrow_times = kill.regrow.times,
      father = c(0, 0), 
      seed = seed,
      depth = 100, 
      min_reads = 5, 
      depth_model = 1,
      sample_strategy = sampling.strategy)
    
    fwrite(sim.data[[1]]$mutation_data,
           paste0(path.to.data, 'vaf_wholetumour_', sim.id, '.csv'))
    
    for (p in 2:7) {
      fwrite(sim.data[[p]]$mutation_data,
             paste0(path.to.data, 'punch_sample_p', p - 1, '_', sim.id, '.csv'))
    }
    
    for (n in 8:9) {
      fwrite(sim.data[[n]]$mutation_data,
             paste0(path.to.data, 'needle_sample_n', n - 1, '_', sim.id, '.csv'))
    }
  }
  
}
