library(data.table)
library(neutralitytestr)
library(ggplot2)
library(gridExtra)
library(ape)
library(ggtree)
library(dplyr)
library(CHESS)

source('scripts/utilities.R')  
source('scripts/figure2.R')  
source('scripts/figure3.R')  
source('scripts/figure4.R')  
source('scripts/figure5.R') 
source('scripts/figure6.R') 

punch.cellcounts <- fread('data/figure_1_5/punch_cellcounts.csv')
needle.cellcounts <- fread('data/figure_1_5/needle_cellcounts.csv')

sim.params <-  fread('data/sim_params.csv')

for (sim.id in 1:4) {
  
  sim.params.curr <- sim.params[sim_id == sim.id]
  
  pdf(paste0('output/plots/screenshot_',sim.id,'.pdf'),
      paper = 'USr', width = 20, height = 20)

  GeneratePlotsFig2(sim.id = sim.id,
                    needle.radius = 4,
                    punch.radius = 27,
                    depth = 100,
                    fmin = .05,
                    fmax = .4,
                    fmin.hist = .03,
                    nclonal.muts = 100,
                    path.to.data = 'data/figure_1_5/',
                    path.to.output = 'output/plots')
  dev.off()

  punch.ind <-  if (sim.id == 4) c(1, 3, 4, 5) else 1:6
  punch.files <- paste0('punch_sample_p',punch.ind,'_',sim.id,'.csv')

  needle.ind <- 7:8
  needle.files <- paste0('needle_sample_n',needle.ind,'_',sim.id,'.csv')

  GeneratePlotsFig3(sim.id = sim.id,
                    bulk.files = punch.files,
                    bulk.ind = punch.ind,
                    axis.name = 'Punch ',
                    bulk.cellcounts = punch.cellcounts,
                    nclonal.muts = 100,
                    depth = 100,
                    path.to.data = 'data/figure_1_5/',
                    path.to.output = 'output/plots/punch_scatters_')

  GeneratePlotsFig3(sim.id = sim.id,
                    bulk.files = needle.files,
                    bulk.ind = needle.ind,
                    axis.name = 'Needle ',
                    bulk.cellcounts = punch.cellcounts,
                    nclonal.muts = 100,
                    depth = 100,
                    path.to.data = 'data/figure_1_5/',
                    path.to.output = 'output/plots/needle_scatters_')

  pdf(paste0('output/plots/trees_',sim.id,'.pdf'),
      paper = 'USr', width = 10, height = 20)

  GeneratePlotsFig4(sim.id = sim.id,
                    path.to.data = 'data/figure_1_5/')
  dev.off()

  GeneratePlotsFig5(sim.id = sim.id,
                    grid = 400,
                    nbulks = 4,
                    radius = 27,
                    mut.rates = sim.params.curr$mutation_rate,
                    birth.rates = sim.params.curr$birth_rate,
                    death.rates = sim.params.curr$death_rate,
                    start.times = sim.params.curr$start_time,
                    kill.times = sim.params.curr$kill_time,
                    push.powers = sim.params.curr$push_power,
                    seed = sim.params.curr$seed[1],
                    path.to.data = 'data/figure_1_5/',
                    path.to.output = 'output/plots/')
}

GeneratePlotsFig6(path.to.data = 'data/figure_6/',
                  path.to.output = 'output/plots/',
                  plot.type = 'error_rate')



 

