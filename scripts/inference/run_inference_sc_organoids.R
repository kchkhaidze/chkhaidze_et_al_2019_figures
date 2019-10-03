library(CHESS)
library(ape)
library(data.table)
library(transport)
library(magrittr)

args <- commandArgs(trailingOnly = TRUE)
#args <- c(1, 200, 20, 8, 200)

patient <- as.numeric(args[1])
N <- as.numeric(args[2])
rounds <- as.numeric(args[3])
ncores <- as.numeric(args[4])
grid <- as.numeric(args[5])

# path.to.data <- 'data/supp/inference_sc_organoids/target_data'
# path.to.output <- 'data/supp/inference_sc_organoids/inferred_params/BT/'

path.to.data <- '../data/target_data/'
path.to.output <- '../data/inferred_params/BL/'

target.tree <- read.nexus(paste0(
  path.to.data,'/pt',patient,'.subs.csv.nex.tre'))
target.tree <- drop.tip(
  target.tree,
  c('normal',
    target.tree$tip.label[substr(target.tree$tip.label, 4, 4) == 'N']))

target.tree.size <- length(target.tree$tip.label)

params <- data.table(
  name = c('mu', 't', 's', 'd', 'a'),
  min = c(1, 4, 1, 0, .1),
  max = c(3000, 50, 4.5, .9, 1),
  fix = c(T, T, F, F, F))

params.to.rec <- 's_d_a'

# vactors of inferred posterior mode values. Set the corresponding fix=T
# in the params table above and rerun the inference to recover 
# the remaining parameters:
mu.rls <- c(1096, 308, 137)
t.rls <- c(12, 19, 27)
d.rls <- c(.2, .4, .5)
s.rls <- c()
a.rls <- c()


ABCSMCwithTreeSamplesBL(
  sim.id = patient,
  N = N,
  rounds = rounds,
  ncores = ncores,
  grid = grid,
  params = params,
  nsamples = target.tree.size,
  target.tree = target.tree,
  target.popsize = .5*3.14*((grid/2)^2),
  params.to.rec = params.to.rec,
  mu.rl = mu.rls[patient],
  s.rl = s.rls[patient],
  t.rl = t.rls[patient],
  d.rl = d.rls[patient],
  a.rl = a.rls[patient],
  output.dir = path.to.output)

