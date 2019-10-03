library(stringr)
library(ape)
library(ggtree)          
library(ggplot2)  
library(dplyr)
library(data.table)

ConvertToNex <- function(mutdata, filename){
  
  mutdata$normal <- 0
  
  mutdata$muts <- NULL
  mutdata <- mutdata[,ncol(mutdata):1]
  
  mutdata[mutdata == 0] <- 'c'
  mutdata[mutdata == 1] <- 't'
  mutdata[] <- lapply(mutdata, factor)
  
  mutdataList <- as.list(as.data.frame(mutdata))
  
  write.nexus.data(mutdataList,
                   paste0(filename, '.nex'), datablock = T)
}

path.to.data <- 'data/supp/inference_sc_organoids/target_data/'

nfiles <- 3
mut.files <- paste0(path.to.data,'pt',1:nfiles,'.subs.csv')
for (i in 1:nfiles) {
  mutdata.org <- data.frame(fread(mut.files[i]))
  mutdata <- 
    mutdata.org[,-c(1:5, which(
      stringr::str_sub(names(mutdata.org), -3) %in% c('DEP')))]
  mutdata.bin <- data.frame(muts = mutdata$ID, 
                            as.data.frame(apply(
                              mutdata[, -1], 2,
                              function(x) replace(x, x > 0, 1))))
  mutdata.bin <- 
    mutdata.bin[, -which(
      substr(colnames(mutdata.bin), 6, 12) %in% c('tissue_', '.tissue'))]
  names(mutdata.bin)[-1] <- substr(names(mutdata.bin)[-1], 1,
                                   nchar(names(mutdata.bin)[-1]) - 4)
  ConvertToNex(mutdata.bin, mut.files[i])
}


trees <- list()
tree.edges <- data.frame()
branching.times <- data.frame()
for (i in 1:nfiles) {
  
  curtree <- read.nexus(paste0(path.to.data,'pt',i,'.subs.csv.nex.tre'))
  trees[[i]] <- curtree
  tree.edges <- rbind(
    tree.edges, data.frame(edge_length = curtree$edge.length,
                           patient = paste0('P', i),
                           drop_normal = 'full_tree'))
  curtree$edge.length[which(curtree$edge.length == 0)] <- 0.0000001
  branching.times <- rbind(
    branching.times, 
    data.frame(branch_times = ape::branching.times(ape::chronopl(curtree,
                                                                 lambda = .1)),
               patient = paste0('P', i),
               drop_normal = 'full_tree'))
  curtree <- drop.tip(phy = curtree,
                      tip = c('normal',
                              curtree$tip.label[substr(
                                curtree$tip.label, 4 ,4) == 'N']))
  tree.edges <- rbind(tree.edges,
                      data.frame(edge_length = curtree$edge.length,
                                 patient = paste0('P', i),
                                 drop_normal = 'drop_normal'))
  branching.times <- rbind(branching.times,
                           data.frame(branch_times = ape::branching.times(
                             ape::chronopl(curtree, lambda = .1)),
                                      patient = paste0('P', i),
                                      drop_normal = 'drop_normal'))
}
class(trees) <- 'multiPhylo'

p <- 
  ggplot(tree.edges, 
            aes(edge_length,
                fill = patient)) + 
  geom_histogram(bins = 30) + 
  facet_wrap(drop_normal~patient, scales = 'free') +
  theme_bw() + 
  scale_fill_brewer(palette = 'Dark2')
ggsave(paste0(path.to.data,'edge_length_distr.pdf'),
       p, width = 9, height = 6)

p <-
ggplot(branching.times,
       aes(branch_times,
           fill = patient)) +
  geom_histogram(bins = 30) + 
  facet_wrap(drop_normal~patient, scales = 'free') +
  theme_bw() + 
  scale_fill_brewer(palette = 'Dark2')
ggsave(paste0(path.to.data,'branching_time_distr.pdf'),
       p, width = 9, height = 6)

p <- ggtree(trees) + facet_wrap(~.id, ncol = 1) + geom_tiplab(size=2)
ggsave(paste0(path.to.data,'organoid_trees.pdf'), p)



trees <- list()
for(patient in 1:3) {
  target.tree <- read.nexus(paste0(
    path.to.data,'/pt',patient,'.subs.csv.nex.tre'))
  target.tree <- drop.tip(
    target.tree,
    c('normal',
      target.tree$tip.label[substr(target.tree$tip.label, 4, 4) == 'N']))
  trees[[patient]] <- target.tree
}
class(trees) <- 'multiPhylo'

p <- ggtree(trees) + facet_wrap(~.id, ncol = 1) + geom_tiplab(size=2)
ggsave(paste0('~/Documents/organoid_trees.pdf'), p)

