library(data.table)
library(dplyr)
library(ggplot2)
library(ggplus)
library(ggsci)
library(reshape2)
library(gridExtra)
library(RColorBrewer)
library(ape)
library(phylobase)
library(adephylo)
library(spdep)
library(neutralitytestr)

source('scripts/utilities.R')

path.to.output <- 'output/plots/supp/'

# Figure S1, S9 ----------------------------------------------------------

path.to.data <- 'data/supp/'
sim.ids <- c(11, 22, 33, 44)

growth.curves <- data.frame()

for (sim.id in sim.ids) {
  
  cell.counts <- fread(paste0(
    path.to.data,'timestep_cellcount_',sim.id,'.csv'))
  
  cell.counts[, c('total', 'time') := list(type1 + type2, round(time,2))]
  cell.counts <- cell.counts[!duplicated(cell.counts$time),]
  cell.counts[, tumour := paste0('T',sim.id)]
  
  growth.curves <- rbind(growth.curves, cell.counts)
}

growth.curves.melt <- reshape2::melt(growth.curves,
                                     id = c('tumour', 'time'))
text.size <- 15
p <- ggplot(growth.curves.melt, 
            aes(time, value, colour = variable)) +
  geom_line() + 
  facet_wrap(~tumour, scales = 'free') +
  scale_colour_manual('Population', 
                      labels = c('WT',
                                 'Mutant',
                                 'Total'),
                      values = c(type1 = 'blue',
                                 type2 = 'red',
                                 total = 'black')) +
  labs(x = 'Gillespie time steps',
       y = 'Number of cells') +
  theme_bw() +
  theme(axis.text = element_text(size = text.size),
        axis.title = element_text(size = text.size),
        legend.text = element_text(size = text.size),
        legend.title = element_text(size = text.size, face = 'bold'),
        strip.text = element_text(size = text.size))

ggsave(paste0(path.to.output,'growth_curves_regrow.pdf'), 
       p, width = 10, height = 8)


# Figure S2, S4, S7, S8 ------------------------------------------------------

source('scripts/figure2.R')
sim.ids <- c(111, 112, 113, 114)

for (sim.id in sim.ids) {
  
  pdf(paste0(path.to.output,'/screenshot_',sim.id,'.pdf'),
      paper = 'USr', width = 20, height = 20)
  GeneratePlotsFig2(sim.id = sim.id,
                    needle.radius = 4,
                    punch.radius = 27,
                    depth = 100,
                    fmin = .05,
                    fmax = .4,
                    fmin.hist = .03,
                    nclonal.muts = 100,
                    path.to.data = path.to.data,
                    path.to.output = path.to.output)
  dev.off()
}


# Figure S3 ---------------------------------------------------------------

radius <- 10
grid <- 400
path.to.data <- 'data/figure_1_5/'

concentric.samples.stack <- data.table()

for (sim.id in 1:4) {
  
  space.clone <- fread(paste0(
    path.to.data,'/space_',sim.id,'.csv'))
  
  space.generation <- fread(paste0(
    path.to.data,'/space_generation_',sim.id,'.csv'))
  
  cell.counts <- fread(paste0(
    path.to.data,'/timestep_cellcount_',sim.id,'.csv'))
  
  max.wt.cells <- max(cell.counts$type1)
  max.mt.cells <- max(cell.counts$type2)
  clone.sizes <- data.table(wt = max.wt.cells, 
                            mutant = max.mt.cells)
  space.wt <- space.mt <- space.generation
  space.wt[space.clone == 2] <- NA
  space.mt[space.clone == 1] <- NA
  
  space.clone[is.na(space.clone)] <- 0
  
  x <- y <- nrow(space.clone)
  space.grid <- matrix(NA, nrow = x, ncol = y)
  
  blk.ind <- c('4W','3W','2W','1C','2E','3E','4E','4S','3S','2S','2N','3N','4N')
  
  base.seq <- seq(50, 350, 50)
  sample.coords.x <- c(base.seq, rep(200, 6))
  sample.coords.y <- c(rep(200, 7), base.seq[base.seq != 200])
  
  cell.type.counts <- data.frame()
  
  for (i in 1:length(blk.ind)) {
    
    x.coords <- (sample.coords.x[i] - radius + 1):(sample.coords.x[i] + radius)
    y.coords <- (sample.coords.y[i] - radius + 1):(sample.coords.y[i] + radius)
    space.grid[x.coords, y.coords] <- 3
    cells.in.sample <- space.clone[x.coords, ..y.coords]
    cell.type.counts <- rbind(
      cell.type.counts,
      data.frame(wt = sum(cells.in.sample == 1),
                 mutant = sum(cells.in.sample == 2)))
    
  }
  
  pdf(paste0('output/plots/supp/screenshot_concentric_',sim.id,'.pdf'),
      paper = 'USr', width = 20, height = 20)
  par(pty = 's')
  image(1:nrow(space.wt),
        1:ncol(space.wt),
        as.matrix(space.wt),
        col = ColorRampAlpha(c('blue4', 'blue'), n = 10, alpha = 1),
        xlab = 'x[cells]',
        ylab = 'y[cells]',
        cex.axis = 1.5, cex.lab = 1.8, cex.main = 1.5)
  image(1:nrow(space.mt),
        1:ncol(space.mt),
        as.matrix(space.mt),
        col = ColorRampAlpha(c('red4', 'red'), n = 10, alpha = 1),
        add = T)
  image(1:nrow(space.grid),
        1:ncol(space.grid),
        as.matrix(space.grid),
        col = scales::alpha('black', 0.4),
        add = T)
  text(sample.coords.x,
       sample.coords.y + 20,
       labels = blk.ind,
       col = scales::alpha('white', 0.5), 
       cex = 1.5)
  dev.off()
  
  sample.files <- paste0(
    path.to.data,'/concentric_samples/T',sim.id,'_sample_',1:13,'.csv')
  samples.stack <- data.frame()
  for (i in 1:length(blk.ind)) {
    sample.mutdata <- fread(sample.files[i])
    if (nrow(sample.mutdata) > 0) {
      sample.mutdata <- CleanMutdata(sample.mutdata,
                                     100,
                                     100,
                                     cell.type.counts[i,])
      sample.mutdata[, c('sample', 'tumour') := 
                       list(blk.ind[i], paste0('T',sim.id))]
      samples.stack <- rbind(samples.stack, sample.mutdata)
    }
  }
  concentric.samples.stack <- rbind(concentric.samples.stack,
                                    samples.stack)
}

mut.clone.type.count.per.sample <- 
  concentric.samples.stack %>% 
  group_by(tumour, sample) %>% 
  summarise(count = n())

models <- c('neutral',
            'selective',
            'neutral peripheral',
            'selective peripheral')

mut.clone.type.count.per.sample$model <- rep(models, c(13,13,13,11))

text.size <- 50
p <- ggplot(mut.clone.type.count.per.sample,
            aes(x = tumour,
                y = count)) +
  geom_point(aes(colour = sample),
             size = 30,
             alpha = .7) ++
  theme_bw() +
  scale_colour_manual(values = colorRampPalette(
    brewer.pal(13, 'Spectral'))(13)) + 
  labs(x = 'Tumour',
       y = 'Number of mutations',
       colour = 'Sample') +
  theme(axis.text = element_text(size = text.size),
        axis.title = element_text(size = text.size),
        legend.text = element_text(size = text.size),
        legend.title = element_text(size = text.size,
                                    face = 'bold'),
        legend.spacing.y = unit(.5, 'cm'),
        strip.text = element_text(size = text.size))

ggsave(paste0(path.to.output, 'mutcount_concentricsamples.pdf'),
       p, width = 45, height = 20)


# Figure S5 ---------------------------------------------------------------

N <- 100
tumours <- paste0('T', c(1:4, 303))
import.vafs <- F

if (import.vafs) {
  
  vafs.stacked <- data.table()
  
  for (tum in tumours) {
    
    path.to.data <- paste0('data/supp/stocheffect/',tum,'/')
    
    for (sim.id in 1:N) {
      
      mutdata.wholet <- fread(paste0(
        path.to.data, 'vaf_wholetumour_',sim.id,'.csv'))
      
      if (tum != 'T303') {
        
        mutdata.wholet <- data.table(
          clone = mutdata.wholet$clone_id,
          alt = mutdata.wholet$alt_count,
          depth = mutdata.wholet$depth,
          id = mutdata.wholet$mutation_id)
      }
      
      mutdata.wholet <- AddClonalMuts(mutdata.wholet, 100, 100)
      mutdata.wholet[, vaf := alt/depth]
      
      vafs.stacked <- rbind(vafs.stacked,
                             data.table(tumour = tum,
                                        sim_id = sim.id,
                                        vaf = mutdata.wholet$vaf,
                                        sample = 'wholet'))
      
      punch.files <- paste0('punch_sample_p',1:6,'_',sim.id,'.csv')
      
      for (i in 1:length(punch.files)) {
        
        mutdata.punch <- fread(paste0(path.to.data, punch.files[i]))
        mutdata.punch <- AddClonalMuts(mutdata.punch, 100, 100)
        mutdata.punch[, vaf := alt/depth]
        
        vafs.stacked <- rbind(vafs.stacked,
                               data.table(tumour = tum,
                                          sim_id = sim.id,
                                          vaf = mutdata.punch$vaf,
                                          sample = paste0('punch', i)))
      }
      
      needle.files <- paste0('needle_sample_n',7:8,'_',sim.id,'.csv')
      
      for(i in 1:length(needle.files)){
        
        mutdata.needle <- fread(paste0(path.to.data, needle.files[i]))
        mutdata.needle <- AddClonalMuts(mutdata.needle, 100, 100)
        mutdata.needle[, vaf := alt/depth]
        
        vafs.stacked <- rbind(vafs.stacked,
                               data.table(tumour = tum,
                                          sim_id = sim.id,
                                          vaf = mutdata.needle$vaf,
                                          sample = paste0('needle', i)))
      }
    }
  }
  
  fwrite(vafs.stacked, 'data/supp/stocheffect/vafs_stacked.csv')
  
} else {
  
  vafs.stacked <- fread('data/supp/stocheffect/vafs_stacked.csv')
}

vafs.stacked$sim_id <- as.factor(vafs.stacked$sim_id)

tumours <- as.character(unique(vafs.stacked$tumour))

text.size <- 25
for (tum in tumours) {
  
  vafs.curr <- vafs.stacked[tumour == tum]
  
  p <- 
    ggplot(vafs.curr, 
           aes(x = vaf, 
               colour = sim_id)) +
    geom_density(alpha = .5, 
                 position = 'identity') +
    facet_wrap(~sample, scales = 'free') +
    labs(x = 'Allele frequency', y = 'Density') +
    theme_bw() + 
    theme(axis.text = element_text(size = text.size - 5),
          axis.title = element_text(size = text.size),
          #legend.position = 'none',
          strip.text = element_text(size = text.size))
  ggsave(paste0(path.to.output, 'S5_',tum,'.pdf'),
         p, width = 20, height = 15)
}

vafs.stacked.t4 <- vafs.stacked[tumour == 'T4']

vafs.stacked.t4$sim_id <- as.character(vafs.stacked.t4$sim_id)

escaped.tum.ids <- fread('data/supp/escaped_tumour_ids_T4.csv')

escaped.fixed <- escaped.tum.ids[fixed == 1]$id
escaped.not.fixed <- escaped.tum.ids[fixed == 0]$id

vafs.stacked.t4$group <- ifelse(
  vafs.stacked.t4$sim_id %in% escaped.fixed,
  'escaped-fixed', 
  ifelse(vafs.stacked.t4$sim_id %in% escaped.not.fixed,
         'escaped-not-fixed',
         'imprisoned'))

vafs.melt <- reshape2::melt(
  vafs.stacked.t4, id.vars = c('tumour', 'sim_id', 'sample', 'group'))

text.size <- 23
p <- 
  ggplot(vafs.melt,
         aes(x = value, 
             colour = sim_id)) +
  geom_density(alpha = .5, position = 'identity') +
  facet_grid(group~sample) +
  labs(x = 'Allele frequency', y = 'Density') +
  theme_bw() + 
  theme(axis.text = element_text(size = text.size - 5),
        axis.title = element_text(size = text.size),
        legend.position = 'none',
        strip.text = element_text(size = text.size))
ggsave(paste0(path.to.output, 'S5_T4_split.pdf'),
       p, width = 25, height = 12)


# Figure S6 ---------------------------------------------------------------

tumours <- paste0('T', 1:4)
fmin <- .05
fmax <- .4

auc.tbl <- data.table()
for (tum in tumours) {
  
  print(tum)
  
  path.to.data <- paste0('data/supp/stocheffect/',tum,'/')
  
  if (tum %in% c('T3', 'T4')) {
    fmin <- .1
    fmax <- .45
  }
  
  for (sim.id in 1:100) {
    
    mutdata <- fread(paste0(
      path.to.data,'/vaf_wholetumour_',sim.id,'.csv'))
    
    mutdata <- data.table(
      clone = mutdata$clone_id,
      alt = mutdata$alt_count,
      depth = mutdata$depth,
      id = mutdata$mutation_id)
    
    mutdata <- AddClonalMuts(mutdata, 100, 100)
    
    mutdata[, vaf := alt/depth]
    
    test.res <- neutralitytest(mutdata$vaf, fmin = fmin, fmax = fmax)
    
    auc <- test.res$area$metric
    auc.pv <- test.res$area$pval
    auc.tbl <- rbind(auc.tbl, data.table(tumour = tum,
                                         sim_id = sim.id,
                                         AUC = auc,
                                         pvalue = auc.pv))
  }
}

auc.tbl$tumour <- factor(auc.tbl$tumour,
                         levels = tumours,
                         labels = c('Homogeneous neutral',
                                    'Homogeneous selective',
                                    'Boundary-driven neutral',
                                    'Boundary-driven selective'))

auc.tbl.melt <- reshape2::melt(
  auc.tbl, id.vars = c('tumour', 'sim_id'))

text.size <- 15
p <- ggplot(auc.tbl.melt[variable == 'pvalue'],
            aes(x = value, fill = tumour, colour = tumour)) +
  geom_histogram(alpha = .5, bins = 50, position = 'identity') +
  scale_fill_brewer(palette = 'Paired') +
  scale_colour_brewer(palette = 'Paired') +
  geom_vline(xintercept = 0.05, colour = 'firebrick3') +
  facet_wrap(~tumour) +
  theme_bw() +
  labs(x = 'AUC p-value') +
  theme(axis.text = element_text(size = text.size),
        axis.title = element_text(size = text.size),
        legend.position = 'none',
        legend.text = element_text(size = text.size),
        legend.title = element_text(size = text.size, face = 'bold'),
        strip.text = element_text(size = text.size))

ggsave(paste0('output/plots/supp/auc_pvalue_distributions.pdf'), 
       p, width = 10, height = 7)

auc.tbl.T4 <- auc.tbl[tumour == 'Boundary-driven selective']

escaped.tum.ids <- fread('data/supp/escaped_tumour_ids_T4.csv')

escaped.fixed <- escaped.tum.ids[fixed == 1]$id
escaped.not.fixed <- escaped.tum.ids[fixed == 0]$id

auc.tbl.T4$group <- ifelse(auc.tbl.T4$sim_id %in% escaped.fixed,
                           'escaped-fixed', 
                           ifelse(auc.tbl.T4$sim_id %in% escaped.not.fixed,
                                  'escaped-not-fixed',
                                  'imprisoned'))

auc.tbl.melt <- reshape2::melt(
  auc.tbl.T4, id.vars = c('tumour', 'sim_id', 'group'))

text.size <- 15
p <- ggplot(auc.tbl.melt[variable == 'pvalue'],
            aes(x = value, fill = group, colour = group)) +
  geom_histogram(alpha = .5, bins = 50, position = 'stack') +
  scale_fill_manual('Mutant subclone',
    values = c('orchid4', 'orchid3', 'lightslategrey')) +
  scale_colour_manual('Mutant subclone',
    values = c('orchid4', 'orchid3',  'lightslategrey')) +
  geom_vline(xintercept = 0.05, colour = 'firebrick3') +
  theme_bw() +
  labs(x = 'AUC p-value') +
  theme(axis.text = element_text(size = text.size),
        axis.title = element_text(size = text.size),
        legend.text = element_text(size = text.size),
        legend.title = element_text(size = text.size, face = 'bold'),
        strip.text = element_text(size = text.size))

ggsave(paste0('output/plots/supp/auc_pvalue_distributions_T4.pdf'), 
       p, width = 10, height = 5)


## Bulk samples
tumours <- paste0('T', 1:4)

f.min <- .05
f.max <- .4

auc.tbl <- data.table()
for (tum in tumours) {
  
  print(tum)
  
  path.to.data <- paste0('data/supp/stocheffect/',tum,'/')
  
  for (sim.id in 1:100) {
    
    mutdata.wholet <- fread(paste0(
      path.to.data, 'vaf_wholetumour_',sim.id,'.csv'))
    
    mutdata.wholet <- data.table(
      clone = mutdata.wholet$clone_id,
      alt = mutdata.wholet$alt_count,
      depth = mutdata.wholet$depth,
      id = mutdata.wholet$mutation_id)
    
    mutdata.wholet <- AddClonalMuts(mutdata.wholet, 100, 100)
    mutdata.wholet[, vaf := alt/depth]
    
    test.res <- neutralitytest(mutdata.wholet$vaf, fmin = f.min, fmax = f.max)
    
    auc.tbl <- rbind(auc.tbl, data.table(tumour = tum,
                                         sim_id = sim.id,
                                         sample = 'Whole tumour',
                                         AUC = test.res$area$metric,
                                         pvalue = test.res$area$pval))
    
    punch.files <- paste0('punch_sample_p',1:6,'_',sim.id,'.csv')
    
    for (i in 1:length(punch.files)) {
      
      mutdata.punch <- fread(paste0(path.to.data, punch.files[i]))
      mutdata.punch <- AddClonalMuts(mutdata.punch, 100, 100)
      mutdata.punch[, vaf := alt/depth]
      
      test.res <- neutralitytest(mutdata.punch$vaf, fmin = f.min, fmax = f.max)
      
      auc.tbl <- rbind(auc.tbl, data.table(tumour = tum,
                                           sim_id = sim.id,
                                           sample = paste0('Punch', i),
                                           AUC = test.res$area$metric,
                                           pvalue = test.res$area$pval))
    }
    
    needle.files <- paste0('needle_sample_n',7:8,'_',sim.id,'.csv')
    
    for(i in 1:length(needle.files)){
      
      mutdata.needle <- fread(paste0(path.to.data, needle.files[i]))
      mutdata.needle <- AddClonalMuts(mutdata.needle, 100, 100)
      mutdata.needle[, vaf := alt/depth]
      
      test.res <- neutralitytest(mutdata.punch$vaf, fmin = f.min, fmax = f.max)
      
      auc.tbl <- rbind(auc.tbl, data.table(tumour = tum,
                                           sim_id = sim.id,
                                           sample = paste0('Needle', i),
                                           AUC = test.res$area$metric,
                                           pvalue = test.res$area$pval))
    }
  }
}

auc.tbl$tumour <- factor(auc.tbl$tumour,
                         levels = tumours,
                         labels = c('Homogeneous neutral',
                                    'Homogeneous selective',
                                    'Boundary-driven neutral',
                                    'Boundary-driven selective'))
auc.pvalues.only <- auc.tbl
auc.pvalues.only$AUC <- NULL
auc.tbl.melt <- reshape2::melt(auc.pvalues.only,
                               id.vars = c('tumour', 'sample', 'sim_id'))

text.size <- 20
p <- ggplot(auc.tbl.melt,
            aes(x = value, 
                fill = tumour, 
                colour = tumour)) +
  geom_histogram(alpha = .5, 
                 bins = 50, 
                 position = 'identity') +
  scale_fill_brewer(palette = 'Paired') +
  scale_colour_brewer(palette = 'Paired') +
  geom_vline(xintercept = 0.05, 
             colour = 'firebrick3') +
  facet_wrap(tumour~sample, ncol = 9) +
  theme_bw() +
  labs(x = 'AUC p-value') +
  theme(axis.text = element_text(size = text.size),
        axis.title = element_text(size = text.size),
        axis.ticks = element_blank(),
        legend.position = 'none',
        strip.text = element_text(size = text.size))

ggsave('output/plots/supp/auc_pvalues_bulksamples.pdf',
       p, width = 35, height = 20)


# Figure S10 ---------------------------------------------------------------

singlecell.vafs <- fread(paste0(
  path.to.data,'singlecell_vafs_pertumour_400x.csv'))

singlecell.vafs[, vaf2 := vaf/2]

p <- ggplot(singlecell.vafs[vaf2 >= .05 & vaf2 <= .4], 
            aes(vaf2,
                fill = tumour)) +
  geom_histogram(bins = 50) +
  scale_fill_brewer(palette = 'Paired') +
  facet_wrap(~tumour, scales = 'free') +
  theme_bw() + 
  labs(x = 'Allele Frequency', 
       y = 'Number of Mutations',
       fill = 'Tumour') +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face = 'bold'),
        legend.spacing.y = unit(.5, 'cm'),
        strip.text = element_text(size = 20))

ggsave(paste0(path.to.output, '/singlecell_vafs_400x.pdf'),
       p, width = 10, height = 8)


# Figure S11 --------------------------------------------------------------

path.to.data <- 'data/supp/entropy_grids/'
rerun.test <- F

if (rerun.test) {
  
  nb <- cell2nb(nrow = 400, ncol = 400, type = 'rook')
  lwb <- nb2listw(nb, style = 'W')
  
  entropy.test.res <- do.call(
    rbind, parallel::mclapply(1:400, function(i) {
      print(i)
      
      space <- fread(paste0(path.to.data,'space_',i,'.csv'))
      space <- as.matrix(space)
      space[is.na(space)] <- 0
      
      clone.types <- as.vector(t(space))
      
      test.res <-
        moran.test(x = clone.types,
                   listw = lwb,
                   randomisation = T,
                   alternative = 'two.sided',
                   na.action = na.pass)
      
      return(
        data.table(sim_id = i,
                   statistic = test.res$statistic,
                   observed = test.res$estimate[1],
                   expected = test.res$estimate[2],
                   variance = test.res$estimate[3],
                   p_value = test.res$p.value)
      )
      
    }, mc.cores = 1))
  
  entropy.test.res[, tumour := rep(paste0('T', 1:4), each = 100)]
  fwrite(entropy.test.res, 'data/supp/moranstest_lattice.csv')
}


entropy.test.res <- fread('data/supp/moranstest_lattice.csv')

entropy.test.res[, obs_median := round(median(observed), 3), by = tumour]

dt.melt <- reshape2::melt(entropy.test.res,
                          id.vars = c('tumour','sim_id'))

obs.medains <- unique(entropy.test.res[, .(tumour, obs_median)])
setnames(obs.medains, 'obs_median', 'observed')

obs.medains.melt <- reshape2::melt(obs.medains)

text.size <- 15
p <- 
  ggplot(dt.melt[variable %in% c('observed')],
         aes(y = value, 
             x = variable,
             fill = tumour,
             colour = tumour)) +
  geom_violin(scale = 'width',
              width = 1,
              position = position_dodge(width = 1),
              alpha = .4) +
  geom_point(aes(fill = tumour),
             position = position_jitterdodge(
               jitter.width = .2, dodge.width = 1),
             size = 3,
             alpha = .8) + 
  geom_text(data = obs.medains.melt,
            aes(x = variable,
                y = .6,
                label = value,
                group = tumour), 
            size = 6,
            fontface = 'bold',
            position = position_dodge(width = .9),
            hjust = .5) + 
  scale_fill_brewer(palette = 'Paired') +
  scale_colour_brewer(palette = 'Paired') +
  theme_bw() +
  labs(x = "Moran's test observed value") +
  theme(axis.text = element_text(size = text.size),
        axis.title = element_text(size = text.size),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = text.size),
        legend.title = element_text(size = text.size, face = 'bold'),
        strip.text = element_text(size = text.size))

ggsave('output/plots/supp/moranstest_lattice.pdf', 
       p, width = 10, height = 6)



# Figure S12 ---------------------------------------------------------------

tumours <- paste0('T', 1:4)
nclonal.muts <- depth <- 100

sfs.stats <- data.table()

for(tum in tumours){
  
  path.to.data <- paste0('data/supp/stocheffect/',tum,'/')
  
  sfs.stats.tumour <- do.call(rbind, parallel::mclapply(1:100, function(i) {
    
    sim.id <- i
    
    fHsub.rAUC.wholet <- Calculate_fHsub_rAUC(
      paste0(path.to.data,'vaf_wholetumour_',sim.id,'.csv'),
      nclonal.muts,
      depth)
    
    fHsub.rAUC.punch <- Calculate_fHsub_rAUC(
      paste0(path.to.data,'punch_sample_p',1:6,'_',sim.id,'.csv'),
      nclonal.muts,
      depth)
    
    fHsub.rAUC.needle <- Calculate_fHsub_rAUC(
      paste0(path.to.data,'needle_sample_n',7:8,'_',sim.id,'.csv'),
      nclonal.muts,
      depth)
    
    fHrs.ksd.punch <- Calculate_fHrs_ksd(
      paste0(path.to.data,'punch_sample_p',1:6,'_',sim.id,'.csv'),
      nclonal.muts,
      depth)
    
    fHrs.ksd.needle <- Calculate_fHrs_ksd(
      paste0(path.to.data,'needle_sample_n',7:8,'_',sim.id,'.csv'),
      nclonal.muts,
      depth)
    
    sfs.stats.wholet <- CaclulateSFSstats(
      paste0(path.to.data,'vaf_wholetumour_',sim.id,'.csv'),
      nclonal.muts,
      depth)
    
    sfs.stats.punch <- CaclulateSFSstats(
      paste0(path.to.data,'punch_sample_p',1:6,'_',sim.id,'.csv'),
      nclonal.muts,
      depth)
    
    sfs.stats.needle <- CaclulateSFSstats(
      paste0(path.to.data,'needle_sample_n',7:8,'_',sim.id,'.csv'),
      nclonal.muts,
      depth)
    
    result <- data.frame(rbind(fHsub.rAUC.wholet,
                               fHsub.rAUC.punch,
                               fHsub.rAUC.needle),
                         rbind(data.frame(fHrs = NA, ksd = NA),
                               fHrs.ksd.punch,
                               fHrs.ksd.needle), 
                         rbind(sfs.stats.wholet,
                               sfs.stats.punch,
                               sfs.stats.needle))
    
    result$sampling <- c('whole tumour',
                         'punch biopsies',
                         'needle biopsies')
    
    return(result)
    
  }, mc.cores=2))
  
  sfs.stats.tumour$tumour <- tum
  sfs.stats <- rbind(sfs.stats, sfs.stats.tumour)
}
  
sfs.stats.melt <- reshape2::melt(sfs.stats, id.vars = c('tumour',
                                                        'sampling'))
sfs.stats.melt$sampling <- factor(sfs.stats.melt$sampling,
                                  levels = c('punch biopsies',
                                             'needle biopsies',
                                             'whole tumour'))
pdf(paste0(path.to.output,'/sfs_stats.pdf'),
    paper = 'USr', width = 50, height = 20)
p <- ggplot(sfs.stats.melt,
            aes(x = value,
                fill = sampling)) +
  geom_density(alpha = .6,
               position = 'identity') +
  theme_bw() + 
  scale_fill_manual(values = c('springgreen4',
                               'darkblue',
                               'deeppink4')) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20,
                                  face = 'bold'),
        strip.text = element_text(size = 15))

facet_multiple(plot = p,
               facets = c('variable',
                          'tumour'), 
               scales = 'free',
               ncol = 4,
               nrow = 2)
dev.off()

tree.stats <- data.frame()
for (tum in tumours) {
  
  print(tum)
  path.to.data <- paste0('data/supp/stocheffect/',tum,'/')
  
  tree.stats.tumour <- do.call(rbind, parallel::mclapply(1:100, function(i) {
    
    sim.id <- i
    
    return(CalculateTreeStats(paste0(path.to.data,'sim_100sc_',sim.id,'.tre')))
    
  }, mc.cores = 1))
  
  tree.stats.tumour$tumour <- tum
  tree.stats <- rbind(tree.stats, tree.stats.tumour)
}

tree.stats.melt <- reshape2::melt(tree.stats, id.vars = 'tumour')

p <- ggplot(tree.stats.melt,
            aes(x = value,
                fill = tumour)) +
  geom_density(alpha = 1, position = 'stack') +
  facet_wrap(~variable, scales = 'free') +
  theme_bw() + 
  scale_fill_brewer(palette = 'Paired')+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20,
                                    face = 'bold'),
        strip.text = element_text(size = 18))
ggsave(paste0(path.to.output,'/tree_stats.pdf'),
       p, width = 15, height = 10)

tum.names <- c('neutral',
               'selective',
               'neutral peripheral',
               'selective peripheral')

metadata <- data.frame()
for(i in 1:length(tumours)) {
  
  print(tumours[i])
  path.to.data <- paste0('data/supp/stocheffect/',tumours[i],'/')
  metadata.tum <- fread(paste0(path.to.data, 'simmetadata.csv'))
   
  pop.size = metadata.tum$wt_cellcount + metadata.tum$mt_cellcount
  
  metadata <- rbind(metadata,
                    data.frame(total_mut_count = metadata.tum$mutcount,
                               population_size = pop.size,
                               f_sub = metadata.tum$mt_cellcount / pop.size,
                               t_end = metadata.tum$tend,
                               tumour = tum.names[i]))
}

metadata.melt <- reshape2::melt(metadata, id.vars = c('tumour'))

p <- ggplot(metadata.melt,
            aes(x = value,
                fill = tumour)) +
  geom_density(alpha = 1, position = 'stack') +
  facet_wrap(~variable, scales = 'free') +
  theme_bw() + 
  scale_fill_brewer(palette = 'Paired') +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.text=element_text(size = 20),
        legend.title=element_text(size = 20, 
                                  face = 'bold'),
        strip.text = element_text(size = 20))

ggsave(paste0(path.to.output,'/tree_stats.pdf'),
       p, width = 12, height = 8)


# Figure S13 ---------------------------------------------------------------

source('scripts/figure6.R')

GeneratePlotsFig6(path.to.data = 'data/figure_6/',
                  path.to.output = 'output/plots/supp/',
                  plot.type = 'posterior')


# Figure S14 ---------------------------------------------------------------

t.range <- 1:20
s.range <- 1:10

file.names <- dir(path = path.to.data, pattern = 'ts_vaf_*')

for (file in file.names) {
  
  ts.vaf <- fread(paste0(path.to.data, file))
  colnames(ts.vaf) <- as.character(s.range)
  ts.vaf[, t := t.range]
  ts.vaf.melt <- reshape2::melt(ts.vaf, id.vars = 't')
  setnames(ts.vaf.melt, 'variable', 's')
  
  p <- ggplot(ts.vaf.melt, aes(t, s, fill = value)) + 
    geom_tile() + 
    theme_bw() +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 20, face = 'bold'))
  
  ggsave(paste0(path.to.output, file, '.pdf'), 
         p, width = 7, height = 5)
}


# Figure S15 --------------------------------------------------------------

path.to.data <- 'data/supp/inference_3D/'
sampling.dirs <- c('needle', 'punch', 'tree', 'wholet')
max.round <- 15

GetParticlesFile <- function(max.round, sampling, rec.param) {
  
  r <- max.round
  repeat {
    particles.file <- paste0(
      path.to.data,sampling,'/rec_',rec.param,'/T9/smc_particles_r',r,'.csv')
    r <- r - 1
    if (file.exists(particles.file)) break
  }
  return(particles.file)
}

particles <- data.table()

for (sampling in sampling.dirs) {
  
  particles.samp <- data.table()
  
  particles.rec <- fread(
    GetParticlesFile(max.round, sampling, 'mu'))
  particles.samp <- rbind(
    particles.samp, 
    data.table(parameter = 'Mutation rate',
               value = particles.rec$mu))
  
  particles.rec <- fread(
    GetParticlesFile(max.round, sampling, 't_s'))
  particles.samp <- rbind(
    particles.samp, 
    data.table(parameter = 'Time to new clone',
               value = particles.rec$t))
  
  particles.rec <- fread(
    GetParticlesFile(max.round, sampling, 's'))
  particles.samp <- rbind(
    particles.samp, 
    data.table(parameter = 'Selective advantage',
               value = particles.rec$s-1))
  
  particles.rec <- fread(
    GetParticlesFile(max.round, sampling, 'd_a'))
  particles.samp <- rbind(
    particles.samp, 
    data.table(parameter = 'Death rate',
               value = particles.rec$d))
  
  particles.samp <- rbind(
    particles.samp, 
    data.table(parameter = 'Aggression',
               value = particles.rec$a))
  
  particles.samp[, sampling := sampling]
  
  particles <- rbind(particles, particles.samp)
  
}

particles.mean.mode <- 
  particles %>% 
  group_by(sampling, parameter) %>%
  summarise(Mode = mean(GetMode(round(value, 1))))

particles.mean.mode$Target <- rep(c(1, 0, 10, 3, 7), 4)

particles.meta.melt <- reshape2::melt(particles.mean.mode)

param.levels <- c('Mutation rate',
                  'Time to new clone',
                  'Selective advantage',
                  'Death rate',
                  'Aggression')

sampling.levels <- c('wholet', 'needle', 'punch', 'tree')

particles$parameter <- factor(
  particles$parameter, levels = param.levels)

particles.meta.melt$parameter <- factor(
  particles.meta.melt$parameter, levels = param.levels)

particles$sampling <- factor(
  particles$sampling, levels = sampling.levels)

particles.meta.melt$sampling <- factor(
  particles.meta.melt$sampling, levels = sampling.levels)

particles.meta.melt <- as.data.table(particles.meta.melt)


setorder(particles, 'parameter')

particles[, true_value := rep(c(10, 7, 3, 0, 1),
                              times = c(table(particles$parameter)))]
particles[, position := rep(c(20, 25, 5, 1.3, 1.3),
                              times = c(table(particles$parameter)))]

text.size <- 35
p <- ggplot(particles, 
            aes(x = sampling,
                y = value, 
                fill = parameter, 
                colour = parameter)) +
  geom_flat_violin(scale = 'width',
                   width = .8,
                   position = position_dodge(width = 1),
                   alpha = .4,
                   trim = F) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'up',
               alpha = .1,
               stackratio = .15,
               dotsize = .6) +
  geom_hline(aes(yintercept = true_value,
                 colour = parameter,
                 group = sampling),
             data = particles,
             linetype = 'dashed') +
  geom_text(data = particles.meta.melt[variable == 'Mode'],
            aes(x = sampling,
                y = value,
                label = value,
                group = parameter), 
            size = 12,
            fontface = 'bold',
            position = position_dodge(width = .9),
            hjust = .5) +
  geom_text(data = particles,
            aes(x = sampling,
                y = position,
                label = true_value,
                group = parameter), 
            size = 12,
            fontface = 'bold',
            position = position_dodge(width = .9),
            hjust = .5) +
  scale_fill_jama() +
  scale_colour_jama() +
  facet_wrap(~parameter, 
             scales = 'free',
             ncol = 2) +
  labs(x = '', y = 'Posterior density') +
  theme_readable() +
  theme(axis.text = element_text(size = text.size),
        axis.title = element_text(size = text.size + 2),
        legend.position = 'none',
        legend.text = element_text(size = text.size),
        legend.title = element_text(size = text.size + 2,
                                    face = 'bold'),
        legend.spacing.y = unit(.5, 'cm'),
        strip.text = element_text(size = text.size,
                                  margin = margin(1,0,1,0, 'cm'))) 

ggsave('output/plots/supp/inference_3D_temp.pdf',
       p, width = 30, height = 25)



# Figure S16 --------------------------------------------------------------

path.to.data <- 'data/supp/inference_roerink2018/'
max.round <- 15
patients <- 1:3

particles <- data.table()

for (patient in patients) {
  
  patient.name <- paste0('Patient ',patient)
  
  r <- max.round
  repeat {
    particles.file <- paste0(
      path.to.data,'/rec_mu/T',patient,'/smc_particles_r',r,'.csv')
    r <- r - 1
    if (file.exists(particles.file)) break
  }
  particles.mu <- fread(particles.file)
  particles.patient <- data.table(patient = patient.name,
                                  parameter = 'Mutation rate',
                                  value = particles.mu$mu)
  
  r <- max.round
  repeat {
    particles.file <- paste0(
      path.to.data,'/rec_t_s_d_a/T',patient,'/smc_particles_r',r,'.csv')
    r <- r - 1
    if (file.exists(particles.file)) break
  }
  
  particles.rest <- fread(particles.file)
  particles.rest <- particles.rest[,.(t,s,d,a)]
  particles.rest$s <- particles.rest$s - 1
  setnames(particles.rest,
           c('t', 's', 'd', 'a'),
           c('Time to new clone',
             'Selective advantage',
             'Death rate',
             'Aggression'))
  
  particles.rest.melt <- reshape2::melt(particles.rest)
  setnames(particles.rest.melt, 'variable', 'parameter')

  particles.patient <- rbind(particles.patient,
                             data.table(patient = patient.name,
                                        particles.rest.melt))
  
  particles <- rbind(particles, particles.patient)
}

particles.meta.set1 <- 
particles[parameter %in% c('Mutation rate',
                           'Time to new clone',
                           'Selective advantage')] %>% 
  group_by(patient, parameter) %>%
  summarise(Mean = mean(value, na.rm = T),
            Mode = mean(GetMode(round(value))))

particles.meta.set2 <- 
  particles[parameter %in% c('Death rate',
                             'Aggression')] %>% 
  group_by(patient, parameter) %>%
  summarise(Mean = mean(value, na.rm = T),
            Mode = mean(GetMode(round(value, 1))))

particles.meta.melt <- as.data.table(reshape2::melt(
  rbind(particles.meta.set1, particles.meta.set2)))


text.size <- 12
p <- ggplot(particles, 
            aes(x = patient,
                y = value, 
                fill = parameter, 
                colour = parameter)) +
  # geom_violin(scale = 'width',
  #             width = .8,
  #             position = position_dodge(width = 1),
  #             alpha = .4, 
  #             trim = F) +
  geom_flat_violin(scale = 'width',
                   width = .8,
                   position = position_dodge(width = 1),
                   alpha = .4,
                   trim = F) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'up',
               alpha = .2,
               stackratio = .5,
               dotsize = .6) +
  geom_text(data = particles.meta.melt[variable == 'Mode'],
            aes(x = patient,
                y = value,
                label = value,
                group = parameter), 
            size = 5,
            fontface = 'bold',
            position = position_dodge(width = .9),
            hjust = .5) +
  scale_fill_jama() +
  scale_colour_jama() +
  facet_wrap(~parameter, 
             scales = 'free',
             ncol = 2) +
  labs(x = '', y = 'Posterior density') +
  theme_readable() +
  theme(axis.text = element_text(size = text.size),
        axis.title = element_text(size = text.size + 3),
        legend.position = 'none',
        legend.text = element_text(size = text.size),
        legend.title = element_text(size = text.size + 3,
                                    face = 'bold'),
        legend.spacing.y = unit(.5, 'cm'),
        strip.text = element_text(size = text.size)) 

ggsave('output/plots/supp/inference_roerink2018_new.pdf',
       p, width = 14, height = 10)


# Figure S17 --------------------------------------------------------------

path.to.data <- 'davros/inferred_params/BL/tree/'
max.round <- 20
patients <- 1:3

particles.all <- data.table()
for (patient in patients) {
  
  patient.name <- paste0('Patient ',patient)
  
  r <- max.round
  repeat {
    particles.file <- paste0(
      path.to.data,'/rec_mu_t_s_d_a/T',patient,'/smc_particles_r',r,'.csv')
    r <- r - 1
    if (file.exists(particles.file)) break
  }
  particles <- fread(particles.file)
  particles.patient <- data.table(patient = patient.name,
                                  parameter = 'Mutation rate',
                                  value = particles$mu)
  
  r <- max.round
  repeat {
    particles.file <- paste0(
      path.to.data,'/rec_t_s_d_a/T',patient,'/smc_particles_r',r,'.csv')
    r <- r - 1
    if (file.exists(particles.file)) break
  }
  particles <- fread(particles.file)
  particles.patient <- rbind(particles.patient,
                             data.table(patient = patient.name,
                                        parameter = 'Time to new clone',
                                        value = particles$t))
  r <- max.round
  repeat {
    particles.file <- paste0(
      path.to.data,'/rec_s_d_a/T',patient,'/smc_particles_r',r,'.csv')
    r <- r - 1
    if (file.exists(particles.file)) break
  }
  particles <- fread(particles.file)
  particles <- particles[,.(s,d,a)]
  particles$s <- particles$s - 1
  setnames(particles,
           c('s', 'd', 'a'),
           c('Selective advantage',
             'Death rate',
             'Aggression'))
  
  particles.melt <- reshape2::melt(particles)
  setnames(particles.melt, 'variable', 'parameter')
  particles.patient <- rbind(particles.patient,
                             data.table(patient = patient.name,
                                        particles.melt))

  
  # r <- max.round
  # repeat {
  #   particles.file <- paste0(
  #     path.to.data,'/rec_t_s_d_a/T',patient,'/smc_particles_r',r,'.csv')
  #   r <- r - 1
  #   if (file.exists(particles.file)) break
  # }
  # 
  # particles <- fread(particles.file)
  # particles <- particles[,.(t,s,d,a)]
  # particles$s <- particles$s - 1
  # setnames(particles,
  #          c('t', 's', 'd', 'a'),
  #          c('Time to new clone',
  #            'Selective advantage',
  #            'Death rate',
  #            'Aggression'))
  # 
  # particles.melt <- reshape2::melt(particles)
  # setnames(particles.melt, 'variable', 'parameter')
  # 
  # particles.patient <- rbind(particles.patient,
  #                            data.table(patient = patient.name,
  #                                       particles.melt))
  
  particles.all <- rbind(particles.all, particles.patient)
}

particles <- particles.all

particles$parameter <- factor(particles$parameter, 
                              levels = c('Mutation rate',
                                         'Time to new clone',
                                         'Selective advantage',
                                         'Death rate',
                                         'Aggression'))

particles.meta.set1 <- 
  particles[parameter %in% c('Mutation rate',
                             'Time to new clone')] %>% 
  group_by(patient, parameter) %>%
  summarise(Mode = round(mean(GetMode(round(value)))))

particles.meta.set2 <- 
  particles[parameter %in% c('Selective advantage',
                             'Death rate',
                             'Aggression')] %>% 
  group_by(patient, parameter) %>%
  summarise(Mode = round(mean(GetMode(round(value, 2))),1))

particles.mode.melt <- as.data.table(reshape2::melt(
  rbind(particles.meta.set1, particles.meta.set2)))


text.size <- 14
p <- ggplot(particles, 
            aes(x = patient,
                y = value, 
                fill = parameter, 
                colour = parameter)) +
  geom_flat_violin(scale = 'width',
              width = .8,
              position = position_dodge(width = 1),
              alpha = .4,
              trim = T) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'up',
               alpha = .2,
               stackratio = .5,
               dotsize = .6) +
  geom_text(data = particles.mode.melt,
            aes(x = patient,
                y = value,
                label = value,
                group = parameter), 
            size = 5,
            fontface = 'bold',
            position = position_dodge(width = .9),
            hjust = .5) +
  scale_fill_jama() +
  scale_colour_jama() +
  facet_wrap(~parameter, 
             scales = 'free',
             ncol = 2) +
  labs(x = '', y = 'Posterior density') +
  theme_readable() +
  theme(axis.text = element_text(size = text.size),
        axis.title = element_text(size = text.size + 3),
        legend.position = 'none',
        legend.text = element_text(size = text.size),
        legend.title = element_text(size = text.size + 3,
                                    face = 'bold'),
        legend.spacing.y = unit(.5, 'cm'),
        strip.text = element_text(size = text.size,
                                  margin = margin(.4,0,.4,0,'cm'))) 

ggsave('output/plots/supp/inference_sc_organoids_grid500.pdf',
       p, width = 20, height = 14)

