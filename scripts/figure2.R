GeneratePlotsFig2 <- function(sim.id,
                              needle.radius,
                              punch.radius, 
                              depth, 
                              fmin,
                              fmax,
                              fmin.hist,
                              nclonal.muts,
                              path.to.data,
                              path.to.output) {
  
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
  
  mt.cells.ratio <- round(
   100*sum(space.clone == 2)/(sum(space.clone > 0)), 2)
  
  cell.type.counts.needle <- data.table()
  xc <- c(0.375, 0.625) * x 
  for (i in 1:length(xc)) {
    
    space.grid[(xc[i] - needle.radius):(xc[i] + needle.radius), 1:y] <- 3
    
    sample.cell.types <- 
      space.clone[(xc[i] - needle.radius):(xc[i] + needle.radius), 1:y]
    
    cell.type.counts.needle <- rbindlist(list(
      cell.type.counts.needle, 
      data.table(wt = sum(sample.cell.types == 1),
                 mutant = sum(sample.cell.types == 2))
      ))
  }
   
  cell.type.counts.punch <- data.table()
  xc <- c(x/4, x/2, 3*x/4)
  yc <- c(x/4, 3*x/4)
  for (i in 1:length(xc)) {
    for (j in 1:length(yc)) {
      
      space.grid[(xc[i] - punch.radius):(xc[i] + punch.radius),
                (yc[j] - punch.radius):(yc[j] + punch.radius)] <- 3
      
      sample.cell.types <- 
        space.clone[(xc[i] - punch.radius):(xc[i] + punch.radius),
                    (yc[j] - punch.radius):(yc[j] + punch.radius)]
      
      cell.type.counts.punch <- rbindlist(list(
        cell.type.counts.punch,
        data.table(wt = sum(sample.cell.types == 1),
                   mutant = sum(sample.cell.types == 2))
      ))
    }
  }
  
  fwrite(cell.type.counts.punch,
         paste0(path.to.data,'punch_cellcounts.csv'))
  fwrite(cell.type.counts.needle,
         paste0(path.to.data,'needle_cellcounts.csv'))
  
  par(pty = 's')
  space.wt.rev <- apply(space.wt, 2, rev)
  space.mt.rev <- apply(space.mt, 2, rev)
  image(1:nrow(space.wt),
        1:ncol(space.wt),
        as.matrix(t(space.wt.rev)),
        col  = ColorRampAlpha(c('blue4', 'blue'), n = 10, alpha = 1),
        xlab = 'x[cells]',
        ylab = 'y[cells]',
        cex.axis = 1.5, cex.lab = 1.8, cex.main = 1.5)
  image(1:nrow(space.mt),
        1:ncol(space.mt),
        as.matrix(t(space.mt.rev)),
        col  = ColorRampAlpha(c('red4', 'red'), n = 10, alpha = 1),
        add = T)
  image(1:nrow(space.grid),
        1:ncol(space.grid),
        as.matrix(space.grid),
        col = scales::alpha('black', 0.4),
        add = T)
  lab <- as.character(1:6)
  xt <- rep(xc, each = 2)
  yt <- rep(yc, 2)
  text(rep(xc, each = 2),
       rep(yc, 2),
       labels = as.character(1:6),
       col = scales::alpha('white', 0.5),
       cex = 3)
  text(c(0.375, 0.625)*x,
       rep(y/2, 2),
       labels = as.character(7:8),
       col=scales::alpha('white', 0.5),
       cex = 3)

  whole.tumour.mutdata <- fread(paste0(
    path.to.data,'/vaf_wholetumour_',sim.id,'.csv'))
  
  whole.tumour.vaf <- GetVAF(
    mutdata = whole.tumour.mutdata,
    nclonal.muts = nclonal.muts, 
    depth = depth, 
    clone.sizes = clone.sizes)
  
  whole.tumour.vaf.hist <- TestNeutrality(
    mutdata = whole.tumour.vaf,
    fmin = fmin,
    fmax = fmax, 
    fmin.hist = fmin.hist, 
    subtitle = paste0('  Tumour ',sim.id))
  
  punch.ind <- 1:6
  needle.ind <- 7:8
  if (sim.id == 4) {
    punch.ind <- c(1,3,4,5)
  }
  
  punch.files <- paste0(
    path.to.data,'punch_sample_p',punch.ind,'_',sim.id,'.csv')
  needle.files <- paste0(
    path.to.data,'needle_sample_n',needle.ind,'_',sim.id,'.csv')
  
  punch.vaf.hists <- list()
  for (i in 1:length(punch.ind)) {
    print(punch.files[i])
    
    punch.mutdata <- fread(punch.files[i])
    
    punch.mutdata.vaf <- GetVAF(
      mutdata = punch.mutdata,
      nclonal.muts = nclonal.muts, 
      depth = depth, 
      clone.sizes = cell.type.counts.punch[i])
    
    punch.vaf.hist <- TestNeutrality(
      mutdata = punch.mutdata.vaf,
      fmin = fmin,
      fmax = fmax, 
      fmin.hist = fmin.hist, 
      subtitle = paste0('  Punch ', punch.ind[[i]]))
    
    punch.vaf.hists[[length(punch.vaf.hists) + 1]] <- punch.vaf.hist
  }
  
  needle.vaf.hists <- list()
  for (i in 1:length(needle.ind)) {
    print(needle.files[i])
    
    needle.mutdata <- fread(needle.files[i])
    
    needle.mutdata.vaf <- GetVAF(
      mutdata = needle.mutdata,
      nclonal.muts = nclonal.muts, 
      depth = depth, 
      clone.sizes = cell.type.counts.needle[i])
    
    needle.vaf.hist <- TestNeutrality(
      mutdata = needle.mutdata.vaf,
      fmin = fmin,
      fmax = fmax, 
      fmin.hist = fmin.hist, 
      subtitle = paste0('  Needle ', needle.ind[[i]]))
    
    needle.vaf.hists[[length(needle.vaf.hists) + 1]] <- needle.vaf.hist
  }
  
  whole.tumour.vaf.hist <- lapply(list(whole.tumour.vaf.hist), ggplotGrob)
  ggsave(paste0(path.to.output,'whole_tumour_vaf_',sim.id,'.pdf'),
         marrangeGrob(whole.tumour.vaf.hist, nrow = 3, ncol = 2),
         width = 25.4, height = 20.2, units = 'cm')
  
  punch.vaf.hists <- lapply(punch.vaf.hists, ggplotGrob)
  ggsave(paste0(path.to.output,'punch_vafs_',sim.id,'.pdf'),
         marrangeGrob(punch.vaf.hists, nrow = 3, ncol = 2),
         width = 25.4, height = 20.2, units = 'cm')
  
  needle.vaf.hists <- lapply(needle.vaf.hists, ggplotGrob)
  ggsave(paste0(path.to.output,'needle_vafs_',sim.id,'.pdf'),
         marrangeGrob(needle.vaf.hists, nrow = 3, ncol = 2),
         width = 25.4, height = 20.2, units = 'cm')
}