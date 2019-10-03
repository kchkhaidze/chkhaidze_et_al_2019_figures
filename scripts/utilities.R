ColorRampAlpha <- function(..., n, alpha) {
  
  colors <- colorRampPalette(...)(n)
  paste(colors, sprintf('%x', ceiling(255*alpha)), sep = '')
}


GetCloneColours <- function() {
  clone.cols <- c('WT clonal' = 'blue',
                  'WT subclonal' = 'steelblue1',
                  'Mutant clonal' = 'red', 
                  'Mutant subclonal' = 'plum',
                  'Clonal' = 'black', 
                  'HH' = 'purple4')
}


SetClonelabels <- function(mutdata, clone.sizes) {
  
 mutdata[, clone := ifelse(clone %in% 1:2, clone, 'HH')]
 
 mutdata[, clone := ifelse(clone == 1, 'WT subclonal', clone)]
 mutdata[, clone := ifelse(clone == 'WT subclonal' &
                             alt == clone.sizes$wt,
                           'WT clonal',
                           clone)]
 mutdata[, clone := ifelse(clone == 2, 'Mutant subclonal', clone)]
 mutdata[, clone := ifelse(clone == 'Mutant subclonal' &
                             alt == clone.sizes$mutant,
                           'Mutant clonal',
                           clone)]
 return(mutdata)
 
} 


AddClonalMuts <- function(mutdata, nclonal.muts, depth, depth.model='Poisson') {
  
  if (depth.model == 'Binomial') {
    
    mutdata <- rbindlist(list(
      mutdata, data.table(clone = 'Clonal',
                          alt = rbinom(nclonal.muts, depth, .5),
                          depth = depth,
                          id = paste0('C', 1:nclonal.muts))))
    
  } else if (depth.model == 'Poisson') {
    
    depth.pois <- rpois(nclonal.muts, depth)
    mutdata <- rbindlist(list(
      mutdata, data.table(clone = 'Clonal',
                          alt = rbinom(nclonal.muts, depth.pois, .5),
                          depth = depth.pois,
                          id = paste0('C', 1:nclonal.muts))))
  } else {
    
    print('Unrecognised depth model')
  }
  return(mutdata)
}
 

GetVAF <- function(mutdata, nclonal.muts, depth, clone.sizes) {
  
  mutdata <- SetClonelabels(mutdata, clone.sizes)
  
  mutdata <- AddClonalMuts(mutdata, nclonal.muts, depth)
  
  mutdata[, vaf := alt/depth]
  
  return(mutdata)
}


GetVAFhistogram <- function(mutdata, fmin.hist, subtitle) {
  
  p <- ggplot(mutdata[vaf >= fmin.hist],
              aes(x = vaf, 
                  fill = clone)) +
    geom_histogram(binwidth = 0.01,
                   alpha = .5) +
    scale_colour_manual(values = GetCloneColours()) +
    scale_fill_manual(values = GetCloneColours()) +
    xlab( 'Allele Frequency' ) +
    ylab('Number of Mutations') +
    xlim(-0.01, .66) + 
    theme_bw() +
    ggtitle(subtitle) +
    theme(plot.title = element_text(size = 15,
                                    face = 'bold'),
          legend.position = 'none',
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 15))
  return(p)
}


TestNeutrality <- function(mutdata, fmin, fmax, fmin.hist, subtitle='') {
  
  test.res <- neutralitytest(mutdata$vaf, fmin, fmax)
  
  hist.subtitle <- paste0('AUC = ', round(test.res$area$metric, 3),
                          '  pv = ', round(test.res$area$pval, 3),
                          '  u = ', round(test.res$mutation.rate, 2))
  
  hist.vaf <- GetVAFhistogram(mutdata,
                              fmin.hist,
                              paste0(hist.subtitle, subtitle))
  return(hist.vaf)
}


SingleCellSampleCoords <- function(x, y, z, nsamples, space) {
  
  radius <- x/2
  sample.strategy <- data.frame() 
  
  while (nrow(sample.strategy) < nsamples) {
    
    x.sample <- sample(1:x, 1)
    y.sample <- sample(1:y, 1)
    
    if (((x.sample - (x - radius))^2 +
        (y.sample - (y - radius))^2 < radius^2) & 
        space[x.sample, y.sample] > 0) {
      sample.strategy <- rbind(sample.strategy, 
                               data.frame(x.sample, y.sample))
    } 
  }
  sample.strategy$z.sample <- z
  
  return(sample.strategy)
}


SingleCellSamplesFromBulk <- function(
  x1, x2, y1, y2, x, y, z, nsamples, space) {
  
  radius <- x/2
  sample.strategy <- data.frame()
  
  while (nrow(sample.strategy) < nsamples) {
    
    sample.x <- sample(x1:x2, 1)
    sample.y <- sample(y1:y2, 1)
    
    if (((sample.x - (x - radius))^2 + (sample.y - (x - radius))^2 < radius^2) & 
        !is.na(space[sample.x, ..sample.y]) &
        space[sample.x, ..sample.y] > 0) {
      
      sample.strategy <- rbind(sample.strategy,
                               data.frame(x = sample.x, y = sample.y))
    } 
  }
  sample.strategy$z <- z
  return(sample.strategy)
}


ConcentricBulkCoords <- function(radius, grid) {
  
  sample.strategy <- list()
  x <- y <- grid
  x.coords <- y.coords <- c()
  gap <- 50
  yc <- grid/2
  xc <- seq(gap, x - gap, gap)
  
  for (i in 1:length(xc)) {
    sample.strategy[[length(sample.strategy) + 1]] <-
      c(xc[i] - radius + 1, yc - radius + 1, 1, xc[i] + radius, yc + radius, 1)
  }
  x.coords <- c(x.coords, xc)
  y.coords <- c(y.coords, rep(yc, length(xc)))
  yc <- xc[xc != x/2]
  xc <- x/2
  for (i in 1:length(yc)) {
    sample.strategy[[length(sample.strategy) + 1]] <-
      c(xc - radius + 1, yc[i] - radius + 1, 1, xc + radius, yc[i] + radius, 1)
  }
  x.coords <- c(x.coords, rep(xc, length(yc)))
  y.coords <- c(y.coords, yc)
  return(list(sample.strategy = sample.strategy,
              x.coords = x.coords,
              y.coords = y.coords))
}


NeedleBiopsyCoords <- function(radius, grid) {
  
  sample.strategy <- list()
  x <- y <- grid
  xc <- yc <- c(0.375, 0.625) * x
  for (i in 1:length(xc)) {
    sample.strategy[[i]] <- c(xc[i] - radius, 1, 1, xc[i] + radius, y, 1)
  } 
  return(sample.strategy)
}


PunchBiopsyCoords <- function(radius, grid) {
  
  sample.strategy <- list()
  x <- y <- grid
  xc <- c(x/4, x/2, 3*x/4)
  yc <- c(.3, .7)*grid
  for (i in 1:length(xc)) {
    for (j in 1:length(yc)) {
      sample.strategy[[length(sample.strategy) + 1]] <- 
        c(xc[i] - radius, yc[j] - radius, 1, xc[i] + radius, yc[j] + radius, 1)
    } 
  }
  return(sample.strategy)
}


GetSamplingStrategy <- function(grid, punch.coords, needle.coords) {
  
  sample.strategy <- list(c(1, 1, 1, grid, grid, 1))
  for (b in 1:length(punch.coords)) {
    sample.strategy[[length(sample.strategy) + 1]] <- punch.coords[[b]]
  } 
  for (n in 1:length(needle.coords)) {
    sample.strategy[[length(sample.strategy) + 1]] <- needle.coords[[n]]
  }
  return(sample.strategy)
}

GetTend <- function(fsub, Nend) {
  log2((1 - fsub)*Nend)
}


GetSfsub <- function(
  tstart, tend, fsub, death.wt = 0, birth.wt = 1) {
  lambda <- birth.wt - death.wt
  (lambda*tstart + log(fsub/(1 - fsub)))/(lambda*(tend - tstart))
}


GetMode <- function(x) { 
  x.tbl <- table(x)
  x.max <- max(x.tbl)
  if (all(x.tbl == x.max)){
    res <- NA
  } else {
    if (is.numeric(x)) {
      res <- as.numeric(names(x.tbl)[x.tbl == x.max])
    } else {
      res <- names(x.tbl)[x.tbl == x.max]
    }
  }
  return(res)
}


MedianNArm <- function(x) median(x, na.rm=T)


Calculate_fHsub_rAUC <- function(sample.files,
                                 nclonal.muts,
                                 depth,
                                 steps=seq(0, 1, .01)) {

  fHsub <- rAUC <- rep(0, length(sample.files))
  
  for (i in 1:length(sample.files)) {
    
    mutdata <- fread(sample.files[i])
    mutdata <- AddClonalMuts(mutdata, nclonal.muts, depth)
    mutdata[, vaf := alt/depth]
    
    auc.merged <- flux::auc(
      steps,sapply(steps, function(x) sum(mutdata$vaf <= x)))
    
    auc.theoretical <- flux::auc(
      steps, sapply(steps, function(x) sum((1:length(mutdata$vaf))^-1 <= x)))
    
    rAUC[i] <- auc.merged/auc.theoretical
    fHsub[i] <- sum(mutdata$vaf > .2)/sum(mutdata$vaf > .08)
  }
  
  fHsub <- sum(fHsub)/length(fHsub)
  rAUC <- sum(rAUC)/length(rAUC)
  
  return(data.frame(fHsub = fHsub, rAUC = rAUC))
}


Calculate_fHrs_ksd <- function(sample.files,
                               nclonal.muts,
                               depth,
                               steps=seq(0, 1, .01)) {
  
  fHrs <- c()
  ksd <- c()
  nfiles <- length(sample.files)
  for (i in 1:nfiles) {
    for (j in 1:nfiles) {
      if (j+i <= nfiles) {
        sample.a <- fread(sample.files[i])
        sample.b <- fread(sample.files[j + i])
        sample.a <- AddClonalMuts(sample.a, nclonal.muts, depth)
        sample.b <- AddClonalMuts(sample.b, nclonal.muts, depth)
        sample.a[, vaf := alt/depth]
        sample.b[, vaf := alt/depth]
        
        fHrs <- append(
          fHrs, sum(sample.a$vaf > .2)/sum(sample.a$vaf > .08) +
                sum(sample.b$vaf  > .2)/sum(sample.b$vaf  > .08))
        
        ksd <- append(
          ksd, max(abs(sapply(steps, function(x) sum(sample.a$vaf <= x)) -
                       sapply(steps, function(x) sum(sample.b$vaf <= x)))))
      }
    }
  }
  
  fHrs <- sum(fHrs)/length(fHrs)
  ksd <- sum(ksd)/length(ksd)
  
  return(data.frame(fHrs = fHrs, ksd = ksd))
}


CaclulateSFSstats <- function(sample.files,
                              nclonal.muts,
                              depth) {
  
  nfiles <- length(sample.files)
  D <- H <- E <- L <- rAUC <- rep(0, nfiles)
  
  for (i in 1:nfiles) {
    
    mutdata <- fread(sample.files[i])
    mutdata <- AddClonalMuts(mutdata, nclonal.muts, depth)
    sfs <-  mapply(function(x) sum(mutdata$alt == x), 0:depth)
    n <- depth - 1
    sfs <- sfs[1:n]/depth
    
    thetaP <- 2 * sum((1:n) * (n - 1:n) * sfs) / (n * (n - 1))
    thetaW <- sum(sfs) / sum((1:n)^-1) 
    thetaH <- 2 * sum(((1:n)^2) * sfs) / (n * (n - 1))
    thetaL <- sum((1:n) * sfs) / (n - 1)
    
    D[i] <- thetaP - thetaW
    H[i] <- thetaP - thetaH
    E[i] <- thetaL - thetaW
    L[i] <- thetaW - thetaH
  }
  return(data.frame(D = mean(D), H = mean(H), E = mean(E), L = mean(L)))
}


CheckRootBinary <- function(tree) {
  
  if (!(is.rooted(tree))) {
    tree <- root(tree, Ntip(tree) + 1, resolve.root = T)
  }
  if (!(is.binary.tree(tree))) {
    tree <- multi2di(tree)
  }
  return(tree)
}


CalculateTreeStats <- function(tree.file) {
  
  tree <- read.tree(tree.file)
  tree <- ape::collapse.singles(CheckRootBinary(tree))
  tree$edge.length[which(tree$edge.length <= 0)] <- 0.0000001
  tree.shape <- apTreeshape::as.treeshape(tree, model = 'yule')
  
  tree.stats <- data.frame(
    sackin = apTreeshape::sackin(tree.shape),
    sackin_pda = apTreeshape::sackin(tree.shape,
                                    norm = 'pda'),
    sackin_yule = apTreeshape::sackin(tree.shape,
                                     norm = 'yule'),
    colless = apTreeshape::colless(tree.shape),
    colless_pda = apTreeshape::colless(tree.shape,
                                      norm = 'pda'),
    colless_yule = apTreeshape::colless(tree.shape,
                                       norm = 'yule'),
    tip_coph = mean(GetMode(
      cophenetic.phylo(tree)[lower.tri(cophenetic.phylo(tree))])),
    node_coph = mean(GetMode(
      dist.nodes(tree)[lower.tri(dist.nodes(tree))])),
    branching_time = mean(ape::branching.times(ape::chronopl(tree,
                                                             lambda = .1))),
    mean_branch_length = mean(tree$edge.length),
    max_branch_length = max(tree$edge.length),
    max_node_depth = max(ape::node.depth(tree)))
  
  return(tree.stats)
}


ConvertMutsBinaryToNexus <- function(muts.binary.file) {
  
  muts.matrix <- as.data.frame(fread(muts.binary.file))
  muts.matrix$normal <- 0
  muts.matrix$muts <- NULL
  
  muts.matrix <- muts.matrix[, ncol(muts.matrix):1]
  
  muts.matrix[muts.matrix == 0] <- 'c'
  muts.matrix[muts.matrix == 1] <- 't'
  
  muts.matrix[] <- lapply(muts.matrix, factor)
  
  muts.matrix.list <- as.list(as.data.frame(muts.matrix))
  
  write.nexus.data(muts.matrix.list,
                   paste0(muts.binary.file, '.nex'),
                   datablock = T)
}



# Function to plot one_side_violin borrowed from stackoverflow ------------

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x - width / 2,
                     xmax = x)
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data, 
                              xmaxv = x,
                              xminv = x + violinwidth * (xmin - x))
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                             plyr::arrange(transform(data, x = xmaxv), -y))
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1,])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                            alpha = NA, linetype = "solid"),
          
          required_aes = c("x", "y")
  )
