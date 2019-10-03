GeneratePlotsFig6 <- function(path.to.data,
                              path.to.output,
                              plot.type) {
  
  sampling.dirs <- c('needle', 'punch', 'tree', 'wholet')
  growth.dirs <- c('exponential', 'death', 'peripheral')
  parambase <- c('mu', 't', 's')
  
  param.sets <- list(exponential = parambase,
                     death = c(parambase, 'd'),
                     peripheral = c(parambase, 'a'))
  tumour.sets <- list(exponential = 1:13,
                      death = 14:24,
                      peripheral = 25:34)
  max.round <- 10
  path.base <- paste0(path.to.data,'/inferred_params/')
  
  params.target <- data.frame(fread(paste0(
    path.to.data,'/params_target.csv')))
  params.target$s  <- params.target$s + 1
  
  param.posteriors <- data.frame()
  param.errors <- data.frame()
  for (sampling in sampling.dirs) {
    for (growth in growth.dirs) {
      
      params <- unlist(param.sets[growth])
      tumours <- as.numeric(unlist(tumour.sets[growth]))
      for (param in params) {
        for (tumour in tumours) {
          
          path.to.data <- paste0(
            path.base, sampling, '/', growth, '/rec_', param,'/T', tumour,'/')
          print(path.to.data)
          
          r <- max.round
          repeat{
            smc.file <- paste0(path.to.data, 'smc_particles_r', r, '.csv')
            r <- r - 1
            if (file.exists(smc.file)) break
          }
          smc.particles <- as.data.frame(fread(smc.file))
          round.by <- ifelse(param %in% c('d','a'), 1, 0)
          
          param.est <- mean(GetMode(
            round(smc.particles[, param], round.by)))
          
          percent.error <- 
            100*(params.target[tumour, param]
                 - param.est)/params.target[tumour, param]
          
          cur.posterior <- data.frame(
            tumour = paste0('T', tumour),
            sampling = sampling,
            growth = growth,
            param = param,
            value = smc.particles[, param])
          
          cur.error <- data.frame(
            tumour = paste0('T', tumour), 
            sampling = sampling,
            growth = growth, 
            param = param, 
            value = percent.error)
          
          if (param == 't' & growth == 'exponential') {
            cur.posterior <- cur.posterior[cur.posterior$value <= 15,]
          } 
          param.posteriors <- rbind(param.posteriors, cur.posterior)
          param.errors <- rbind(param.errors, cur.error)
        }
      }
    }
  }
  
  param.errors$param <- factor(
    param.errors$param,
    levels = c('mu', 't', 's', 'd', 'a'),
    labels = c('Mutation rate',
               'Time to new clone',
               'Selective advantage',
               'Death rate',
               'Aggression'))
    
  if (plot.type == 'error_rate') {
    
    median.errors <- param.errors %>% 
      group_by(param, growth, sampling) %>% 
      mutate_at(vars(value), funs(MedianNArm))
    
    median.errors <- as.data.frame(median.errors)
    median.errors$value <- paste0(round(median.errors$value, 1), '%')
    
    median.errors <- median.errors[order(median.errors$param), ]
    median.errors.uniq <- data.frame(unique(median.errors[, -1]))
    median.errors.uniq$pos <- c(rep(rep(c(rep(130, 3),rep(100, 3)), 2), 2), 
                                rep(c(rep(90, 3),rep(78, 3)), 2),
                                rep(100, 4),
                                rep(130, 4))
    
    colours <- c('springgreen4', 'darkblue', 'salmon3', 'deeppink4')
    
    p <- ggplot(param.errors,
                aes(y = value,
                    x = growth,
                    fill = sampling,
                    colour = sampling)) +
      geom_violin(scale = 'width',
                  width=.8,
                  position = position_dodge(width = .8),
                  alpha = .5) +
      geom_point(aes(fill = sampling),
                 position = position_jitterdodge(dodge.width = 0.8),
                 size = 2,
                 alpha = .7) + 
      scale_fill_manual(values = colours) +
      scale_colour_manual(values = colours) +
      geom_text(data = median.errors.uniq,
                aes(x = growth,
                    y = pos,
                    label = value,
                    group = sampling), 
                size = 3,
                fontface = 'bold',
                position = position_dodge(width = .9),
                hjust = .5) + 
      facet_wrap(~param, scales = 'free', ncol = 2) +
      theme_bw() + 
      labs(x = '', y = '% error', fill = 'sampling') +
      theme(axis.text = element_text(size = 15),
            axis.title = element_text(size = 20),
            legend.text = element_text(size = 20),
            legend.title = element_text(size = 20, face = 'bold'),
            legend.position='none',
            legend.spacing.y = unit(.5, 'cm'),
            strip.text = element_text(size = 16, face = 'bold')) +
      geom_hline(yintercept = 0, linetype = 'dotted')
    
    ggsave(paste0(path.to.output, '/abcsmc_percentage_error.pdf'),
           p, width = 12, height = 10)
    
  } else if (plot.type == 'posterior') {
    
    observed.data <- as.data.table(params.target)
    setnames(observed.data, 'id', 'tumour')
    observed.data[, c('tumour', 'growth', 's') := 
                    list(paste0('T', tumour),
                         c(rep('exponential', 13),
                           rep('death', 11),
                           rep('peripheral', 10)),
                         s - 1)]
    observed.data.melt <- reshape2::melt(observed.data)
    setnames(observed.data.melt, 'variable', 'param')
    
    param.posteriors <- as.data.table(param.posteriors)
    
    for (growth.model in growth.dirs) {
      
      param.posteriors.sub <- param.posteriors[growth == growth.model]
      observed.data.melt.sub <- observed.data.melt[growth == growth.model]
      observed.data.melt.sub <- 
        observed.data.melt.sub[param %in% unlist(param.sets[growth.model]),]
      
      colours <- c('springgreen4', 'darkblue', 'salmon3','deeppink4')
      
      p <- ggplot(param.posteriors.sub,
                  aes(y = value,
                      x = tumour,
                      fill = sampling,
                      colour = sampling)) +
        geom_violin(trim = T,
                    scale = 'width',
                    width = .8,
                    position = position_dodge(width = .8),
                    alpha = .5) +
        scale_fill_manual(values = colours) +
        scale_colour_manual(values = colours) +
        geom_text(data = observed.data.melt.sub,
                  aes(x = tumour,
                      y = value,
                      label = value),
                  position = position_dodge(width = .9),
                  hjust = .5,
                  size = 5,
                  col = 'black',
                  alpha = 1) + 
        facet_wrap(~param,
                   scales = 'free',
                   ncol = 1) +
        theme_bw() + 
        labs(x = '', y = '% error', fill = 'sampling') +
        theme(axis.text = element_text(size = 15),
              axis.title = element_text(size = 20),
              legend.text = element_text(size = 20),
              legend.title = element_text(size = 20,
                                          face = 'bold'),
              legend.spacing.y = unit(.5, 'cm'),
              strip.text = element_text(size = 20)) +
        geom_hline(yintercept = 0, linetype = 'dotted')
      
      ggsave(paste0(path.to.output,'param_posteriors_',growth.model,'.pdf'),
             p, width = 12, height = 10)
    }
  }
}