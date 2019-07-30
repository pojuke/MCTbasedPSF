#################################################################################
#### Effects of soil microbes on plant competition: a perspective from modern coexistence theory
#### Ke & Wan (accepted July 2019) Ecological Monographs
####
#### This R script creates Figure 3 in the main text
#### Figure 3: Examples of how plant???soil microbe interactions and plant???plant competition together determine competition outcome
#################################################################################



######################################
#### Load packages
######################################
library('ggplot2')
library('grid')
library('gridExtra')
library('magrittr')
library('dplyr')
library('RColorBrewer')
library('cowplot')



######################################
#### 1. Plot setting and supporting functions
######################################
#### Define functions to calculate niche overlap. 
#### Input: paramters 
#### Output: niche overlap based on MCT (Chesson, 2000)
#### Note that we used the notation where competition is represented by negative c_{ij}
get.rho <- function(cAA = NA, cAB = NA, cBA = NA, cBB = NA,
                    sAA = NA, sAB = NA, sBA = NA, sBB = NA,
                    phiA = NA, phiB = NA){
  result <- sqrt((cBA + sBA*phiA) * (cAB + sAB*phiB) / ((cAA + sAA*phiA) * (cBB + sBB*phiB)))
  return(result)
}

#### Define functions to calculate niche overlap. 
#### Input: paramters
#### Output: fitness ratio based on MCT (Chesson, 2000)
#### Note that the function calculates f_{A}/f_{B} if "A.over.B"=T, but flips the parameters to calculate f_{B}/f_{A} if "A.over.B"=F
get.fit <- function(cAA = NA, cAB = NA, cBA = NA, cBB = NA,
                    sAA = NA, sAB = NA, sBA = NA, sBB = NA,
                    phiA = NA, phiB = NA, A.over.B = T){
  if (!A.over.B) return(get.fit(cAA = cBB, cAB = cBA, cBA = cAB, cBB = cAA,
                                sAA = sBB, sAB = sBA, sBA = sAB, sBB = sAA,
                                phiA = phiB, phiB = phiA))
  result <- sqrt((cBB + sBB*phiB) * (cBA + sBA*phiA) / ((cAA + sAA*phiA) * (cAB + sAB*phiB)))
  return(result)
}

#### Function to generate data frame for plotting 
#### Input: (1) "params": default parameters used for the calculations
####        (2) "vary.variable": the c_{ij} parameter that would be varied
####        (3) "vary.values": the range that the c_{ij} parameter would be varied
####        (4) "strengths": the range that the sigma_{ij} parameter would be varied (i.e., the strength of the PSF scenario)
####        (5) "get.params": a function for each scenario that modifies the sigma_{ij} parameter for that scenario
####        (6) "A.over.B": determining whether to calculate f_{A}/f_{B} ("A.over.B"=T) or f_{B}/f_{A} ("A.over.B"=F)
#### Output: a data.frame with four columns: (i) strength: the sigma_{ij} spectrum; (ii) vary: the c_{ij} spectrum; (iii) niche overlap; (iv) fitness ratio
plot.frame <- function(params = NA, 
                       vary.variable = NA, 
                       vary.values = NA, 
                       strengths = NA, 
                       get.params = NA, 
                       A.over.B = T){
  
  # Generate a grid of varying sigma_{ij} (i.e., set by "strength") and c_{ij} coefficients (i.e., set by "vary.values")
  grid <- expand.grid(strength=strengths, vary=vary.values)
  
  # Modify the sigma_{ij} parameters that we are interested in (i.e., set by each "get.params" function) based on the PSF strength spectrum
  new.params <- get.params(params, grid$strength)
  
  # Modify the c_{ij} parameters that we are interested in (i.e., set by "vary.variable") based on the PSF strength spectrum
  new.params[[vary.variable]] <- grid$vary
  
  # Calculate niche overlap and fitness ratio
  rhos <- do.call(get.rho, new.params)
  new.params[['A.over.B']] <- A.over.B
  fits <- do.call(get.fit, new.params)
  
  # Put results and parameters together in a data frame for plotting
  result <- cbind(grid, data.frame(rho = rhos, fit = fits))
  return(result)
}

#### Generate a function that streamlines the calculation for each c_{ij}
#### Input: (1) name: the name vector of c_{ij}s that would be varied 
####        (2) fn: the function that calls the "plot.frame" function for each c_{ij}
#### Output: a list, names are the c_{ij} that were varied, each element is a data.frame created by "plot.frame"
name.apply <- function(name, fn){
  result <- lapply(name, fn)
  names(result) <- name
  return(result)
}

#### Supporting functions of the main plot -- make c_{ij} a factor
abs.order.factor <- function(x) {factor(x, unique(x) %>% .[order(abs(.))])}

#### Supporting functions of the main plot -- plot fitness ratio ~ sigma_{ij}
#### Input: a data.frame created by "plot.frame"
#### Output: a ggplot of fitness ratio ~ sigma_{ij}
fit.plot <- function(frame, A.over.B = T){
  result <- 
    ggplot(frame, 
           aes(x = strength, 
               y = fit, 
               group = abs.order.factor(vary), 
               color = abs.order.factor(vary))) +
    geom_abline(intercept = 1.0, slope = 0, linetype = 2) +
    geom_line(size = 1) + 
    ylab(ifelse(A.over.B,
                expression(paste('fitness ratio, ', italic('f'['A']/'f'['B']))),
                expression(paste('fitness ratio, ', italic('f'['B']/'f'['A'])))))
  return(result)
}

#### Supporting functions to generate the main plot -- plot niche overlap ~ sigma_{ij}
#### Input: a data.frame created by "plot.frame"
#### Output: a ggplot of niche overlap ~ sigma_{ij}
rho.plot <- function(frame){
  result <- ggplot(frame, 
                   aes(x = strength, 
                       y = rho, 
                       group = abs.order.factor(vary), 
                       color = abs.order.factor(vary))) + 
    geom_line(size = 1) +
    ylab(expression(paste('niche overlap, ', italic(rho))))
  return(result)
}

#### Supporting functions to generate the main plot -- plot niche/fitness difference on its parameter space
#### Input: a data.frame created by "plot.frame"
#### Output: a ggplot of niche/fitness difference on its parameter space
battleaxe.plot <- function(frame, 
                           xlim = c(-0.5, 1), ylim = c(0, 2), 
                           start.strength = -2, end.strength = 0, cross.linetype = 2,
                           line.avoid.radius = 0.05, arrow.avoid.radius = 0.04,
                           A.over.B = T){
  
  # Create data used to plot the competitive outcome boundaries
  my.seq <- seq(1 - min(xlim), 1 - max(xlim), by = -0.001)
  battleaxe <- data.frame(rho = my.seq, bottom = my.seq, top = 1 / my.seq)
  
  # Clean up/rearrange the data frame based on plot range
  filtered.data <- filter(frame, 
                          strength >= min(start.strength, end.strength), 
                          strength <= max(start.strength, end.strength))
  if (start.strength < end.strength) {
    filtered.data <- arrange(filtered.data, strength, abs(vary))
  } else {
    filtered.data <- arrange(filtered.data, desc(strength), abs(vary))
  }
  
  # Find starting point for the plot
  starts <- 
    filtered.data %>% arrange(abs(strength - start.strength)) %>% 
    group_by(vary) %>% slice(1) %>% ungroup()
  
  # Find end point for the plot
  ends <- 
    filtered.data %>% arrange(abs(strength - end.strength)) %>%
    group_by(vary) %>% slice(1) %>% ungroup()
  
  # Combine all date
  enriched.data <- filtered.data %>% 
    left_join(transmute(starts, vary = vary, start.rho = rho, start.fit = fit), by = 'vary') %>%
    left_join(transmute(ends, vary = vary, end.rho = rho, end.fit = fit), by = 'vary') %>%
    mutate(start.dist = sqrt((rho - start.rho)^2 + (fit - start.fit)^2),
           end.dist = sqrt((rho - end.rho)^2 + (fit - end.fit)^2)) %>%
    mutate(min.dist = pmin(start.dist, end.dist))
  
  # Plot
  result <- 
    ggplot() +
    
    # Make axe-shaped plot
    geom_ribbon(data = battleaxe, 
                aes(x = 1 - rho, 
                    ymin = bottom, 
                    ymax = top), 
                alpha=0.1) + 
  
    # Add cross for referencing values
    geom_hline(yintercept = 1.0, linetype = cross.linetype) +
    geom_vline(xintercept = 0, linetype = cross.linetype) + 
  
    # Add line between points
    geom_line(data = filter(enriched.data, min.dist>line.avoid.radius), 
              aes(x = 1 - rho, 
                  y = fit, 
                  group = abs.order.factor(vary), 
                  color = abs.order.factor(vary)), 
              size=1
              ) + 
  
    # Add arrow
    geom_path(data = filter(enriched.data, end.dist <= 2*arrow.avoid.radius, end.dist >= arrow.avoid.radius), 
              aes(x = 1 - rho, 
                  y = fit, 
                  group = abs.order.factor(vary), 
                  color = abs.order.factor(vary)), 
              size = 0, 
              arrow = arrow(type = 'closed', length = unit(6, 'pt')), 
              show.legend = F) + 
  
    # Add start point
    geom_point(data = starts, 
               aes(x = 1 - rho, 
                   y = fit, 
                   color = abs.order.factor(vary)), 
               show.legend = F,
               size = 2, stroke = 1, shape = 21, fill = 'transparent') + 
  
    # Add end point
    geom_point(data = ends, 
               aes(x = 1 - rho, 
                   y = fit, 
                   color = abs.order.factor(vary)), 
               show.legend = F, 
               size = 2) +
    coord_cartesian(xlim = xlim, ylim = ylim, expand = F) + 
  
    # Make theme nice
    cowplot::theme_cowplot() + 
    theme(legend.position = 'bottom', legend.title.align = 0.5) +
    ylab(ifelse(A.over.B,
                expression(paste('fitness ratio, ', italic('f'['A']/'f'['B']))),
                expression(paste('fitness ratio, ', italic('f'['B']/'f'['A']))))) +
    xlab(expression(paste('niche difference, ', 1-italic(''~rho))))
  
  # Return plot
  return(result)
}

#### Complex composite plot recipes 
full.plot <- function(
  # All data frames for each c_{ij}
  all.frames, 
  
  # Component plot parameters
  xlim = NA, rho.ylim = NA, fit.ylim = NA, xlabel = 'strength',
  
  # Battleaxe plot parameters
  battleaxe.xlim = c(-0.5, 1), battleaxe.ylim = c(0, 2), 
  start.strength = -1, end.strength = 0,
  line.avoid.radius = 0.05, arrow.avoid.radius = 0.04, cross.linetype = 2,
  
  # Color
  palettes = NA, 
  panel.labels = c('A', 'B', 'C', 'D'), 
  panel.font = gpar(fontsize = 16, 
                    fontfamily = 'Helvetica', 
                    fontface = 'bold'),
  horiz.margin = 0.25, vert.margin = 0, units = 'in', 
  A.over.B = T,
  row.labels = list('cAA' = 'A', 'cAB' = 'B', 'cBA' = 'C', 'cBB' = 'D'),
  my.theme = cowplot::theme_cowplot(),
  display = F){

  # Calculate ranges for x and y
  xlim <- c(start.strength, end.strength)
  if (any(is.na(xlim))) xlim <- range(sapply(all.frames, function(x) x$strength))
  if (any(is.na(rho.ylim))) rho.ylim <- range(sapply(all.frames, function(x) x$rho))
  if (any(is.na(fit.ylim))) fit.ylim <- range(sapply(all.frames, function(x) x$fit))
  
  # Generate base panel ggplots for niche overlap and fitness ratio
  fit.coord <- coord_cartesian(xlim = xlim, ylim = fit.ylim, expand = F)
  rho.coord <- coord_cartesian(xlim = xlim, ylim = rho.ylim, expand = F)
  rho.plots <- name.apply(names(all.frames), 
                          function(coef){rho.plot(all.frames[[coef]]) + 
                              rho.coord + 
                              palettes[[coef]] +
                              xlab(xlabel)}
                          )
  
  fit.plots <- name.apply(names(all.frames), 
               function(coef){
                 fit.plot(all.frames[[coef]], A.over.B=A.over.B) + 
                   fit.coord + 
                   palettes[[coef]] + 
                   xlab(xlabel)}
               )
  
  if(xlim[1] > xlim[2]) {
    fit.plots <- plyr::llply(fit.plots, function(x){x + scale_x_reverse()})
    rho.plots <- plyr::llply(rho.plots, function(x){x + scale_x_reverse()})
  }
  
  # Battleaxes
  battleaxes <- name.apply(names(all.frames), 
                           function(coef){battleaxe.plot(all.frames[[coef]], 
                                                         xlim = battleaxe.xlim, ylim = battleaxe.ylim,
                                                         start.strength = start.strength, end.strength = end.strength, cross.linetype=cross.linetype, 
                                                         line.avoid.radius = line.avoid.radius, arrow.avoid.radius = arrow.avoid.radius, 
                                                         A.over.B = A.over.B) + 
                                          palettes[[coef]]}
                           )
 
  # Adjust margins and x title spacing
  mar <- margin(l = horiz.margin, r = horiz.margin, t = vert.margin, b = vert.margin, unit = units)
  thm <- my.theme +
    theme(plot.margin = mar,
          axis.title.x = element_text(margin = margin(t = 8, unit = 'pt')),
          legend.position = 'bottom', 
          legend.title.align = 0.5)
  
  # Create multipanel figure
  arrange.func <- ifelse(display, grid.arrange, arrangeGrob)
  
  # Hack to hide legend without affecting layout
  hide.legend <- guides(colour = guide_legend(override.aes = list(alpha = 0), 
                                              label.theme = element_blank(), title.theme = element_blank()))
  
  # Gather rho plots
  rho.plots <- name.apply(names(rho.plots), 
                          function(coef){rho.plots[[coef]] + 
                                         hide.legend + 
                                         ggtitle(row.labels[[coef]])}
                          )
  
  # Gather fitness ratio plots
  fit.plots <- plyr::llply(fit.plots, 
                           function(x){x + ggtitle(' ')})
  
  # Gather niche/fitness difference parameter space plots
  battleaxes <- plyr::llply(battleaxes, 
                            function(x){x + hide.legend + ggtitle(' ')})
  
  # Combine all plots
  all.plots <- name.apply(names(row.labels), 
                          function(coef){list(rho.plots[[coef]], 
                                              fit.plots[[coef]], 
                                              battleaxes[[coef]])}
    ) %>% unlist(recursive=F)
  all.plots <- lapply(all.plots, function(p){p + thm})
  result <- do.call(arrange.func, c(all.plots, list(nrow=4)))
  
  # Return result if we are not displaying figure
  if (!display){return(result)}
}

#### General parameters and plot setting for all plots
spacer <- '  '
row.labels <- list(
  'cAA' = expression(paste(bold('(a) '), '  Changing ', italic('c'['AA']), 
                         ' (effect of intraspecific competition on plant A)')),
  'cAB' = expression(paste(bold('(b) '), '  Changing ', italic('c'['AB']),
                         ' (effect of interspecific competition on plant A)')),
  'cBB' = expression(paste(bold('(c) '), '  Changing ', italic('c'['BB']),
                         ' (effect of intraspecific competition on plant B)')),
  'cBA' = expression(paste(bold('(d) '), '  Changing ', italic('c'['BA']),
                         ' (effect of interspecific competition on plant B)')))
plot.theme <- 
  cowplot::theme_cowplot(font_size = 10) + 
  theme(legend.box.spacing = unit(0.0625 / 2, 'in'), 
        plot.title = element_text(hjust = 0.08, vjust = 2, size = 14))
palette.names <- c('cAA' = 'Reds', 
                   'cAB' = 'Blues',
                   'cBA' = 'Greens', 
                   'cBB' = 'Purples')
scale.labels <- c('cAA' = expression(paste(italic('c'['AA']), '   ')),
                  'cAB' = expression(paste(italic('c'['AB']), '   ')),
                  'cBA' = expression(paste(italic('c'['BA']), '   ')),
                  'cBB' = expression(paste(italic('c'['BB']), '   ')))
A.over.B <- F

#### Color setting for all plots
get.colors <- function(coef, total.colors = 7, skip = 2) {
  labels <- paste(as.character(coef.values[[coef]]), spacer)
  breaks <- coef.values[[coef]]
  indices <- seq(total.colors - (length(breaks) - 1) * skip, total.colors, by = skip)
  colors <- brewer.pal(total.colors, palette.names[[coef]])[indices]
  result <- scale_color_manual(
    scale.labels[[coef]], values = colors,
    breaks = breaks)
  return(result)
}



######################################
#### 2. Janzen-Connell scenario
######################################
#### "get.params" function for Janzen-Connell: sigma_{AA} & sigma_{BB}
get.JC.params <- function(params, strengths){
  new.params <- params
  new.params[['sAA']] <- strengths
  new.params[['sBB']] <- strengths
  return(new.params)
}

#### Define default parameters within a list
default.params <- list(
  cAA = -1, cAB = -1, cBA = -0.8, cBB = -0.8,
  sAA = -0.4, sAB = -0.4, sBA = -0.5, sBB = -0.5,
  phiA = 1.6, phiB = 2.0)

#### Generate plotting range of the Janzen-Connell parameters
min.str.JC <- -0.32
max.str.JC <- -6
strengths.JC <- seq(min.str.JC, max.str.JC, by = -0.001)

#### Define values for competition coefficient (varied for each panel)
my.values <- c(-0.5, -1.0, -1.5)
coef.values <- list('cAA' = my.values, 'cAB' = my.values,
                    'cBA' = my.values, 'cBB' = my.values)

#### Generate plotting data
JC.frames <- name.apply(names(coef.values), 
                        function(coef){plot.frame(params = default.params, 
                                                  vary.variable = coef, 
                                                  vary.values = coef.values[[coef]], 
                                                  get.params = get.JC.params,
                                                  strengths = strengths.JC, 
                                                  A.over.B = A.over.B)}
                        )

#### Generate palettes for each parameter
JC.palettes <- name.apply(names(palette.names), get.colors)

#### Generate x-axis label
JC.label <- expression(
  paste('Janzen-Connell effect, ', italic(sigma['AA']) == italic(''~sigma['BB'])))

#### Plot in full
JC.full <- full.plot(JC.frames, 
                     xlim = c(min.str, max.str), xlabel = JC.label,
                     battleaxe.xlim = c(-0.5, 1), battleaxe.ylim = c(0.73, 1.27),
                     start.strength = min.str.JC, end.strength = max.str.JC,
                     cross.linetype = 2, 
                     palettes = JC.palettes,
                     horiz.margin = 0.125, vert.margin = 0.0625, 
                     A.over.B = A.over.B,
                     row.labels = row.labels, 
                     my.theme = plot.theme)

#### Save plots
ggsave('MCTbasedPSF_FigS1.pdf', JC.full, width = 7.5, height = 10, units = 'in')



######################################
#### 3. Enemy release scenario
######################################
#### "get.params" function for enemy release: sigma_{BA} & sigma_{BB}
get.ER.params <- function(params, strengths){
  new.params <- params
  new.params[['sBA']] <- strengths
  new.params[['sBB']] <- strengths
  return(new.params)
}

#### Define default parameters within a list
default.params <- list(
  cAA = -1.5, cAB = -1, cBA = -1, cBB = -1.5,
  sAA = -0.5, sAB = -0.5, sBA = -2, sBB = -2,
  phiA = 0.5, phiB = 0.5)

#### Generate plotting range of the enemy release parameters
min.str.ER <- -2
max.str.ER <- 0
strengths.ER <- seq(min.str.ER, max.str.ER, by = 0.001)

#### Define values for competition coefficient (varied for each panel)
my.values <- c(-0.5, -1, -1.5)
coef.values <- list('cAA' = my.values, 'cAB' = my.values,
                    'cBA' = my.values, 'cBB' = my.values)

#### Generate plotting data
ER.frames <- name.apply(names(coef.values), 
                        function(coef){plot.frame(params = default.params, 
                                                  vary.variable = coef, 
                                                  vary.values = coef.values[[coef]], 
                                                  get.params = get.ER.params,
                                                  strengths = strengths.ER, 
                                                  A.over.B = A.over.B)}
                        )

#### Generate palettes for each parameter
ER.palettes <- name.apply(names(palette.names), get.colors)

#### Generate x-axis label
ER.label <- expression(
  paste('enemy impact on B, ', italic(sigma['BA']) == italic(''~sigma['BB'])))

#### Plot in full
ER.full <- full.plot(ER.frames, 
                     xlim = c(min.str, max.str), xlabel = ER.label,
                     battleaxe.xlim = c(-0.225, 0.55), battleaxe.ylim = c(0,2.25), 
                     start.strength = min.str.ER, end.strength = max.str.ER, 
                     cross.linetype = 2, 
                     palettes = ER.palettes,
                     horiz.margin = 0.125, vert.margin = 0.0625,
                     A.over.B = A.over.B, 
                     row.labels = row.labels, 
                     my.theme = plot.theme)

#### Save plots
ggsave('MCTbasedPSF_FigS2.pdf', ER.full, width = 7.5, height = 10, units = 'in')



######################################
#### 4. Mutual facilitation scenario
######################################
#### "get.params" function for mutual facilitation: sigma_{AB} & sigma_{BA}
get.MF.params <- function(params, strengths) {
  new.params <- params
  new.params[['sAB']] <- strengths
  new.params[['sBA']] <- strengths
  return(new.params)
}

#### Define default parameters within a list
default.params <- list(
  cAA = -1, cAB = -1, cBA = -1, cBB = -1,
  sAA = 0.5, sAB = 0.5, sBA = 0.5, sBB = 0.5,
  phiA = 0.15, phiB = 0.2)

#### Generate plotting range of the mutual facilitation parameters
min.str.MF <- 0
max.str.MF <- 2
strengths.MF <- seq(min.str.MF, max.str.MF, by = 0.001)

#### Define values for competition coefficient (varied for each panel)
my.values <- c(-0.5, -1, -1.5)
coef.values <- list('cAA' = my.values, 'cAB' = my.values,
                    'cBA' = my.values, 'cBB' = my.values)

#### Generate plotting data
MF.frames <- name.apply(names(coef.values), 
                        function(coef){plot.frame(params = default.params, 
                                                  vary.variable = coef, 
                                                  vary.values = coef.values[[coef]],
                                                  get.params = get.MF.params,
                                                  strengths = strengths.MF, 
                                                  A.over.B = A.over.B)}
                        )

#### Generate palettes for each parameter
MF.palettes <- name.apply(names(palette.names), get.colors)

#### Generate x-axis label
MF.label <- expression(
  paste('mutualist faciliation, ', italic(sigma['AB']) == italic(''~sigma['BA'])))

#### Plot in full
MF.full <- full.plot(MF.frames, 
                     xlim = c(min.str, max.str), xlabel = MF.label,
                     battleaxe.xlim = c(-0.75, 1), battleaxe.ylim = c(0.3, 1.8),
                     start.strength = min.str.MF, end.strength = max.str.MF, 
                     cross.linetype = 2,
                     palettes = MF.palettes,
                     horiz.margin = 0.125, vert.margin = 0.0625,
                     A.over.B = A.over.B, 
                     row.labels = row.labels, 
                     my.theme = plot.theme)

#### Save plots
ggsave('MCTbasedPSF_FigS3.pdf', MF.full, width = 7.5, height = 10, units = 'in')



######################################
#### 5. Soil conditioning scenario
######################################
#### "get.params" function for soil conditioning: phi_{B}
get.SC.params <- function(params, strengths) {
  new.params <- params
  new.params[['phiB']] <- strengths
  return(new.params)
}

#### Define default parameters within a list
default.params <- list(
  cAA = -1, cAB = -1, cBA = -0.9, cBB = -0.9,
  sAA = 0.5, sAB = 0.25, sBA = 0.10, sBB = 0.20,
  phiA = 0.025, phiB = 0.025)

#### Generate plotting range of the soil conditioning parameters
min.str.SC <- 0.025
max.str.SC <- 2.5
by.str <- 0.0001
strengths.SC <- seq(min.str.SC, max.str.SC, by = by.str)

#### Define values for competition coefficient (varied for each panel)
my.values <- c(-0.8, -1, -1.2)
coef.values <- list('cAA' = my.values, 'cAB' = my.values,
                    'cBA' = my.values, 'cBB' = my.values)

#### Generate plotting data
SC.frames <- name.apply(names(coef.values), 
                        function(coef){plot.frame(params = default.params,
                                                  vary.variable = coef, 
                                                  vary.values = coef.values[[coef]], 
                                                  get.params = get.SC.params,
                                                  strengths = strengths.SC,
                                                  A.over.B = A.over.B)}
                        )

#### Generate palettes for each parameter
SC.palettes <- name.apply(names(palette.names), get.colors)

#### Generate x-axis label
SC.label <- expression(
  paste('soil conditioning of B, ', italic(phi1['B'])))

#### Plot in full
SC.full <- full.plot(SC.frames, 
                     xlim = c(min.str, max.str), xlabel = SC.label,
                     battleaxe.xlim = c(-0.4, 0.4), battleaxe.ylim = c(0.6, 1.4),
                     start.strength = min.str.SC, end.strength = max.str.SC, 
                     line.avoid.radius = 0.02, arrow.avoid.radius = 0.01,
                     cross.linetype = 2,
                     palettes = SC.palettes,
                     horiz.margin = 0.125, vert.margin = 0.0625,
                     A.over.B = A.over.B, 
                     row.labels = row.labels, 
                     my.theme = plot.theme)

#### Save plots
ggsave('MCTbasedPSF_FigS4.pdf', SC.full, width = 7.5, height = 10, units = 'in')



######################################
#### 6. PLot Figure 3
######################################
#### Get individual niche/fitness difference plots
battleaxe.JC.cAA <- battleaxe.plot(JC.frames[['cAA']],
                                   xlim = c(-0.5, 1), ylim = c(0.73, 1.27),
                                   start.strength = min.str.JC, end.strength = max.str.JC, 
                                   A.over.B = F)
battleaxe.ER.cBA <- battleaxe.plot(ER.frames[['cBA']],
                                   xlim = c(-0.5, 1), ylim = c(0, 2),
                                   start.strength = min.str.ER, end.strength = max.str.ER, 
                                   A.over.B = F)
battleaxe.MF.cAB <- battleaxe.plot(MF.frames[['cAB']],
                                   xlim = c(-0.75, 0.75), ylim = c(0.2, 1.4),
                                   start.strength = min.str.MF, end.strength = max.str.MF, 
                                   A.over.B = F)
battleaxe.SC.cAB <- battleaxe.plot(SC.frames[['cAB']],
                                   xlim = c(-0.201, 0.401), ylim = c(0.65, 1.35),
                                   start.strength = min.str.SC, end.strength = max.str.SC, 
                                   line.avoid.radius = 0.02, arrow.avoid.radius = 0.01,
                                   A.over.B = F)

#### Adjust margins and x title spacing
horiz.margin <- 0.25; vert.margin <- 0.125; units <- 'in'
mar <- margin(l = horiz.margin, r = horiz.margin, t = vert.margin, b = vert.margin, unit=units)
thm <- cowplot::theme_cowplot() + 
       theme(plot.margin = mar, 
             axis.title.x = element_text(margin = margin(t = 8, unit = 'pt')),
             legend.position = 'bottom', 
             legend.title.align = 0.5)

#### Throwaway function to add labels to plot panels
panel.font <- gpar(fontsize = 16, fontfamily = 'Helvetica')
add.label <- function(p, label) arrangeGrob(p, 
                                            top = textGrob(label, just = c("left", "top"),
                                                           x = unit(6, 'pt') + unit(horiz.margin,units),
                                                           y = unit(24, 'pt') - unit(vert.margin,units)))

#### Assemble the plot
all.together <- arrangeGrob(
  add.label(battleaxe.JC.cAA + thm + JC.palettes[['cAA']], expression(paste(bold('(a)'), '   Janzen-Connell'))), 
  add.label(battleaxe.ER.cBA + thm + ER.palettes[['cBA']], expression(paste(bold('(b)'), '   Enemy release'))), 
  add.label(battleaxe.MF.cAB + thm + MF.palettes[['cAB']], expression(paste(bold('(c)'), '   Mutual facilitation'))), 
  add.label(battleaxe.SC.cAB + thm + SC.palettes[['cAB']], expression(paste(bold('(d)'), '   Differential soil conditioning'))),
  ncol = 2)

#### Save the final plot!
ggsave('MCTbasedPSF_Fig3.pdf', plot = all.together, width = 7.5, height = 7.5)
