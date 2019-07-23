#################################################################################
#### Effects of soil microbes on plant competition: a perspective from modern coexistence theory
#### Ke & Wan (accepted July 2019) Ecological Monographs
####
#### This R script creates Figure 2 in the main text
#### Figure 2: Potential effects of soil microbes on the outcome of plant competition, 
####           visualized on the parameter space of niche difference (1-\rho, x-axis) and fitness ratio (fB/fA, y axis).
#################################################################################



######################################
#### Load packages
######################################
library(tidyr)
library(ggplot2)
library(cowplot)
library(Cairo)
library(grid)
library(gridExtra)



######################################
#### 1. Function setting for Chesson's plot with two-color ribbon
######################################
ggplot.MCT.stablizing.logscale <- function(YMIN){

  #### Create data for ribbon
  x.1 <- seq(1, 0, by=-0.001)
  x.2 <- seq(0, -1, by=-0.001)
  y1.1 <- 1/((1-x.1))
  y2.1 <- 1-x.1
  y1.2 <- 1/((1-x.2))
  y2.2 <- 1-x.2

  #### Coexistence ribbon
  cf.df.1 <- data.frame(my.x = c(x.1, x.1),
                        my.y = c(y1.1, y2.1),
                        whichy = c(rep("y1", times = length(seq(1, 0, by=-0.001))),
                                    rep("y2", times = length(seq(1, 0, by=-0.001)))))
  rib.dims.1 <- data.frame(min.dim = y1.1, 
                           max.dim = y2.1, 
                           x.dim = x.1)

  #### Priority effect ribbon
  cf.df.2 <- data.frame(my.x = c(x.2, x.2),
                        my.y = c(y1.2, y2.2),
                        whichy = c(rep("y1", times = length(seq(0, -1, by=-0.001))),
                                    rep("y2", times = length(seq(0, -1, by=-0.001)))))
  rib.dims.2 <- data.frame(min.dim = y1.2, 
                           max.dim = y2.2, 
                           x.dim = x.2)
  
  #### Plot
  ggplot() +

    # Shaded region and boundaries
    geom_ribbon(data = rib.dims.1,
                aes(x = x.dim,
                    ymin = min.dim,
                    ymax = max.dim),
                alpha=0.2) +
    geom_line(data = cf.df.1,
              aes(x = my.x,
                  y = my.y,
                  linetype = whichy,
                  group = whichy),
              col = "black") +
    geom_ribbon(data = rib.dims.2,
                aes(x = x.dim,
                    ymin = min.dim,
                    ymax = max.dim),
                alpha=0.2) +
    geom_line(data = cf.df.2,
              aes(x = my.x,
                  y = my.y,
                  linetype = whichy,
                  group = whichy),
              col = "black") +
    scale_linetype_manual(values = c("y2" = "solid",
                                     "y1" = "dotted")) +

    # Cross reference lines
    geom_hline(yintercept = 1.0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +

    # Scaling and other settings
    scale_y_log10() +
    coord_cartesian(expand = c(0,0), 
                    ylim=c(1/(YMIN*0.9), (YMIN)*0.9)) +
    theme(plot.margin = unit(c(1, 1, 1, 1), "lines")) +
    xlab(expression(paste('niche difference, ', 1-italic(rho)))) +
    ylab(expression(paste('fitness ratio, ', italic('f'['B']/'f'['A'])))) +
    theme(legend.position = "none")
  
}



######################################
#### 2. Plot Figure 2
######################################
#### Plot stylistic setting
mylims.x.Chesson <- 0.5
mylims.y.Chesson <- 1.6
cross.linetype <- 2
line.avoid.radius <- 0.015
arrow.avoid.radius <- 0.020

#### Point position setting
point.data <- data.frame(start.x = c(0.3, 0.12, -0.12),
                         start.y = c(1 /(1-0.3), 1/(1-0.15), (1+0.15)) + 0.08,
                         end.x = c(0.45, 0.15, -0.20),
                         end.y = c(1.45, 0.9, 0.7))

#### Arrow position setting
length.data <- sqrt((point.data$start.x - point.data$end.x)^2 + (point.data$start.y - point.data$end.y)^2)
arrow.data <- data.frame(start.x = point.data$start.x + c(1, 1, -1) * abs(point.data$start.x - point.data$end.x) * (line.avoid.radius/length.data),
                         start.y = point.data$start.y + c(-1, -1, -1) * abs(point.data$start.y - point.data$end.y) * (line.avoid.radius/length.data),
                         end.x = point.data$end.x + c(-1, -1, 1) * abs(point.data$start.x - point.data$end.x) * (arrow.avoid.radius/length.data),
                         end.y = point.data$end.y + c(1, 1, 1) * abs(point.data$start.y - point.data$end.y) * (arrow.avoid.radius/length.data))

#### Text position setting
text.data <- data.frame(x = point.data$end.x + c(0.025, 0.04, 0.04),
                        y = point.data$end.y + c(0.04, 0.03, 0.025),
                        text = c('(a)', '(b)', '(c)'))

#### Plot
Chesson.conceptual <-

  # Plot Chesson framework using previous function
  ggplot.MCT.stablizing.logscale(2) +

  # Add text for coexistence/priority effect region:
  geom_text(aes(x = 0.0,
                y = exp(log(1.7)/2)),
            label = "Plant B wins",
            size = 5.0) +
  geom_text(aes(x = 0.0,
                y = exp(-log(1.7)/2)),
            label = "Plant A wins",
            size = 5.0) +
  geom_text(aes(x = 0.3,
                y = exp(log(1.7)/30)),
            label = "Coexistence",
            size = 5.0) +
  geom_text(aes(x = -0.3,
                y = exp(log(1.7)/30)),
            label = "Priority effect",
            size = 5.0) +
  
  # Add starting points
  geom_point(data = point.data,
             aes(x = start.x,
                 y = start.y),
             color = "#de2d26",
             show.legend = F,
             size = 2, stroke = 1,
             shape = 21, fill = 'transparent') +

  # Add end points
  geom_point(data = point.data,
             aes(x = end.x,
                 y = end.y),
             color = "#de2d26",
             size = 2,
             show.legend = F) +

  # Add arrows
  geom_segment(data = arrow.data,
               aes(x = start.x,
                   y = start.y,
                   xend = end.x,
                   yend = end.y),
               color = "#de2d26",
            size = 1,
            linejoin = 'mitre',
            arrow = arrow(type = 'closed', length = unit(7, 'pt'))) +

  # Add text for arrows
  geom_text(data = text.data,
            aes(x = x,
                y = y,
                label = text)) +

  # Coordinate and axes settings
  scale_y_continuous(breaks = c(0.7, 1, 1.4),
                     labels = c("0.7", "1", "1.4")) +
  scale_x_continuous(breaks = c(-0.5, -0.25, 0, 0.25, 0.5),
                     labels = c("-0.5", "-0.25", "0", "0.25", "0.5")) +
  coord_cartesian(expand = c(0, 0),
                  xlim = c(-mylims.x.Chesson, mylims.x.Chesson),
                  ylim = c(1/mylims.y.Chesson, mylims.y.Chesson)) +
  xlab(expression(paste('niche difference, ', 1-italic(rho)))) +
  ylab(expression(paste('fitness ratio, ', italic('f'['B']/'f'['A'])))) +

  # Other settings
  cowplot::theme_cowplot() +
  theme(legend.position = "none",
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))




######################################
#### 3. Save plots
######################################
dev.off()
CairoPDF(file = "MCTbasedPSF_Fig2.pdf", width = 6.5, height = 6.5)
Chesson.conceptual
dev.off()
