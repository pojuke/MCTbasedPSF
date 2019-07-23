#################################################################################
#### Effects of soil microbes on plant competition: a perspective from modern coexistence theory
#### Ke & Wan (accepted July 2019) Ecological Monographs
####
#### This R script creates Figure 5b and 5c in the main text
#### Figure 5: Applying modern coexistence theory to understand the effects of soil microbes on plant competitive outcome using data from Aguilera et al. (2017) as a case study
#################################################################################



######################################
#### Load packages
######################################
library('ggplot2')
library('directlabels')
library('RColorBrewer')
library('viridis')
library('scales')
library('cowplot')
library('dplyr')
library('tibble')
library('ggrepel')



######################################
#### 1. Coefficients calculated from Aguilera 2017 
######################################
#### \alpha_{ij} calculated from sterile soil treatments
aguilera.sterile <- list(
  A_sat_sat = -0.2417684,
  A_sat_ser = -0.1379671,
  A_ser_sat = -0.236798,
  A_ser_ser = -0.1736271  
)

#### \alpha_{ij} calculated from live soil treatments
aguilera.live <- list(
  A_sat_sat = -0.2417684,
  A_sat_ser = -0.1942,
  A_ser_sat = -0.2429021,
  A_ser_ser = -0.2314036  
)

#### Set up some hypothetical scenarios, mixing coefficients from treatments
aguilera.ser.only <- aguilera.sterile
aguilera.ser.only['A_sat_ser'] <- aguilera.live['A_sat_ser']
aguilera.ser.only['A_ser_ser'] <- aguilera.live['A_ser_ser']

aguilera.sat.only <- aguilera.sterile
aguilera.sat.only['A_ser_sat'] <- aguilera.live['A_ser_sat']
aguilera.sat.only['A_sat_sat'] <- aguilera.live['A_sat_sat']

aguilera.no.spillover <- aguilera.sterile
aguilera.no.spillover['A_ser_ser'] <- aguilera.live['A_ser_ser']
aguilera.no.spillover['A_sat_sat'] <- aguilera.live['A_sat_sat']

aguilera.only.spillover <- aguilera.sterile
aguilera.only.spillover['A_ser_sat'] <- aguilera.live['A_ser_sat']
aguilera.only.spillover['A_sat_ser'] <- aguilera.live['A_sat_ser']

### A named list of scenarios of interest
scenarios <- list(
  'all' = aguilera.live,
  'none' = aguilera.sterile,
  'ser soil only' = aguilera.ser.only,
  'sat soil only' = aguilera.sat.only
)



######################################
#### 2. Calculate MCT components for scenarios
######################################
#### Helper function to calculate components
calculate.components <- function(coefficients) {
  result <- with(coefficients, {
    rho = sqrt(((A_sat_ser) * (A_ser_sat)) / ((A_sat_sat) * (A_ser_ser)))
    f2f1 = sqrt(((A_ser_ser) * (A_ser_sat)) / ((A_sat_sat) * (A_sat_ser)))
    ND = 1 - rho
    FR = f2f1
    c(ND = ND, FR = FR)
  })
  return(result)
}

#### Calculate for all scenarios and convert to data frame
components.all <- sapply(scenarios, calculate.components) %>%
  t() %>% as.data.frame() %>%
  rownames_to_column('scenario')

#### Separate out sterile soil
components.sterile <- filter(components.all, scenario == 'none')
components <- filter(components.all, scenario != 'none')

#### Plot setting: Distance to leave from head and tail of the arrow
line.avoid.radius <- 0.0025   
arrow.avoid.radius <- 0.0025 

#### Calculate arrow coordinates (slightly offset from start/end coordinates)
scenario.arrows <- components %>%
  
  # All scenarios start at the sterile point
  mutate(ND0 = components.sterile$ND, FR0 = components.sterile$FR) %>%
  
  # Calculate change in ND and FR, as well as length of line
  mutate(dND = ND - ND0, dFR = FR - FR0, length = sqrt(dND^2+dFR^2)) %>%
  
  # Offset arrow start/end by the required length
  mutate(ND1 = ND0 + (line.avoid.radius / length) * dND,
         FR1 = FR0 + (line.avoid.radius / length) * dFR,
         ND2 = ND - (arrow.avoid.radius / length) * dND,
         FR2 = FR - (arrow.avoid.radius / length) * dFR)

#### Put start and end points in tidy format
scenario.points <- components %>%
  mutate(ND = components.sterile$ND, FR = components.sterile$FR, point.type = 'start') %>%
  bind_rows(mutate(components, point.type = 'end'))



######################################
#### 3. Plot real and hypothetical scenarios
######################################
#### Plot Figure 5b
real.panel <- 
  ggplot() +
  
  # Plot Chesson parameter space
  geom_ribbon(data = data.frame(x = seq(0, 0.3, 0.001)),
              aes(x = x, 
                  ymin = 1-x, 
                  ymax = 1/(1-x)),
              fill = "grey", alpha = 0.3) +
  
  # Add starting and end points
  geom_point(data = filter(components.all, scenario %in% c('all', 'none')),
             aes(x = ND, 
                 y = FR, 
                 color = scenario, 
                 shape = scenario),
             size = 2, stroke = 1) +
  
  # Add arrows
  geom_segment(data = filter(scenario.arrows, scenario == 'all'), 
               aes(x = ND1, 
                   y = FR1, 
                   xend = ND2, 
                   yend = FR2, 
                   color = scenario),
               size = 1, linejoin = 'mitre',
               arrow = arrow(type = 'closed', length = unit(7, 'pt')),
               show.legend = F) +
  
  # Add text
  geom_text(aes(x = components.all[components.all$scenario == "none", ]$ND,
                y = components.all[components.all$scenario == "none", ]$FR + 0.01,
                label = "without\nmicrobes"),
            size = 4.0) +
  geom_text(aes(x = components.all[components.all$scenario == "all", ]$ND,
                y = components.all[components.all$scenario == "all", ]$FR - 0.01,
            label = "with\nmicrobes"),
            size = 4.0) +
  
  # Color and text setting
  scale_color_manual('microbial effects', 
                     values = c("all" = "black", 
                                "none" = "black",
                                "ser soil only" = "#b20016", 
                                "sat soil only" = "#9ecae1"),
                     breaks = c('none', 'sat soil only', 'ser soil only', 'all'),
                     guide = F) +
  scale_shape_manual('microbial effects', 
                     values = c("all" = 16, 
                                "none" = 1,
                                "ser soil only" = 16, 
                                "sat soil only" = 16),
                     breaks = c('none', 'sat soil only', 'ser soil only', 'all'),
                     guide = F) +
  
  # Coordinate and axes settings
  coord_cartesian(xlim = c(0.07, 0.125), ylim = c(1.07, 1.14)) +
  ylab(expression(paste("fitness ratio, ", italic(f[sat]), "/", italic(f[ser])))) +
  xlab(expression(paste("niche difference, ", 1-rho))) +
  ggtitle(bold('(b) ')~'observed (total) effects') + 
  theme(plot.title = element_text(face = 'plain', hjust=0))


#### Plot Figure 5c
hypothetical.panel <- 
  ggplot() +
  
  # Plot Chesson parameter space
  geom_ribbon(data = data.frame(x = seq(0, 0.3, 0.001)), 
              aes(x = x, 
                  ymin = 1-x,
                  ymax = 1/(1-x)),
              fill = "grey", alpha = 0.3) +
  
  # Add starting and end points
  geom_point(data = filter(components.all, scenario != 'all'), 
             aes(x = ND, 
                 y = FR, 
                 color = scenario,  
                 shape = scenario),
             size = 2, stroke = 1) +
  
  # Add arrows
  geom_segment(data = filter(scenario.arrows, scenario != 'all'), 
               aes(x = ND2 - 0.05 * dND, 
                   y = FR2 - 0.05 * dFR, 
                   xend = ND2, 
                   yend = FR2, 
                   color = scenario),
               size = 1, linejoin = 'mitre', 
               arrow = arrow(type = 'closed', length = unit(7, 'pt')),
               show.legend = F) +
  geom_segment(data = filter(scenario.arrows, scenario != 'all'), 
               aes(x = ND1, 
                   y = FR1, 
                   xend = ND2 - 0.05 * dND, 
                   yend = FR2 - 0.05 * dFR, 
                   color = scenario),
               size = 1, linejoin = 'mitre', linetype = '32',
               show.legend = F) +
  
  # Add text
  geom_text(data = filter(components, scenario != 'all'),
            aes(x = ND,
                y = FR + ifelse(FR > 1.12, 0.009, -0.005),
                label = ifelse(scenario == 'ser soil only', "italic('L. serriola')~soil", "italic('L. sativa')~soil")),
            parse = T) +  
  geom_text(data = filter(components, scenario != 'all'),
            aes(x = ND,
                y = FR + ifelse(FR > 1.12, 0.004, -0.010),
                label = "conditioning only")) +
  
  # Color and text setting
  scale_color_manual('microbial effects', 
                     values = c("all" = "black", 
                                "none" = "black",
                                "ser soil only" = "#b20016", 
                                "sat soil only" = "#9ecae1"),
                     breaks = c('none', 'sat soil only', 'ser soil only', 'all'),
                     guide = F) +
  scale_shape_manual('microbial effects', 
                     values = c("all" = 16, 
                                "none" = 1,
                                "ser soil only" = 16, 
                                "sat soil only" = 16),
                     breaks = c('none', 'sat soil only', 'ser soil only', 'all'),
                     guide = F) +
  
  # Coordinate and axes settings
  coord_cartesian(xlim = c(0.07, 0.125), ylim = c(1.07, 1.14)) +
  ylab(expression(paste("fitness ratio, ", italic(f[sat]), "/", italic(f[ser])))) +
  xlab(expression(paste("niche difference, ", 1-rho))) +
  ggtitle(bold('(c) ')~'soil conditioning effects') + 
  theme(plot.title=element_text(face = 'plain', hjust=0))



######################################
#### 4. Save plots
######################################
plot_grid(real.panel, hypothetical.panel)
ggsave("MCTbasedPSF_Fig5.pdf", width=7.5, height=4)
