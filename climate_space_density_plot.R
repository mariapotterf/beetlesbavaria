# climate space: get density plots


# get climatic characteristics
# temp
# prec
# from ERA-NET data


#library(terra)
#library(rnaturalearth)
#library(rnaturalearthdata)
#library(sf)
library(tidyr)
library(dplyr)
library(data.table)


library(PCICt)
#library(zoo)                # for as.Date() specification

library(ggplot2)
library(ggpubr)
library(stringr)

# get variables
spring.months         = 3:5
veg.months            = 3:10
study.period.extended = 2012:2021
study.period          = 2015:2021



# get data from Franconia ------------------------------------------------
df_clim <- fread('outTable/xy_clim_DWD.csv')
df_spei_months <- fread('outTable/xy_spei_all_DWD.csv')


# classify warm years
df_mean_spei <- df_spei_months %>% 
  dplyr::filter(scale == 3) %>% 
  dplyr::filter(month %in% veg.months) %>% 
  mutate(cat = case_when(year %in% study.period ~ 'study_period',
                         .default = "ref")) %>% 
  group_by(falsto_name, cat, year) %>%
  summarize(spei = mean(spei)) #%>%  # sum annual precip, convert from meters to mm




# classify warm years
df_mean <- df_clim %>% 
  dplyr::filter(month %in% veg.months) %>% 
  mutate(cat = case_when(year %in% study.period ~ 'study_period',
                         .default = "ref")) %>% 
  group_by(falsto_name, cat, year) %>%
  summarize(tmp = mean(tmp),# mean annual temp
         prcp = sum(prcp)) %>%
  right_join(df_mean_spei)

# time series plot: ----------------------------------------

# keep variation between locations per year
# classify warm years
df_mean_spei <- df_spei_months %>% 
  dplyr::filter(scale == 3) %>% 
  dplyr::filter(month %in% veg.months) %>% 
  mutate(cat = case_when(year %in% study.period ~ 'study_period',
                         .default = "ref")) %>% 
  group_by(falsto_name, cat, year) %>%
  summarize(spei = mean(spei)) #%>%  # sum annual precip, convert from meters to mm




# classify warm years
df <- df_clim %>% 
  dplyr::filter(month %in% veg.months) %>% 
  mutate(cat = case_when(year %in% study.period ~ 'study_period',
                         .default = "ref")) %>% 
  group_by(falsto_name, cat, year) %>%
  summarize(tmp = mean(tmp),# mean annual temp
            prcp = sum(prcp)) %>%
  right_join(df_mean_spei)


# get reference values for tmp and spei
tmp_ref  = df %>% 
  ungroup() %>% 
  filter(year %in% 1980:2010) %>% 
  summarise(tmp = mean(tmp)) %>% 
  pull()

spei_ref  = df %>% 
  ungroup() %>% 
  filter(year %in% 1980:2010) %>% 
  summarise(spei = mean(spei)) %>% 
  pull()


p.time.series.spei <- df %>% 
  ggplot(aes(x = year,
             y = spei,
             col = cat)) + 
  scale_color_manual(values = c('black', 'red')) +
  geom_hline(yintercept = spei_ref, lty = 'dashed', col = 'grey50' ) +
  stat_summary(fun.data = "mean_sdl") +
  theme_classic() +
  theme(legend.position = 'none') +
  labs(y =  "SPEI [dim.]", 
       x = '') +
  scale_x_continuous(n.breaks = 5)



p.time.series.tmp <- df %>% 
  ggplot(aes(x = year,
             y = tmp,
             col = cat)) + 
  
  stat_summary(fun.data = "mean_sdl") +
  geom_hline(yintercept = tmp_ref, lty = 'dashed', col = 'grey50' ) +
  
  scale_color_manual(values = c('black', 'red')) +
  theme_classic() +
  theme(legend.position = 'none') +
  labs(y =  bquote('Temperature [' * degree * 'C]'), 
       x = '') +
  scale_x_continuous(n.breaks = 5)
  

p.time <- ggarrange(p.time.series.tmp,
          p.time.series.spei, nrow = 2)

windows(4, 4)
p.time

ggsave(filename = 'outFigs/time_series.png', 
       plot = p.time, width = 4, height = 4, dpi = 300)




# get data fr reference
df_ref <- df_mean %>%
  dplyr::filter(cat == 'ref')  %>%
  dplyr::filter(year %in% 1980:2010) 
  


# get data for warm years
df_dist <- df_mean %>%
  dplyr::filter(cat == 'study_period') %>% 
  dplyr::filter(year %in% 2018:2020) 


rescale_density <- function(x) {
  scales::rescale(x, to = c(0.001, 1))
}

# get summary table for franconia
df_summary <- df_dist %>%
  mutate(pairID =  substr(falsto_name, 1, nchar(falsto_name) - 2)) %>% 
  group_by(pairID) %>%
  summarize(
    spei_mean = mean(spei),
    spei_sd = sd(spei),
    spei_min = mean(spei) - sd(spei),#min(tmp),
    spei_max = mean(spei) + sd(spei),#max(tmp),
    tmp_mean = mean(tmp),
    tmp_sd = sd(tmp),
    tmp_min = mean(tmp) - sd(tmp),#min(tmp),
    tmp_max = mean(tmp) + sd(tmp),#max(tmp),
    prec_mean = mean(prcp),
    prec_sd = sd(prcp),
    prec_min = mean(prcp) - sd(prcp),#min(prcp),
    prec_max = mean(prcp) + sd(prcp),#max(prcp),
    .groups = 'drop'
  )

df_dist
df_summary



# get final plot PRCP vs TMP --------------------------------------------------------------

library(MASS)

# Function to calculate levels
calculate_levels <- function(density, percentages) {
  sorted_density <- sort(density, decreasing = TRUE)
  cum_density <- cumsum(sorted_density) / sum(sorted_density)
  sapply(percentages, function(p) {
    sorted_density[max(which(cum_density <= p))]
  })
}

# Calculate 2D density
d <- kde2d(df_dist$tmp, df_dist$prcp, n = 500)
density_values <- d$z

# Set densities outside the range of your data to zero
density_values[density_values < 1e-6] <- 0

# Calculate the levels for specified percentages
levels <- calculate_levels(as.vector(density_values), c(0.50, 0.75, 0.90, 0.99, 1))

# Prepare data for ggplot
plot_data <- expand.grid(tmp = d$x, prec = d$y)
plot_data$density <- as.vector(density_values)

# Use cut to create factor levels, including one for zero density
plot_data$level <- cut(plot_data$density, breaks = c(-Inf, levels, Inf), labels = FALSE, include.lowest = TRUE)

# Define colors (from red to yellow, plus white for zero density)

library(RColorBrewer)
blue_colors <- brewer.pal(5, "Blues")

# Add 'white' at the beginning
my_colors <- c("white", blue_colors)


# Create the plot
p<-
  ggplot(plot_data) +
  geom_raster(aes(x = tmp, y = prec, fill = factor(level)), alpha = 0.8) +
  scale_fill_manual(values = my_colors,
                    labels = c("", rev(c("50%", "75%", "90%", "99%", "100%"))),
                    name = "") +
   geom_errorbar(data = df_summary,
                aes(x = tmp_mean,

                    ymin = prec_min,
                    ymax = prec_max, group = pairID),
                width = 0.1,
                linewidth = 0.2,
                color = 'grey50') +
  geom_errorbarh(data = df_summary,
                 aes(
                   y = prec_mean,
                   xmin = tmp_min,
                   xmax = tmp_max, group = pairID),
                 height =20,
                 linewidth = 0.2,
                 color = 'grey50') +
    geom_point(data = df_summary, 
               aes(x = tmp_mean, y = prec_mean, group = pairID), 
               color = 'white', size = 1.5) +
    geom_point(data = df_summary, 
               aes(x = tmp_mean, y = prec_mean, group = pairID), 
               color = 'black', size = 0.8) +
    
  ylab('Mean seasonal precipitation [mm]') +
  xlab(bquote('Mean seasonal temperature [' * degree * 'C]')) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_rect(fill = "white", colour = NA),
        aspect.ratio = 1,
        axis.title.x = element_text(size = 10),  # X-axis label
        axis.title.y = element_text(size = 10),  # Y-axis label
        legend.title = element_blank(),  # Legend title
        legend.text = element_text(size = 10),   #   # Legend item text)
        legend.key.size = unit(0.3, "cm"))  

windows(3,3)
(p)

ggsave("franc_in_clim_space.png", 
       plot = p, 
       device = "png", 
       units = "inch",
       dpi = 300,
       width = 12, height = 10)






# Get final plot: SPEI vs TMP ---------------------------------------------



# Calculate 2D density
d <- kde2d(df_dist$tmp, df_dist$spei, n = 500)
density_values <- d$z

# Set densities outside the range of your data to zero
density_values[density_values < 1e-6] <- 0

# Calculate the levels for specified percentages
levels <- calculate_levels(as.vector(density_values), c(0.50, 0.75, 0.90, 0.99, 1))

# Prepare data for ggplot
plot_data <- expand.grid(tmp = d$x, spei = d$y)
plot_data$density <- as.vector(density_values)

# Use cut to create factor levels, including one for zero density
plot_data$level <- cut(plot_data$density, breaks = c(-Inf, levels, Inf), labels = FALSE, include.lowest = TRUE)

# Define colors (from red to yellow, plus white for zero density)

library(RColorBrewer)
blue_colors <- brewer.pal(5, "Blues")

# Add 'white' at the beginning
my_colors <- c("white", blue_colors)


# Create the plot
#p<-
ggplot(plot_data) +
  geom_raster(aes(x = tmp, y = spei, fill = factor(level)), alpha = 0.6) +
  scale_fill_manual(values = my_colors,
                    labels = c("", rev(c("50%", "75%", "90%", "99%", "100%"))),
                    name = "") +
  geom_errorbar(data = df_summary,
                aes(x = tmp_mean,
                    ymin = spei_min,
                    ymax = spei_max, group = pairID),
                #width = 0.1,
                #linewidth = 0.2,
                color = 'grey50') #+
  # geom_errorbarh(data = df_summary, 
  #                aes(
  #                  y = spei_mean, 
  #                  xmin = spei_min, 
  #                  xmax = spei_max, group = pairID),
  #                height =20, 
  #                linewidth = 0.2,
  #                color = 'grey50') +
  geom_point(data = df_summary, 
             aes(x = tmp_mean, y = spei_mean, group = pairID), 
             color = 'white', size = 1.5) +
  geom_point(data = df_summary, 
             aes(x = tmp_mean, y = spei_mean, group = pairID), 
             color = 'black', size = 0.8) +
  ylab('Mean seasonal SPEI [dim.]') +
  xlab(bquote('Mean seasonal temperature [' * degree * 'C]')) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_rect(fill = "white", colour = NA),
        aspect.ratio = 1,
        axis.title.x = element_text(size = 10),  # X-axis label
        axis.title.y = element_text(size = 10),  # Y-axis label
        legend.title = element_blank(),  # Legend title
        legend.text = element_text(size = 10),   #   # Legend item text)
        legend.key.size = unit(0.3, "cm"))  

(p)

