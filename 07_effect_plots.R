# Make effect plots

library(ggplot2)
library(ggeffects)
library(ggpubr)

load('outData/fin_models.Rdata')


# define plotiing functions


# Create effect plot function with an additional argument to control y-axis labels
create_effect_plot <- function(data, line_color = "blue", x_title = "X-axis", y_title = "Y-axis", y_lim = c(100,80000), show_y_axis = TRUE) {
  p <- ggplot(data, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high)) +
    geom_line(color = line_color) +
    geom_ribbon(alpha = 0.1, fill = line_color) +
    labs(x = x_title) +
    ylim(y_lim) +
    theme_minimal(base_size = 10) +
    theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white", colour = "black"))
  
  if (show_y_axis) {
    p <- p + labs(y = y_title)
  } else {
    p <- p + theme(axis.title.y = element_blank(), axis.text.y = element_blank())
  }
  
  return(p)
}

# change y-axis into days (from 0-1)
create_effect_Y_trans <- function(data, line_color = "blue", x_title = "X-axis", y_title = "Y-axis", y_lim = c(0, 1), show_y_axis = TRUE) {
  data$predicted <- (data$predicted * (304 - 60)) + 60  # Reverse transformation for y-axis
  data$conf.low  <- (data$conf.low * (304 - 60)) + 60
  data$conf.high <- (data$conf.high * (304 - 60)) + 60
  
  p<- ggplot(data, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high)) +
    geom_line(color = line_color) +
    geom_ribbon(alpha = 0.1, fill = line_color) +
    labs(x = x_title, y = y_title) +
    ylim(y_lim) +
    theme_minimal(base_size = 10) +
    theme(aspect.ratio = 1, 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "white", colour = "black"))
  
  if (show_y_axis) {
    p <- p + labs(y = y_title)
  } else {
    p <- p + theme(axis.title.y = element_blank(), 
                   axis.text.y = element_blank())
  }
}



# Beetle counts -----------------------------------------------------------


# Assuming 'model' is your glm.nb model
p1 <- ggpredict(fin.m.counts, terms = "veg_tmp [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.counts, terms = "previous_spei3_2 [all]", allow.new.levels = TRUE)

p1.counts <- create_effect_plot(p1, line_color = "red", x_title = "Temperature [z-score]", y_title = "Beetle counts", show_y_axis = TRUE)
p2.counts <- create_effect_plot(p2, line_color = "blue", x_title = "SPEI [z-score]", y_title = "", show_y_axis = FALSE)

# Arrange the plots side by side with the same size
p.effect.counts <- ggarrange(p1.counts, p2.counts, ncol = 2, align = "v")
print(p.effect.counts)





# DOY aggregation ---------------------------------------------------------

# Assuming 'model' is your glm.nb model
p1 <- ggpredict(fin.m.agg, terms = "veg_tmp [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.agg, terms = "previous_spei3_2 [all]", allow.new.levels = TRUE)


p1.agg <- create_effect_Y_trans(p1, line_color = "red", x_title = "Temperature [z-score]", y_title = "Aggregation day", y_lim = c(80,250),show_y_axis = TRUE  )
p2.agg <- create_effect_Y_trans(p2, line_color = "blue", x_title = "SPEI [z-score]", y_title = "Aggregation day" , y_lim = c(80,250), show_y_axis = FALSE)



p.effect.agg <- ggarrange(p1.agg,p2.agg, ncol = 2, align = "v")

(p.effect.agg)

###### PEak Effect plots ------------------------------------------------------------

# Assuming 'model' is your glm.nb model
p1 <- ggpredict(fin.m.peak, terms = "veg_tmp [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.peak, terms = "previous_spei6 [all]", allow.new.levels = TRUE)


p1.peak <- create_effect_Y_trans(p1, line_color = "red", x_title = "Temperature [z-score]", y_title = "Peak day", y_lim = c(130,200), show_y_axis = TRUE )
p2.peak <-create_effect_Y_trans(p2, line_color = "blue", x_title = "SPEI [z-score]", y_title = "Peak day", y_lim = c(130,200), show_y_axis = F)

# effect plots  ----------------------------------------------------
p.effect.peak <- ggarrange(p1.peak,p2.peak, ncol = 2, align =  'v')

p.effect.peak



# effect plot peak dif ----------------------------------------------------
# Assuming 'model' is your glm.nb model
p1 <- ggpredict(fin.m.peak.diff, terms = "spring_tmp [all]", allow.new.levels = TRUE)
p2 <- ggpredict(fin.m.peak.diff, terms = "previous_spei3_2 [all]", allow.new.levels = TRUE)



p1.peak.diff <- create_effect_plot(p1, line_color = "red", x_title = "Spring Temperature [z-score]", y_title = "Peak difference", y_lim = c(220,650), show_y_axis = TRUE)
p2.peak.diff <-create_effect_plot(p2, line_color = "blue", x_title = "SPEI [z-score]", y_title = "Peak difference", y_lim = c(220,650), show_y_axis = FALSE)

# effect plots Sum IPS ----------------------------------------------------
p.effect.peak.diff <- ggarrange(p1.peak.diff,p2.peak.diff, ncol = 2, align = 'v')

p.effect.peak.diff


# all Effect plots --------------------------------------------------------
windows(5,9)
ggarrange(p.effect.counts, p.effect.agg, p.effect.peak, p.effect.peak.diff, nrow=4 , align = 'hv')




