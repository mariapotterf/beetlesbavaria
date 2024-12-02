


# plot labels: 
lab_popul_level       = "Population level [#]"
lab_colonization_time = "Aggregation timing [DOY]"
lab_peak_time         = "Peak sw. timing [DOY]"
lab_peak_growth       = "Peak sw. intensity [#]"


# extract p value from the model and keep it as label

format_p_value_label <- function(model_summary, term) {
  # Extract p-value based on whether the term is parametric or smooth
  if (term %in% rownames(model_summary$p.table)) {
    # Parametric term
    p_value <- model_summary$p.table[term, "Pr(>|t|)"]
  } else if (term %in% rownames(model_summary$s.table)) {
    # Smooth term
    p_value <- model_summary$s.table[term, "p-value"]
  } else {
    # Term not found
    return("Not found")
  }
  
  # Format the p-value
  if (p_value < 0.001) {
    formatted_p <- "< 0.001"
  } else {
    formatted_p <- sprintf("%.3f", p_value)  # Round to 3 decimal places
  }
  
  # Return formatted p-value
  return(formatted_p)
}





# Create effects plots and scatter points below

my_theme_square <- function() {
  theme_minimal(base_size = 8) +
    theme(
      #aspect.ratio = 1, 
      axis.ticks.y = element_line(),
      axis.ticks.x = element_line(),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.background = element_rect(fill = "white", colour = "black"),
      legend.position = "bottom",
      axis.title.y = element_blank()
    ) 
}






# Create effect plot function with an additional argument to control y-axis labels
create_effect_previous_year <- function(data, avg_data, line_color = "grey40", 
                                        x_col = "sum_ips_lag1", 
                                        y_col = "sum_ips", 
                                        x_title = "X-axis", 
                                        y_title = "Y-axis", 
                                        my_title = '',
                                        x_annotate = 0, lab_annotate = "lab ann") {
  x_col <- ensym(x_col)
  y_col <- ensym(y_col)
  
  p <- ggplot() + 
    geom_jitter(data = avg_data, aes(x = !!x_col, y = !!y_col, group = pairID), 
                size = 1, alpha = 0.2, color = 'grey') +
    
    geom_ribbon(data = data, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.25, fill = line_color) +
    geom_line(data = data, aes(x = x, y = predicted), color = line_color,linewidth = 1) +
    
    labs(x = x_title,
         title = my_title,
         y = y_title) +
    my_theme_square() + 
    annotate("text", x = x_annotate, y = Inf, label = lab_annotate, hjust = 0.5, vjust = 1.5)
  
  return(p)
}


# Create effect plot function with additional arguments to select columns from avg_data
create_effect_plot <- function(data, 
                               avg_data, 
                               x_col = "tmp_z_lag1", 
                               y_col = "sum_ips", 
                               line_color = "blue", 
                               x_title = "X-axis", 
                               y_title = "Y-axis", my_title = '',
                               x_annotate = 0, lab_annotate = "lab ann") {
  
  x_col <- ensym(x_col)
  y_col <- ensym(y_col)
  
  p <- ggplot() +
    geom_jitter(data = avg_data, aes(x = !!x_col, y = !!y_col, group = pairID), 
                size = 1, alpha = 0.2, color = 'grey') +
    geom_line(data = data, aes(x = x, y = predicted), color = line_color,linewidth = 1) +
    geom_ribbon(data = data, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.25, fill = line_color) +
    labs(x = x_title,
         title = my_title,
         y = y_title) +
    my_theme_square() + 
    annotate("text", x = x_annotate, y = Inf, label = lab_annotate, hjust = 0.5, vjust = 1.5)
  
  return(p)
}



# plot interactions
plot_effect_interactions <- function(data, 
                                     avg_data, 
                                     x_col = "tmp_z_lag1", 
                                     y_col = "sum_ips", 
                                     temp_label = 'Temp',
                                     y_title = 'y_title',
                                     x_annotate = 0, 
                                     lab_annotate = "lab ann") {
  #library(ggplot2)
  x_col <- ensym(x_col)
  y_col <- ensym(y_col)
  
  p<- ggplot() +
    geom_jitter(data = avg_data, aes(x = !!x_col, y = !!y_col, group = pairID), 
                size = 1, alpha = 0.2, color = 'grey') +
    geom_ribbon(data=data, aes(x = x,
                               ymin = conf.low, 
                               ymax = conf.high, 
                               fill = group), alpha = 0.25) +
    geom_line(data = data, 
              aes(x = x, y = predicted,
                  color = group), linewidth = 1) +#, linetype = group
    labs(x = temp_label,
         y = y_title,
         fill = "SPEI levels",
         color = "SPEI levels",
         linetype = "SPEI levels") +  # Fixed "y_title" to "y" for correct y-axis label argument
    
    guides(color = guide_legend(nrow = 1), 
           fill = guide_legend(nrow = 1),
           linetype = guide_legend(nrow = 1)) +
    annotate("text", x = x_annotate, y = Inf, 
             label = lab_annotate, hjust = 0.5, vjust = 1.5) +
    my_theme_square() +
    theme(legend.position = "bottom")
  
  return(p)
}


plot_data_with_average <- function(data, y_var, y_label, my_title) {
  # Convert the y_var argument to a symbol to use in aes()
  y_var_sym <- rlang::sym(y_var)
  
  data %>%
    ungroup() %>%
    filter(year %in% 2015:2021) %>%
    ggplot(aes(x = year, y = !!y_var_sym, group = pairID)) +
    labs(x = 'Year', 
         y = y_label, 
         title = my_title) +
    geom_line(alpha = 0.2, linewidth = 0.2, color = 'grey70') +  
    stat_summary(
      aes(x = year, y = !!y_var_sym, group = 1), 
      fun = mean,  # Calculate the mean for each year
      geom = "line", 
      color = "#E69F00",  # Ensure the average line is also red
      linewidth = 1  # Make the average line slightly thicker than individual lines
    ) + my_theme_square()
  
}
