


# analyse buffer beetle data dn ips counts


load(file =  "outData/buffers.Rdata")

# df_fin,            # final table with coordinates, counts and disturbance and present forest per year
# xy_sf_expand,      # all xy trap coordinates (shifted over years)



names(df_fin)


# [1] "falsto_name"     "globalid"        "year"            "id"              "wind_beetle"     "harvest"         "forest_freq"     "all_dist_sum"    "cum_removed"    
#[10] "remained_forest" "sum_ips"         "lag1_dist_sum"   "lag1_harvest"    "lag2_harvest"    "lag1_beetles"    "lag2_beetles"   
  
# LM only for high beetle data
df_fin %>% 
  #filter(sum_ips > 50000) %>% 
  ggplot(aes(x = harvest*0.09,
             y =  sum_ips)) +
    geom_point(alpha = 0.5) +
  geom_smooth() +
 # facet_wrap(year ~.,scales = 'free') +
  #ggtitle('beetles') + 
  theme_bw()


df_fin %>% 
  #filter(sum_ips > 50000) %>% 
  ggplot(aes(x = sum_ips,
             y = lag1_beetles*0.09  )) +
  # stat_poly_line() +
  #stat_poly_eq(use_label(c( "R2", 'p'))) + #"eq",
  geom_point(alpha = 0.5) +
  facet_wrap(year ~.,) +
  ggtitle('Wind+ beetles') + 
  theme_bw()



df_fin %>% 
  #filter(sum_ips > 50000) %>% 
  ggplot(aes(x = sum_ips,
             y = wind*0.09 + harvest*0.09  )) +
  # stat_poly_line() +
  #stat_poly_eq(use_label(c( "R2", 'p'))) + #"eq",
  geom_point(alpha = 0.5) +
  facet_wrap(year ~.,) +
  ggtitle('Harvest + Wind+ beetles') + 
  theme_bw()




df_fin %>% 
  ggplot(aes(harvest*0.09))+ 
  geom_histogram(binwidth = 0.1) + 
  facet_wrap(year~.)

df_fin %>% 
  ggplot(aes(sum_ips))+ 
  geom_histogram() + 
  facet_wrap(year~.)


library(ggpmisc)
df_fin %>% 
  ggplot(aes(x = sum_ips/10000,
             y = all_dist_sum )) +
  stat_poly_line() +
  stat_poly_eq(use_label(c( "R2", 'p'))) + #"eq",
  geom_point(alpha = 0.5) +
  facet_wrap(year ~. ) +
  ggtitle('All disturbances')



df_fin %>% 
  ggplot(aes(x = sum_ips/10000,
             y = harvest*0.09  )) +
  # stat_poly_line() +
  #stat_poly_eq(use_label(c( "R2", 'p'))) + #"eq",
  geom_point(alpha = 0.5) +
  facet_wrap(year ~.,) +
  ggtitle('Harvest') + 
  theme_bw()



# need to lag values by one-two year! eg high number of beetles in year n
# triggers high mortality in next one to two years
# export all disturbances = sum, maybe that will work?


df_fin %>%
  as.data.frame() %>% 
  ggplot(aes(x = sum_ips/10000,
             y = lag1_dist_sum  )) +
  stat_poly_line() +
  stat_poly_eq(use_label(c( "R2", 'p'))) + #"eq",
  geom_point(alpha = 0.5) +
  # geom_smooth('lm') +
  facet_wrap(year ~., scales = 'free') +
  ggtitle('Lag1 all disturb')


# buffer area = 780 ha
df_fin %>% 
  ggplot(aes(x = sum_ips/10000,
             y = lag1_harvest*0.09  )) +
  #stat_poly_line() +
  # stat_poly_eq(use_label(c( "R2", 'p'))) + #"eq",
  geom_point(alpha = 0.5) +
  facet_wrap(year ~.) + # , scales = 'free'
  ggtitle('Harvest +1year') + 
  theme_bw()



df_fin %>% 
  ggplot(aes(x = sum_ips,
             y = lag2_harvest  )) +
  stat_poly_line() +
  stat_poly_eq(use_label(c( "R2", 'p'))) + #"eq",
  geom_point(alpha = 0.5) +
  # geom_smooth('lm') +
  facet_wrap(year ~., scales = 'free') +
  ggtitle('Lag2 Harvest')






