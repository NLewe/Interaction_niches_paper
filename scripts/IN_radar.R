### Radar of species 


library (ggradar)

library (scales)

All_metrics_M1_8Sp <-  read_rds("data/ALL_metrics_M1_8Sp_df.rds") %>% 
  as_tibble () %>%  add_column (Exp = "E2")

## radar plots per experiment ####
plot_radar_E1<- 
  All_metrics_M0_8Sp %>%  
  select (-PlantFamily) %>% 
  mutate_at (vars (-PlantSpeciesfull, -Exp), rescale) %>% 
  select (-Exp) %>% 
  ggradar(legend.title = "Plant species", plot.title = "Experiment 1")

plot_radar_E2 <- 
All_metrics_M1_8Sp %>%  
  mutate_at (vars (-PlantSpeciesfull,-PlantFamily,  -Exp), rescale) %>% 
  select (-PlantFamily, -Exp) %>% 
  relocate (PlantSpeciesfull) %>% 
  ggradar(legend.title = "Plant species", plot.title = "Experiment 2") 


ggarrange (plot_radar_E1, plot_radar_E2, common.legend = T, legend = "bottom", label = "AUTO")

# radar plot for each plant species ######
radar_plots <- 
  All_metrics_M0_8Sp %>%  
  mutate_at (vars (-PlantSpeciesfull, -Exp, -PlantFamily), rescale)  %>% 
  bind_rows (All_metrics_M1_8Sp %>%     
               mutate_at (vars (-PlantSpeciesfull,-PlantFamily,  -Exp), rescale)) %>%  
  select (-PlantFamily) %>% 
  group_split(PlantSpeciesfull) %>% 
  map (~select (., -PlantSpeciesfull) %>% relocate (Exp)) %>% 
  map (~ggradar(.,)) 


ggarrange (radar_plots[[1]], radar_plots[[2]], radar_plots[[3]], radar_plots[[4]], 
           radar_plots[[5]], radar_plots[[6]], radar_plots[[7]], radar_plots[[8]], 
           labels = c("A. capillaris", "A. millefolium", "B. willdenowii", "C. intybus", 
                      "H. lanatus", "P. cita", "P. lanceolata", "S. arundinaceus"),
           common.legend = T, nrow = 2, ncol=4)


### 
LMER5B_DW <- lmer (SqrtDW_abovePlaSpe ~  SqrtDW_above_focal + (1|focal.sp)  + Observed.y + (1|GSPlaSpe), data = BC_richness_LMdata) 

BC_richness_LMdata %>%  ggplot (aes (x = Observed.y, y = Observed.x, color = PlantSpeciesfull)) + geom_point () + geom_smooth(method = "lm")





# Interaction properties experiment 1 ####
Metrics_exp1  <-
adiv_richness %>% select (sampleID, Observed, Shannon, Chao1) %>% 
  left_join(metaM0 %>%  select (sampleID, PlantSpeciesfull))  %>%  
  group_by(PlantSpeciesfull)  %>% 
  left_join (stand_pd_Glo_all %>% select (sampleID, pd.obs)) %>% 
  left_join(ses.MPD_Glo %>%  select (sampleID, mpd.obs))  %>%  
  drop_na () %>% 
  add_column (exp = "Exp1")
  
  





# interaction properties experiment 2 ####
adiv_richness_M1 %>% select (sampleID, Observed, Shannon, Chao) %>% 


  
## Pro 
  PCA_all_metrics <- PCA(ALL_metricsCh2_8, quali.sup = 9,   scale.unit = T, graph = F)




