


ALL_m_scaled <- scale.default(All_metrics[-20,1:7])


maxs <- apply(All_metrics[-20,1:7], 2, max)
mins <- apply(All_metrics[-20,1:7], 2, min)
scale(All_metrics[,1:7], center = mins, scale = maxs - mins)

ALL_m_scaled <- scale(All_metrics[-20,1:7], center=(maxs+mins)/2, scale=(maxs-mins)/2)


ALL_m_scaled %>%  as_tibble (rownames = "PlantSpeciesfull") %>%  
  mutate (`ß(core)` = -1* `ß(core)`) %>% 
  pivot_longer ( !PlantSpeciesfull, names_to = "metric", values_to = "scaled_m") %>% 
  left_join(meannumberASVsper_species)  %>% 
    ggplot (aes (x= metric, y = reorder (PlantSpeciesfull, -mean_n_ASV_per_species), fill = scaled_m)) +
  geom_tile (color = "black") +
  scale_fill_gradient2(low="red", high="green", mid = "white")
  

