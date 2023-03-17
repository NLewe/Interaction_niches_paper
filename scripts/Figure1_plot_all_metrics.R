### test scrpt ##


adiv_richness %>% 
  left_join(metaM0)  %>%  
  group_by(PlantSpeciesfull)  %>% 
  #mutate (Chao1 = mean (Chao1)) %>% 
  mutate (Shannon = mean (Shannon)) %>%  
  select (sampleID, Shannon) %>% ## chao1
  select (!sampleID) %>% 
  unique()  %>%  
  left_join (comp_units %>% select (!repl) %>% select (!"1-CU")) %>%   #includes CU and unique and richness
  left_join(cores) %>% 
  select (!c(PlantType,uniqueASVsperPlSpe, PlantFamily))  %>% 
  left_join (stand_pd_Glo_all %>% 
               left_join(ses.MPD_Glo, by ="sampleID")  %>%  
               left_join (metaM0)  %>% 
               drop_na ()  %>% 
               group_by (PlantSpeciesfull)  %>%   
               mutate (mean_pd = mean (pd.obs), mean_mpd = mean (mpd.obs)) %>%  
               select (PlantSpeciesfull,   mean_pd, mean_mpd,PlantFamily)   %>% 
               unique ()  )  %>% 
  dplyr::rename("Richness"= "mean_n_ASV_per_species", "CU" = "CUnits", "β(core)"  = "perc_core",
                "PD" = "mean_pd", "MPD" = "mean_mpd" ) %>% 
  add_column (Ex = "E1")

All_metrics_M1  <-
  adiv_richness_M1 %>% 
  group_by(PlantFamily, PlantSpeciesfull)  %>% 
  #mutate (Chao1 = mean (Chao1)) %>% 
  mutate (Shannon = mean (Shannon)) %>%  
  select (sampleID, Shannon) %>% 
  select (!sampleID) %>% 
  unique()  %>%  
  left_join (comp_units %>%  select (!'1-CU'))    %>%   #includes CU and unique and richness
  left_join(cores_roots_M1) %>% 
  left_join (PD_MPD_M1) %>% 
  dplyr::rename("Richness"= "meanASV", "CU" = "CUnits", "β(core)"  = "perc_core",
                "PD" = "mean_pd", "MPD" = "mean_mpd" , "uniqueASV" = "uniqueASVsperPlSpe") %>% 
  add_column (Exp = "E2")



All_metrics %>%  add_column (Exp = "E1")
  
  mutate_at(c(1:8),scale)  %>%    # z score standardisation
  group_split (focal.sp) %>%  # split into each plant species
  map (~ dist (., method = "euclidean"),.id= focal.sp) %>%  ## calculate distances
  map (  ~melt (as.matrix(.), varnames = c("row", "col"))) %>%   # put together into tibble again
  map (~filter (., col =="1" & row != "1")) %>%   ## remove unimportant distances
  map (~add_column(., TM = c("AC", "BC"))) %>% 
  map_dfr (~select (.,value, TM)) %>%  
  add_column (focal.sp= c ("G1", "G1", "G2", "G2",
                           "G3", "G3", "G4", "G4", 
                           "S1", "S1", "S2", "S2", 
                           "S3", "S3", "S4", "S4")
  ) %>% 
  mutate (treatment =str_c ( TM, focal.sp, sep = "-")) %>%  
  left_join (meta_plants) %>% 
  mutate (x= str_c (TM, PlantSpeciesfull, sep =" ")) %>% 
  arrange (M1_M2_order) %>% 
  add_column (dist_order = c(1:16))
  
  
  # all Exp1 plots barplot together ######
  ### save all 450 x 800 

  #Plot richness ####
  
  p1  <-   adiv_richness %>% 
    left_join(metaM0, by="sampleID") %>% 
    filter (Observed!=0) %>% 
    group_by (PlantSpeciesfull) %>% 
    mutate(meannumberASVs=mean(Observed)) %>% 
    #select (sampleID,Observed, Shannon, PlantSpeciesfull, meannumberASVs, PlantFamily)  %>% 
    #pivot_longer(cols=c( Observed, Shannon), names_to ="diversity_measure", values_to = "div", values_drop_na = TRUE) %>% 
    arrange(meannumberASVs)   %>% 
    # mutate (diversity_measure = str_replace(diversity_measure, "Observed","Richness S")) %>% 
    #   mutate (diversity_measure = str_replace(diversity_measure, "Shannon","Shannon's H'")) %>% 
    mutate (PlantFamily = str_replace(PlantFamily, "Soil","Soil control")) %>% 
    ggplot(aes(y=reorder(PlantSpeciesfull, -meannumberASVs), x=Observed, fill=PlantFamily))  +
    geom_boxplot()  #+
   ## facet_wrap(vars(diversity_measure), scales = "free_x")
  
  plot1_richness  <- 
    p1 +  
    stat_summary (fun=mean, geom="point", shape = 23)  +
    xlab("Richness") +
    ylab ("Plant species and soil control") +
    theme_bw () +
    theme (axis.text.y= element_text(family= "sans", face= "italic", size = 11)) +
    theme (axis.text.x=element_text(family= "sans", size = 11), legend.position = "bottom")  +
    theme (axis.title = element_text(family = "sans", size = 12 ),  
           legend.title=element_text(size=12), 
           legend.text=element_text(size=11)) +
    scale_fill_manual(values=c("Asteraceae"= "#edc948" ,"Cyperaceae" ="#5A6351" ,
                               "Fabaceae" = "#f28e2b" , 
                               "Plantaginaceae"= "#4e79a7","Soil control"= "#9c755f",
                               "Poaceae" = "#7F9A65" ), name = "Plant family and soil control") 
  
 
  #Plot beta diversity among replicates
  #β-diversity among AMF communities associated with individual plants.species 
  #A) numeric β-diversity illustrating the relative compositional heterogeneity of the AMF communities associated with the plant species. 
  
  # comp_units   <- meannumberASVsper_species  %>% 
  #   left_join(uniqueASVsperPlaSpe)  %>% 
  #   left_join (replicates_samples_after_Glo) %>% 
  #   group_by (PlantSpeciesfull) %>% 
  #   mutate (unique = mean(uniqueASVsperPlSpe)) %>% 
  #   mutate (CUnits = 5* unique/(mean_n_ASV_per_species*repl)) %>% # repl adjusted for , times 5 to get back to number of 5 replicates 
  #   mutate ("1-CU" = 1 - CUnits)
  # 
  # 
# Plot2 Shannon ####  
  p1  <-   adiv_richness %>% 
    left_join(metaM0, by="sampleID") %>% 
    filter (Observed!=0) %>% 
    group_by (PlantSpeciesfull) %>% 
    mutate(meannumberASVs=mean(Observed)) %>% 
    #select (sampleID,Observed, Shannon, PlantSpeciesfull, meannumberASVs, PlantFamily)  %>% 
    #pivot_longer(cols=c( Observed, Shannon), names_to ="diversity_measure", values_to = "div", values_drop_na = TRUE) %>% 
    arrange(meannumberASVs)   %>% 
    # mutate (diversity_measure = str_replace(diversity_measure, "Observed","Richness S")) %>% 
    #   mutate (diversity_measure = str_replace(diversity_measure, "Shannon","Shannon's H'")) %>% 
    mutate (PlantFamily = str_replace(PlantFamily, "Soil","Soil control")) %>% 
    ggplot(aes(y=reorder(PlantSpeciesfull, -meannumberASVs), x=Shannon, fill=PlantFamily))  +
    geom_boxplot()  #+
  ## facet_wrap(vars(diversity_measure), scales = "free_x")

  plot2_Shannon  <- 
    p1 +  
    stat_summary (fun=mean, geom="point", shape = 23)  +
    xlab("Shannon's H") +
   # ylab ("Plant species and soil control") +
    theme_bw () +
    theme(axis.text.y= element_blank()) +
    theme (axis.text.x=element_text(family= "sans", size = 11), legend.position = "bottom")  +
    theme (axis.title.x  = element_text(family = "sans", size = 12 ),  
           legend.title=element_text(size=12), 
           legend.text=element_text(size=11), 
           axis.title.y = element_blank(), 
           axis.ticks = element_blank()) +
    scale_fill_manual(values=c("Asteraceae"= "#edc948" ,"Cyperaceae" ="#5A6351" ,
                               "Fabaceae" = "#f28e2b" , 
                               "Plantaginaceae"= "#4e79a7","Soil control"= "#9c755f",
                               "Poaceae" = "#7F9A65" ), name = "Plant family and soil control") 
# Plot 3 PD ####  
  p1  <-
    stand_pd_Glo_all  %>% 
    left_join (metaM0, by= "sampleID")  %>% 
    #filter (!is.na(pd.obs)) %>% 
    group_by (PlaSpe) %>% 
    mutate(mean_pd.obs =mean(pd.obs)) %>% 
    arrange(mean_pd.obs)   %>% 
    left_join ((meannumberASVsper_species %>% ungroup () %>% select (-PlantFamily, PlantSpeciesfull, mean_n_ASV_per_species)), by = "PlantSpeciesfull")  %>% 
    mutate (PlantFamily = str_replace(PlantFamily, "Soil","Soil control")) %>% 
    ggplot (aes(x=reorder (PlantSpeciesfull,-mean_n_ASV_per_species),
                y=pd.obs, fill=PlantFamily))  +
    geom_boxplot() #  +   geom_point (data = (stand_pd_Glo_PlaSpe %>%  left_join(metaM0, by= c("sampleID"= "PlantSpeciesfull"))), aes (x= sampleID, y =pd.obs ), shape = 15)
  #that are single points to see 
  plot3_PD   <- p1 +coord_flip () +
    stat_summary (fun=mean, geom="point", shape = 23)  +
    ylab("Faith's phylogenetic diversity") +
    theme_bw () +
    theme(axis.text.y= element_blank()) +
    theme (axis.text.x=element_text(family= "sans", size = 11))  +
    theme (axis.title.x = element_text(family = "sans", size = 12 ),  
           legend.title=element_text(size=12), 
           legend.text=element_text(size=11), 
           axis.title.y = element_blank(), 
           legend.position = "bottom",
           axis.ticks = element_blank()) +
    scale_fill_manual(values=c("Asteraceae"= "#edc948" ,"Cyperaceae" ="#5A6351" ,
                               "Fabaceae" = "#f28e2b" , 
                               "Plantaginaceae"= "#4e79a7","Soil control"= "#9c755f",
                               "Poaceae" = "#7F9A65" ), name = "Plant family and soil control") 
  #Plot Faiths PD-mean effect size####
  p1  <-  stand_pd_Glo_all %>% 
    left_join (metaM0, by= "sampleID")  %>% 
   # filter (!is.na(pd.obs.z)) %>% 
    group_by (PlaSpe) %>% 
    mutate(mean_pd.obs.z =mean(pd.obs.z)) %>% 
    #arrange(mean_pd.obs.z)   %>% 
    left_join ((meannumberASVsper_species %>% ungroup () %>% 
                  select (-PlantFamily, PlantSpeciesfull, mean_n_ASV_per_species)), by = "PlantSpeciesfull")  %>% #to sort like in richness plot
    mutate (PlantFamily = str_replace(PlantFamily, "Soil","Soil control")) %>% 
    ggplot (aes(x=reorder (PlantSpeciesfull,-mean_n_ASV_per_species), 
                y=pd.obs.z, fill=PlantFamily))  +  
    geom_boxplot() +
    geom_point (data = (stand_pd_Glo_PlaSpe %>% 
                          left_join (metaM0, by= c("sampleID" = "PlantSpeciesfull")) %>% 
                          filter (!between (pd.obs.p, 0.05, 0.95))%>% 
                          select (sampleID, PlantFamily, pd.obs.z)  %>% unique ()),
                aes(x=sampleID, y= pd.obs.z), shape = 15, size =3) 
  
  
  plot4_PDEff <- p1  +
    coord_flip () +
    stat_summary (fun=mean, geom="point", shape = 23)  +
    ylab("Faith's PD effect size (z-score)") +
    theme_bw () +
    theme(axis.text.y= element_blank()) +
    theme (axis.text.x =element_text(family= "sans", size = 11))  +
    theme (axis.title.x = element_text(family = "sans", size = 12 ),  
           legend.title=element_text(size=12), 
           legend.text=element_text(size=11), 
           axis.title.y = element_blank(), 
           legend.position = "bottom", 
           axis.ticks = element_blank()) +
    scale_fill_manual(values=c("Asteraceae"= "#edc948" ,"Cyperaceae" ="#5A6351" ,
                               "Fabaceae" = "#f28e2b" , 
                               "Plantaginaceae"= "#4e79a7","Soil control"= "#9c755f",
                               "Poaceae" = "#7F9A65" ), name = "Plant family and soil control") 
 

  # Plot MPD ####

  p1   <- 
    ses.MPD_Glo  %>% 
    left_join (metaM0, by= "sampleID")  %>% 
   # filter (!is.na(mpd.obs)) %>% 
    group_by (PlaSpe) %>% 
    mutate(mean_mpd.obs =mean(mpd.obs)) %>% 
    arrange(mean_mpd.obs)   %>% 
    left_join ((meannumberASVsper_species %>% ungroup () %>%
                  select (-PlantFamily, PlantSpeciesfull, mean_n_ASV_per_species)), by = "PlantSpeciesfull")  %>% 
    mutate (PlantFamily = str_replace(PlantFamily, "Soil","Soil control")) %>% 
    ggplot (aes(x=reorder (PlantSpeciesfull,-mean_n_ASV_per_species),
                y=mpd.obs, fill=PlantFamily))  +
    geom_boxplot() 
  #+geom_point (data = (ses.MPD_Glo_PlaSpe %>% left_join (metaM0, by= c("sampleID" = "PlantSpeciesfull"))), aes(x=sampleID, y= mpd.obs), shape = 15)
  #small numbers are generalists here!
  
  plot5_mpd   <- p1 +coord_flip () +
    stat_summary (fun=mean, geom="point", shape = 23)  +
    ylab("Mean phylogenetic distance (MPD)") +
  theme_bw () +
    theme(axis.text.y= element_blank()) +
    theme (axis.text.x =element_text(family= "sans", size = 11))  +
    theme (axis.title.x = element_text(family = "sans", size = 12 ),  
           legend.title=element_text(size=12), 
           legend.text=element_text(size=11), 
           axis.title.y = element_blank(), 
           legend.position = "bottom", 
           axis.ticks = element_blank()) +
    scale_fill_manual(values=c("Asteraceae"= "#edc948" ,"Cyperaceae" ="#5A6351" ,
                               "Fabaceae" = "#f28e2b" , 
                               "Plantaginaceae"= "#4e79a7","Soil control"= "#9c755f",
                               "Poaceae" = "#7F9A65" ), name = "Plant family and soil control") +
    scale_x_discrete(labels = NULL)

  
  # Plot MPD - effect size  ####
  p2   <- 
    ses.MPD_Glo  %>% 
    left_join (metaM0, by= "sampleID")  %>% 
   # filter (!is.na(mpd.obs)) %>% 
    group_by (PlaSpe) %>% 
    mutate(mean_mpd.obs.z =mean(mpd.obs.z)) %>% 
    arrange(mean_mpd.obs.z)   %>% 
    left_join ((meannumberASVsper_species %>% ungroup () %>%
                  select (-PlantFamily, PlantSpeciesfull, mean_n_ASV_per_species)), by = "PlantSpeciesfull")  %>% 
    mutate (PlantFamily = str_replace(PlantFamily, "Soil","Soil control")) %>% 
    ggplot (aes(x=reorder (PlantSpeciesfull,-mean_n_ASV_per_species),
                y=mpd.obs.z, fill=PlantFamily))  +
    geom_boxplot()  +
    geom_point (data = (ses.MPD_Glo_PlaSpe %>% 
                          left_join (metaM0, by= c("sampleID" = "PlantSpeciesfull")) %>% 
                          filter (! between(mpd.obs.p, 0.05, 0.95))), 
                aes(x=sampleID, y= mpd.obs.z), shape = 15, size =3)
  
  #small numbers are generalists here!
  
  plot_mpdeff  <- p2 +coord_flip () +
    stat_summary (fun=mean, geom="point", shape = 23)  +
    ylab("MPD effect size (z-score)") +
    # xlab ("Plant Species") +
    theme_bw () +
    theme(axis.text.y= element_blank()) +
    theme (axis.text.x =element_text(family= "sans", size = 11))  +
    theme (axis.title.x = element_text(family = "sans", size = 12 ),  
           legend.title=element_text(size=12), 
           legend.text=element_text(size=11), 
           axis.title.y = element_blank(), 
           legend.position = "bottom", 
           axis.ticks = element_blank()) +
    scale_fill_manual(values=c("Asteraceae"= "#edc948" ,"Cyperaceae" ="#5A6351" ,
                               "Fabaceae" = "#f28e2b" , 
                               "Plantaginaceae"= "#4e79a7","Soil control"= "#9c755f",
                               "Poaceae" = "#7F9A65" ), name = "Plant family and soil control") +
    scale_x_discrete(labels = NULL)
  
  
  
  
  # Plot 6
  # Plot unique  ####
  p1  <-  ASV_table_Glo  %>%  
    select(-c(Kingdom, Phylum, Class, Order, Species))  %>% 
    gather (-c(ASV_ID, Genus, Family), key=sampleID, value = ASV_counts)  %>%
    left_join(metaM0, by="sampleID") %>% 
    filter (ASV_counts!=0)  %>% 
    group_by (ASV_ID, PlantSpeciesfull, Genus) %>% 
    tally () %>%
    group_by(Genus, PlantSpeciesfull )  %>% 
    tally () %>% 
    replace_na(list (Genus="unidentified"))   %>% 
    left_join (meannumberASVsper_species, by = "PlantSpeciesfull")  %>% 
    ggplot (aes(x=reorder (PlantSpeciesfull,-mean_n_ASV_per_species), y= n, fill=Genus))  +
    geom_bar(position="stack", stat ="identity") 
  
  
  plot_unique  <-
    p1+ 
    coord_flip () +
    ylab("Unique AMF ASVs per plant species") +
    theme_bw () +
    theme(axis.text.y= element_blank()) +
    theme (axis.text.x=element_text(family= "sans", size = 11), legend.position = "bottom")  +
    theme (axis.title.x  = element_text(family = "sans", size = 12 ),  
           legend.title=element_text(size=12), 
           legend.text=element_text(size=11), 
           axis.title.y = element_blank(), 
           axis.ticks = element_blank()) +
    scale_fill_manual(values = c( "unidentified" ="#C7EAE5", "Archaeospora"  = "#287D8EFF" , 
                                  "Acaulospora" =  "#C1A363" , 
                                  "Funneliformis" = "#F6E8C3", "Cetraspora"  = "#DFC27D", 
                                  "Claroideoglomus" =  "#20A386FF", 
                                  "Glomus" = "#80CDC1" , "Paraglomus" = "#35978F", 
                                  "Scutellospora"= "#01665E" ), name = "AMF genus"    )  
  
  
  
  
  
  

 ggarrange (plot1_richness, plot2_Shannon, plot3_PD,  plot4_PDEff,  
            plot5_mpd,plot_mpdeff, plot_unique, widths =  c(2,1,1,1,1,1,1),
            common.legend = T, nrow =1)
##height 777

All_metrics_df <-   All_metrics %>%  
  filter (PlantFamily != "Asteraceae") %>%  
  filter (PlantFamily != "Fabaceae")  %>% 
  filter (PlantFamily != "Cyperaceae") %>% 
  #select (!Richness) %>% 
  data.frame (row.names = "PlantSpeciesfull")

df <- scale (All_metrics_df [, 1:7])
  # Correlation-based distance method
res.dist <- get_dist(df, method = "euclidean")


# Compute hierarchical clustering
res.hc <- hclust(res.dist, method = "ward.D2")
# Visualize
plot(res.hc, cex = 0.5)

library("factoextra")

# Enhanced k-means clustering
res.km <- eclust(df, "kmeans", k=5)

fviz_gap_stat(res.km$gap_stat)
# Silhouette plot
fviz_silhouette(res.km)
  