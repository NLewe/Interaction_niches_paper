# clusterin ###

library (factoextra)
library (FactoMineR)

library (ggrepel)


All_metrics_gen_df <-   All_metrics_E1 %>%  
  filter (PlantFamily != "Asteraceae") %>%  
  filter (PlantFamily != "Fabaceae")  %>% 
  filter (PlantFamily != "Cyperaceae") %>% 
  #select (!Richness) %>% 
  data.frame (row.names = "PlantSpeciesfull")

df <- scale (All_metrics_gen_df [, 1:7])
# Correlation-based distance method
res.dist <- get_dist(df, method = "euclidean")


# Compute hierarchical clustering
res.hc <- hclust(res.dist, method = "ward.D2")
# Visualize
plot(res.hc, cex = 0.5)


# Enhanced k-means clustering
res.km <- eclust(df, "kmeans", k=5)

fviz_gap_stat(res.km$gap_stat)
# Silhouette plot
fviz_silhouette(res.km)



# Compute hierarchical k-means clustering
res.hk <-hkmeans(df, 5)
# Elements returned by hkmeans()
names(res.hk)

# Visualize the hkmeans final clusters
fviz_cluster(res.hk, palette = "jco", repel = TRUE,
             ggtheme = theme_classic())


# PCA #####
#PCA_all_metrics  <- vegan::rda (All_metrics_generalists, scale = T)
PCA_all_metrics_generalists <- PCA(All_metrics_gen_df[,-c(1,3)], quali.sup = c(7,6),   scale.unit = T, graph = F)


# Plot PCA #### 

PCA_arrows_metrics_gen<-
  PCA_all_metrics_generalists$var$coord %>%  
  as_tibble (rownames = "metric") %>% 
  dplyr::rename("D1end" = "Dim.1", "D2end"= "Dim.2", "D3end"="Dim.3")

PCA_metric_data_gen   <- PCA_all_metrics_generalists$ind$coord  %>%  
  as_tibble(rownames = "PlantSpeciesfull") %>% 
  left_join(metaM0 %>%  select (PlantSpeciesfull, PlantFamily) %>%  unique ())

PCA_eig_m_Dim1_gen <- round (PCA_all_metrics_generalists$eig[1,2],1)
PCA_eig_m_Dim2_gen <- round (PCA_all_metrics_generalists$eig[2,2],1)
PCA_eig_m_Dim3_gen <- round (PCA_all_metrics_generalists$eig[3,2],1)

plotPca_gen1  <- 
  PCA_metric_data_gen %>% 
  mutate (PlantFamily = str_replace(PlantFamily, "Soil","Soil control")) %>% 
  ggplot (aes (x = Dim.1, y = Dim.2, color = PlantFamily)) + 
  geom_point (size =3) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  #stat_ellipse(aes (x= Dim.1, y = Dim.2, color = group))  +   # thisis 95%confidence
  geom_segment(data = PCA_arrows_metrics_gen, aes (x=0, xend= D1end*1, y = 0, yend = D2end*1), 
               arrow = arrow(length = unit(0.3, "picas")), color = "darkblue", inherit.aes = F)  +
  geom_text_repel ( data = PCA_arrows_metrics_gen, aes (x  = D1end*1, y = D2end*1, label = metric), 
                    color = "darkblue", inherit.aes = F , force = 0.6) + 
  geom_text_repel(data = PCA_metric_data_gen, aes (x = Dim.1, y = Dim.2, label = PlantSpeciesfull), fontface = "italic", inherit.aes = F) +
  theme_minimal() + 
  xlab(label = "PC1 (61.5%)") +
  ylab ("PC2 (22.7 %)") + 
  guides (color= guide_legend( "Plant family")) +
  theme (legend.position = "bottom") +
  scale_color_manual(values=c(
    "Plantaginaceae"= "#4e79a7","Soil control"= "#9c755f",
    "Poaceae" = "#7F9A65" )) 

plotPca_gen2  <- 
  PCA_metric_data_gen %>%
  mutate (PlantFamily = str_replace(PlantFamily, "Soil","Soil control")) %>% 
  ggplot (aes (x = Dim.2, y = Dim.3, color = PlantFamily)) + 
  geom_point (size =3) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  #stat_ellipse(aes (x= Dim.1, y = Dim.2, color = group))  +   # thisis 95%confidence
  geom_segment(data = PCA_arrows_metrics_gen, aes (x=0, xend= D2end*1, y = 0, yend = D3end*1), 
               arrow = arrow(length = unit(0.3, "picas")), color = "darkblue", inherit.aes = F)  +
  geom_text_repel ( data = PCA_arrows_metrics_gen, aes (x  = D2end*1, y = D3end*1, label = metric), 
                    color = "darkblue", inherit.aes = F , force = 0.6) + 
  geom_text_repel(data = PCA_metric_data_gen, aes (x = Dim.2, y = Dim.3, label = PlantSpeciesfull), fontface = "italic", inherit.aes = F) +
  theme_minimal() + 
  xlab(label = "PC1 (21.1 %)") +
  ylab ("PC2 (5.9 %)") + 
  guides (color= guide_legend( "Plant family")) +
  theme (legend.position = "bottom") +
  scale_color_manual(values=c(
    "Plantaginaceae"= "#4e79a7","Soil control"= "#9c755f",
    "Poaceae" = "#7F9A65" )) 


# plotPca  <- fviz_pca_biplot(PCA_all_metrics,axes = c(1,2), col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
#                 repel =T)
# plotPca2  <- fviz_pca_biplot(PCA_all_metrics,axes = c(2,3), col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
#                 repel =T)

plot1 <- fviz_contrib(PCA_all_metrics_generalists, choice = "var", axes = 1,
                      font.main = c(size =12 ), title= "Contribution of variables to PC1")
plot2 <- fviz_contrib(PCA_all_metrics_generalists, choice = "var", axes = 2, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC2")
plot3 <- fviz_contrib(PCA_all_metrics_generalists, choice = "var", axes = 3, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC3")

ggarrange (ggarrange (plotPca_gen1,plotPca_gen2, ncol = 2, nrow =1, labels= c("A", "B"), common.legend = T, legend = "bottom"), 
           ggarrange (plot1, plot2, plot3, nrow = 1, ncol= 3,
                      widths = 2, labels= c( "C", "D", "E")), 
           nrow=2, ncol=1, heights= c(2,1) )