# ## Examiners report 
# 
# 
# # Tree  chapter 2#####
# 
# 
# # #get files#### ---- later, when rarefy needed, we have to use ALL sequences... and then select the Glomeromycota
# # #tree####
# #Tree for Glo sequences only#### tree for full dataset is in dada2 script
# #Fasta file with all sequences
# ASV_fasta  <- seqinr::read.fasta(file = "resultsCh2/ASVs.fa", seqtype = "DNA", as.string = TRUE)
# #Subset fasta file for Glo sequences
# ASV_fasta_Glo  <- ASV_fasta[names(ASV_fasta) %in% ASV_table_Glo$ASV_ID]
# seqinr::write.fasta(ASV_fasta_Glo, names = names(ASV_fasta_Glo), file.out= "resultsCh2/ASV_fasta_Glo.fa")
# #Import fasta file as dna sequences for later alignment
# dna_Glo   <- readDNAStringSet("resultsCh2/ASV_fasta_Glo.fa")
# #Alignment with ClustalW###
# ASVs_align_Glo  <-msa(dna_Glo, "ClustalW")
# #Convert for different package, here seqinr
# ASVs_align_Glo_forseqinr  <- msaConvert(ASVs_align_Glo, type= "seqinr::alignment")
# 
# #Build distance matrix as prep for tree###
# d_Glo <- seqinr::dist.alignment(ASVs_align_Glo_forseqinr, "identity")  #seqinr package
# #phylogenetic tree
# #neighbor joining
# ASV_tree_Glo  <- nj(d_Glo)  #ape package
# 
# # #check if nJ was appropriate
# # # x  <- as.vector(d_Glo)
# # # y  <- as.vector(as.dist(cophenetic(ASV_tree_Glo)))
# # #
# # # plot (x, y, xlab = "original distance", ylab = "distance in the tree",
# # #       main = "Is NJ appropriate?", pch = 20, col = transp("black", 0.1), cex = 3)
# # # abline(lm(y ~ x))
# # # cor(x, y)^2
# # ##when points on line, then fine. looked good enough
# #
# 
# 
# 
# ## Bootstrapping 
# # https://fuzzyatelin.github.io/bioanth-stats/module-24/module-24.html
# dna <- fasta2DNAbin (file = "resultsCh2/ASV_fasta_Glo.fa")
# myBoots <- boot.phylo(treeNJ, dna, function(e) njs(dist.dna(e)))

#
#
# ##Matrix =otu-table herstellen für das Package####
# matrix  <-  ASV_table_Glo  %>%
#   select(-c(Kingdom, Phylum, Class, Family, Order, Genus, Species))  %>%
#   as.data.frame()
# rownames(matrix)  <- matrix$ASV_ID
# matrix  <- matrix[,-1]
#
# #OTU-table prep
# #here ASV abundances, ASVs in rows, samples in columns, nur GLO
# OTU  <- otu_table(as.matrix(matrix), taxa_are_rows = TRUE) # makes otu_table for phyloseq
#
# #Sample data
# # is my meta table
# sampledata  <- sample_data (data.frame (metaM0, row.names = "sampleID")) ##sampledata for phyloseq
#
# #taxa
# # Add taxa data - maybe do as explained for their trait data#
# ##data frame with rows = ASVs, colums = traits (here: taxa names)
# taxa_df_Glo  <-  taxa_table %>%
#   filter(Phylum=="Glomeromycota")  %>%
#   select (!Kingdom) %>%
#   as.data.frame()
# rownames(taxa_df_Glo)  <- taxa_df_Glo$ASV_ID
# TAX  <- tax_table(as.matrix(taxa_df_Glo[,-1])) #taxa table for phyloseq
#
# #Tree
# #class(ASV_tree_Glo) #is it phylo? Yes, is okay already
#
# ##Put together into phyloseq object
# ps_Glo  <- phyloseq(OTU, TAX, ASV_tree_Glo, sampledata)


# #Phyloseq_all####
# #OTUtable
# OTU_ALL  <- otu_table ( ASV_counts, taxa_are_rows =  TRUE) # change to OTU_table object
#
# #taxa table
# taxa_test <- as.matrix(read.table("results/ASVs_taxonomy_simple.tsv", header=T,row.names=1, check.names=F, sep="\t")) # load
# TAX_ALL  <- tax_table(gsub(taxa_test[, colnames(taxa_test)],   pattern = "[a-z]__", replacement = ""))
#                 # remove name strings p__ etc with gsub. change to phyloseq taxa table with tax_table ()
#
# #Sample data
# # is my meta table
# sampledata  <- sample_data (data.frame (metaM0, row.names = "sampleID")) ##sampledata for phyloseq

# ##Phyloseqtree for the full dataset here for reload if needed:####
# #
# dna<- readDNAStringSet ("results/ASVs.fa")
# #
# #
# ASVs_align  <-msa(dna, "ClustalW")
# ASVs_align_forseqinr  <- msaConvert(ASVs_align, type= "seqinr::alignment")
# #
### Citing this package
# If you use this package for research that is published later,
#you are kindly asked to cite it as follows:
#   U. Bodenhofer, E. Bonatesta, C. Horejš-Kainrath, and S. Hochreiter (2015).
#msa: an R package for multiple sequence alignment. Bioinformatics 31(24):3997-3999. DOI: 10.1093/bioinformatics/btv494.
 # Moreover, we insist that, any time you use/cite the package,
#you also cite the original paper in which the algorithm/method/package
#that you have been using has been introduced:
#   ClustalW:
#   J. D. Thompson, D. G. Higgins, and T. J. Gibson (1994).CLUSTAL W: improving the sensitivity of progressive multiple sequence alignment through sequence weighting, position-specific gap
#penalties and weight matrix choice. Nucleic Acids Res., 22(22):4673 4680. DOI: 10.1093/nar/22.22.4673.

# #Build distance matrix as prep for tree
# d <- dist.alignment(ASVs_align_forseqinr, "identity")  #seqinr package
# #tree
# #neighbor joining
# ASV_tree  <- nj(d)  #ape package
# ps_ALL  <- phyloseq(OTU_ALL, TAX_ALL, sampledata, ASV_tree)
# write_rds(ps_ALL, "data/ps_ALL.rds")




phang.align <- phyDat(as(ASVs_align_Glo_forseqinr, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
### oder bootstrap 
fit2 = bootstrap.pml(fit, bs = 100)


tree_w_BS <- plotBS(fit$tree, fit2, type = "none")



## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

#### Exchange old tree for new tree!!
phy_tree(ps_Glo)  <- phy_tree(fitGTR$tree)

#  Run Ch2 analysis (PD, MPD, core etc ahagian and PCA!)








# #Exp1 - all metrics in one df #####
# All_metrics  <-
#   adiv_richness %>% 
#   left_join(metaM0)  %>%  
#   group_by(PlantSpeciesfull)  %>% 
#   mutate (Chao1 = mean (Chao1)) %>% 
#   mutate (Shannon = mean (Shannon)) %>%  
#   select (sampleID, Chao1, Shannon) %>% 
#   select (!sampleID) %>% 
#   unique()  %>%  
#   left_join (comp_units %>% select (!repl) %>% select (!"1-CU")) %>%   #includes CU and unique and richness
#   left_join(cores) %>% 
#   select (!c(PlantType,uniqueASVsperPlSpe, PlantFamily))  %>% 
#   left_join (stand_pd_Glo_all %>% 
#                left_join(ses.MPD_Glo, by ="sampleID")  %>%  
#                left_join (metaM0)  %>% 
#                drop_na ()  %>% 
#                group_by (PlantSpeciesfull)  %>%   
#                mutate (mean_pd = mean (pd.obs), mean_mpd = mean (mpd.obs)) %>%  
#                select (PlantSpeciesfull,   mean_pd, mean_mpd,PlantFamily)   %>% 
#                unique ()  )  %>% 
#   dplyr::rename("Richness"= "mean_n_ASV_per_species", "CU" = "CUnits", "β(core)"  = "perc_core",
#                 "PD" = "mean_pd", "MPD" = "mean_mpd" )  %>% 
#   as.data.frame(row.names = NULL)
# 
# rownames(All_metrics)  <- All_metrics$PlantSpeciesfull
# All_metrics  <- All_metrics[,-1]

#  Who are the generalists?######
# scaled_metrics  <-
#   All_metrics %>%  
#   as_tibble (rownames =  "tm") %>%  
#   mutate_at (c(2:9), scale) %>% 
#   mutate (b_core = `ß(core)`[,1] * (-1)) %>% 
#   select (!(`ß(core)`)) %>% 
#   pivot_longer(cols = c(2:8,10), values_to = "valueINP", names_to = "metric") %>% 
#   ggplot(aes(x= tm, y=metric))  +
#   geom_point (aes(size=valueINP, color = valueINP)) +
#   coord_flip() +
#   scale_color_gradient2(low = "#f90025", mid= "#d0bbb8", high = "#689178", na.value="#89e1ff") +
#   ylab("Diversity metric") +
#   xlab ("Plant Species") +
#   theme_light () +
#   theme(axis.text.y= element_text(family= "sans", face= "italic", size = 11)) +
#   theme (axis.text.x=element_text(family= "sans", size = 11))  +
#   theme (axis.title = element_text(family = "sans", size = 13 ),  
#          legend.title=element_text(size=12), 
#          legend.text=element_text(size=11))



#PCA Exp 1 #####

# Ch2 only 8 species #####


#Get values for all metrics in one tibble#####
ALL_metricsCh2_8  <-
  adiv_richness %>% 
  left_join(metaM0)  %>%  
  group_by(PlantSpeciesfull)  %>% 
  mutate (Chao1 = mean (Chao1)) %>% 
  mutate (Shannon = mean (Shannon)) %>%  
  select (sampleID, Chao1, Shannon) %>% 
  select (!sampleID) %>% 
  unique()  %>%  
  left_join (comp_units %>% select (!repl) %>% select (!"1-CU")) %>%   #includes CU and unique and richness
  left_join(cores) %>% 
  select (!c(PlantType, Chao1, unique, PlantFamily))  %>% 
  left_join (stand_pd_Glo_all %>% 
               left_join(ses.MPD_Glo, by ="sampleID")  %>%  
               left_join (metaM0)  %>% 
               drop_na ()  %>% 
               group_by (PlantSpeciesfull)  %>%   
               mutate (mean_pd = mean (pd.obs), mean_mpd = mean (mpd.obs)) %>%  
               select (PlantSpeciesfull,   mean_pd, mean_mpd,PlantFamily)   %>% 
               unique ()  )  %>% 
  dplyr::rename("Richness"= "mean_n_ASV_per_species", "CU" = "CUnits", "β(core)"  = "perc_core",
                "PD" = "mean_pd", "MPD" = "mean_mpd" )  %>% 
  right_join(meta_plants %>%  filter (PlaSpe != "soil") %>% select (PlantSpeciesfull) , by = "PlantSpeciesfull") %>% 
  as.data.frame(row.names = NULL)

rownames(ALL_metricsCh2_8)  <- ALL_metricsCh2_8$PlantSpeciesfull
ALL_metricsCh2_8  <- ALL_metricsCh2_8[,-1]

# #  Who are the generalists?
# scaled_metrics  <-
#   ALL_metricsCh2_8 %>%  
#   as_tibble (rownames =  "tm") %>%  
#   mutate_at (c(2:9), scale) %>% 
#   mutate (b_core = `ß(core)`[,1] * (-1)) %>% 
#   select (!(`ß(core)`)) %>% 
#   pivot_longer(cols = c(2:8,10), values_to = "valueINP", names_to = "metric") %>% 
#   ggplot(aes(x= tm, y=metric))  +
#   geom_point (aes(size=valueINP, color = valueINP)) +
#   coord_flip() +
#   scale_color_gradient2(low = "#f90025", mid= "#d0bbb8", high = "#689178", na.value="#89e1ff") +
#   ylab("Diversity metric") +
#   xlab ("Plant Species") +
#   theme_light () +
#   theme(axis.text.y= element_text(family= "sans", face= "italic", size = 11)) +
#   theme (axis.text.x=element_text(family= "sans", size = 11))  +
#   theme (axis.title = element_text(family = "sans", size = 13 ),  
#          legend.title=element_text(size=12), 
#          legend.text=element_text(size=11))
# 

#PCA Exp 1 ####
#PCA_all_metrics  <- vegan::rda (All_metrics, scale = T)




PCA_all_metrics_E1 <- PCA(ALL_metricsCh2_8, quali.sup = 8,   scale.unit = T, graph = F)

PCA_arrows_metrics<-
  PCA_all_metrics_E1$var$coord %>%  
  as_tibble (rownames = "metric") %>% 
  dplyr::rename("D1end" = "Dim.1", "D2end"= "Dim.2")

PCA_metric_data   <- PCA_all_metrics_E1$ind$coord  %>%  as_tibble(rownames = "PlantSpeciesfull") %>% 
  left_join(metaM0 %>%  select (PlantSpeciesfull, PlantFamily) %>%  unique ())

PCA_eig_m_Dim1 <- round (PCA_all_metrics_E1$eig[1,2],1)

PCA_eig_m_Dim2 <- round (PCA_all_metrics_E1$eig[2,2],1)
#Plot 

plotPca  <- 
  PCA_metric_data %>%  
  mutate (PlantFamily = str_replace(PlantFamily, "Soil","Soil control")) %>% 
  ggplot (aes (x = Dim.1, y = Dim.2, color = PlantFamily)) + 
  geom_point (size =3) +
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) + 
  #stat_ellipse(aes (x= Dim.1, y = Dim.2, color = group))  +   # thisis 95%confidence
  geom_segment(data = PCA_arrows_metrics, aes (x=0, xend= D1end*4.5, y = 0, yend = D2end*4.5), 
               arrow = arrow(length = unit(0.3, "picas")), color = "blue", inherit.aes = F)  +
  geom_text_repel ( data = PCA_arrows_metrics, aes (x  = D1end*4.5, y = D2end*4.5, label = metric), 
                    color = "blue", inherit.aes = F , force = 0.6) + 
  geom_text_repel(data = PCA_metric_data, aes (x = Dim.1, y = Dim.2, label = PlantSpeciesfull), fontface = "italic", inherit.aes = F) +
  theme_minimal() + 
  xlab(label = "PC1 (62.1.9 %)") +
  ylab ("PC2 (27.8.8 %)") + 
  guides (color= guide_legend( "Plant family")) +
  theme (legend.position = "bottom") +
  scale_color_manual(values=c("Asteraceae"= "#edc948" ,"Cyperaceae" ="#5A6351" ,
                              "Fabaceae" = "#f28e2b" , 
                              "Plantaginaceae"= "#4e79a7","Soil control"= "#9c755f",
                              "Poaceae" = "#7F9A65" )) 



# plotPca  <- fviz_pca_biplot(PCA_all_metrics,axes = c(1,2), 
#                             col.var = "contrib", 
#                             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
#                             repel =T)

plot1 <- fviz_contrib(PCA_all_metrics_E1, choice = "var", axes = 1, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC1")
plot2 <- fviz_contrib(PCA_all_metrics_E1, choice = "var", axes = 2, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC2")
plot3 <- fviz_contrib(PCA_all_metrics_E1, choice = "var", axes = 3, 
                      font.main = c(size =12 ), title= "Contribution of variables to PC3")

ggarrange (plotPca,
           ggarrange (plot1, plot2, plot3, nrow = 3, ncol= 1, labels= c("B", "C", "D")), 
           nrow=1, ncol=2, widths =  c(2, 1), labels= "A") 



# Exp 2 - all metrics ####


# together PD, MPD
PD_MPD_M1 <-  
  stand_pd_Glo_all_M1 %>%  
  select (sampleID, pd.obs.z) %>% 
  left_join(ses.MPD_Glo_M1, by = "sampleID") %>%  
  left_join (meta_M1) %>% 
  select (sampleID, pd.obs, mpd.obs, PlaSpe, PlantSpeciesfull) %>% 
  drop_na() %>% 
  group_by (PlantSpeciesfull) %>% 
  mutate (mean_pd = mean (pd.obs.z), mean_mpd = mean (mpd.obs.z)) %>%  
  select (PlantSpeciesfull,  mean_pd, mean_mpd)   %>% 
  unique ()  





#Get values for all metrics in one tibble
All_metrics_M1  <-
  adiv_richness_M1 %>% 
  group_by(PlantFamily, PlantSpeciesfull)  %>% 
  #mutate (Chao1 = mean (Chao1)) %>% 
  mutate (Shannon = mean (Shannon)) %>%  
  select (sampleID, Shannon) %>% 
  select (!sampleID) %>% 
  unique()  %>%  
  left_join (comp_units_M1 %>%  select (!'1-CU'))    %>% select (!"uniqueASVsperPlSpe") %>%   #includes CU and unique and richness
  left_join(cores_roots_M1) %>% 
  left_join (PD_MPD_M1) %>% 
  dplyr::rename("Richness"= "meanASV", "CU" = "CUnits", "β(core)"  = "perc_core",
                "PD" = "mean_pd", "MPD" = "mean_mpd" )  %>% 
  as.data.frame(row.names = NULL) %>% 
  relocate (where (is.numeric))
rownames(All_metrics_M1)  <- All_metrics_M1$PlantSpeciesfull

##calc of eigenvalues for PC1 for text
PCA_all_metrics_M1 <- PCA(All_metrics_M1, quali.sup = c(7,8),   scale.unit = T, graph = F, ncp=3)




All_metrics_M1_8Sp_df<- All_metrics_M1_8Sp  %>%  data.frame(row.names = "PlantSpeciesfull")

# PCA Exp 2 ####
PCA_all_metrics_E2 <- PCA(All_metrics_M1_8Sp_df , quali.sup = c(8,9),   scale.unit = T, graph = F, ncp = 3)


##Procrustes analysis between PCAs#####


## Needs 
# symmetry = T - both solutions are scaled to#
# unit variance
# results in procrustes m2



procr <- procrustes (PCA_all_metrics_E1$ind$coord, PCA_all_metrics_E2$ind$coord, symmetric = T)
#choice = c(1,2)  )  # with choice the number of dimensions is chosen
summary (procr)


# Plotting procrustes results##

plot (procr, kind =2)
#Kind 2 plots show the residuals for each sample. 
#This allows identification of samples with the worst fit. 
#The horizontal lines, from bottom to top, are 
#the 25% (dashed), 50% (solid), and 75% (dashed) quantiles of the residuals.

##Test of Significance

#Function protest is a permutational test of the significance of the procrustes result.
#It is based on the correlation from a symmetric Procrustes analysis.


procrtest <- protest(X = PCA_all_metrics_E1$ind$coord , Y = PCA_all_metrics_E2$ind$coord , 
                     scores = "sites", permutations = 99999)



## Procrustes protest raw data  #####
procr_test<- 
All_metrics_M0_8Sp %>%  select (-Exp, -PlantFamily) %>% 
  #mutate_at (vars (-PlantSpeciesfull), rescale) %>%   ## in case scaled beforehand?? Like in the radar plots , would be measured relative to other plants
  pivot_longer(!PlantSpeciesfull, names_to = "metric", values_to = "Mvalue") %>%  
  left_join(All_metrics_M1_8Sp %>%  select (-Exp, -PlantFamily) %>% 
              #mutate_at (vars (-PlantSpeciesfull), rescale) %>% 
              pivot_longer(!PlantSpeciesfull, names_to = "metric", values_to = "valuesE2")) %>%
  split(~PlantSpeciesfull) %>% 
  map (~protest(X= .$Mvalue, Y = .$valuesE2, scale = F, scores = "sites", permutations = 9999, symmetric = F))


# randomisation (permutational) test estimate the significance of an observed statistic
# relative to a large number of values for the same statistic that are generataed 
# by permuting the original data #
# get results in table #
map_dfr(procr_test, ~unlist (c(.$ss, .$svd$d, .$signif)) ) %>% 
  add_column ("procr" = c("M2", "correl", "sig")) %>%  
  mutate (round (is.numeric(), 2)) %>% 
write.csv("results/Procrustes_test.csv")


## get single results ###
procr <- procrustes (test[[1]]$Mvalue, test[[1]]$valuesE2, symmetric = T)
#If symmetric is set to TRUE, both solutions are first scaled to unit variance, 
#giving a more scale-independent and symmetric statistic, often known as Procrustes m2.
# look at differences##
plot (procr, kind = 2, choices = 1   )
