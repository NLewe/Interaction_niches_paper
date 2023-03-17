## Traits for the plants 


# installed package TR8

# packages ####


library (tidyverse)

library(nlme)
library (vegan)
library (lme4)
library (car)
library(stats)
#library(ade4)
#library(ape)
#library(adegenet)
#library(phangorn)
#library(geiger)
#library (concatipede)
#library(phytools)





#Pylogenetic analyses LM#####
## Build tree ####
### Use of different genes to build a phylogenetic tree:
#

## The package concatipede is used to merge the genes into one file after they were aligned per gene.
# https://github.com/tardipede/concatipede
#Please cite this package as: Vecchi M. & Bruneaux M. (2021). concatipede: 
#an R package to concatenate fasta sequences easily. http://doi.org/10.5281/zenodo.5130603

# get package
# devtools::install_github("tardipede/concatipede")

dnaTRNL<-read.dna("data/trnl_fastas/trnlFasta8.fasta",format="fasta")


dnaTRNL
labels(dnaTRNL) #checken
rownames(dnaTRNL)<-c("AchMil" ,"CicInt","PlaLan", "HolLan" ,"PoaCit", "BroWil","SchAru", "AgrCap")

##erst distance matrix machen
##welche nehmen?
#test mit verschiedenen distnace berechnungen
distTN93<-dist.dna(dnaTRNL, model="TN93")
#TN93: Tamura and Nei (1993) developed a model which assumes distinct rates for both kinds of transition 
#(A <-> G versus C <-> T), and transversions. The base frequencies are not assumed to be equal and are estimated from the data. 
#A gamma correction of the inter-site variation in substitution rates is possible.

distTNK81<-dist.dna(dnaTRNL, model="K81")

class(distTN93)  #check dass es auch distance matrix ist




## test tree
tree<-nj(distTN93)
class(tree)
tree

plot(tree,cex=0.6)

library (phangorn)


## trooted tree
# 
 tree_jalView<-read.tree("data/trnl_fastas/NJ_tree_8")
tree_jalView$tip.label # check
plot(tree_jalView)
# 
# # rename 
tree_jalView$tip.label <- c("PlaLan","CicInt","AchMil", "BroWil","AgrCap","SchAru", "PoaCit","HolLan"  ) 
## optimise tree using maximum likelihood
dna2 <- as.phyDat(dnaTRNL) 
class(dna2)




tre.ini<-nj(dist.dna(dnaTRNL, model="TN93"))
tre.ini



#To  initialize  the  optimization  procedure,  we  need  an  initial  fit  for  the model chosen.  
#This is computed using pml
fit.ini<-pml(tre.ini, dna2,k=4)



fit.ini<-pml(tree_jalView, dna2,k=4)


#Optimize tree
fit<-optim.pml(fit.ini, optNni = TRUE, optBf = TRUE, optQ = TRUE,optGamma = TRUE)
fit
class(fit)
names(fit) ##all the useful parameters of the model and the optimal tree
#optimal tree is stored in fit$tree

##test if optimized tree is better than original
anova(fit.ini,fit)
AIC(fit.ini)
AIC(fit)

treeML<-root(fit$tree,1)
plot(treeML)


##PGLS ####
###And now the PGLS
#http://www.phytools.org/Cordoba2017/ex/4/PGLS.html
# Symonds 2014 A primer --
# PGLS == phylogenetic generalised least squares 
# also called :
# PGLS = phylogenetic regression or phylogenetic general linear models

# Background  ####
# because of shared evolutionary history, species do not provide independent data points for analysis
 # respnse variable can't be discrete!! predictor can  binary discrete coding with 0 1 for specific trait can be dealt with
# PGLS likewise identifies from phylogeny the amount of expected correlation between species based on their shared
#evolutionary history, and weights for this in the generalised least squares regression calculation.

#it allows one to incorporate information on the extent of phylogenetic signal in the data
# If there is no phylogenetic signal in the data, then PGLS will return estimates identical to an 
#ordinary least squares regression analysis.

#The simplest way to think of PGLS is as a weighted regression. In a standard regression, each independent data point 
# contributes equally to the estimation of the regression line. 
#By contrast, PGLS ‘downweights’ points that derive from species with shared phylogenetic history.

#important issue is the estimation of effect sizes and associated confidence
#intervals from GLS models. The parameter estimates of slopes (for continuous 
#predictors) and the intercept and differences between means (for
# categorical predictors) are the most important results of the analyses.


# do we have a phylogenetic signal??
# calculate pagels lambda (pagel 1999)
# To test if a phylogenetic signal exists in the data, pagels lambda was calculated in R package phytools 
# based on the tree
# A k value of 0 is consistent with no phylogenetic signal in the trait, whereas a value of 1 is consistent with strong
# phylogenetic signal. Intermediate values of k indicate intermediate phylogenetic
# signal. # I think the problem might be that this calculates the phy signal in the "trait" that I downloaded as 

# Pagel's lambda  ####
x<-fastBM(treeML)
phylosig(treeML,x,method="lambda",test=TRUE)


#define covariance structure
bm<-corBrownian(1,tree_jalView, form = ~ PlaSpe)
# lit Brownian motion model (see earlier), as described by Felsenstein (1985).
#Felsenstein J (1985) Phylogenies and the comparative method. Am Nat 125:1–15


library(tidyverse) # the most important package in the universe
library(readxl)
library(ggplot2)
library (flextable)
library(ggpubr)
library(phyloseq)
library(phylosmith)
library (broom) #needed for tidy linear modelling
library (lmerTest)
library (picante)
library(rstatix)
library(FactoMineR) # simple ordinations
library(factoextra)
library (rmarkdown)
library (ggsignif) # adds signif to ggplot
library (lme4)
library (nlme)
library (broom) # for pretty lm tables
library(broom.mixed) # for pretty mixed model tables
#library (wPerm) #permutation ANOVA etc
#library(paletteer) 


library (stats)
library (ggtree)
library(ggrepel)
library (car)



## AMF biomass in the soil - correlated to PC of plant species ' interaction niches?


# Useful data #####
dataNLrootstosoil <- read_rds ("data/dataNLrootstosoil.rds")
dataNLrootstosoil_gen <- read_rds ("data/dataNLrootstosoil_gen.rds")
metaM0 <- read_xlsx("dataCh2/meta/M0_meta.xlsx")
meta_plants <- read_rds ("data/meta_plants.rds")
meta_M1Wcontrol <-  read_rds ("data/meta_M1Wcontrol.rds")
meta_M1Wcontrol <-  read_rds ("data/meta_M1Wcontrol.rds")
PL_FA_conc_Soil <- read_rds("data/PL_FA_conc_soil.rds")
All_metrics_M1_8Sp <- read_rds ("data/All_metrics_M1_8Sp.rds")
PC_data_M1roots  <-  read_rds ("data/PC_data_M1roots.rds")


Model1_AMF  <-
dataNLrootstosoil %>%  # run   dataNLrootstosoil_gen
  select (sampleID, PlaSpe,AMF, Dim.1, Dim.2, Dim.3, DW_roots, DW_above)  %>% 
  dplyr::filter(!is.na (AMF))%>% 
  ungroup() %>% 
  select (!sampleID) %>% 
  pivot_longer(cols =  !c(AMF, PlaSpe, DW_roots, DW_above), names_to = "Dimx", values_to = "dimValue" ) %>% 
  split (.$Dimx) %>% 
  purrr::map ( ~ nlme::gls (AMF ~ dimValue*DW_roots* DW_above, data = ., correlation = bm)) %>% 
  purrr::map (~ summary ())


## Add dat richness etc into table ###


Metrics_exp2  <-
  adiv_richness_M1 %>% select (sampleID, Observed, Shannon) %>% 
  left_join(meta_M1 %>%  select (sampleID, PlaSpe))  %>%  
  left_join (stand_pd_Glo_all_M1 %>% select (sampleID, pd.obs)) %>% 
  left_join(ses.MPD_Glo_M1 %>%  select (sampleID, mpd.obs))  %>%  
  drop_na () %>% 
  add_column (exp = "Exp2")

dataNLrootstosoil  <-
dataNLrootstosoil %>%  left_join(Metrics_exp2)

### AMFroot ~ ####


#investigate models
model1<-nlme::gls(AMFroot ~ Dim.1  * DW_roots *DW_above , data=dataNLrootstosoil, correlation =bm, method = "ML")
summary(model1)
Anova (model1, type = 3)

model1LM <- update (model1, . ~ ., method = "ML")
# library (MASS)
# 
# fit <- stepAIC (model1LM,direction="both",trace=0)
# detach("package:MASS", unload = TRUE)
# does not work
# modell is overfitted, AIC cannot be calculated , is that because we have a tree with only 8 tips, are we 
# I was not able to fit a rooted tree using a moss, tree was catastrophe
#

model1B<-nlme::gls(AMFroot ~DW_roots * Dim.1 , data=dataNLrootstosoil, correlation =bm)
summary(model1B)
Anova (model1B, type = 3)
model1BLM <- update (model1B, . ~ ., method = "ML")
Anova (model1BLM, type = 3)


model1C<-nlme::gls(AMFroot ~ DW_roots  , data=dataNLrootstosoil, correlation =bm)
summary(model1C)
Anova (model1C, type = 3)
# all other combinations tested, too!

model1D<-nlme::gls(AMFroot ~ DW_roots * DW_above, data=dataNLrootstosoil, correlation =bm)
summary(model1D)
Anova (model1D, type = 3)

model1E<-nlme::gls(AMFroot ~ Dim.1  + DW_roots + DW_above , data=dataNLrootstosoil, correlation =bm, method = "ML")
summary(model1E)
Anova (model1E, type = 3)

# Best Model 1 ####
model1C
model1C<-nlme::gls(AMFroot ~ DW_roots  , data=dataNLrootstosoil, correlation =bm)
model1C.LM <- lm (AMFroot ~ DW_roots  , data=dataNLrootstosoil)
summary (model1C.LM)
anova (model1C, model1C.LM)
 # Change form REML to ML
# model1LM  <- update (model1, . ~ ., method = "ML")
#  
# model1BLM  <- update (model1B, . ~ ., method = "ML")
# model1CLM  <- update (model1C, . ~ ., method = "ML")
# 
# 
# Anova( model1BLM, type =3)
# anova (model1BLM, model1CLM)

# geht nicht##

# careful, in pgls it is not the residuals that are tested but 
#check if residuals are okay
plot (model1C, select = c(1))
# p_mdl1 <- round (glance (selectedmdl1)$p.value[[1]],3)
# 
# #Check assumptions
#   #normality
resid(model1C) %>%  shapiro.test()
# non normal 


### AMF (soil) ~ ####
modelAMF1  <-nlme::gls(AMF ~  mpd.obs + DW_roots + DW_above, data=dataNLrootstosoil, correlation =bm, 
            na.action = na.omit, method = "ML" )

#fit<-stepAIC(modelAMF1,direction="both",trace=0)  ## ovrfitted ##
summary (modelAMF1)
Anova (modelAMF1, type =3)


model1LM  <- update (modelAMF1, . ~ ., method = "ML")

modelAMF1B  <-nlme::gls(AMF ~    mpd.obs, data=dataNLrootstosoil, correlation =bm, 
                        na.action = na.omit )
model1BLM  <- update (modelAMF1B, . ~ ., method = "ML")

Anova (modelAMF1B, type =3)


#AIC (model1BLM, model1LM)

modelAMF2  <-gls(AMF ~ DW_roots + DW_above, data=dataNLrootstosoil, correlation =bm, 
                  na.action = na.omit )

modelAMF3  <-gls(AMF ~ DW_roots , data=dataNLrootstosoil, correlation =bm, 
                 na.action = na.omit )

modelAMF4  <-gls(AMF ~  DW_above, data=dataNLrootstosoil, correlation =bm, 
                 na.action = na.omit )
### No relevant variables


modeltotAMF1  <-nlme::gls(totalAMF ~   DW_roots * DW_above * Dim.1, data=dataNLrootstosoil, correlation =bm, 
                       na.action = na.omit )
summary (modeltotAMF1)
Anova (modeltotAMF1, type =3)

modeltotAMF2  <-nlme::gls(totalAMF ~ DW_roots * DW_above, data=dataNLrootstosoil, correlation =bm, 
                 na.action = na.omit )

modeltotAMF3  <-nlme::gls(totalAMF ~ DW_roots , data=dataNLrootstosoil, correlation =bm, 
                 na.action = na.omit )

modeltotAMF4  <-nlme::gls(totalAMF ~  DW_above, data=dataNLrootstosoil, correlation =bm, 
                 na.action = na.omit )


## no relevant model ###
### PL AMF biomarker ####



# presenting PGLS results ##
# estimates and errors
# if appropriate r, t, or F p values
# the estimate of the phylogenetic signal associated with the regressions along with its confidence intervals
# (indicates as to the extent the phylogeny affects the error structure of the data)

# Zuur 2009 p 75 ####


 M.lm <- nlme::gls(AMF ~ Dim.1 * DW_roots, data=dataNLrootstosoil, na.action = na.omit)


anova(M.lm, model1C)


# test assumptions #####

#check if residuals are okay
plot (model1C, select = c(1))
# p_mdl1 <- round (glance (selectedmdl1)$p.value[[1]],3)
# 
# #Check assumptions
#   #normality
resid_selectedmdl1  <-resid(model1C)
shapiro.resid_mdl1  <- shapiro.test(resid_selectedmdl1) # if p> 0.05, data is normal # NOT NORMAL ##qqnorm (resid_selectedmdl1)

model1C<-nlme::gls(AMFroot ~ DW_roots  , data=dataNLrootstosoil, correlation = corBrownian(1,tree))
bartlett.test()

lillie.test (chol(solve(vcv(tree)))
             


#The concatipede package allows to concatenate sequences from different Multiple Sequence Alignments (MSAs) based on a correspondence table that indicates how sequences from different MSAs are linked to each other.

# For this package to work, all the MSAs in fasta format should be placed in the same folder (that can be the working directory).
#Only the fasta files to be concatenated should be present in the directory.
# If other fasta files are present and you don´t want to include them, see the concatipede_prepare(exclude = "") option.

# Save the path to the initial directory for later clean-up
old_dir <- getwd()
# Create a directory to put the fasta files for this example
dir.create("concatipede_test")
# Set it as the working directory
setwd("concatipede_test")
# Copy the example fasta files shipped with the package into that directory
example_files = list.files(system.file("extdata", package="concatipede"), full.names = TRUE)
file.copy(from = example_files, to = getwd())


library(concatipede)
library(tidyverse)
find_fasta()
concatipede_prepare(out = "seqnames")

# Clean-up: go back to initial directory
setwd(old_dir)
