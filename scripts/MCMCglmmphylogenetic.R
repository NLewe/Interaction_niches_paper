## Phylogenetic informed analysis####


# packages #####
library(tidyverse)
#library (ape)
library (broom.mixed)
library (brms)
library (shinystan)
library (tidybayes)
# data  ####

dataNLrootstosoil <- read_rds ("data/dataNLrootstosoil.rds")
metaM0 <- read_xlsx("dataCh2/meta/M0_meta.xlsx")
meta_plants <- read_rds ("data/meta_plants.rds")
meta_M1Wcontrol <-  read_rds ("data/meta_M1Wcontrol.rds")
PL_FA_conc_Soil <- read_rds("data/PL_FA_conc_soil.rds")




All_metrics_E2_repl <-
  adiv_richness_M1 %>% 
  left_join (comp_units_M1 %>%  select (!'1-CU'))    %>%   #includes CU and unique and richness
  left_join(cores_roots_M1) %>% 
  left_join (stand_pd_Glo_all_M1 %>%  select (!PlantSpeciesfull)) %>% 
  left_join (ses.MPD_Glo_M1 %>%  select (sampleID, mpd.obs)) %>% 
  dplyr::rename("Richness"= "meanASV", "CU" = "CUnits", "β(core)"  = "perc_core",
                 "uniqueASV" = "uniqueASVsperPlSpe") %>% 
  add_column (Exp = "E2")




dataNL <- 
  dataNLrootstosoil %>%  
  select (sampleID, AMF, AMFroot, PlaSpe, 
          PlantSpeciesfull, DW_above, DW_roots, totalAMF) %>% 
  left_join (All_metrics_E2_repl) %>% 
  left_join (PCA_all_metrics_E2$ind$coord %>%  as_tibble (rownames = "PlantSpeciesfull")) %>% 
  mutate (AMF = 100 * AMF, AMFroot = 100 * AMFroot, totalAMF = 100 * totalAMF)


saveRDS (dataNL, "data/dataNL.rds")
# dataNL2 <- 
#   dataNL %>%  filter (PlaSpe != "AchMil") %>% 
#   filter (PlaSpe != "CicInt")

# check out data ###


(plot_hist <- ggplot(dataNL, aes(x = AMF)) +
    geom_histogram(binwidth = 0.1) +
    theme_classic())


## slightly right skewed

# Poisson?? 

# prepare for models #
options(mc.cores = parallel::detectCores())
# to use multiple cors in parallel#

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


# Background  ###
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






#----------------------------------------------------------------#


# BRMS package MCMCglmm #####

## Gutes tutoria: ###
#   https://ourcodingclub.github.io/tutorials/brms/   #

# data ##
head (dataNL)

#To enable parallel computing, 
#you can run this line of code and 
#then later on in the model code, you can specify how many cores you want to use.
options(mc.cores = parallel::detectCores())

# get variance structure for model ###
A <- ape::vcv.phylo(treeML) 

# test first model #
model_simple <- brm(
  AMF ~ Dim.1* Dim.2 * Dim.3 +DW_roots + DW_above + (1|gr(PlaSpe, cov = A)),
  data = dataNL,
  family = gaussian(),
  data2 = list(A = A)
)
# change from 0.8 to higher (0.8 to 1) if divergent transitions ##
control = list(adapt_delta = 0.9)

summary(model_simple)

#On the top of the output, some general information on the model is given, such as family, formula, number of iterations and chains. 

#Next, group-level effects are displayed separately for each grouping factor in terms of standard deviations and 
#(in case of more than one group-level effect per grouping factor; not displayed here) correlations between group-level effects. 

#On the bottom of the output, population-level effects (i.e. regression coefficients) are displayed.
#If incorporated, autocorrelation effects and family specific parameters (e.g., the residual standard deviation ‘sigma’ in normal models) are also given.

# If Rhat is considerably greater than 1 (i.e., > 1.1), the chains have not yet converged and it is
 #necessary to run more iterations and/or set stronger priors
#in general only fully trust the sample if R-hat is less than 1.01. In early workflow, R-hat below 1.1 is often sufficien

# EFF more than 1000 is a good sign


# One way to assess model convergence is by visually examining the trace plots. 
# They should be fuzzy with no big gaps, breaks or gigantic spikes.

plot(model_simple, N = 2, ask = FALSE)

# Effects of population-level predictors can also be visualized with the
 #conditional_effects method
plot(conditional_effects(model_simple), points = TRUE)


hyp <- "sd_PlaSpe__Intercept^2 / (sd_PlaSpe__Intercept^2 + sigma^2) = 0"
hyp <- hypothesis(model_simple, hyp, class = NULL)

#Phylogenetic Model with Repeated Measurements ####

# add the means to the table
dataNL$mean_DW_roots <-
  with(dataNL, sapply(split(DW_roots, PlaSpe), mean)[PlaSpe])

dataNL$mean_DW_above <-
  with(dataNL, sapply(split(DW_above, PlaSpe), mean)[PlaSpe])
#The variable mean_xx just contains the mean of the cofactor for each species. 
head(dataNL)
#The code for the repeated measurement phylogenetic model looks as follows:


#The most important reason to use control is to decrease (or eliminate at best) the number of divergent transitions
#that cause a bias in the obtained posterior samples. Whenever you see the warning
#"There were x divergent transitions after warmup.", you should really think about
#increasing adapt_delta. To do this, write control = list(adapt_delta = <x>), where
#<x> should usually be a value between 0.8 (current default) and 1. Increasing adapt_delta
#will slow down the sampler but will decrease the number of divergent transitions threatening
#the validity of your posterior samples.

# Model AMF ~  Dim.1 ####

## Set appropriate priors using get_prior function #### 
prior <- get_prior ( AMF ~ Dim.1 + mean_DW_above + mean_DW_roots + (1|gr(PlaSpe, cov = A)) + (1|PlantSpeciesfull), data = dataNL, data2 = list(A = A))

# model 1
 ##
model_repeat1 <- brm(
  AMF ~ Dim.1 + mean_DW_above + mean_DW_roots + (1|gr(PlaSpe, cov = A)) + (1|PlantSpeciesfull),
  data = dataNL,
  family = gaussian(),
  data2 = list(A = A),
  prior = prior,  sample_prior = TRUE, chains = 2, cores = 2,
  iter = 4000, warmup = 1000
)

control = list(adapt_delta = 0.9) ## then rerun!!

#The variables PlaSpe and PlantSpeciesfull are identical as they are both identifiers of the species. 
#However, we model the phylogenetic covariance only for PlaSpe and thus the PlantSpeciesfull variable accounts for any specific effect 
#that would be independent of the phylogenetic relationship between species 
#(e.g., environmental or niche effects). 

#obtain model summaries as well as estimates of the phylogenetic signal.
summary(model_repeat1)

# estimate of phylogenetic signal!!
hyp <- paste(
  "sd_PlaSpe__Intercept^2 /",
  "(sd_PlaSpe__Intercept^2 + sd_PlantSpeciesfull__Intercept^2 + sigma^2) = 0"
)
hyp <- hypothesis(model_repeat1, hyp, class = NULL)

hyp
### means we ahave a phylogenetc signal 
plot (hyp)


#So far, we have completely ignored the variability of the cofactor within species. 
#To incorporate this into the model, we define

dataNL$within_spec_DW_roots <- dataNL$DW_roots - dataNL$mean_DW_roots

dataNL$within_spec_DW_above <- dataNL$DW_above - dataNL$mean_DW_above

#and then fit it again using within_spec_xx as an additional predictor.

model_repeat2 <- update(
  model_repeat1, formula = ~ . + within_spec_DW_roots + within_spec_DW_above,
  newdata = dataNL, chains = 2, cores = 2,
  iter = 4000, warmup = 1000
)

summary (model_repeat2)
# compare with first model
summary (model_repeat1)


# estimate phylogenetic signal ###
hyp <- paste(
  "sd_PlaSpe__Intercept^2 /",
  "(sd_PlaSpe__Intercept^2 + sd_PlantSpeciesfull__Intercept^2 + sigma^2) = 0"
)
hyp <- hypothesis(model_repeat1, hyp, class = NULL)

hyp
plot (hyp)




# fit new model, 

model_repeat3 <- update(model_repeat2, formula. = ~ . - mean_DW_above + Dim.1)

model_repeat4 <- update(model_repeat3, formula. = ~ . - within_spec_DW_roots)


model_repeat5 <- update(model_repeat4, formula. = ~ . - (1|PlantSpeciesfull))
model_repeat6 <- update(model_repeat5, formula. = ~ . - within_spec_DW_above)
control = list(adapt_delta = 0.85) ## then rerun!!

loo (model_repeat5, model_repeat6)

plot (conditional_effects(model_repeat6))


# Model AMF ~  mpd.obs ####

## get default  priors  #### 
prior <- get_prior ( AMF ~ mpd.obs + mean_DW_above + mean_DW_roots + (1|gr(PlaSpe, cov = A)) 
                     #+ (1+ mean_DW_roots|PlantSpeciesfull) ### hi is only needed to account for diff between the species 
                     #+ OTHEr than phylogenetic (environmental factors, niches)
                     , data = dataNL, data2 = list(A = A))

# model


mpd_fit1 <- brm(
  AMF ~ mpd.obs + mean_DW_above + mean_DW_roots + (1|gr(PlaSpe, cov = A)) + (1+ mean_DW_roots|PlantSpeciesfull),
  data = dataNL,
  family = gaussian(),
  data2 = list(A = A),
  prior = prior,  sample_prior = TRUE, chains = 4, cores = 8,
  iter = 4000, warmup = 1000
)

control = list(adapt_delta = 0.9) ## then rerun!!

#The variables PlaSpe and PlantSpeciesfull are identical as they are both identifiers of the species. 
#However, we model the phylogenetic covariance only for PlaSpe and thus the PlantSpeciesfull variable accounts for any specific effect 
#that would be independent of the phylogenetic relationship between species 
#(e.g., environmental or niche effects). 

#obtain model summaries as well as estimates of the phylogenetic signal.
summary(mdp_fit1)

# estimate of phylogenetic signal!!
hyp <- paste(
  "sd_PlaSpe__Intercept^2 /",
  "(sd_PlaSpe__Intercept^2 + sd_PlantSpeciesfull__Intercept^2 + sigma^2) = 0"
)
hyp <- hypothesis(mpd_fit1, hyp, class = NULL)

hyp
### means we ahave a phylogenetc signal 
plot (hyp)


#So far, we have completely ignored the variability of the cofactor within species. 
#To incorporate this into the model, we define

## why?? isn't that what the random part would do if 
# 1 + DW_roots | covariancestructure  

dataNL$within_spec_DW_roots <- dataNL$DW_roots - dataNL$mean_DW_roots

dataNL$within_spec_DW_above <- dataNL$DW_above - dataNL$mean_DW_above

#and then fit it again using within_spec_xx as an additional predictor.

mpd_fit2 <- update(
  mdp_fit1, formula = ~ . + within_spec_DW_roots + within_spec_DW_above,
  newdata = dataNL, chains = 4, cores = 8,
  iter = 5000, warmup = 1000
)

summary (mpd_fit2)
# compare with first model
summary (mpd_fit1)


# estimate phylogenetic signal ###
hyp <- paste(
  "sd_PlaSpe__Intercept^2 /",
  "(sd_PlaSpe__Intercept^2 + sd_PlantSpeciesfull__Intercept^2 + sigma^2) = 0"
)
hyp <- hypothesis(mpd_fit2, hyp, class = NULL)

hyp
plot (hyp)




# fit new model, 

mpd_fit3 <- update(mpd_fit2, formula. = ~ . - mean_DW_above)

mpd_fit4 <- update(mpd_fit3, formula. = ~ . - within_spec_DW_roots)

mpd_fit5 <- update(mpd_fit4, formula. = ~ . - within_spec_DW_above)
loo (mpd_fit4, mpd_fit5)

plot (conditional_effects(mpd_fit5))


model_repeat5 <- update(model_repeat4, formula. = ~ . - (1|PlantSpeciesfull))
model_repeat6 <- update(model_repeat5, formula. = ~ . - within_spec_DW_roots)
control = list(adapt_delta = 0.85) ## then rerun!!


VarCorr (mpd_fit5)

head (dataNL)
## Set appropriate priors using get_prior function #### 
prior <- get_prior ( AMF ~ Observed* Shannon * mpd.obs * pd.obs + 
                       mean_DW_above + mean_DW_roots + (1|gr(PlaSpe, cov = A)) + (1|PlantSpeciesfull), data = dataNL, data2 = list(A = A))

# model 1
##
model_all1 <- brm(
  AMF ~ Observed* Shannon * mpd.obs * pd.obs + 
  mean_DW_above + mean_DW_roots + (1|gr(PlaSpe, cov = A)) + (1|PlantSpeciesfull),
  data = dataNL,
  family = gaussian(),
  data2 = list(A = A),
  prior = prior,  sample_prior = TRUE, chains = 2, cores = 2,
  iter = 4000, warmup = 1000
)

control = list(adapt_delta = 0.9) ## then rerun!!

#The variables PlaSpe and PlantSpeciesfull are identical as they are both identifiers of the species. 
#However, we model the phylogenetic covariance only for PlaSpe and thus the PlantSpeciesfull variable accounts for any specific effect 
#that would be independent of the phylogenetic relationship between species 
#(e.g., environmental or niche effects). 

#obtain model summaries as well as estimates of the phylogenetic signal.
summary(model_all1)

# estimate of phylogenetic signal!!
hyp <- paste(
  "sd_PlaSpe__Intercept^2 /",
  "(sd_PlaSpe__Intercept^2 + sd_PlantSpeciesfull__Intercept^2 + sigma^2) = 0"
)
hyp <- hypothesis(model_all1, hyp, class = NULL)

hyp
### means we ahave a phylogenetc signal 
plot (hyp)


#So far, we have completely ignored the variability of the cofactor within species. 
#To incorporate this into the model, we define

dataNL$within_spec_DW_roots <- dataNL$DW_roots - dataNL$mean_DW_roots

dataNL$within_spec_DW_above <- dataNL$DW_above - dataNL$mean_DW_above

#and then fit it again using within_spec_xx as an additional predictor.

model_all22 <- update(
  model_all1, formula = ~ . + within_spec_DW_roots + within_spec_DW_above,
  newdata = dataNL, chains = 2, cores = 2,
  iter = 4000, warmup = 1000
)

summary (model_all22)
# compare with first model




# fit new model, 

model_all3 <- update(model_all22, formula. = ~ . - mean_DW_above - Observed - Shannon:Observed:pd.obs   -Shannon:Observed   -Observed:pd.obs        )

summary (model_all3) 

model_all4 <- update(model_all3, formula. = ~ . - Shannon - pd.obs -within_spec_DW_roots)

summary (model_all4)

model_repeat5 <- update(model_repeat4, formula. = ~ . - (1|PlantSpeciesfull))
model_repeat6 <- update(model_repeat5, formula. = ~ . - within_spec_DW_above)
control = list(adapt_delta = 0.85) ## then rerun!!

loo (model_repeat5, model_repeat6)

plot (conditional_effects(model_repeat6))



