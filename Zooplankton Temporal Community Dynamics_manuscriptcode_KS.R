###Installing Packages####
install.packages("tidyverse")
install.packages("rstudioapi")
install.packages("hillR")
install.packages("vegan")
install.packages("ggplot2")
install.packages("lme4")
install.packages("lmerTest")
install.packages("emmeans")
install.packages("dplyr")
install.packages("ggsignif")
install.packages("ggrepel")
install.packages("tibble")
install.packages("indicspecies")
install.packages("permute")
install.packages("igraph")
install.packages("fBasics")
install.packages("gtools")
install.packages("ggeffects")
install.packages("car")
install.packages("DHARMa")
install.packages("glmmTMB")
install.packages("TITAN2")

###Loading Libraries####
library(vegan)
library(dplyr)
library(tibble)
library(tidyverse)
library(ggplot2)
library(ggsignif)
library(ggeffects)
library(ggrepel)
library(gtools)
library(rstudioapi)
library(hillR)
library(lme4)
library(lmerTest)
library(emmeans)
library(DHARMa)
library(glmmTMB)
library(car)
library(indicspecies)
library(TITAN2)
library(permute)

# library(igraph)
# library(fBasics)

###Loading Data####
####Set working directory####
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

####Read in cleaned datasets####
taxa.density.dataNEW = read.csv("taxa.density.dataNEW.csv", header = T, sep = ",")
chem.hach.dataWide = read.csv("chem.hach.dataWide.csv", header = T, sep = ",")
diversity_nutrient_data = read.csv("diversity_nutrient_data.csv", header = T, sep = ",")
taxa.biomass.dataNEW = read.csv("taxa.biomass.dataNEW.csv", header = T, sep = ",")
Hellenger_biomass = read.csv("Hellenger_biomass_site_andRound.csv", header = T, sep = ",")
PA_ILMS_June_July_site_resid = read.csv("PA_ILMS_June_July_site_resid.csv", header = T, sep = ",")
PA_ILMS_June_July_species_resid = read.csv("PA_ILMS_June_July_species_resid.csv", header = T, sep = ",")
PA_ILMS_July_August_site_resid = read.csv("PA_ILMS_July_August_site_resid.csv", header = T, sep = ",")
PA_ILMS_July_August_species_resid = read.csv("PA_ILMS_July_August_species_resid.csv", header = T, sep = ",")
PA_ILMS_June_August_site_resid = read.csv("PA_ILMS_June_August_site_resid.csv", header = T, sep = ",")
PA_ILMS_June_August_species_resid = read.csv("PA_ILMS_June_August_species_resid.csv", header = T, sep = ",")


###Zooplankton Taxa Summary#####################################################
table(taxa.biomass.dataNEW$Taxon)

taxa_by_round <- taxa.biomass.dataNEW %>%
  filter(TotalBiomassPerSample > 0) %>%
  distinct(Taxon, Round) %>%
  group_by(Taxon) %>%
  summarise(
    n_months = n_distinct(Round),
    months_present = paste(sort(Round), collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(n_months == 1)

taxa_count_by_round <- taxa.biomass.dataNEW %>%
  group_by(Round) %>%
  summarise(
    n_taxa = n_distinct(Taxon),
    .groups = "drop"
  )
taxa_count_by_round


###N:P Site Characterization####################################################
##Biplot of TN & TP with Site & Round
chem.hach.dataWide$Round <-
  factor(chem.hach.dataWide$Round)

chem.hach.dataWide <- chem.hach.dataWide %>%
  mutate(
    Sampling_Month = case_when(
      Round %in% c("1") ~ "June",
      Round %in% c("2") ~ "July",
      Round %in% c("3") ~ "August"
    )
  )
chem.hach.dataWide$Sampling_Month <- factor(chem.hach.dataWide$Sampling_Month,
                                                             levels = c("June", "July", "August"))

ggplot(chem.hach.dataWide,
       aes(x = Total.Nitrogen,
           y = Total.Phosphorus,
           color = Sampling_Month,
           shape = WetlandType,
           label = Site)) +
  geom_point(size = 2.25, alpha = 1) +
  #Eutrophic line
  geom_hline(
    yintercept = 0.024,
    linetype = "dashed",
    color = "red",
    alpha = 0.20
  ) +
  #hypereutrophic line
  geom_hline(
    yintercept = 0.096,
    linetype = "dashed",
    color = "red",
    alpha = 0.20
  ) +
  geom_text_repel(size = 2.5) +
  labs(
    x = "Total Nitrogen (mg/L)",
    y = "Total Phosphorus (mg/L)",
    color = "Sampling_Month"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "inside",
    legend.position.inside = c(0.90, 0.84),
    legend.title = element_blank(),
    legend.text = element_text(size = 6.5), #face = "bold"
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
#ggsave("TN to TP.png", width = 7, height = 5, dpi = 300)
#ggsave("TN to TPsmall.png", width = 6, height = 4, dpi = 800)


#OBJECTIVE 1####################################################################
#Diversity state, seasonal change, and nutrient context

###Calculate Hill Diversity#####################################################
####Site AND ROUND
community_diversity_bysite_andRound <- taxa.density.dataNEW %>%
  group_by(Site, Round, Taxa) %>%
  summarise(TotalCount = sum(SubsampleCount, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(
    names_from = Taxa,
    values_from = TotalCount,
    values_fill = 0
  )

invert_hill_results_bysite_andRound <- data.frame(
  Site = community_diversity_bysite_andRound$Site,
  Round = community_diversity_bysite_andRound$Round,
  q0.invert = hill_taxa(select(community_diversity_bysite_andRound, -Site, -Round), q = 0),
  q1.invert = hill_taxa(select(community_diversity_bysite_andRound, -Site, -Round), q = 1),
  q2.invert = hill_taxa(select(community_diversity_bysite_andRound, -Site, -Round), q = 2)
)

###1A.Seasonal change in Hill Diversity#########################################
###Models testing seasonal change in diversity

##Data setup
#create sampling month column & define as factor
invert_hill_results_bysite_andRound <- invert_hill_results_bysite_andRound %>%
  mutate(
    Sampling_Month = case_when(
      Round %in% c("1") ~ "June",
      Round %in% c("2") ~ "July",
      Round %in% c("3") ~ "August"
    )
  )
invert_hill_results_bysite_andRound$Sampling_Month <- factor(invert_hill_results_bysite_andRound$Sampling_Month,
                                                             levels = c("June", "July", "August"))

invert_hill_results_bysite_andRound$Round <- factor(invert_hill_results_bysite_andRound$Round)

#Add wetland classifications to diversity dataset
invert_hill_results_bysite_andRound <- invert_hill_results_bysite_andRound %>%
  mutate(
    WetlandType = case_when(
      Site %in% c("DD15N", "SS", "DD8", "DD65", "LP", "WW", "FCC2",
                  "KS", "RS") ~ "CREP",
      Site %in% c("UH1", "UH2") ~ "Restored for habitat",
      Site %in% c("Bjork") ~ "Remnant"
    )
  )

####Richness######################################################################
##fit model
#generalized poisson
q0_model_GP<- glmmTMB(
  q0.invert ~ Sampling_Month + (1 | Site),
  family = genpois(link = "log"),
  data = invert_hill_results_bysite_andRound
) 
summary(q0_model_GP)

##Inspect model fit
sim_GP <- simulateResiduals(q0_model_GP)
plot(sim_GP)
testDispersion(sim_GP)

##Pairwise comparisons
pairwise_q0 <- emmeans(q0_model_GP, ~ Sampling_Month)
pairs(pairwise_q0)

##Plot of observed + predicted Hill diversity through time
#Generate prediction dataset
pred <- ggpredict(q0_model_GP, terms = "Sampling_Month")
pred <- pred %>%
  rename(Sampling_Month = "x")

#Richness plot
TemporalHillRichness <- ggplot(invert_hill_results_bysite_andRound, aes(Sampling_Month, q0.invert)) +
  geom_line(aes(group = Site,
                color = WetlandType),
            linewidth = 0.5,
            alpha = 1) +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  geom_point(aes(color = WetlandType),
             size = 1.5) +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  scale_x_discrete(expand = expansion(add = c(0.05, 0.05))) +
  geom_line(data = pred,
            aes(x = Sampling_Month,
                y = predicted,
                group = 1),
            inherit.aes = FALSE,
            linewidth = 1.25,
            color = "black") +
  geom_point(data = pred,
             aes(x = Sampling_Month,
                 y = predicted),
             inherit.aes = FALSE,
             size = 2.25,
             color = "black") +
  theme_minimal() +
  labs(
    x = "Sampling Month",
    y = "Richness (q=0)"
  ) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    #legend.text = element_text(size = 8),
    #legend.title = element_text(size = 10),
    panel.grid.minor = element_blank()
  )
TemporalHillRichness

#ggsave("TemporalHillDiversity_ObservedandMean_Richness.png", width = 7, height = 5, units = "in", dpi = 300)


####Shannons######################################################################
##fit model 
hist(invert_hill_results_bysite_andRound$q1.invert)

q1_model <- glmmTMB(
  q1.invert ~ Sampling_Month + (1 | Site),
  family = gaussian(link = "log"),
  data = invert_hill_results_bysite_andRound
)
summary(q1_model)

#Inspect model fit
sim <- simulateResiduals(q1_model)
plot(sim)
testDispersion(sim)

#Pairwise comparisons
pairwise <- emmeans(q1_model, ~ Sampling_Month, bias_correction = TRUE)
pairs(pairwise)

##Plot of observed + predicted Hill diversity through time
#Generate prediction dataset
pred <- ggpredict(q1_model, terms = "Sampling_Month")
pred <- pred %>%
  rename(Sampling_Month = "x")

#Exponential Shannons plot
TemporalHillShannons <- ggplot(invert_hill_results_bysite_andRound, aes(Sampling_Month, q1.invert)) +
  geom_line(aes(group = Site,
                color = WetlandType),
            linewidth = 0.5,
            alpha = 1) +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  geom_point(aes(color = WetlandType),
             size = 1.5) +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  scale_x_discrete(expand = expansion(add = c(0.05, 0.05))) +
  geom_line(data = pred,
            aes(x = Sampling_Month,
                y = predicted,
                group = 1),
            inherit.aes = FALSE,
            linewidth = 1.25,
            color = "black") +
  geom_point(data = pred,
             aes(x = Sampling_Month,
                 y = predicted),
             inherit.aes = FALSE,
             size = 2.25,
             color = "black") +
  theme_minimal() +
  labs(
    x = "Sampling Month",
    y = "Exponential Shannon (q=1)"
  ) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    #legend.text = element_text(size = 8),
    #legend.title = element_text(size = 10),
    panel.grid.minor = element_blank()
  )
TemporalHillShannons

#ggsave("TemporalHillDiversity_ObservedandMean_Shannon.png", width = 7, height = 5, units = "in", dpi = 300)


####Inverse Simpsons##############################################################
#fit model
hist(invert_hill_results_bysite_andRound$q2.invert)

q2_model <- glmmTMB(
  q2.invert ~ Sampling_Month + (1 | Site),
  family = gaussian(link = "log"),
  data = invert_hill_results_bysite_andRound
)
summary(q2_model)

##Inspect model fit
sim <- simulateResiduals(q2_model)
plot(sim)
testDispersion(sim)

##Pairwise comparisons
pairwise_q2 <- emmeans(q2_model, ~ Sampling_Month)
pairs(pairwise_q2)

##Plot of observed + predicted Hill diversity through time
#Generate prediction dataset
pred <- ggpredict(q2_model, terms = "Sampling_Month")
pred <- pred %>%
  rename(Sampling_Month = "x")

#Inverse Simpsons Plot
TemporalHillSimpsons <- ggplot(invert_hill_results_bysite_andRound, aes(Sampling_Month, q2.invert)) +
  geom_line(aes(group = Site,
                color = WetlandType),
            linewidth = 0.5,
            alpha = 1) +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  geom_point(aes(color = WetlandType),
             size = 1.5) +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  scale_x_discrete(expand = expansion(add = c(0.05, 0.05))) +
  geom_line(data = pred,
            aes(x = Sampling_Month,
                y = predicted,
                group = 1),
            inherit.aes = FALSE,
            linewidth = 1.25,
            color = "black") +
  geom_point(data = pred,
             aes(x = Sampling_Month,
                 y = predicted),
             inherit.aes = FALSE,
             size = 2.25,
             color = "black") +
  theme_minimal() +
  labs(
    x = "Sampling Month",
    y = "Inverse Simpson (q=2)",
    color = "Wetland Type"
  ) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    #legend.text = element_text(size = 8),
    #legend.title = element_text(size = 10),
    panel.grid.minor = element_blank()
  )
TemporalHillSimpsons

#ggsave("TemporalHillDiversity_ObservedandMean_Simpson.png", width = 7, height = 6, units = "in", dpi = 300)


####Biomass#####################################################################
##Data setup
#Total biomass to round/sampling month
biomass_round <- taxa.biomass.dataNEW %>%
  group_by(Round, Site) %>%
  summarise(
    total_biomass = sum(TotalBiomassPerSample, na.rm = TRUE),
    .groups = "drop"
  )
biomass_round <- biomass_round %>%
  filter(Round !=0)

biomass_round<- biomass_round %>%
  mutate(
    Round = as.factor(Round)
  )

invert_hill_results_bysite_andRound<- invert_hill_results_bysite_andRound %>%
  mutate(
    Round = as.factor(Round)
  )

#add to invert hill dataset
invert_hill_results_bysite_andRound <- invert_hill_results_bysite_andRound %>%
  left_join(
    biomass_round,
    by = c("Round", "Site")
  )


##fit model
hist(invert_hill_results_bysite_andRound$total_biomass) #right skewed so starting with gamma

biomass_gamm <- glmmTMB(
  total_biomass ~ Sampling_Month + (1 | Site),
  family = Gamma(link = "log"),
  data = invert_hill_results_bysite_andRound
)
summary(biomass_gamm)

##Inspect model fit
sim <- simulateResiduals(biomass_gamm)
plot(sim)
testDispersion(sim)

##Pairwise comparisons
pairwise_biomass <- emmeans(biomass_gamm, ~ Sampling_Month)
pairs(pairwise_biomass)

##Plot of observed + predicted Hill diversity through time
#Generate prediction dataset
pred <- ggpredict(biomass_gamm, terms = "Sampling_Month")
pred <- pred %>%
  rename(Sampling_Month = "x")

#Biomass Plot
TemporalBiomass <- ggplot(invert_hill_results_bysite_andRound, aes(Sampling_Month, total_biomass)) +
  geom_line(aes(group = Site,
                color = WetlandType),
            linewidth = 0.5,
            alpha = 1) +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  geom_point(aes(color = WetlandType),
             size = 1.5) +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  scale_x_discrete(expand = expansion(add = c(0.05, 0.05))) +
  geom_line(data = pred,
            aes(x = Sampling_Month,
                y = predicted,
                group = 1),
            inherit.aes = FALSE,
            linewidth = 1.25,
            color = "black") +
  geom_point(data = pred,
             aes(x = Sampling_Month,
                 y = predicted),
             inherit.aes = FALSE,
             size = 2.25,
             color = "black") +
  theme_minimal() +
  labs(
    x = "Sampling Month",
    y = "Biomass (ug/L)",
    color = "Wetland Type"
  ) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    #legend.text = element_text(size = 8),
    #legend.title = element_text(size = 10),
    panel.grid.minor = element_blank()
  )
TemporalBiomass

#ggsave("TemporalBiomass_ObservedandMean.png", width = 7, height = 6, units = "in", dpi = 300)


####Four-panel Temporal Diversity Figure##########################################
#install.packages("patchwork")
#install.packages("cowplot")
library(patchwork)
library(cowplot)

#make legend
legend <- cowplot::get_legend(
  TemporalBiomass +
    theme(
      legend.position = "bottom",
      legend.justification = "center"
    )
)

#plot
Temporal_Patterns <- (TemporalHillRichness | TemporalHillShannons) /
  (TemporalHillSimpsons | TemporalBiomass) +
  plot_annotation(tag_levels = "A")

#combine plot with legend
Temporal_Patterns <- cowplot::plot_grid(
  Temporal_Patterns,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.10),
  align = "v"
)
Temporal_Patterns

#ggsave("Temporal_DiversityPanel.png", width = 12, height = 8, units = "in", dpi = 300, bg = "white")
#ggsave("Temporal_DiversityPanelsmall.png", width = 8, height = 6, units = "in", dpi = 300, bg = "white")


###1B.Seasonal change in nutrients##############################################
###Models testing seasonal change in nutrients

##Data setup
#create sampling month column & define as factor
chem.hach.dataWide <- chem.hach.dataWide %>%
  mutate(
    Sampling_Month = case_when(
      Round %in% c("1") ~ "June",
      Round %in% c("2") ~ "July",
      Round %in% c("3") ~ "August"
    )
  )
chem.hach.dataWide$Sampling_Month <- factor(chem.hach.dataWide$Sampling_Month,
                                            levels = c("June", "July", "August"))

#Wetland type
chem.hach.dataWide <- chem.hach.dataWide %>%
  mutate(
    WetlandType = case_when(
      Site %in% c("DD15N", "SS", "DD8", "DD65", "LP", "WW", "FCC2",
                  "KS", "RS") ~ "CREP",
      Site %in% c("UH1", "UH2") ~ "Restored for habitat",
      Site %in% c("Bjork") ~ "Remnant"
    )
  )

####Ammonia######################################################################
#fit model
hist(chem.hach.dataWide$Ammonia)

Ammonia_model <- glmmTMB(
  Ammonia ~ Sampling_Month + (1 | Site),
  family = Gamma(link= "log"),
  data = chem.hach.dataWide
)
summary(Ammonia_model)

#Inspect model fit
sim <- simulateResiduals(Ammonia_model)
plot(sim)
testDispersion(sim)

##Pairwise comparison
pairwise_ammonia <- emmeans(Ammonia_model, ~ Sampling_Month)
pairs(pairwise_ammonia)

##Plot of observed + predicted Hill diversity through time
#Generate prediction dataset
pred <- ggpredict(Ammonia_model, terms = "Sampling_Month")
pred <- pred %>%
  rename(Sampling_Month = "x")

#Ammonia Plot
TemporalNH3_N <- ggplot(chem.hach.dataWide, aes(Sampling_Month, Ammonia)) +
  geom_line(aes(group = Site,
                color = WetlandType),
            linewidth = 0.5,
            alpha = 1) +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  geom_point(aes(color = WetlandType),
             size = 1.5) +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  scale_x_discrete(expand = expansion(add = c(0.05, 0.05))) +
  geom_line(data = pred,
            aes(x = Sampling_Month,
                y = predicted,
                group = 1),
            inherit.aes = FALSE,
            linewidth = 1.25,
            color = "black") +
  geom_point(data = pred,
             aes(x = Sampling_Month,
                 y = predicted),
             inherit.aes = FALSE,
             size = 2.25,
             color = "black") +
  theme_minimal() +
  labs(
    x = "Sampling Month",
    y = "NH3-N (mg/L)"
  ) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    #legend.text = element_text(size = 8),
    #legend.title = element_text(size = 10),
    panel.grid.minor = element_blank()
  )
TemporalNH3_N

#ggsave("TemporalAmmonia.png", width = 7, height = 5, units = "in", dpi = 300)


####Nitrate-Nitrite###############################################################
#fit model
hist(chem.hach.dataWide$Nitrate.Nitrite)

Nitrate.Nitrite_model <- glmmTMB(
  Nitrate.Nitrite ~ Sampling_Month + (1 | Site),
  family = Gamma(link= "log"),
  data = chem.hach.dataWide
)
summary(Nitrate.Nitrite_model)

##Inspect model fit
sim <- simulateResiduals(Nitrate.Nitrite_model)
plot(sim)
testDispersion(sim)

#Pairwise comparisons
pairwise_NN <- emmeans(Nitrate.Nitrite_model, ~ Sampling_Month)
pairs(pairwise_NN)

##Plot of observed + predicted nutrients through time
#create prediction plot
pred <- ggpredict(Nitrate.Nitrite_model, terms = "Sampling_Month", bias_correction = TRUE)
pred$x <- factor(pred$x)
pred <- pred %>%
  rename(Sampling_Month = "x")

#Nitrate-Nitrite Plot
TemporalNN <- ggplot(chem.hach.dataWide, aes(Sampling_Month, Nitrate.Nitrite)) +
  geom_line(aes(group = Site,
                color = WetlandType),
            linewidth = 0.5,
            alpha = 0.5) +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  geom_point(aes(color = WetlandType),
             size = 1.5) +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  scale_x_discrete(expand = expansion(add = c(0.05, 0.05))) +
  geom_line(data = pred,
            aes(x = Sampling_Month,
                y = predicted,
                group = 1),
            inherit.aes = FALSE,
            linewidth = 1.25,
            color = "black") +
  geom_point(data = pred,
             aes(x = Sampling_Month,
                 y = predicted),
             inherit.aes = FALSE,
             size = 2.25,
             color = "black") +
  theme_minimal() +
  labs(
    x = "Sampling Month",
    y = "NO3-N + NO2-N (mg/L)"
  )+
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    #legend.text = element_text(size = 8),
    #legend.title = element_text(size = 10),
    panel.grid.minor = element_blank()
  )
TemporalNN

#ggsave("TemporalNitrate.Nitrite.png", width = 7, height = 5, units = "in", dpi = 300)


####TKN#########################################################################
#fit model
hist(chem.hach.dataWide$TKN)

TKN_model <- glmmTMB(
  TKN ~ Sampling_Month + (1| Site),
  dispformula = ~ Sampling_Month, #Dispersion model; variability changes w/ predictor
  family = Gamma(link= "log"),
  data = chem.hach.dataWide
)
summary(TKN_model)

##Model Residuals
sim <- simulateResiduals(TKN_model)
plot(sim) 
testDispersion(sim)

##Pairwise comparisons
pairwise_TKN <- emmeans(TKN_model, ~ Sampling_Month)
pairs(pairwise_TKN)

##Plot of observed + predicted Hill diversity through time
#generate predictions manually due to dispersion model
TKN_model_mean <- mean(chem.hach.dataWide$TKN, na.rm = TRUE)
TKN_model_sd <- sd(chem.hach.dataWide$TKN, na.rm = TRUE)

#prediction dataframe
predict_data <- data.frame(
  Sampling_Month = factor(
    c("June", "July", "August"),
    levels = levels(chem.hach.dataWide$Sampling_Month)
  )
)

#predict from gamma+dispersion glmm
predict_data$predicted <- predict(
  TKN_model,
  newdata = predict_data,
  type = "response",
  re.form = NA
)

#TKN plot
TemporalTKN <- ggplot(chem.hach.dataWide, aes(Sampling_Month, TKN)) +
  geom_line(aes(group = Site,
                color = WetlandType),
            linewidth = 0.5,
            alpha = 1) +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  geom_point(aes(color = WetlandType),
             size = 1.5) +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  scale_x_discrete(expand = expansion(add = c(0.05, 0.05))) +
  geom_line(data = predict_data,
            aes(x = Sampling_Month,
                y = predicted,
                group = 1),
            inherit.aes = FALSE,
            linewidth = 1.25,
            color = "black") +
  geom_point(data = predict_data,
             aes(x = Sampling_Month,
                 y = predicted),
             inherit.aes = FALSE,
             size = 2.25,
             color = "black") +
  theme_minimal() +
  labs(
    x = "Sampling Month",
    y = "Total Kjeldahl Nitrogen (mg/L)"
  )+
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    #legend.text = element_text(size = 8),
    #legend.title = element_text(size = 10),
    panel.grid.minor = element_blank()
  )
TemporalTKN 

#ggsave("TemporalTKN.png", width = 7, height = 4, units = "in", dpi = 300)


####TN###########################################################################
#fit model
hist(chem.hach.dataWide$Total.Nitrogen)

TN_model <- glmmTMB(
  Total.Nitrogen ~ Sampling_Month + (1 | Site),
  family = Gamma(link= "log"), 
  data = chem.hach.dataWide
)
summary(TN_model)

##Inspect model fit
sim <- simulateResiduals(TN_model)
plot(sim)
testDispersion(sim)

##Pairwise comparisons
pairwise_TN <- emmeans(TN_model, ~ Sampling_Month)
pairs(pairwise_TN)

##Plot of observed + predicted Hill diversity through time
#prediction dataset
pred <- ggpredict(TN_model, terms = "Sampling_Month", bias_correction = TRUE)
pred <- pred %>%
  rename(Sampling_Month = "x")

#TN plot
TemporalTN <- ggplot(chem.hach.dataWide, aes(Sampling_Month, Total.Nitrogen)) +
  geom_line(aes(group = Site,
                color = WetlandType),
            linewidth = 0.5,
            alpha = 1) +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  geom_point(aes(color = WetlandType),
             size = 1.5) +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  scale_x_discrete(expand = expansion(add = c(0.05, 0.05))) +
  geom_line(data = pred,
            aes(x = Sampling_Month,
                y = predicted,
                group = 1),
            inherit.aes = FALSE,
            linewidth = 1.25,
            color = "black") +
  geom_point(data = pred,
             aes(x = Sampling_Month,
                 y = predicted),
             inherit.aes = FALSE,
             size = 2.25,
             color = "black") +
  theme_minimal() +
  labs(
    x = "Sampling Month",
    y = "Total Nitrogen (mg/L)"
  )+
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    #legend.text = element_text(size = 8),
    #legend.title = element_text(size = 10)
  )
TemporalTN

#ggsave("TemporalTN.png", width = 7, height = 4, units = "in", dpi = 300)


####TP###########################################################################
#fit model
hist(chem.hach.dataWide$Total.Phosphorus)

TP_model <- glmmTMB(
  Total.Phosphorus ~ Sampling_Month + (1 | Site),
  family = Gamma(link= "log"), 
  data = chem.hach.dataWide
)
summary(TP_model)

#Inspect model fit
sim <- simulateResiduals(TP_model)
plot(sim)
testDispersion(sim)

##Pairwise comparisons
pairwise_TP <- emmeans(TP_model, ~ Sampling_Month)
pairs(pairwise_TP)

##Plot of observed + predicted Hill diversity through time
#prediction dataset
pred <- ggpredict(TP_model, terms = "Sampling_Month", bias_correction = TRUE)
pred <- pred %>%
  rename(Sampling_Month = "x")

#TP plot
TemporalTP <- ggplot(chem.hach.dataWide, aes(Sampling_Month, Total.Phosphorus)) +
  geom_line(aes(group = Site,
                color = WetlandType),
            linewidth = 0.5,
            alpha = 1) +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  geom_point(aes(color = WetlandType),
             size = 1.5) +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  scale_x_discrete(expand = expansion(add = c(0.05, 0.05))) +
  geom_line(data = pred,
            aes(x = Sampling_Month,
                y = predicted,
                group = 1),
            inherit.aes = FALSE,
            linewidth = 1.25,
            color = "black") +
  geom_point(data = pred,
             aes(x = Sampling_Month,
                 y = predicted),
             inherit.aes = FALSE,
             size = 2.25,
             color = "black") +
  theme_minimal() +
  labs(
    x = "Sampling Month",
    y = "Total Phosphorus (mg/L)"
  )+
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    #legend.text = element_text(size = 8),
    #legend.title = element_text(size = 10)
  )
TemporalTP

#ggsave("TemporalTP.png", width = 7, height = 4, units = "in", dpi = 300)


####N:P###########################################################################
#fit model
hist(chem.hach.dataWide$NP.ratio)

NP_model <- glmmTMB(
  NP.ratio ~ Sampling_Month + (1 | Site),
  family = Gamma(link= "log"), 
  data = chem.hach.dataWide
)
summary(NP_model)

##Inspect model fit
sim <- simulateResiduals(NP_model)
plot(sim)
testDispersion(sim)

##Pairwise comparisons
pairwise_NP <- emmeans(NP_model, ~ Sampling_Month)
pairs(pairwise_NP)

##Plot of observed + predicted Hill diversity through time
#prediction dataset
pred <- ggpredict(NP_model, terms = "Sampling_Month", bias_correction = TRUE)
pred <- pred %>%
  rename(Sampling_Month = "x")

#TN:TP plot
TemporalN_P <- ggplot(chem.hach.dataWide, aes(Sampling_Month, NP.ratio)) +
  geom_line(aes(group = Site,
                color = WetlandType),
            linewidth = 0.5,
            alpha = 1) +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  geom_point(aes(color = WetlandType),
             size = 1.5) +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  scale_x_discrete(expand = expansion(add = c(0.05, 0.05))) +
  geom_line(data = pred,
            aes(x = Sampling_Month,
                y = predicted,
                group = 1),
            inherit.aes = FALSE,
            linewidth = 1.25,
            color = "black") +
  geom_point(data = pred,
             aes(x = Sampling_Month,
                 y = predicted),
             inherit.aes = FALSE,
             size = 2.25,
             color = "black") +
  theme_minimal() +
  labs(
    x = "Sampling Month",
    y = "TN: TP (mg/L)"
  )+
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  )
TemporalN_P

#ggsave("TemporalNP.png", width = 7, height = 4, units = "in", dpi = 300)


####Six-panel Temporal Nutrient Figure##########################################
#install.packages("patchwork")
#install.packages("cowplot")
library(patchwork)
library(cowplot)

#create legend
legend <- cowplot::get_legend(
  TemporalN_P +
    theme(
      legend.position = "bottom",
      legend.justification = "center"
    )
)

#create plots
plots <- list(
  TemporalNH3_N + theme(legend.position = "none"),
  TemporalNN + theme(legend.position = "none"),
  TemporalTKN + theme(legend.position = "none"),
  TemporalTN + theme(legend.position = "none"),
  TemporalTP + theme(legend.position = "none"),
  TemporalN_P + theme(legend.position = "none")
)

row1 <- cowplot::plot_grid(
  plots[[1]], plots[[2]],
  labels = c("A", "B"),
  ncol = 2,
  label_size = 14,
  label_fontface = "bold",
  label_x = 0.02,
  label_y = 0.98,
  hjust = 0,
  vjust = 1
)

row2 <- cowplot::plot_grid(
  plots[[3]], plots[[4]],
  labels = c("C", "D"),
  ncol = 2,
  label_size = 14,
  label_fontface = "bold",
  label_x = 0.02,
  label_y = 0.98,
  hjust = 0,
  vjust = 1
)

row3 <- cowplot::plot_grid(
  plots[[5]], plots[[6]],
  labels = c("E", "F"),
  ncol = 2,
  label_size = 14,
  label_fontface = "bold",
  label_x = 0.02,
  label_y = 0.98,
  hjust = 0,
  vjust = 1
)

#combine legend with plots
Temporal_Nutrient_Patterns <- cowplot::plot_grid(
  row1,
  row2,
  row3,
  legend,
  ncol = 1,
  rel_heights = c(1, 1, 1, 0.12)
)
Temporal_Nutrient_Patterns

#ggsave("Temporal_NutrientsPanel.png", width = 9, height = 10, units = "in", dpi = 300, bg = "white")


###1C. Hill Diversity vs nutrient, accounting for season########################
##Models testing temporal diversity association to nutrients

##Data setup
#combine diversity & nutrient datasets
diversity_nutrient_data <- invert_hill_results_bysite_andRound %>%
  left_join(
    chem.hach.dataWide %>%
      select(Site, Round,
             Ammonia,
             Nitrate.Nitrite,
             TKN,
             Total.Nitrogen,
             Total.Phosphorus,
             NP.ratio),
    by = c("Site", "Round")
  )

diversity_nutrient_data <- diversity_nutrient_data %>%
  mutate(
    Sampling_Month = case_when(
      Round %in% c("1") ~ "June",
      Round %in% c("2") ~ "July",
      Round %in% c("3") ~ "August"
    )
  )
diversity_nutrient_data$Sampling_Month <- factor(diversity_nutrient_data$Sampling_Month,
                                            levels = c("June", "July", "August"))


#Log nutrients before scaling
diversity_nutrient_data$TN_logged <- as.numeric(log(diversity_nutrient_data$Total.Nitrogen))
diversity_nutrient_data$TP_logged <- as.numeric(log(diversity_nutrient_data$Total.Phosphorus))
diversity_nutrient_data$NP_logged <- as.numeric(log(diversity_nutrient_data$NP.ratio))

#Scale nutrients for all models
diversity_nutrient_data$TN_scaled <- as.numeric(scale(diversity_nutrient_data$TN_logged))
diversity_nutrient_data$TP_scaled <- as.numeric(scale(diversity_nutrient_data$TP_logged))
diversity_nutrient_data$NP_scaled <- as.numeric(scale(diversity_nutrient_data$NP_logged))

#Save as dataframe
#write.csv(diversity_nutrient_data, "diversity_nutrient_data.csv")

#####Richness###################################################################
######TN########################################################################
#fit model
q0_TN_model_quad_GP <- glmmTMB(
  q0.invert ~ Sampling_Month + TN_scaled + I(TN_scaled^2) + (1 | Site),
  dispformula = ~ Sampling_Month,
  family = genpois(link = "log"),
  data = diversity_nutrient_data
)
summary(q0_TN_model_quad_GP)

#Inspect model fit
sim_quad_GP <- simulateResiduals(q0_TN_model_quad_GP)
plot(sim_quad_GP)
testDispersion(sim_quad_GP)

##Predicted curve plot
q0_TN_mean <- mean(diversity_nutrient_data$Total.Nitrogen, na.rm = TRUE)
q0_TN_sd <- sd(diversity_nutrient_data$Total.Nitrogen, na.rm = TRUE)

#prediction dataframe
predict_data <- data.frame(
  TN_mg_l = seq(
    from = max(0, min(diversity_nutrient_data$Total.Nitrogen, na.rm = TRUE)),
    to = max(diversity_nutrient_data$Total.Nitrogen, na.rm = TRUE),
    length.out = 200
  ),
  Sampling_Month = factor(
    "June",
    levels = levels(diversity_nutrient_data$Sampling_Month)
  )
)

#recalculate tn_scaled
predict_data$TN_scaled <- (predict_data$TN_mg_l - q0_TN_mean) / q0_TN_sd

#predict from Poisson glmm
pred <- predict(
  q0_TN_model_quad_GP, #change model used in plot here
  newdata = predict_data,
  type = "link",
  se.fit = TRUE,
  re.form = NA
)

predict_data$predicted <- exp(pred$fit)
predict_data$conf.low <- exp(pred$fit - 1.96 * pred$se.fit)
predict_data$conf.high <- exp(pred$fit + 1.96 * pred$se.fit)

#Richness ~ TN plot
Rich_TN_plot <- ggplot() +
  geom_point(
    data = diversity_nutrient_data,
    aes(x = Total.Nitrogen,
        y = q0.invert),
    alpha = 0.3
  ) +
  geom_ribbon(
    data = predict_data,
    aes(x = TN_mg_l,
        ymin = conf.low,
        ymax = conf.high),
    alpha = 0.2,
    fill = "steelblue"
  ) +
  geom_line(
    data = predict_data,
    aes(x = TN_mg_l,
        y = predicted),
    linewidth = 1.2,
    color = "steelblue"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.03))
  ) +
  coord_cartesian(
    xlim = c(0, max(diversity_nutrient_data$Total.Nitrogen, na.rm = TRUE))
  ) +
  labs(
    x = "Total nitrogen (mg/L)",
    y = "q=0"
  ) +
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )
Rich_TN_plot

#ggsave("Temporalq0_TN.png", width = 8, height = 6, units = "in", dpi = 300)
#ggsave("Temporalq0_TN_quad.poisson.png", width = 8, height = 6, units = "in", dpi = 300)


######TP########################################################################
#fit model
q0_TP_model_LGP <- glmmTMB(
  q0.invert ~ Sampling_Month + TP_scaled + (1 | Site),
  family = genpois(link = "log"),
  data = diversity_nutrient_data
)
summary(q0_TP_model_LGP)

#Inspect model fit
sim_q0_TP_LGP <- simulateResiduals(q0_TP_model_LGP)
plot(sim_q0_TP_LGP)
testDispersion(sim_q0_TP_LGP)

##Predicted curve plot
q0_TP_mean <- mean(diversity_nutrient_data$Total.Phosphorus, na.rm = TRUE)
q0_TP_sd <- sd(diversity_nutrient_data$Total.Phosphorus, na.rm = TRUE)

#prediction dataframe
predict_data <- data.frame(
  TP_mg_l = seq(
    from = max(0, min(diversity_nutrient_data$Total.Phosphorus, na.rm = TRUE)),
    to = max(diversity_nutrient_data$Total.Phosphorus, na.rm = TRUE),
    length.out = 200
  ),
  Sampling_Month = factor(
  "June",
  levels = levels(diversity_nutrient_data$Sampling_Month)
  )
)

#recalculate tp_scaled
predict_data$TP_scaled <- (predict_data$TP_mg_l - q0_TP_mean) / q0_TP_sd

#predict from model
pred <- predict(
  q0_TP_model_LGP,
  newdata = predict_data,
  type = "link",
  se.fit = TRUE,
  re.form = NA
)

predict_data$predicted <- exp(pred$fit)
predict_data$conf.low <- exp(pred$fit - 1.96 * pred$se.fit)
predict_data$conf.high <- exp(pred$fit + 1.96 * pred$se.fit)

#Richness ~ TP plot
Rich_TP_plot <- ggplot() +
  geom_point(
    data = diversity_nutrient_data,
    aes(x = Total.Phosphorus,
        y = q0.invert),
    alpha = 0.3
  ) +
  geom_ribbon(
    data = predict_data,
    aes(x = TP_mg_l,
        ymin = conf.low,
        ymax = conf.high),
    alpha = 0.2,
    fill = "steelblue"
  ) +
  geom_line(
    data = predict_data,
    aes(x = TP_mg_l,
        y = predicted),
    linewidth = 1.2,
    color = "steelblue"
  ) +
  labs(
    x = "Total phosphorus (mg/L)",
    y = "q=0"
  ) +
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )
Rich_TP_plot

#ggsave("Temporalq0_TP.genP.png", width = 8, height = 6, units = "in", dpi = 300) #gen poisson


######TN:TP#####################################################################
#fit model
##needs observation-level random effect (OLRE)
diversity_nutrient_data$obs_id <- factor(seq_len(nrow(diversity_nutrient_data)))

q0_NP_model_LGP <- glmmTMB(
  q0.invert ~ Sampling_Month + NP_scaled + (1 | Site) + (1 | obs_id),
  dispformula = ~ Sampling_Month,
  family = genpois(link = "log"),
  data = diversity_nutrient_data
)
summary(q0_NP_model_LGP)

#Inspect Model Fit
sim_q0_NP_LGP <- simulateResiduals(q0_NP_model_LGP)
plot(sim_q0_NP_LGP)
testDispersion(sim_q0_NP_LGP)

##Predicted curve plot
q0_NP_mean <- mean(diversity_nutrient_data$NP.ratio, na.rm = TRUE)
q0_NP_sd <- sd(diversity_nutrient_data$NP.ratio, na.rm = TRUE)

#prediction dataframe
predict_data <- data.frame(
  NP_mg_l = seq(
    from = max(0, min(diversity_nutrient_data$NP.ratio, na.rm = TRUE)),
    to = max(diversity_nutrient_data$NP.ratio, na.rm = TRUE),
    length.out = 200
  ),
  Sampling_Month = factor(
  "June",
  levels = levels(diversity_nutrient_data$Sampling_Month)
  )
)

#recalculate np_scaled
predict_data$NP_scaled <- (predict_data$NP_mg_l - q0_NP_mean) / q0_NP_sd


#predict from general Poisson glmm
pred <- predict(
  q0_NP_model_LGP,
  newdata = predict_data,
  type = "link",
  se.fit = TRUE,
  re.form = NA
)

predict_data$predicted <- exp(pred$fit)
predict_data$conf.low <- exp(pred$fit - 1.96 * pred$se.fit)
predict_data$conf.high <- exp(pred$fit + 1.96 * pred$se.fit)

#Richness ~ NP plot
Rich_NP_plot <- ggplot() +
  geom_point(
    data = diversity_nutrient_data,
    aes(x = NP.ratio,
        y = q0.invert),
    alpha = 0.3
  ) +
  geom_ribbon(
    data = predict_data,
    aes(x = NP_mg_l,
        ymin = conf.low,
        ymax = conf.high),
    alpha = 0.2,
    fill = "steelblue"
  ) +
  geom_line(
    data = predict_data,
    aes(x = NP_mg_l,
        y = predicted),
    linewidth = 1.2,
    color = "steelblue"
  ) +
  labs(
    x = "Total nitrogen:Total phosphorus (mg/L)",
    y = "q=0"
  ) +
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )
Rich_NP_plot

#ggsave("Temporalq0_NP.png", width = 8, height = 6, units = "in", dpi = 300)


#####Shannons###################################################################
######TN########################################################################
ggplot(diversity_nutrient_data, aes(x = Total.Nitrogen, y = q1.invert)) +
  geom_point(aes(color = Sampling_Month), alpha = 0.7)

##fit model 
q1_TN_model <- glmmTMB(
  q1.invert ~ Sampling_Month + TN_scaled + (1 | Site),
  family = gaussian(link= "log"),
  data = diversity_nutrient_data
)
summary(q1_TN_model)

#Inspect model fit
sim <- simulateResiduals(q1_TN_model)
plot(sim)
testDispersion(sim)

##Predicted curve plot
q1_TN_mean <- mean(diversity_nutrient_data$Total.Nitrogen, na.rm = TRUE)
q1_TN_sd <- sd(diversity_nutrient_data$Total.Nitrogen, na.rm = TRUE)

#prediction dataframe
predict_data <- data.frame(
  TN_mg_l = seq(
    from = max(0, min(diversity_nutrient_data$Total.Nitrogen, na.rm = TRUE)),
    to = max(diversity_nutrient_data$Total.Nitrogen, na.rm = TRUE),
    length.out = 200
  ),
  Sampling_Month = factor(
    "June",
    levels = levels(diversity_nutrient_data$Sampling_Month)
  )
)

#recalculate tn_scaled
predict_data$TN_scaled <- (predict_data$TN_mg_l - q1_TN_mean) / q1_TN_sd


#predict from model
pred <- predict(
  q1_TN_model,
  newdata = predict_data,
  type = "link",
  se.fit = TRUE,
  re.form = NA
)

predict_data$predicted <- exp(pred$fit)
predict_data$conf.low <- exp(pred$fit - 1.96 * pred$se.fit)
predict_data$conf.high <- exp(pred$fit + 1.96 * pred$se.fit)

#Exponential Shannon ~ TN plot
Shan_TN_plot <- ggplot() +
  geom_point(
    data = diversity_nutrient_data,
    aes(x = Total.Nitrogen,
        y = q1.invert),
    alpha = 0.3
  ) +
  geom_ribbon(
    data = predict_data,
    aes(x = TN_mg_l,
        ymin = conf.low,
        ymax = conf.high),
    alpha = 0.2,
    fill = "steelblue"
  ) +
  geom_line(
    data = predict_data,
    aes(x = TN_mg_l,
        y = predicted),
    linewidth = 1.2,
    color = "steelblue"
  ) +
  labs(
    x = "Total nitrogen (mg/L)",
    y = "q=1"
  ) +
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )
Shan_TN_plot

#ggsave("Temporalq1_TN.png", width = 7, height = 4, units = "in", dpi = 300)


######TP########################################################################
ggplot(diversity_nutrient_data, aes(x = Total.Phosphorus, y = q1.invert)) +
  geom_point(aes(color = Sampling_Month), alpha = 0.7)

#fit model 
q1_TP_model <- glmmTMB(
  q1.invert ~ Sampling_Month + TP_scaled + (1 | Site),
  family = gaussian(link= "log"),
  data = diversity_nutrient_data
)
summary(q1_TP_model)

##Verify model fit
sim <- simulateResiduals(q1_TP_model)
plot(sim)
testDispersion(sim)

##Predicted curve plot
q1_TP_mean <- mean(diversity_nutrient_data$Total.Phosphorus, na.rm = TRUE)
q1_TP_sd <- sd(diversity_nutrient_data$Total.Phosphorus, na.rm = TRUE)

#prediction dataframe
predict_data <- data.frame(
  TP_mg_l = seq(
    from = max(0, min(diversity_nutrient_data$Total.Phosphorus, na.rm = TRUE)),
    to = max(diversity_nutrient_data$Total.Phosphorus, na.rm = TRUE),
    length.out = 200
  ),
  Sampling_Month = factor(
    "June",
    levels = levels(diversity_nutrient_data$Sampling_Month)
  )
)

#recalculate tp_scaled
predict_data$TP_scaled <- (predict_data$TP_mg_l - q1_TP_mean) / q1_TP_sd

#predict from model
pred <- predict(
  q1_TP_model,
  newdata = predict_data,
  type = "link",
  se.fit = TRUE,
  re.form = NA
)

predict_data$predicted <- exp(pred$fit)
predict_data$conf.low <- exp(pred$fit - 1.96 * pred$se.fit)
predict_data$conf.high <- exp(pred$fit + 1.96 * pred$se.fit)

#Exponential shannon ~ TP plot
Shan_TP_plot <- ggplot() +
  geom_point(
    data = diversity_nutrient_data,
    aes(x = Total.Phosphorus,
        y = q1.invert),
    alpha = 0.3
  ) +
  geom_ribbon(
    data = predict_data,
    aes(x = TP_mg_l,
        ymin = conf.low,
        ymax = conf.high),
    alpha = 0.2,
    fill = "steelblue"
  ) +
  geom_line(
    data = predict_data,
    aes(x = TP_mg_l,
        y = predicted),
    linewidth = 1.2,
    color = "steelblue"
  ) +
  labs(
    x = "Total phosphorus (mg/L)",
    y = "q=1"
  ) +
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )
Shan_TP_plot

#ggsave("Temporalq1_TP.png", width = 7, height = 4, units = "in", dpi = 300)


######TN:TP#####################################################################
ggplot(diversity_nutrient_data, aes(x = NP.ratio, y = q1.invert)) +
  geom_point(aes(color = Sampling_Month), alpha = 0.7)

#fit model 
q1_NP_model <- glmmTMB(
  q1.invert ~ Sampling_Month + NP_scaled + (1 | Site),
  family = gaussian(link= "log"),
  data = diversity_nutrient_data
)
summary(q1_NP_model)

#verify model fit
sim <- simulateResiduals(q1_NP_model)
plot(sim)
testDispersion(sim)

##Predicted curve plot
q1_NP_mean <- mean(diversity_nutrient_data$NP.ratio, na.rm = TRUE)
q1_NP_sd <- sd(diversity_nutrient_data$NP.ratio, na.rm = TRUE)

#prediction dataframe
predict_data <- data.frame(
  NP_mg_l = seq(
    from = max(0, min(diversity_nutrient_data$NP.ratio, na.rm = TRUE)),
    to = max(diversity_nutrient_data$NP.ratio, na.rm = TRUE),
    length.out = 200
  ),
  Sampling_Month = factor(
    "June",
    levels = levels(diversity_nutrient_data$Sampling_Month)
  )
)

#recalculate np_scaled
predict_data$NP_scaled <- (predict_data$NP_mg_l - q1_NP_mean) / q1_NP_sd

#predict from model
pred <- predict(
  q1_NP_model,
  newdata = predict_data,
  type = "link",
  se.fit = TRUE,
  re.form = NA
)

predict_data$predicted <- exp(pred$fit)
predict_data$conf.low <- exp(pred$fit - 1.96 * pred$se.fit)
predict_data$conf.high <- exp(pred$fit + 1.96 * pred$se.fit)

#Exponential shannon ~ TN:TP plot
Shan_NP_plot <- ggplot() +
  geom_point(
    data = diversity_nutrient_data,
    aes(x = NP.ratio,
        y = q1.invert),
    alpha = 0.3
  ) +
  geom_ribbon(
    data = predict_data,
    aes(x = NP_mg_l,
        ymin = conf.low,
        ymax = conf.high),
    alpha = 0.2,
    fill = "steelblue"
  ) +
  geom_line(
    data = predict_data,
    aes(x = NP_mg_l,
        y = predicted),
    linewidth = 1.2,
    color = "steelblue"
  ) +
  labs(
    x = "Total nitrogen:Total phosphorus (mg/L)",
    y = "q=1"
  ) +
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )
Shan_NP_plot

#ggsave("Temporalq1_NP.png", width = 7, height = 4, units = "in", dpi = 300)


#####Inverse Simpson############################################################
######TN########################################################################
##fit model 
q2_TN_model <- glmmTMB(
  q2.invert ~ Sampling_Month + TN_scaled + (1 | Site),
  family = gaussian(link= "log"),
  data = diversity_nutrient_data
)
summary(q2_TN_model)

#verify model fit
sim <- simulateResiduals(q2_TN_model)
plot(sim)
testDispersion(sim)

##Predicted curve plot
q2_TN_mean <- mean(diversity_nutrient_data$Total.Nitrogen, na.rm = TRUE)
q2_TN_sd <- sd(diversity_nutrient_data$Total.Nitrogen, na.rm = TRUE)

#prediction dataframe
predict_data <- data.frame(
  TN_mg_l = seq(
    from = max(0, min(diversity_nutrient_data$Total.Nitrogen, na.rm = TRUE)),
    to = max(diversity_nutrient_data$Total.Nitrogen, na.rm = TRUE),
    length.out = 200
  ),
  Sampling_Month = factor(
    "June",
    levels = levels(diversity_nutrient_data$Sampling_Month)
  )
)

#recalculate tn_scaled
predict_data$TN_scaled <- (predict_data$TN_mg_l - q2_TN_mean) / q2_TN_sd

#predict from model
pred <- predict(
  q2_TN_model,
  newdata = predict_data,
  type = "link",
  se.fit = TRUE,
  re.form = NA
)

predict_data$predicted <- exp(pred$fit)
predict_data$conf.low <- exp(pred$fit - 1.96 * pred$se.fit)
predict_data$conf.high <- exp(pred$fit + 1.96 * pred$se.fit)

#Inverse Simpson ~ TN plot
Sim_TN_Plot <- ggplot() +
  geom_point(
    data = diversity_nutrient_data,
    aes(x = Total.Nitrogen,
        y = q2.invert),
    alpha = 0.3
  ) +
  geom_ribbon(
    data = predict_data,
    aes(x = TN_mg_l,
        ymin = conf.low,
        ymax = conf.high),
    alpha = 0.2,
    fill = "steelblue"
  ) +
  geom_line(
    data = predict_data,
    aes(x = TN_mg_l,
        y = predicted),
    linewidth = 1.2,
    color = "steelblue"
  ) +
  labs(
    x = "Total nitrogen (mg/L)",
    y = "q=2"
  ) +
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )
Sim_TN_Plot

#ggsave("Temporalq2_TN.png", width = 7, height = 4, units = "in", dpi = 300)


######TP########################################################################
#fit model 
q2_TP_model <- glmmTMB(
  q2.invert ~ Sampling_Month + TP_scaled + (1 | Site),
  family = gaussian(link= "log"),
  data = diversity_nutrient_data
)
summary(q2_TP_model)

#verify model fit
sim <- simulateResiduals(q2_TP_model)
plot(sim)
testDispersion(sim)

##Predicted curve plot
q2_TP_mean <- mean(diversity_nutrient_data$Total.Phosphorus, na.rm = TRUE)
q2_TP_sd <- sd(diversity_nutrient_data$Total.Phosphorus, na.rm = TRUE)

#prediction dataframe
predict_data <- data.frame(
  TP_mg_l = seq(
    from = max(0, min(diversity_nutrient_data$Total.Phosphorus, na.rm = TRUE)),
    to = max(diversity_nutrient_data$Total.Phosphorus, na.rm = TRUE),
    length.out = 200
  ),
  Sampling_Month = factor(
    "June",
    levels = levels(diversity_nutrient_data$Sampling_Month)
  )
)

#recalculate tp_scaled
predict_data$TP_scaled <- (predict_data$TP_mg_l - q2_TP_mean) / q2_TP_sd

#predict from model
pred <- predict(
  q2_TP_model,
  newdata = predict_data,
  type = "link",
  se.fit = TRUE,
  re.form = NA
)

predict_data$predicted <- exp(pred$fit)
predict_data$conf.low <- exp(pred$fit - 1.96 * pred$se.fit)
predict_data$conf.high <- exp(pred$fit + 1.96 * pred$se.fit)

#Inverse Simpson ~ TP plot
Sim_TP_Plot <- ggplot() +
  geom_point(
    data = diversity_nutrient_data,
    aes(x = Total.Phosphorus,
        y = q2.invert),
    alpha = 0.3
  ) +
  geom_ribbon(
    data = predict_data,
    aes(x = TP_mg_l,
        ymin = conf.low,
        ymax = conf.high),
    alpha = 0.2,
    fill = "steelblue"
  ) +
  geom_line(
    data = predict_data,
    aes(x = TP_mg_l,
        y = predicted),
    linewidth = 1.2,
    color = "steelblue"
  ) +
  labs(
    x = "Total phosphorus (mg/L)",
    y = "q=2"
  ) +
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )
Sim_TP_Plot 

#ggsave("Temporalq2_TP.png", width = 7, height = 4, units = "in", dpi = 300)


######TN:TP#####################################################################
#fit model 
q2_NP_model <- glmmTMB(
  q2.invert ~ Sampling_Month + NP_scaled + (1 | Site),
  family = gaussian(link= "log"),
  data = diversity_nutrient_data
)
summary(q2_NP_model)

#verify model fit
sim <- simulateResiduals(q2_NP_model)
plot(sim)
testDispersion(sim)

##Predicted curve plot
q2_NP_mean <- mean(diversity_nutrient_data$NP.ratio, na.rm = TRUE)
q2_NP_sd <- sd(diversity_nutrient_data$NP.ratio, na.rm = TRUE)

#prediction dataframe
predict_data <- data.frame(
  NP_mg_l = seq(
    from = max(0, min(diversity_nutrient_data$NP.ratio, na.rm = TRUE)),
    to = max(diversity_nutrient_data$NP.ratio, na.rm = TRUE),
    length.out = 200
  ),
  Sampling_Month = factor(
    "June",
    levels = levels(diversity_nutrient_data$Sampling_Month)
  )
)

#recalculate np_scaled
predict_data$NP_scaled <- (predict_data$NP_mg_l - q2_NP_mean) / q2_NP_sd

#predict from model
pred <- predict(
  q2_NP_model,
  newdata = predict_data,
  type = "link",
  se.fit = TRUE,
  re.form = NA
)

predict_data$predicted <- exp(pred$fit)
predict_data$conf.low <- exp(pred$fit - 1.96 * pred$se.fit)
predict_data$conf.high <- exp(pred$fit + 1.96 * pred$se.fit)

#Inverse Simpson ~ TN:TP plot
Sim_NP_Plot <- ggplot() +
  geom_point(
    data = diversity_nutrient_data,
    aes(x = NP.ratio,
        y = q2.invert),
    alpha = 0.3
  ) +
  geom_ribbon(
    data = predict_data,
    aes(x = NP_mg_l,
        ymin = conf.low,
        ymax = conf.high),
    alpha = 0.2,
    fill = "steelblue"
  ) +
  geom_line(
    data = predict_data,
    aes(x = NP_mg_l,
        y = predicted),
    linewidth = 1.2,
    color = "steelblue"
  ) +
  labs(
    x = "Total nitrogen:Total phosphorus (mg/L)",
    y = "q=2"
  ) +
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )
Sim_NP_Plot

#ggsave("Temporalq2_NP.png", width = 7, height = 4, units = "in", dpi = 300)


####Biomass#####################################################################
##add total biomass to dataset
biomass_round <- taxa.biomass.dataNEW %>%
  group_by(Round, Site) %>%
  summarise(
    total_biomass = sum(TotalBiomassPerSample, na.rm = TRUE),
    .groups = "drop"
  )
biomass_round <- biomass_round %>%
  filter(Round !=0)

biomass_round<- biomass_round %>%
  mutate(
    Round = as.factor(Round)
  )

diversity_nutrient_data<- diversity_nutrient_data %>%
  mutate(
    Round = as.factor(Round)
  )

diversity_nutrient_data <- diversity_nutrient_data%>%
  left_join(
    biomass_round,
    by = c("Round", "Site")
  )

######TN########################################################################
#fit model
biomass_TN_model_quad <- glmmTMB(
  total_biomass ~ Sampling_Month + TN_scaled + I(TN_scaled^2) + (1 | Site),
  family = Gamma(link = "log"),
  data = diversity_nutrient_data
)
summary(biomass_TN_model_quad)

#Verify Model Fit
sim_biomass_TN_quad <- simulateResiduals(biomass_TN_model_quad)
plot(sim_biomass_TN_quad)
testDispersion(sim_biomass_TN_quad)

##Predicted curve plot
bio_TN_mean <- mean(diversity_nutrient_data$Total.Nitrogen, na.rm = TRUE)
bio_TN_sd <- sd(diversity_nutrient_data$Total.Nitrogen, na.rm = TRUE)

#prediction dataframe
predict_data <- data.frame(
  TN_mg_l = seq(
    from = max(0, min(diversity_nutrient_data$Total.Nitrogen, na.rm = TRUE)),
    to = max(diversity_nutrient_data$Total.Nitrogen, na.rm = TRUE),
    length.out = 200
  ),
  Sampling_Month = factor(
    "June",
    levels = levels(diversity_nutrient_data$Sampling_Month)
  )
)

#recalculate tn_scaled
predict_data$TN_scaled <- (predict_data$TN_mg_l - bio_TN_mean) / bio_TN_sd

#predict from Poisson glmm
pred <- predict(
  biomass_TN_model_quad, #change model used in plot here
  newdata = predict_data,
  type = "link",
  se.fit = TRUE,
  re.form = NA
)

predict_data$predicted <- exp(pred$fit)
predict_data$conf.low <- exp(pred$fit - 1.96 * pred$se.fit)
predict_data$conf.high <- exp(pred$fit + 1.96 * pred$se.fit)

#Biomass ~ TN plot
Biomass_TN_plot <- ggplot() +
  geom_point(
    data = diversity_nutrient_data,
    aes(x = Total.Nitrogen,
        y = total_biomass),
    alpha = 0.3
  ) +
  geom_ribbon(
    data = predict_data,
    aes(x = TN_mg_l,
        ymin = conf.low,
        ymax = conf.high),
    alpha = 0.2,
    fill = "steelblue"
  ) +
  geom_line(
    data = predict_data,
    aes(x = TN_mg_l,
        y = predicted),
    linewidth = 1.2,
    color = "steelblue"
  ) +
  labs(
    x = "Total nitrogen (mg/L)",
    y = "Total biomass (ug/L)"
  ) +
  annotate("text",
           x = -Inf,
           y = Inf,
           label = "*",
           hjust = -0.5,
           vjust = 1.0,
           size = 10) +
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )
Biomass_TN_plot

#ggsave("Biomass_TN.png", width = 8, height = 6, units = "in", dpi = 300)


######TP########################################################################
#fit model
biomass_TP_model <- glmmTMB(
  total_biomass ~ Sampling_Month + TP_scaled + (1 | Site),
  dispformula = ~ Sampling_Month,
  family = Gamma(link = "log"),
  data = diversity_nutrient_data
)
summary(biomass_TP_model)

#Verify Model Fit
sim_biomass_TP <- simulateResiduals(biomass_TP_model)
plot(sim_biomass_TP)
testDispersion(sim_biomass_TP)

##Predicted curve plot
bio_TP_mean <- mean(diversity_nutrient_data$Total.Phosphorus, na.rm = TRUE)
bio_TP_sd <- sd(diversity_nutrient_data$Total.Phosphorus, na.rm = TRUE)

#prediction dataframe
predict_data <- data.frame(
  TP_mg_l = seq(
    from = max(0, min(diversity_nutrient_data$Total.Phosphorus, na.rm = TRUE)),
    to = max(diversity_nutrient_data$Total.Phosphorus, na.rm = TRUE),
    length.out = 200
  ),
  Sampling_Month = factor(
    "June",
    levels = levels(diversity_nutrient_data$Sampling_Month)
  )
)

#recalculate tp_scaled
predict_data$TP_scaled <- (predict_data$TP_mg_l - bio_TP_mean) / bio_TP_sd

#predict from model
pred <- predict(
  biomass_TP_model,
  newdata = predict_data,
  type = "link",
  se.fit = TRUE,
  re.form = NA
)

predict_data$predicted <- exp(pred$fit)
predict_data$conf.low <- exp(pred$fit - 1.96 * pred$se.fit)
predict_data$conf.high <- exp(pred$fit + 1.96 * pred$se.fit)

#Biomass ~ TP plot
Biomass_TP_plot <- ggplot() +
  geom_point(
    data = diversity_nutrient_data,
    aes(x = Total.Phosphorus,
        y = total_biomass),
    alpha = 0.3
  ) +
  geom_ribbon(
    data = predict_data,
    aes(x = TP_mg_l,
        ymin = conf.low,
        ymax = conf.high),
    alpha = 0.2,
    fill = "steelblue"
  ) +
  geom_line(
    data = predict_data,
    aes(x = TP_mg_l,
        y = predicted),
    linewidth = 1.2,
    color = "steelblue"
  ) +
  labs(
    x = "Total phosphorus (mg/L)",
    y = "Total biomass (ug/L)"
  ) +
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )
Biomass_TP_plot

#ggsave("Biomass_TP.png", width = 7, height = 4, units = "in", dpi = 300)


######TN:TP#####################################################################
#fit model
biomass_NP_model <- glmmTMB(
  total_biomass ~ Sampling_Month + NP_scaled + (1 | Site),
  family = Gamma(link = "log"),
  data = diversity_nutrient_data
)
summary(biomass_NP_model)

#Verify Model Fit
sim_biomass_NP <- simulateResiduals(biomass_NP_model)
plot(sim_biomass_NP)
testDispersion(sim_biomass_NP)

##Predicted curve plot
bio_NP_mean <- mean(diversity_nutrient_data$NP.ratio, na.rm = TRUE)
bio_NP_sd <- sd(diversity_nutrient_data$NP.ratio, na.rm = TRUE)

#prediction dataframe
predict_data <- data.frame(
  NP_mg_l = seq(
    from = max(0, min(diversity_nutrient_data$NP.ratio, na.rm = TRUE)),
    to = max(diversity_nutrient_data$NP.ratio, na.rm = TRUE),
    length.out = 200
  ),
  Sampling_Month = factor(
    "June",
    levels = levels(diversity_nutrient_data$Sampling_Month)
  )
)

#recalculate np_scaled
predict_data$NP_scaled <- (predict_data$NP_mg_l - bio_NP_mean) / bio_NP_sd

#predict from model
pred <- predict(
  biomass_NP_model,
  newdata = predict_data,
  type = "link",
  se.fit = TRUE,
  re.form = NA
)

predict_data$predicted <- exp(pred$fit)
predict_data$conf.low <- exp(pred$fit - 1.96 * pred$se.fit)
predict_data$conf.high <- exp(pred$fit + 1.96 * pred$se.fit)

#Biomass ~ TN:TP plot
Biomass_NP_plot <- ggplot() +
  geom_point(
    data = diversity_nutrient_data,
    aes(x = NP.ratio,
        y = total_biomass),
    alpha = 0.3
  ) +
  geom_ribbon(
    data = predict_data,
    aes(x = NP_mg_l,
        ymin = conf.low,
        ymax = conf.high),
    alpha = 0.2,
    fill = "steelblue"
  ) +
  geom_line(
    data = predict_data,
    aes(x = NP_mg_l,
        y = predicted),
    linewidth = 1.2,
    color = "steelblue"
  ) +
  labs(
    x = "Total nitrogen:Total phosphorus (mg/L)",
    y = "Total biomass (ug/L)"
  ) +
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8)
  )
Biomass_NP_plot

#ggsave("Biomass_NP.png", width = 7, height = 4, units = "in", dpi = 300)


####Twelve-panel Temporal Nutrient/Diversity Figure#############################
#install.packages("patchwork")
library(patchwork)

#create plot
Temporal_Diversity_Nutrient_Models <- (Rich_TN_plot | Rich_TP_plot | Rich_NP_plot) /
  (Shan_TN_plot | Shan_TP_plot | Shan_NP_plot) + (Sim_TN_Plot | Sim_TP_Plot | Sim_NP_Plot) + 
  (Biomass_TN_plot | Biomass_TP_plot | Biomass_NP_plot) +
  plot_annotation(tag_levels = "A")
Temporal_Diversity_Nutrient_Models

#ggsave("Temporal_Diversity_Nutrient_Models.png", width = 11, height = 10, units = "in", dpi = 300)
#ggsave("Temporal_Diversity_Nutrient_Models_med.png", width = 9, height = 10, units = "in", dpi = 300)
#ggsave("Temporal_Diversity_Nutrient_Models_small.png", width = 8, height = 9, units = "in", dpi = 300)



#OBJECTIVE 2####################################################################
#Whole-community seasonal trajectories and temporal instability

###Hellinger Transformation of biomass data####################################

##Create Taxon column that gives genus, otherwise family, or infraclass
taxa.biomass.dataNEW <- taxa.biomass.dataNEW %>%
  mutate(
    Taxon = case_when(
      !is.na(Genus) & Genus != "" ~ Genus,
      !is.na(Family) & Family != "" ~ Family,
      !is.na(Order) & Order != "" ~ Order,
      !is.na(Infraclass) & Infraclass != "" ~ Infraclass,
      TRUE ~ NA_character_
    )
  )
#make smaller dataset to be transformed
Hellenger_biomass <- taxa.biomass.dataNEW %>%
  select(Round, Site, Taxon, TotalBiomassPerSample)

Hellenger_biomass$Round <- factor(Hellenger_biomass$Round)

#filter out round 0
Hellenger_biomass <- Hellenger_biomass %>%
  filter(Round !=0)

#combine duplicate taxa
#sum total biomass for duplicate taxa
Hellenger_biomass <- Hellenger_biomass %>%
  group_by(Round, Site, Taxon) %>%
  summarise(
    TotalBiomassPerSample = sum(TotalBiomassPerSample, na.rm = TRUE),
    .groups = "drop"
  )

#create sample id column to convert to rownames
Hellenger_biomass <- Hellenger_biomass %>%
  mutate(
    sample_id = paste(Site, Round, sep = "_")
  )

#metadata table
meta <- Hellenger_biomass %>%
  distinct(sample_id, Site, Round)

#switch from long dataset to wide
Hellenger_biomass_wide <- Hellenger_biomass %>%
  select (sample_id, Taxon, TotalBiomassPerSample) %>%
  pivot_wider(
    names_from = Taxon,
    values_from = TotalBiomassPerSample,
    values_fill = 0
  )

#sample id to rownames
Hellenger_biomass_wide <- Hellenger_biomass_wide %>%
  column_to_rownames("sample_id")

#Hellenger Transform
Trans_Hellenger_biomass <- decostand(
  Hellenger_biomass_wide,
  method = "hellinger"
)

#Save as dataframe
#write.csv(Trans_Hellenger_biomass, "Hellenger_biomass_per_site.csv")


#2A.Seasonal Composition########################################################
###Unconstrained PCA############################################################
##PCA with all wetland sites
#Run PCA
pca <- rda(Trans_Hellenger_biomass)

#create scores dataframe & add site & round columns
pca_scores <- as.data.frame(
  scores(pca, display = "sites", choices = 1:2)
  ) %>%
  rownames_to_column("sample_id") %>%
  left_join(meta, by = "sample_id")

#create sampling month column & define as factor
pca_scores<- pca_scores %>%
  mutate(
    Sampling_Month = case_when(
      Round %in% c("1") ~ "June",
      Round %in% c("2") ~ "July",
      Round %in% c("3") ~ "August"
    )
  )
pca_scores$Sampling_Month <- factor(pca_scores$Sampling_Month,
                                                levels = c("June", "July", "August"))

####Extract Eigenvalues#########################################################
#extract eigenvalues
ev <- eigenvals(pca)

#calc proportion & cumulative proportion of variance explained
prop_var <- ev / sum(ev)
cumu_var <- cumsum(prop_var)

#combine into df
pca_eigen_table <- data.frame(
  Axis = names(ev),
  Eigenvalue = round(as.numeric(ev), 4),
  Variance_explained = paste0(round(as.numeric(prop_var) * 100, 2), "%"),
  Cumulative_variance = paste0(round(as.numeric(cumu_var) * 100, 2), "%")
)

print(pca_eigen_table)


##PCA plot
#calculate convex hull for each Sampling Month
hulldata <- pca_scores %>%
  group_by(Sampling_Month) %>%
  slice(chull(PC1, PC2))

#calculate variance explained
var_explained <- pca$CA$eig / sum(pca$CA$eig) * 100
var_explained[1:2]

##Create vectors
#Small dataset with desired nutrients
Nutrients <- chem.hach.dataWide %>%
  select(Round, Site, Ammonia, Nitrate.Nitrite, TKN, Total.Nitrogen, Total.Phosphorus, NP.ratio)

#add same sample id column
Nutrients <- Nutrients %>%
  mutate(
    sample_id = paste(Site, Round, sep = "_")
  )

#select env variables for envfit
env <- Nutrients %>%
  select(
    Ammonia,
    Nitrate.Nitrite,
    TKN,
    Total.Nitrogen,
    Total.Phosphorus,
    NP.ratio
  )

ef <- envfit(
  pca,
  env,
  permutations = 999
)
ef

#extract vector coordinates for plotting
env_vectors <- scores(
  ef,
  display = "vectors"
) %>%
  as.data.frame() %>%
  rownames_to_column("Variable")

#pvalues
env_vectors <- env_vectors %>%
  mutate(
    r2 = ef$vectors$r,
    p = ef$vectors$pvals
  )

#scale arrows for plotting
env_vectors <- env_vectors %>%
  mutate(
    scale = case_when(
      Variable == "Nitrate.Nitrite" ~ 1.5,
      Variable == "Total.Phosphorus" ~ 3.5,
      Variable == "Total.Nitrogen" ~ 3.3,
      Variable == "Ammonia" ~ 1.75,
      TRUE ~ 2.05
    ),
    PC1 = PC1 * scale,
    PC2 = PC2 * scale,
    label = paste0(Variable, " (r^2 = ", round(r2, 2), ")")
  )

#add wetland type to pca scores dataset
pca_scores <- pca_scores %>%
  mutate(
    WetlandType = case_when(
      Site %in% c("DD15N", "SS", "DD8", "DD65", "LP", "WW", "FCC2",
                  "KS", "RS") ~ "CREP",
      Site %in% c("UH1", "UH2") ~ "Restored for habitat",
      Site %in% c("Bjork") ~ "Remnant"
    )
  )

####Sites PCA plot with vectors#################################################
PCA_sites <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  geom_point(
    aes(color = Sampling_Month, shape = WetlandType),
    size = 2,
    alpha = 0.8
  ) +
  geom_polygon(
    data = hulldata,
    aes(
      fill = Sampling_Month,
      color = Sampling_Month,
      group = Sampling_Month
    ),
    alpha = 0.15,
    linewidth = 0.7
  ) +
  geom_path(
    aes(group = Site),
    color = "grey60",
    alpha = 0.4,
    linewidth = 0.3,
    arrow = arrow(length = unit(0.12, "cm"))
  ) +
  # Environmental vectors
  geom_segment(
    data = env_vectors,
    aes(
      x = 0,
      y = 0,
      xend = PC1,
      yend = PC2
    ),
    inherit.aes = FALSE,
    linewidth = 0.55,
    color = "black",
    alpha = 0.5,
    arrow = arrow(length = unit(0.18, "cm"))
  ) +
  # Environmental vector labels
  geom_text(
    data = env_vectors,
    aes(
      x = PC1,
      y = PC2,
      label = Variable
    ),
    inherit.aes = FALSE,
    size = 3.55,
    fontface = "bold",
    color = "black",
    hjust = -0.1,
    vjust = -0.2
  ) + 
  scale_x_continuous(
    breaks = seq(-2, 2, by = 0.25)
  ) +
  scale_y_continuous(
    breaks = seq(-2, 2, by = 0.25)
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "inside",
    legend.position.inside = c(0.01, 0.99),
    #legend.box = "horizontal",
    legend.justification = c(0, 1),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
    color = "Sampling_Month",
    fill = "Sampling_Month"
  )
PCA_sites

#ggsave("WetlandTrajectoriesPanelB.png", width = 6, height = 6, units = "in", dpi = 300)


####Species PCA plot with vectors###############################################
#create species scores dataframe & add site & round columns
pca_scores_species <- as.data.frame(
  scores(pca, display = "species", choices = 1:2)
) %>%
  rownames_to_column("sample_id") %>%
  left_join(meta, by = "sample_id")

PCA_species <- ggplot(pca_scores_species, aes(x = PC1, y = PC2, label = abbreviate(sample_id, minlength = 3))) +
  geom_text_repel(
    max.overlaps = Inf,  
    size = 3.5,
    box.padding = 0.15,
    point.padding = 0.1,
    force = 2,
    force_pull = 0.5,
    max.time = 2,
    max.iter = 20000,
    min.segment.length = 0.2,
    segment.size = 0.25,
    segment.color = "grey65",
    segment.alpha = 0.45
  ) +
  # Environmental vectors
  geom_segment(
    data = env_vectors,
    aes(
      x = 0,
      y = 0,
      xend = PC1,
      yend = PC2
    ),
    inherit.aes = FALSE,
    linewidth = 0.55,
    color = "black",
    alpha = 0.4,
    arrow = arrow(length = unit(0.18, "cm"))
  ) +
  # Environmental vector labels
  geom_text(
    data = env_vectors,
    aes(
      x = PC1,
      y = PC2,
      label = Variable
    ),
    inherit.aes = FALSE,
    size = 3.55,
    fontface = "bold",
    color = "black",
    hjust = -0.1,
    vjust = -0.2
  ) + 
  scale_x_continuous(
    breaks = seq(-2, 2, by = 0.25)
  ) +
  scale_y_continuous(
    breaks = seq(-2, 2, by = 0.25)
  ) +
  #coord_cartesian(xlim = c(-0.75, 0.50), ylim = c(-0.75, 0.75)) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    # legend.position = "inside",
    # legend.position.inside = c(0.03, 0.97),
    # #legend.box = "horizontal",
    # legend.justification = c(0, 1),
    # legend.title = element_blank(),
    # legend.text = element_text(size = 8),
    # legend.key.size = unit(0.3, "cm"),
    # legend.spacing.y = unit(0.1, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)")
  )
PCA_species

#ggsave("PCAGenera.png", width = 7, height = 7, units = "in", dpi = 300)

####PCA Panel plot##############################################################
PCA_Panel <- (PCA_sites)/(PCA_species) +
   plot_annotation(tag_levels = "A")
PCA_Panel

#ggsave("PCA_Panel.png", width = 8, height = 9, units = "in", dpi = 300)


###2B.Wetland Trajectory Metrics################################################
##Wetland path, length, displacement, & directionality Calculations

##Full transformed community matrix (not just PC1 & PC2)
#calculate euclidean distances from transformed community
D_hel <- as.matrix(dist(Trans_Hellenger_biomass, method = "euclidean"))

#Nutrient dataset
Nutrients <- Nutrients %>%
  mutate(sample_id = as.character(sample_id))

#Site-level trajectory 
traj_metrics <- Nutrients %>%
  group_by(Site) %>%
  summarise(
    June = sample_id[Round == "1"],
    July = sample_id[Round == "2"],
    August = sample_id[Round == "3"],
    
    TN_mean = mean(Total.Nitrogen, na.rm = TRUE),
    TP_mean = mean(Total.Phosphorus, na.rm = TRUE),
    NP_mean = mean(NP.ratio, na.rm = TRUE),
    
    TN_max = max(Total.Nitrogen, na.rm = TRUE),
    TP_max = max(Total.Phosphorus, na.rm = TRUE),
    
    TN_cv = sd(Total.Nitrogen, na.rm = TRUE) / mean(Total.Nitrogen, na.rm = TRUE),
    TP_cv = sd(Total.Phosphorus, na.rm = TRUE) / mean(Total.Phosphorus, na.rm = TRUE),
    
    .groups = "drop"
  )

####Calculate Trajectory metrics################################################
traj_metrics <- traj_metrics %>%
  rowwise() %>%
  mutate(
    #distance between consecutive sampling periods
    d_June_July = D_hel[June, July],
    d_July_August = D_hel[July, August],
    #direct distance from beginning to end
    d_June_August = D_hel[June, August],
    #total trajectory length
    path_length = d_June_July + d_July_August,
    #net displacement
    net_displacement = d_June_August,
    #directionality
    directionality = net_displacement / path_length
  ) %>%
  ungroup()

#Table
traj_metrics %>%
  select(
    Site,
    path_length,
    net_displacement,
    directionality
  )

##Plot trajectory metrics against nutrients
####TP Trajcetory Plots#########################################################
#add wetland Type to traj_metrics
traj_metrics <- traj_metrics %>%
  mutate(
    WetlandType = case_when(
      Site %in% c("DD15N", "SS", "DD8", "DD65", "LP", "WW", "FCC2",
                  "KS", "RS") ~ "CREP",
      Site %in% c("UH1", "UH2") ~ "Restored for habitat",
      Site %in% c("Bjork") ~ "Remnant"
    )
  )


#Path length vs mean total phosphorus
path_TP_plot <- ggplot(traj_metrics, aes(x = TP_mean, y = path_length)) +
  geom_point(
    aes(color = WetlandType),
    size = 1.5 #3
  ) +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  scale_x_continuous(
    breaks = seq(0, 5, by = 0.25)
  ) +
  scale_y_continuous(
    breaks = seq(0, 5, by = 0.25)
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black"
  ) +
  theme_bw() +
  theme(
   axis.title.x = element_text(size = 10),
   axis.title.y = element_text(size = 10),
   axis.text.x = element_text(size = 8),
   axis.text.y = element_text(size = 8),
   legend.position = "none",
   #legend.position.inside = c(0.90, 0.90), #0.92, 0.92
   #legend.title = element_blank(),
   #legend.text = element_text(size = 4), #6
   #legend.key.size = unit(0.1, "cm"), #0.3
   panel.grid.major = element_blank(),
   panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Mean total phosphorus (mg/L)",
    y = "Seasonal path length"
    #color = "Wetland Type"
  )
path_TP_plot

#ggsave("SeasonalPathLength_TP.png", width = 7, height = 4, units = "in", dpi = 300)
#ggsave("SeasonalPathLength_TP_small.png", width = 4, height = 2, units = "in", dpi = 300)


#Directionality vs. mean total phosphorus
direction_TP_plot <- ggplot(traj_metrics, aes(x = TP_mean, y = directionality)) +
  geom_point(
    aes(color = WetlandType),
    size = 1.5 #3
  ) +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  scale_x_continuous(
    breaks = seq(0, 5, by = 0.25)
  ) +
  scale_y_continuous(
    breaks = seq(0, 5, by = 0.25)
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black"
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "none",
    #legend.position.inside = c(0.92, 0.92),
    #legend.title = element_blank(),
    #legend.text = element_text(size = 6),
    #legend.key.size = unit(0.3, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Mean total phosphorus (mg/L)",
    y = "Trajectory directionality"
    #color = "Wetland Type"
  )
direction_TP_plot

#ggsave("TrajectoryDirectionality_TP.png", width = 7, height = 4, units = "in", dpi = 300)
#ggsave("TrajectoryDirectionality_TP_small.png", width = 4, height = 2, units = "in", dpi = 300)


####TN Trajectory Plots#########################################################
#Path length vs mean total nitrogen
path_TN_plot <- ggplot(traj_metrics, aes(x = TN_mean, y = path_length)) +
  geom_point(
    aes(color = WetlandType),
    size = 1.5 #3
  ) +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  scale_y_continuous(
    breaks = seq(0.4, 4, by = 0.2)
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black"
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "none",
    #legend.position.inside = c(0.92, 0.91),
    #legend.title = element_blank(),
    #legend.text = element_text(size = 6),
    #legend.key.size = unit(0.3, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Mean total nitrogen (mg/L)",
    y = "Seasonal path length"
    #color = "Wetland Type"
  )
path_TN_plot

#ggsave("SeasonalPathLength_TN.png", width = 7, height = 4, units = "in", dpi = 300)
#ggsave("SeasonalPathLength_TN_small.png", width = 4, height = 2, units = "in", dpi = 300)


#Directionality vs. mean total nitrogen
direction_TN_plot <- ggplot(traj_metrics, aes(x = TN_mean, y = directionality)) +
  geom_point(
    aes(color = WetlandType),
    size = 1.5 #3
  ) +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  scale_y_continuous(
    breaks = seq(0.2, 2, by = 0.2)
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "black"
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "none",
    #legend.position.inside = c(0.10, 0.92),
    #legend.title = element_blank(),
    #legend.text = element_text(size = 6),
    #legend.key.size = unit(0.3, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Mean total nitrogen (mg/L)",
    y = "Trajectory directionality",
    color = "Wetland Type"
  )
direction_TN_plot

#ggsave("TrajectoryDirectionality_TN.png", width = 7, height = 4, units = "in", dpi = 300)
#ggsave("TrajectoryDirectionality_TN_small.png", width = 4, height = 2, units = "in", dpi = 300)


####Four-panel Trajectory Figure################################################
#create legend
metrics_legend <- cowplot::get_legend(
  direction_TN_plot +
    theme(
      legend.position = "bottom",
      legend.justification = "center"
    )
)

#create plot
Wetland_Trajectory_Nutrient_Patterns <- (path_TP_plot| direction_TP_plot) /
  (path_TN_plot | direction_TN_plot) +
  plot_annotation(tag_levels = "A")

#combine plot with legend
Wetland_Trajectory_Nutrient_Patterns <- cowplot::plot_grid(
  Wetland_Trajectory_Nutrient_Patterns,
  metrics_legend,
  ncol = 1,
  rel_heights = c(1, 0.10),
  align = "v"
)
Wetland_Trajectory_Nutrient_Patterns

#ggsave("Wetland_Trajectory_Nutrient_PatternsPanel.png", width = 12, height = 8, units = "in", dpi = 300, bg = "white")
#ggsave("Wetland_Trajectory_Nutrient_PatternsPanelsmall.png", width = 7, height = 5, units = "in", dpi = 300, bg = "white")


####Simple Regression & Spearman Correlation####################################
#####Total Phosphorus###########################################################
#Path Length
summary(lm(path_length ~ TP_mean, data = traj_metrics))
cor.test(traj_metrics$path_length, traj_metrics$TP_mean, method = "spearman")

#Directionality
summary(lm(directionality ~ TP_mean, data = traj_metrics))
cor.test(traj_metrics$directionality, traj_metrics$TP_mean, method = "spearman")


#####Total Nitrogen#############################################################
#Path Length
summary(lm(path_length ~ TN_mean, data = traj_metrics))
cor.test(traj_metrics$path_length, traj_metrics$TN_mean, method = "spearman")

#Directionality
summary(lm(directionality ~ TN_mean, data = traj_metrics))
cor.test(traj_metrics$directionality, traj_metrics$TN_mean, method = "spearman")



###2c.ILMS Procrustes residuals for discordant seasonal change##################
#Add Round & Site columns back into Trans Hellenger Biomass dataframe
Trans_Hellenger_biomass <- Trans_Hellenger_biomass %>%
  rownames_to_column("sample_id") %>%
  separate(sample_id, into = c("Site", "Round"), sep = "_")

#Order sites & species the same way for all round datasets
species_order <- c(
  "Polyarthra", "Plationus", "Platyias", "Lepadella", "Lecane", "Mytilina", "Brachionus",
  "Trichocerca", "Testudinella", "Euchlanis", "Rotaria", "Neocopepoda", "Cyclopidae",
  "Tropocyclops", "Acanthocyclops", "Microcyclops", "Macrocyclops", "Eucyclops", "Picripleuroxus", "Chydorus",
  "Ceriodaphnia", "Alona", "Collotheca", "Trichotria", "Simocephalus", "Ascomorpha",
  "Asplanchna", "Pleuroxus", "Scapholeberis", "Cyclops", "Anuraeopsis", "Bosmina",
  "Pompholyx", "Synchaeta", "Daphnia", "Asplanchnopus", "Hexarthra", "Kurzia", "Filinia",
  "Keratella", "Calanoida", "Diaphanosoma", "Camptocercus", "Monommata", "Ectocyclops", "Ilyocryptus",
  "Notholca", "Orthocyclops", "Harpacticoida", "Paracyclops", "Philodina", "Macrothrix", "Mesocyclops",
  "Skistodiaptomus", "Dunhevedia", "Leptodiaptomus", "Eurycercus", "Sida", "Attheyella"
)

site_order <- c(
  "UH1", "DD15N", "DD8", "RS", "Bjork", "DD65", "KS", "SS", "WW", "FCC2", "UH2", "LP"
)

#Create separate community matrices for each round
#will need repeated for each comparison of sampling periods
R1 <- Trans_Hellenger_biomass %>%
  filter(Round == 1) %>%
  mutate(Site = factor(Site, levels = site_order)) %>%
  arrange(Site) %>%
  as.data.frame() %>%
  select(Site, all_of(species_order), -Round)

R2 <- Trans_Hellenger_biomass %>%
  filter(Round == 2) %>%
  mutate(Site = factor(Site, levels = site_order)) %>%
  arrange(Site) %>%
  as.data.frame() %>%
  select(Site, all_of(species_order), -Round)

R3 <- Trans_Hellenger_biomass %>%
  filter(Round == 3) %>%
  mutate(Site = factor(Site, levels = site_order)) %>%
  arrange(Site) %>%
  as.data.frame() %>%
  select(Site, all_of(species_order), -Round)

#packages
library(gllvm)
library(vegan)
library(FactoMineR)
library(factoextra)
library(fBasics)#take trace
library(combinat)#permutation
library(readxl)
library(gtools) #for combination()
library(ggplot2)
library(ggrepel)  # for avoiding overlapping labels
library(patchwork)

#PA & ILMS code sources (code via Tang & Jackson [DOI:10.1002/env.70073])
source("resistant_codes_LMS_LTS.txt")
source("test of permutation.txt")
#when running sites ILMS_permute 'result1= LMS_LTS_estimator(P_mat_new,Q_mat,0,s_l_LMS,s_u_LMS) line needs 0
#when running species ILMS_permute 'result1= LMS_LTS_estimator(P_mat_new,Q_mat,1000,s_l_LMS,s_u_LMS) line needs 1000
#must repeat for each comparison of sampling periods

####Round 1 vs Round 2##########################################################
#re-run R1,2,3 dataframe creation after sites before species
dim(R1)
dim(R2)

rownames(R1) = R1[,1]
R1 = R1[,-1]
dim(R1)

rownames(R2) = R2[,1]
R2 = R2[,-1]
dim(R2)

#PCA1
PCA_R1 = rda(R1) 
summary(PCA_R1)
species_R1 = scores(PCA_R1, display = "species", choices = 1:5)
sites_R1 = scores(PCA_R1, display = "sites", choices = 1:5)

#PCA2
PCA_R2 = rda(R2)
summary(PCA_R2)
species_R2 = scores(PCA_R2, display = "species", choices = 1:2)
sites_R2 = scores(PCA_R2, display = "sites", choices = 1:2)


# #sites comparison
P_mat = sites_R1
Q_mat = sites_R2

#zooplankton species comparison
#P_mat = species_R1
#Q_mat = species_R2


#####PA
result_proc = Procrustes_fun(P_mat,Q_mat)
procrustes_obj = result_proc$s
res_procrustes = Procrustes_permute(P_mat,Q_mat,999,procrustes_obj)
sum(res_procrustes$converge_vec)
res_procrustes$p_value #P-value tests for concordance (low value = high concordance, high value = low/no concordance)

vector_1 = rep(1,dim(P_mat)[1])
result_procrustes_residual = 
  Q_mat - result_proc$scal*(P_mat)%*%(result_proc$rot)-vector_1%*%t(result_proc$trans_P)

#calculate Euclidean distance of each row
procrustes_residual_dis = sqrt(rowSums(result_procrustes_residual*result_procrustes_residual))
procrustes_residual_dis

# #save SITE resids for June-July
# June_July_site_resid <- data.frame(
#   Site = rownames(P_mat),
#   Procrustes_Residual = procrustes_residual_dis
# )

# #save SPECIES resids for June-July
# June_July_species_resid <- data.frame(
#   Taxon = rownames(P_mat),
#   Procrustes_Residual = procrustes_residual_dis
# )
# 
# ##Save residual Barplots (options for base R barplot and ggplot)
# #Sites baseR
# # png("PA_June_July_SiteResids.png", width = 7, height = 4, units = "in", res = 300)
# # PA_June_July_SiteResids <- barplot(procrustes_residual_dis,las = 2,main="PA: June-July",cex.names = 0.65)
# # dev.off()
# 
# #Sites ggplot
# procrustes_df <- data.frame(
#   Site = names(procrustes_residual_dis),
#   Residual = as.numeric(procrustes_residual_dis)
# )
# ggplot(procrustes_df, aes(x = Site, y = Residual)) +
#   geom_col(
#     fill = "grey70",
#     color = "black",
#     linewidth = 0.25
#   ) +
#   labs(
#     title = "Procrustes Analysis: Sites June-July",
#     x = NULL,
#     y = "Procrustes residual"
#   ) +
#   theme_bw() +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.text.x = element_text(
#       angle = 80,
#       hjust = 1,
#       size = 8
#     )
#   )
# ggsave("PA_June_July_SiteResids.png", width = 7, height = 4, units = "in", dpi = 300)
# 
#
# #Species BaseR
# # png("PA_June_July_SpeciesResids.png", width = 7, height = 4, units = "in", res = 300)
# # PA_June_July_SpeciesResids <- barplot(procrustes_residual_dis,las = 2,main="PA: June-July",cex.names = 0.65)
# # dev.off()
# 
# #Species ggplot
# procrustes_df <- data.frame(
#   Site = names(procrustes_residual_dis),
#   Residual = as.numeric(procrustes_residual_dis)
# )
# ggplot(procrustes_df, aes(x = Site, y = Residual)) +
#   geom_col(
#     fill = "grey70",
#     color = "black",
#     linewidth = 0.25
#     ) +
#   labs(
#     title = "Procrustes Analysis: Species June-July",
#     x = NULL,
#     y = "Procrustes residual"
#   ) +
#   theme_bw() +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.text.x = element_text(
#       angle = 80,
#       hjust = 1,
#       size = 8
#     )
#   )
# ggsave("PA_June_July_SpeciesResids.png", width = 7, height = 4, units = "in", dpi = 300)


#####ILMS#######################################################################
#result_ILMS= LMS_LTS_estimator(P_mat,Q_mat,1000,s_l_LMS,s_u_LMS)#for species as there are lots of combinations
result_ILMS= LMS_LTS_estimator(P_mat,Q_mat,0,s_l_LMS,s_u_LMS) # for ponds as there are only 220 combinations when retain two dimensions
ILMS_obj = result_ILMS$s
#res_ILMS = ILMS_permute(P_mat,Q_mat,999,ILMS_obj) #run readRDS to load saved datasets

##save res_ILMS for species to reload anytime and not have to re-calculate in R
#saveRDS(res_ILMS, file = "res_ILMS_sites_June_July.L")
#saveRDS(res_ILMS, file = "res_ILMS_species_June_July.L")
res_ILMS <- readRDS(file = "res_ILMS_sites_June_July.L") #run this line now whenever re-running sites for June-July
#res_ILMS <- readRDS(file = "res_ILMS_species_June_July.L") #run this line now whenever re-running species for June-July

ILMS_p_value = res_ILMS$p_value
ILMS_p_value #P-value tests for concordance (low value = high concordance, high value = low/no concordance)
#save(ILMS_p_value, file = "ILMS_p_value.RData")

vector_1 = rep(1,dim(P_mat)[1])
result_ILMS_residual =
  Q_mat - result_ILMS$scal*(P_mat)%*%(result_ILMS$rot)-vector_1%*%t(result_ILMS$trans_P)

#calculate Euclidean distance of each row
result_ILMS_residual_dis = sqrt(rowSums(result_ILMS_residual*result_ILMS_residual))
result_ILMS_residual_dis

# #save SITE resids for June-July
# June_July_site_resid <- data.frame(
#   Site = rownames(P_mat),
#   Procrustes_Residual = procrustes_residual_dis,
#   ILMS_Residual = result_ILMS_residual_dis
# )
# #Save as dataframe
# write.csv(June_July_site_resid, "PA_ILMS_June_July_site_resid.csv")

# #save species resids for June-July
# June_July_species_resid <- data.frame(
#   Taxon = rownames(P_mat),
#   Procrustes_Residual = procrustes_residual_dis,
#   ILMS_Residual = result_ILMS_residual_dis
# )
# #Save as dataframe
# write.csv(June_July_species_resid, "PA_ILMS_June_July_species_resid.csv")


##Save residual Barplots (options for base R barplot and ggplot)
#Sites baseR
# png("ILMS_June_July_SiteResids.png", width = 7, height = 4, units = "in", res = 300)
# ILMS_June_July_SiteResids <- barplot(procrustes_residual_dis,las = 2,main="ILMS: June-July",cex.names = 0.65)
# dev.off()
#
# #Sites ggplot
# ILMS_df <- data.frame(
#   Site = names(result_ILMS_residual_dis),
#   Residual = as.numeric(result_ILMS_residual_dis)
# )
# 
# Jun_Jul_ILMSresid <- ggplot(ILMS_df, aes(x = Site, y = Residual)) +
#   geom_col(
#     fill = "grey70",
#     color = "black",
#     linewidth = 0.25
#   ) +
#   labs(
#     #title = "Improved Least Median of Squares: Sites June-July",
#     x = NULL,
#     y = "June-July"
#   ) +
#   theme_bw() +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.text.x = element_text(
#       angle = 80,
#       hjust = 1,
#       size = 8
#     )
#   )
# Jun_Jul_ILMSresid
# 
# ggsave("ILMS_June_July_SiteResids.png", width = 7, height = 4, units = "in", dpi = 300)
# 
# #Species BaseR
# # png("ILMS_June_July_SpeciesResids.png", width = 7, height = 4, units = "in", res = 300)
# # ILMS_June_July_SpeciesResids <- barplot(result_ILMS_residual_dis,las = 2,main="ILMS: June-July",cex.names = 0.65)
# # dev.off()
# 
# #Species ggplot
# ILMS_df <- data.frame(
#   Site = names(result_ILMS_residual_dis),
#   Residual = as.numeric(result_ILMS_residual_dis)
# )
# ggplot(ILMS_df, aes(x = Site, y = Residual)) +
#   geom_col(
#     fill = "grey70",
#     color = "black",
#     linewidth = 0.25
#   ) +
#   labs(
#     title = "Improved Least Median of Squares: Species June-July",
#     x = NULL,
#     y = "ILMS residual"
#   ) +
#   theme_bw() +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.text.x = element_text(
#       angle = 80,
#       hjust = 1,
#       size = 8
#     )
#   )
# ggsave("ILMS_June_July_SpeciesResids.png", width = 7, height = 4, units = "in", dpi = 300)


##Plot residuals against nutrients
######TP Plots against nutrients################################################
#Create dataset combining all site residuals
resids <- PA_ILMS_June_July_site_resid %>%
  select(Site, Procrustes_Residual, ILMS_Residual) %>%
  rename(
    Proc_jn_jl = Procrustes_Residual,
    ILMS_jn_jl = ILMS_Residual
  ) %>%
  full_join(
    PA_ILMS_July_August_site_resid %>%
      select(Site, Procrustes_Residual, ILMS_Residual) %>%
      rename(
        Proc_jl_au = Procrustes_Residual,
        ILMS_jl_au = ILMS_Residual
      ),
    by = "Site"
  ) %>%
  full_join(
    PA_ILMS_June_August_site_resid %>%
      select(Site, Procrustes_Residual, ILMS_Residual) %>%
      rename(
        Proc_jn_au = Procrustes_Residual,
        ILMS_jn_au = ILMS_Residual
      ),
    by = "Site"
  )

#Join to traj_metrics dataframe
traj_metrics <- traj_metrics %>%
  left_join(resids, by = "Site")

traj_metrics <- traj_metrics %>%
  select(-June, -July, -August)

##Plot
#ILMS + TP (June-July)
ILMS_TP_JUN_JULY <- ggplot(traj_metrics,
       aes(x = TP_mean,
           y = ILMS_jn_jl)) +
  geom_point(
    aes(color = WetlandType),
    size = 1.5) + #3
  geom_smooth(method = "lm",
              se = TRUE,
              color = "black") +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "none",
    #legend.position.inside = c(0.09, 0.92),
    #legend.title = element_blank(),
    #legend.text = element_text(size = 6),
    #legend.key.size = unit(0.3, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Mean total phosphorus (mg/L)",
    y = "June-July ILMS residual",
    color = "Wetland Type"
  ) +
  annotate("text",
           x = -Inf,
           y = Inf,
           label = "p=0.036*",
           hjust = -0.5,
           vjust = 1.4,
           size = 3)
ILMS_TP_JUN_JULY

#ggsave("ILMS_TP_Jn_Jl.png", width = 7, height = 4, units = "in", dpi = 300)
#ggsave("ILMS_TP_Jn_Jl_small.png", width = 4, height = 2, units = "in", dpi = 300)

#Run linear model
MeanTP_ILMS_jn_jl <- lm(
  ILMS_jn_jl ~ TP_mean,
  data = traj_metrics
)
summary(MeanTP_ILMS_jn_jl)


######TN Plots against nutrients################################################
#ILMS + TN (June-July)
ILMS_TN_JUN_JULY <- ggplot(traj_metrics,
       aes(x = TN_mean,
           y = ILMS_jn_jl)) +
  geom_point(
    aes(color = WetlandType),
    size = 1.5) + #3
  geom_smooth(method = "lm",
              se = TRUE,
              color = "black") +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "none",
    #legend.position.inside = c(0.92, 0.92),
    #legend.title = element_blank(),
    #legend.text = element_text(size = 6),
    #legend.key.size = unit(0.3, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Mean total nitrogen (mg/L)",
    y = "June-July ILMS residual",
    color = "Wetland Type"
  )
ILMS_TN_JUN_JULY

ggsave("ILMS_TN_Jn_Jl.png", width = 7, height = 4, units = "in", dpi = 300)
ggsave("ILMS_TN_Jn_Jl_small.png", width = 4, height = 2, units = "in", dpi = 300)

#run linear model
MeanTN_ILMS_jn_jl <- lm(
  ILMS_jn_jl ~ TN_mean,
  data = traj_metrics
)
summary(MeanTN_ILMS_jn_jl)



####Round 2 vs Round 3##########################################################
R1 <- Trans_Hellenger_biomass %>%
  filter(Round == 1) %>%
  mutate(Site = factor(Site, levels = site_order)) %>%
  arrange(Site) %>%
  as.data.frame() %>%
  select(Site, all_of(species_order), -Round)

R2 <- Trans_Hellenger_biomass %>%
  filter(Round == 2) %>%
  mutate(Site = factor(Site, levels = site_order)) %>%
  arrange(Site) %>%
  as.data.frame() %>%
  select(Site, all_of(species_order), -Round)

R3 <- Trans_Hellenger_biomass %>%
  filter(Round == 3) %>%
  mutate(Site = factor(Site, levels = site_order)) %>%
  arrange(Site) %>%
  as.data.frame() %>%
  select(Site, all_of(species_order), -Round)

source("resistant_codes_LMS_LTS.txt")
#source("resistant_codes_LMS_LTS_save_memory.txt")
source("test of permutation.txt")
#when running sites ILMS_permute 'result1= LMS_LTS_estimator(P_mat_new,Q_mat,0,s_l_LMS,s_u_LMS) line needs 0
#when running species ILMS_permute 'result1= LMS_LTS_estimator(P_mat_new,Q_mat,1000,s_l_LMS,s_u_LMS) line needs 1000

dim(R2)
dim(R3)

rownames(R2) = R2[,1]
R2 = R2[,-1]
dim(R2)

rownames(R3) = R3[,1]
R3 = R3[,-1]
dim(R3)

#PCA1
PCA_R2 = rda(R2)
summary(PCA_R2)
species_R2 = scores(PCA_R2, display = "species", choices = 1:2)
sites_R2 = scores(PCA_R2, display = "sites", choices = 1:2)

#PCA2
PCA_R3 = rda(R3) 
summary(PCA_R3)
species_R3 = scores(PCA_R3, display = "species", choices = 1:5)
sites_R3 = scores(PCA_R3, display = "sites", choices = 1:5)


#sites comparison
P_mat = sites_R2
Q_mat = sites_R3

# #zooplankton species comparison
# P_mat = species_R2
# Q_mat = species_R3


#####PA
result_proc = Procrustes_fun(P_mat,Q_mat)
procrustes_obj = result_proc$s
res_procrustes = Procrustes_permute(P_mat,Q_mat,999,procrustes_obj)
sum(res_procrustes$converge_vec)
res_procrustes$p_value #P-value tests for concordance (low value = high concordance, high value = low/no concordance)

vector_1 = rep(1,dim(P_mat)[1])
result_procrustes_residual = 
  Q_mat - result_proc$scal*(P_mat)%*%(result_proc$rot)-vector_1%*%t(result_proc$trans_P)

#calculate Euclidean distance of each row
procrustes_residual_dis = sqrt(rowSums(result_procrustes_residual*result_procrustes_residual))
procrustes_residual_dis

# #save SITE resids for June-July
# July_August_site_resid <- data.frame(
#   Site = rownames(P_mat),
#   Procrustes_Residual = procrustes_residual_dis
# )
# 
# #save SPECIES resids for June-July
# July_August_species_resid <- data.frame(
#   Taxon = rownames(P_mat),
#   Procrustes_Residual = procrustes_residual_dis
# )


##Save residual Barplots (options for base R barplot and ggplot)
#Sites baseR
# png("PA_July_August_SiteResids.png", width = 7, height = 4, units = "in", res = 300)
# PA_July_August_SiteResids <- barplot(procrustes_residual_dis,las = 2,main="PA: July-August",cex.names = 0.65)
# dev.off()

# #Sites ggplot
# procrustes_df <- data.frame(
#   Site = names(procrustes_residual_dis),
#   Residual = as.numeric(procrustes_residual_dis)
# )
# ggplot(procrustes_df, aes(x = Site, y = Residual)) +
#   geom_col(
#     fill = "grey70",
#     color = "black",
#     linewidth = 0.25
#   ) +
#   labs(
#     title = "Procrustes Analysis: Sites July-August",
#     x = NULL,
#     y = "Procrustes residual"
#   ) +
#   theme_bw() +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.text.x = element_text(
#       angle = 80,
#       hjust = 1,
#       size = 8
#     )
#   )
# ggsave("PA_July_August_SiteResids.png", width = 7, height = 4, units = "in", dpi = 300)
# 
# #Species BaseR
# # png("PA_July_August_SpeciesResids.png", width = 7, height = 4, units = "in", res = 300)
# # PA_July_August_SpeciesResids <- barplot(procrustes_residual_dis,las = 2,main="PA: July-August",cex.names = 0.65)
# # dev.off()
# 
# #Species ggplot
# procrustes_df <- data.frame(
#   Site = names(procrustes_residual_dis),
#   Residual = as.numeric(procrustes_residual_dis)
# )
# ggplot(procrustes_df, aes(x = Site, y = Residual)) +
#   geom_col(
#     fill = "grey70",
#     color = "black",
#     linewidth = 0.25
#   ) +
#   labs(
#     title = "Procrustes Analysis: Species July-August",
#     x = NULL,
#     y = "Procrustes residual"
#   ) +
#   theme_bw() +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.text.x = element_text(
#       angle = 80,
#       hjust = 1,
#       size = 8
#     )
#   )
# ggsave("PA_July_August_SpeciesResids.png", width = 7, height = 4, units = "in", dpi = 300)


#####ILMS#######################################################################
#result_ILMS= LMS_LTS_estimator(P_mat,Q_mat,1000,s_l_LMS,s_u_LMS)#for species as there are lots of combinations
result_ILMS= LMS_LTS_estimator(P_mat,Q_mat,0,s_l_LMS,s_u_LMS) # for ponds as there are only 220 combinations when retain two dimensions
ILMS_obj = result_ILMS$s
#res_ILMS = ILMS_permute(P_mat,Q_mat,999,ILMS_obj) #run readRDS to load saved datasets

#save res_ILMS for species to reload anytime and not have to re-calculate in R
#saveRDS(res_ILMS, file = "res_ILMS_sites_July_August.L")
#saveRDS(res_ILMS, file = "res_ILMS_species_July_August.L")
res_ILMS <- readRDS(file = "res_ILMS_sites_July_August.L") #run this line now whenever re-running sites for June-July
#res_ILMS <- readRDS(file = "res_ILMS_species_July_August.L") #run this line now whenever re-running species for June-July

ILMS_p_value = res_ILMS$p_value
ILMS_p_value #P-value tests for concordance (low value = high concordance, high value = low/no concordance)
save(ILMS_p_value, file = "ILMS_p_value_July_August_species.RData") 

vector_1 = rep(1,dim(P_mat)[1])
result_ILMS_residual =
  Q_mat - result_ILMS$scal*(P_mat)%*%(result_ILMS$rot)-vector_1%*%t(result_ILMS$trans_P)

#calculate Euclidean distance of each row
result_ILMS_residual_dis = sqrt(rowSums(result_ILMS_residual*result_ILMS_residual))
result_ILMS_residual_dis

# #save SITE resids for July-August
# July_August_site_resid <- data.frame(
#   Site = rownames(P_mat),
#   Procrustes_Residual = procrustes_residual_dis,
#   ILMS_Residual = result_ILMS_residual_dis
# )
# #Save as dataframe
# write.csv(July_August_site_resid, "PA_ILMS_July_August_site_resid.csv")
# 
# #save species resids for July-August
# July_August_species_resid <- data.frame(
#   Taxon = rownames(P_mat),
#   Procrustes_Residual = procrustes_residual_dis,
#   ILMS_Residual = result_ILMS_residual_dis
# )
# #Save as dataframe
# write.csv(July_August_species_resid, "PA_ILMS_July_August_species_resid.csv")
# 
# 
# ##Save residual Barplots (options for base R barplot and ggplot)
# #Sites baseR
# # png("ILMS_July_August_SiteResids.png", width = 7, height = 4, units = "in", res = 300)
# # ILMS_July_August_SiteResids <- barplot(procrustes_residual_dis,las = 2,main="ILMS: July-August",cex.names = 0.65)
# # dev.off()
# 
# #Sites ggplot
# ILMS_df <- data.frame(
#   Site = names(result_ILMS_residual_dis),
#   Residual = as.numeric(result_ILMS_residual_dis)
# )
# 
# Jul_Aug_ILMSresid <- ggplot(ILMS_df, aes(x = Site, y = Residual)) +
#   geom_col(
#     fill = "grey70",
#     color = "black",
#     linewidth = 0.25
#   ) +
#   labs(
#     #title = "Improved Least Median of Squares: Sites July-August",
#     x = NULL,
#     y = "July-August"
#   ) +
#   theme_bw() +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.text.x = element_text(
#       angle = 80,
#       hjust = 1,
#       size = 8
#     )
#   )
# Jul_Aug_ILMSresid
# 
# ggsave("ILMS_July_August_SiteResids.png", width = 7, height = 4, units = "in", dpi = 300)
# 
# #Species BaseR
# # png("ILMS_July_August_SpeciesResids.png", width = 7, height = 4, units = "in", res = 300)
# # ILMS_July_August_SpeciesResids <- barplot(result_ILMS_residual_dis,las = 2,main="ILMS: July-August",cex.names = 0.65)
# # dev.off()
# 
# #Species ggplot
# ILMS_df <- data.frame(
#   Site = names(result_ILMS_residual_dis),
#   Residual = as.numeric(result_ILMS_residual_dis)
# )
# ggplot(ILMS_df, aes(x = Site, y = Residual)) +
#   geom_col(
#     fill = "grey70",
#     color = "black",
#     linewidth = 0.25
#   ) +
#   labs(
#     title = "Improved Least Median of Squares: Species July-August",
#     x = NULL,
#     y = "ILMS residual"
#   ) +
#   theme_bw() +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.text.x = element_text(
#       angle = 80,
#       hjust = 1,
#       size = 8
#     )
#   )
# ggsave("ILMS_July_August_SpeciesResids.png", width = 7, height = 4, units = "in", dpi = 300)


######TP Plots against nutrients################################################
#ILMS + TP (July-August)
ILMS_TP_JULY_AUG <- ggplot(traj_metrics,
       aes(x = TP_mean,
           y = ILMS_jl_au)) +
  geom_point(
    aes(color = WetlandType),
    size = 1.5) + #3
  geom_smooth(method = "lm",
              se = TRUE,
              color = "black") +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "none",
    #legend.position.inside = c(0.92, 0.92),
    #legend.title = element_blank(),
    #legend.text = element_text(size = 6),
    #legend.key.size = unit(0.3, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Mean total phosphorus (mg/L)",
    y = "July-August ILMS residual",
    color = "Wetland Type"
  )
ILMS_TP_JULY_AUG

ggsave("ILMS_TP_Jl_Au.png", width = 7, height = 4, units = "in", dpi = 300)
ggsave("ILMS_TP_Jl_Au_small.png", width = 4, height = 2, units = "in", dpi = 300)

#linear regression model
MeanTP_ILMS_jl_au <- lm(
  ILMS_jl_au ~ TP_mean,
  data = traj_metrics
)
summary(MeanTP_ILMS_jl_au)


######TN Plots against nutrients################################################
#ILMS + TN (July-August)
ILMS_TN_JULY_AUG <- ggplot(traj_metrics,
       aes(x = TN_mean,
           y = ILMS_jl_au)) +
  geom_point(
    aes(color = WetlandType),
    size = 1.5) + #3
  geom_smooth(method = "lm",
              se = TRUE,
              color = "black") +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "none",
    #legend.position.inside = c(0.11, 0.92),
    #legend.title = element_blank(),
    #legend.text = element_text(size = 6),
    #legend.key.size = unit(0.3, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Mean total nitrogen (mg/L)",
    y = "July-August ILMS residual",
    color = "Wetland Type"
  )
ILMS_TN_JULY_AUG

ggsave("ILMS_TN_Jl_Au.png", width = 7, height = 4, units = "in", dpi = 300)
ggsave("ILMS_TN_Jl_Au_small.png", width = 4, height = 2, units = "in", dpi = 300)

#linear model
MeanTN_ILMS_jl_au <- lm(
  ILMS_jl_au ~ TN_mean,
  data = traj_metrics
)
summary(MeanTN_ILMS_jl_au)



####Round 1 vs Round 3##########################################################
R1 <- Trans_Hellenger_biomass %>%
  filter(Round == 1) %>%
  mutate(Site = factor(Site, levels = site_order)) %>%
  arrange(Site) %>%
  as.data.frame() %>%
  select(Site, all_of(species_order), -Round)

R2 <- Trans_Hellenger_biomass %>%
  filter(Round == 2) %>%
  mutate(Site = factor(Site, levels = site_order)) %>%
  arrange(Site) %>%
  as.data.frame() %>%
  select(Site, all_of(species_order), -Round)

R3 <- Trans_Hellenger_biomass %>%
  filter(Round == 3) %>%
  mutate(Site = factor(Site, levels = site_order)) %>%
  arrange(Site) %>%
  as.data.frame() %>%
  select(Site, all_of(species_order), -Round)

source("resistant_codes_LMS_LTS.txt")
#source("resistant_codes_LMS_LTS_save_memory.txt")
source("test of permutation.txt")
#when running sites ILMS_permute 'result1= LMS_LTS_estimator(P_mat_new,Q_mat,0,s_l_LMS,s_u_LMS) line needs 0
#when running species ILMS_permute 'result1= LMS_LTS_estimator(P_mat_new,Q_mat,1000,s_l_LMS,s_u_LMS) line needs 1000

dim(R1)
dim(R3)

rownames(R1) = R1[,1]
R1 = R1[,-1]
dim(R1)

rownames(R3) = R3[,1]
R3 = R3[,-1]
dim(R3)

#PCA1
PCA_R1 = rda(R1)
summary(PCA_R1)
species_R1 = scores(PCA_R1, display = "species", choices = 1:2)
sites_R1 = scores(PCA_R1, display = "sites", choices = 1:2)

#PCA2
PCA_R3 = rda(R3) #conduct PCA here
summary(PCA_R3)
species_R3 = scores(PCA_R3, display = "species", choices = 1:5)
sites_R3 = scores(PCA_R3, display = "sites", choices = 1:5)


#sites comparison
P_mat = sites_R1
Q_mat = sites_R3

# #zooplankton species comparison
# P_mat = species_R1
# Q_mat = species_R3


#####PA
result_proc = Procrustes_fun(P_mat,Q_mat)
procrustes_obj = result_proc$s
res_procrustes = Procrustes_permute(P_mat,Q_mat,999,procrustes_obj)
sum(res_procrustes$converge_vec)
res_procrustes$p_value #P-value tests for concordance (low value = high concordance, high value = low/no concordance)

vector_1 = rep(1,dim(P_mat)[1])
result_procrustes_residual = 
  Q_mat - result_proc$scal*(P_mat)%*%(result_proc$rot)-vector_1%*%t(result_proc$trans_P)

#calculate Euclidean distance of each row
procrustes_residual_dis = sqrt(rowSums(result_procrustes_residual*result_procrustes_residual))
procrustes_residual_dis


# #save SITE resids for June-July
# June_August_site_resid <- data.frame(
#   Site = rownames(P_mat),
#   Procrustes_Residual = procrustes_residual_dis
# )
# 
# #save SPECIES resids for June-July
# June_August_species_resid <- data.frame(
#   Taxon = rownames(P_mat),
#   Procrustes_Residual = procrustes_residual_dis
# )


##Save residual Barplots (options for base R barplot and ggplot)
#Sites baseR
# png("PA_June_August_SiteResids.png", width = 7, height = 4, units = "in", res = 300)
# PA_June_August_SiteResids <- barplot(procrustes_residual_dis,las = 2,main="PA: June-August",cex.names = 0.65)
# dev.off()

# #Sites ggplot
# procrustes_df <- data.frame(
#   Site = names(procrustes_residual_dis),
#   Residual = as.numeric(procrustes_residual_dis)
# )
# ggplot(procrustes_df, aes(x = Site, y = Residual)) +
#   geom_col(
#     fill = "grey70",
#     color = "black",
#     linewidth = 0.25
#   ) +
#   labs(
#     title = "Procrustes Analysis: Sites June-August",
#     x = NULL,
#     y = "Procrustes residual"
#   ) +
#   theme_bw() +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.text.x = element_text(
#       angle = 80,
#       hjust = 1,
#       size = 8
#     )
#   )
# ggsave("PA_June_August_SiteResids.png", width = 7, height = 4, units = "in", dpi = 300)

# #Species BaseR
# # png("PA_June_August_SpeciesResids.png", width = 7, height = 4, units = "in", res = 300)
# # PA_June_August_SpeciesResids <- barplot(procrustes_residual_dis,las = 2,main="PA: June-August",cex.names = 0.65)
# # dev.off()
# 
# #Species ggplot
# procrustes_df <- data.frame(
#   Site = names(procrustes_residual_dis),
#   Residual = as.numeric(procrustes_residual_dis)
# )
# ggplot(procrustes_df, aes(x = Site, y = Residual)) +
#   geom_col(
#     fill = "grey70",
#     color = "black",
#     linewidth = 0.25
#   ) +
#   labs(
#     title = "Procrustes Analysis: Species June-August",
#     x = NULL,
#     y = "Procrustes residual"
#   ) +
#   theme_bw() +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.text.x = element_text(
#       angle = 80,
#       hjust = 1,
#       size = 8
#     )
#   )
# ggsave("PA_June_August_SpeciesResids.png", width = 7, height = 4, units = "in", dpi = 300)


#####ILMS#######################################################################
#result_ILMS= LMS_LTS_estimator(P_mat,Q_mat,1000,s_l_LMS,s_u_LMS)#for species as there are lots of combinations
result_ILMS= LMS_LTS_estimator(P_mat,Q_mat,0,s_l_LMS,s_u_LMS) # for ponds as there are only 220 combinations when retain two dimensions
ILMS_obj = result_ILMS$s
#res_ILMS = ILMS_permute(P_mat,Q_mat,999,ILMS_obj) #run readRDS to load saved datasets

#save res_ILMS for species to reload anytime and not have to re-calculate in R
#saveRDS(res_ILMS, file = "res_ILMS_sites_June_August.L")
#saveRDS(res_ILMS, file = "res_ILMS_species_June_August.L")
res_ILMS <- readRDS(file = "res_ILMS_sites_June_August.L") #run this line now whenever re-running sites for June-August
#res_ILMS <- readRDS(file = "res_ILMS_species_June_August.L") #run this line now whenever re-running species for June-August

ILMS_p_value = res_ILMS$p_value
ILMS_p_value #P-value tests for concordance (low value = high concordance, high value = low/no concordance)
#save(ILMS_p_value, file = "ILMS_p_value.RData")

vector_1 = rep(1,dim(P_mat)[1])
result_ILMS_residual =
  Q_mat - result_ILMS$scal*(P_mat)%*%(result_ILMS$rot)-vector_1%*%t(result_ILMS$trans_P)

#calculate Euclidean distance of each row
result_ILMS_residual_dis = sqrt(rowSums(result_ILMS_residual*result_ILMS_residual))
result_ILMS_residual_dis

# #save SITE resids for July-August
# June_August_site_resid <- data.frame(
#   Site = rownames(P_mat),
#   Procrustes_Residual = procrustes_residual_dis,
#   ILMS_Residual = result_ILMS_residual_dis
# )
# #Save as dataframe
# write.csv(June_August_site_resid, "PA_ILMS_June_August_site_resid.csv")
# 
# #save species resids for June-August
# June_August_species_resid <- data.frame(
#   Taxon = rownames(P_mat),
#   Procrustes_Residual = procrustes_residual_dis,
#   ILMS_Residual = result_ILMS_residual_dis
# )
# #Save as dataframe
# write.csv(June_August_species_resid, "PA_ILMS_June_August_species_resid.csv")


##Save residual Barplots (options for base R barplot and ggplot)
#Sites baseR
# png("ILMS_June_August_SiteResids.png", width = 7, height = 4, units = "in", res = 300)
# ILMS_June_August_SiteResids <- barplot(procrustes_residual_dis,las = 2,main="ILMS: June-August",cex.names = 0.65)
# dev.off()

# #Sites ggplot
# ILMS_df <- data.frame(
#   Site = names(result_ILMS_residual_dis),
#   Residual = as.numeric(result_ILMS_residual_dis)
# )
# 
# Jun_Aug_ILMSresid <- ggplot(ILMS_df, aes(x = Site, y = Residual)) +
#   geom_col(
#     fill = "grey70",
#     color = "black",
#     linewidth = 0.25
#   ) +
#   labs(
#     #title = "Improved Least Median of Squares: Sites June-August",
#     x = NULL,
#     y = "June-August"
#   ) +
#   theme_bw() +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.text.x = element_text(
#       angle = 80,
#       hjust = 1,
#       size = 8
#     )
#   )
# Jun_Aug_ILMSresid
# 
# ggsave("ILMS_June_August_SiteResids.png", width = 7, height = 4, units = "in", dpi = 300)
# 
# #Species BaseR
# # png("ILMS_June_August_SpeciesResids.png", width = 7, height = 4, units = "in", res = 300)
# # ILMS_June_August_SpeciesResids <- barplot(result_ILMS_residual_dis,las = 2,main="ILMS: June-August",cex.names = 0.65)
# # dev.off()
# 
# #Species ggplot
# ILMS_df <- data.frame(
#   Site = names(result_ILMS_residual_dis),
#   Residual = as.numeric(result_ILMS_residual_dis)
# )
# ggplot(ILMS_df, aes(x = Site, y = Residual)) +
#   geom_col(
#     fill = "grey70",
#     color = "black",
#     linewidth = 0.25
#   ) +
#   labs(
#     title = "Improved Least Median of Squares: Species June-August",
#     x = NULL,
#     y = "ILMS residual"
#   ) +
#   theme_bw() +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.text.x = element_text(
#       angle = 80,
#       hjust = 1,
#       size = 8
#     )
#   )
# ggsave("ILMS_June_August_SpeciesResids.png", width = 7, height = 4, units = "in", dpi = 300)


######TP Plots against nutrients################################################
#ILMS + TP (June-August)
ILMS_TP_JUN_AUG <- ggplot(traj_metrics,
       aes(x = TP_mean,
           y = ILMS_jn_au)) +
  geom_point(
    aes(color = WetlandType),
    size = 1.5) + #3
  geom_smooth(method = "lm",
              se = TRUE,
              color = "black") +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "none",
    #legend.position.inside = c(0.88, 0.075),
    #legend.title = element_blank(),
    #legend.text = element_text(size = 6),
    #legend.key.size = unit(0.3, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Mean total phosphorus (mg/L)",
    y = "June-August ILMS residual",
    color = "Wetland Type"
  )
ILMS_TP_JUN_AUG

#ggsave("ILMS_TP_Jn_Au.png", width = 7, height = 4, units = "in", dpi = 300)
#ggsave("ILMS_TP_Jn_Au_small.png", width = 4, height = 2, units = "in", dpi = 300)

#linear model
MeanTP_ILMS_jn_au <- lm(
  ILMS_jn_au ~ TP_mean,
  data = traj_metrics
)
summary(MeanTP_ILMS_jn_au)


######TN Plots against nutrients################################################
#ILMS + TN (June-August)
ILMS_TN_JUN_AUG <- ggplot(traj_metrics,
       aes(x = TN_mean,
           y = ILMS_jn_au)) +
  geom_point(
    aes(color = WetlandType),
    size = 1.5) + #3
  geom_smooth(method = "lm",
              se = TRUE,
              color = "black") +
  scale_color_manual(
    values = c(
      "CREP" = "#9E9AC8",
      "Restored for habitat" = "#A6D854",
      "Remnant" = "#D55E00"
    )
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "none",
    #legend.position.inside = c(0.92, 0.92),
    #legend.title = element_blank(),
    #legend.text = element_text(size = 6),
    #legend.key.size = unit(0.3, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Mean total nitrogen (mg/L)",
    y = "June-August ILMS residual",
    color = "Wetland Type"
  )
ILMS_TN_JUN_AUG

#ggsave("ILMS_TN_Jn_Au.png", width = 7, height = 4, units = "in", dpi = 300)
#ggsave("ILMS_TN_Jn_Au_small.png", width = 4, height = 2, units = "in", dpi = 300)

#linear model
MeanTN_ILMS_jn_au <- lm(
  ILMS_jn_au ~ TN_mean,
  data = traj_metrics
)
summary(MeanTN_ILMS_jn_au)

####Six-panel ILMS Resid + Nutrients Figure#####################################
#install.packages("patchwork")
#install.packages("cowplot")
library(patchwork)
library(cowplot)

ILMS_legend <- cowplot::get_legend(
  ILMS_TN_JUN_AUG +
    theme(
      legend.position = "bottom",
      legend.justification = "center"
    )
)

plots <- list(
  ILMS_TP_JUN_JULY + theme(legend.position = "none"),
  ILMS_TN_JUN_JULY + theme(legend.position = "none"),
  ILMS_TP_JULY_AUG + theme(legend.position = "none"),
  ILMS_TN_JULY_AUG + theme(legend.position = "none"),
  ILMS_TP_JUN_AUG + theme(legend.position = "none"),
  ILMS_TN_JUN_AUG + theme(legend.position = "none")
)

row1 <- cowplot::plot_grid(
  plots[[1]], plots[[2]],
  labels = c("A", "B"),
  ncol = 2,
  label_size = 14,
  label_fontface = "bold",
  label_x = 0.02,
  label_y = 0.98,
  hjust = 0,
  vjust = 1
)

row2 <- cowplot::plot_grid(
  plots[[3]], plots[[4]],
  labels = c("C", "D"),
  ncol = 2,
  label_size = 14,
  label_fontface = "bold",
  label_x = 0.02,
  label_y = 0.98,
  hjust = 0,
  vjust = 1
)

row3 <- cowplot::plot_grid(
  plots[[5]], plots[[6]],
  labels = c("E", "F"),
  ncol = 2,
  label_size = 14,
  label_fontface = "bold",
  label_x = 0.02,
  label_y = 0.98,
  hjust = 0,
  vjust = 1
)

ILMS_Nutrient_Patterns <- cowplot::plot_grid(
  row1,
  row2,
  row3,
  ILMS_legend,
  ncol = 1,
  rel_heights = c(1, 1, 1, 0.12)
)
ILMS_Nutrient_Patterns

#ggsave("ILMS_NutrientsPanel.png", width = 9, height = 10, units = "in", dpi = 300)


####Three-panel ILMS Residual Figure############################################
ILMS_Residuals <- (Jun_Jul_ILMSresid) /
  (Jul_Aug_ILMSresid) / (Jun_Aug_ILMSresid) +
  plot_annotation(tag_levels = "A")
ILMS_Residuals

#ggsave("ILMS_residual_plot_panel.png", width = 12, height = 8, units = "in", dpi = 300)
#ggsave("ILMS_residual_plot_panel_small.png", width = 6, height = 5, units = "in", dpi = 300)



###OBJECTIVE 3##################################################################
##Community-level environmental association

####3A.PERMANOVA################################################################
#supporting evidence for trajectory analysis

#Add sampling month & nutrient regime to Nutrients dataset
Nutrients <- Nutrients %>%
  mutate(
    Sampling_Month = case_when(
      Round %in% c("1") ~ "June",
      Round %in% c("2") ~ "July",
      Round %in% c("3") ~ "August"
    )
  )
Nutrients$Sampling_Month <- factor(Nutrients$Sampling_Month,
                                            levels = c("June", "July", "August"))

#Add molar N:P nutrient regime classifications to trajectory dataset
Nutrients <- Nutrients %>%
  left_join(
    NP_Regimes %>%
      select(Site, Nutrient_group_molL_ratio),
    by = "Site"
  )

#Permanova
perm_site <- how(nperm = 999)
setBlocks(perm_site) <- Nutrients$Site

D <- dist(Trans_Hellenger_biomass, method = "euclidean")

adonis_season <- adonis2(
  D ~ Sampling_Month,
  data = Nutrients,
  permutations = perm_site
)
adonis_season

adonis_TP <- adonis2(
  D ~ Sampling_Month + Total.Phosphorus,
  data = Nutrients,
  permutations = perm_site
)
adonis_TP #(0.227)

#testing out logging predictor vs response
#response
# adonis_TP2 <- adonis2(
#   D ~ Sampling_Month + Total.Phosphorus,
#   family = gaussian(link = log),
#   data = Nutrients,
#   permutations = perm_site
# )
# adonis_TP2 #lower p (0.203 compared to 0.227)
# #predictor
# adonis_TP3 <- adonis2(
#   D ~ Sampling_Month + log(Total.Phosphorus),
#   data = Nutrients,
#   permutations = perm_site
# )
# adonis_TP3 #higher p (0.248)
# #response & predictor
# adonis_TP4 <- adonis2(
#   D ~ Sampling_Month + log(Total.Phosphorus),
#   family = gaussian(link = log),
#   data = Nutrients,
#   permutations = perm_site
# )
# adonis_TP4 #higher p (0.229)


adonis_TN <- adonis2(
  D ~ Sampling_Month + Total.Nitrogen,
  data = Nutrients,
  permutations = perm_site
)
adonis_TN

#testing out logging predictor vs response
#response
# adonis_TN2 <- adonis2(
#   D ~ Sampling_Month + Total.Nitrogen,
#   family = gaussian(link = log),
#   data = Nutrients,
#   permutations = perm_site
# )
# adonis_TN2 #lower p (0.087 compared to 0.089)
# #predictor
adonis_TN3 <- adonis2(
  D ~ Sampling_Month + log(Total.Nitrogen),
  data = Nutrients,
  permutations = perm_site
)
adonis_TN3 #lower p (0.046* compared to 0.089)
#response & predictor
# adonis_TN4 <- adonis2(
#   D ~ Sampling_Month + log(Total.Nitrogen),
#   family = gaussian(link = log),
#   data = Nutrients,
#   permutations = perm_site
# )
# adonis_TN4 #higher p (0.054)


adonis_TNTP <- adonis2(
  D ~ Sampling_Month + NP.ratio,
  data = Nutrients,
  permutations = perm_site
)
adonis_TNTP #p = 0.03
(adonis_TNTP)

#testing out logging predictor vs response
#response
# adonis_TNTP2 <- adonis2(
#   D ~ Sampling_Month + NP.ratio,
#   family = gaussian(link = log),
#   data = Nutrients,
#   permutations = perm_site
# )
# adonis_TNTP2 #0.03
# #predictor
# adonis_TNTP3 <- adonis2(
#   D ~ Sampling_Month + log(NP.ratio),
#   data = Nutrients,
#   permutations = perm_site
# )
# adonis_TNTP3 #0.136
# #response & predictor
# adonis_TNTP4 <- adonis2(
#   D ~ Sampling_Month + log(NP.ratio),
#   family = gaussian(link = log),
#   data = Nutrients,
#   permutations = perm_site
# )
# adonis_TNTP4 #0.128


adonis_TP_int <- adonis2(
  D ~ Sampling_Month * Total.Phosphorus,
  data = Nutrients,
  permutations = perm_site
)
adonis_TP_int #0.556

#testing out logging predictor vs response
#response
# adonis_TP_int2 <- adonis2(
#   D ~ Sampling_Month * Total.Phosphorus,
#   family = gaussian(link = log),
#   data = Nutrients,
#   permutations = perm_site
# )
# adonis_TP_int2 #0.544
#predictor
# adonis_TP_int3 <- adonis2(
#   D ~ Sampling_Month * log(Total.Phosphorus),
#   data = Nutrients,
#   permutations = perm_site
# )
# adonis_TP_int3 #0.606
# #response & predictor
# adonis_TP_int4 <- adonis2(
#   D ~ Sampling_Month * log(Total.Phosphorus),
#   family = gaussian(link = log),
#   data = Nutrients,
#   permutations = perm_site
# )
# adonis_TP_int4 #0.604

adonis_TN_int <- adonis2(
  D ~ Sampling_Month * Total.Nitrogen,
  data = Nutrients,
  permutations = perm_site
)
adonis_TN_int #0.016

#testing out logging predictor vs response
# #response
# adonis_TN_int2 <- adonis2(
#   D ~ Sampling_Month * Total.Nitrogen,
#   family = gaussian(link = log),
#   data = Nutrients,
#   permutations = perm_site
# )
# adonis_TN_int2 #0.026
# #predictor
# adonis_TN_int3 <- adonis2(
#   D ~ Sampling_Month * log(Total.Nitrogen),
#   data = Nutrients,
#   permutations = perm_site
# )
# adonis_TN_int3 #0.049
# #response & predictor
# adonis_TN_int4 <- adonis2(
#   D ~ Sampling_Month * log(Total.Nitrogen),
#   family = gaussian(link = log),
#   data = Nutrients,
#   permutations = perm_site
# )
# adonis_TN_int4 #0.041


adonis_TNTP_int <- adonis2(
  D ~ Sampling_Month * NP.ratio,
  data = Nutrients, 
  permutations = perm_site
)
adonis_TNTP_int #0.083

#testing out logging predictor vs response
#response
# adonis_TNTP_int2 <- adonis2(
#   D ~ Sampling_Month * NP.ratio,
#   family = gaussian(link = log),
#   data = Nutrients, 
#   permutations = perm_site
# )
# adonis_TNTP_int2 #0.09
# #predictor
# adonis_TNTP_int3 <- adonis2(
#   D ~ Sampling_Month * log(NP.ratio),
#   data = Nutrients, 
#   permutations = perm_site
# )
# adonis_TNTP_int3 #0.325
# #response & predictor
# adonis_TNTP_int4 <- adonis2(
#   D ~ Sampling_Month * log(NP.ratio),
#   family = gaussian(link = log),
#   data = Nutrients, 
#   permutations = perm_site
# )
# adonis_TNTP_int4 #0.332


#check for dispersion
bd_season <- betadisper(D, Nutrients$Sampling_Month)
anova(bd_season)
permutest(bd_season)

bd_regime <- betadisper(D, Nutrients$Nutrient_group_molL_ratio)
anova(bd_regime)
permutest(bd_regime)

####3B.Complementary RDA########################################################
rda_TP <- rda(Trans_Hellenger_biomass ~ Sampling_Month + Total.Phosphorus, data = Nutrients)
anova(rda_TP, permutations = perm_site)
anova(rda_TP, by = "term", permutations = perm_site)

rda_TN <- rda(Trans_Hellenger_biomass ~ Sampling_Month + Total.Nitrogen, data = Nutrients)
anova(rda_TN, permutations = perm_site)
anova(rda_TN, by = "term", permutations = perm_site)

rda_NP <- rda(Trans_Hellenger_biomass ~ Sampling_Month + NP.ratio, data = Nutrients)
anova(rda_NP, permutations = perm_site)
anova(rda_NP, by = "term", permutations = perm_site)



###OBJECTIVE 4##################################################################
#Taxon indicators of extreme nutrient regimes and continuous gradients

###create categorical molar nutrient regime dataframe
NP_Regimes <- chem.hach.dataWide %>%
  group_by(Site) %>%
  summarise(
    Total.Nitrogen = mean(Total.Nitrogen, na.rm = TRUE),
    Total.Phosphorus = mean(Total.Phosphorus, na.rm = TRUE),
    .groups = "drop"
  )

##Nutrient regime categories based on molar ratio
NP_Regimes <- NP_Regimes %>%
  mutate(
    # Convert TP reported as PO4 to elemental P
    P_mgL = Total.Phosphorus * (30.97 / 94.97),
    # Convert TN and P from mg/L to mol/L
    N_molL = Total.Nitrogen / 14010,
    P_molL = P_mgL / 30970,
    # Molar N:P ratio
    NP_ratio_molL = N_molL / P_molL
  )

#use molar N:P ratio tertiles to define nutrient regimes
q <- quantile(
  NP_Regimes$NP_ratio_molL,
  probs = c(1/3, 2/3),
  na.rm = TRUE
)

NP_Regimes <- NP_Regimes %>%
  mutate(
    Nutrient_group_molL_ratio = cut(
      NP_ratio_molL,
      breaks = c(-Inf, q[1], q[2], Inf),
      labels = c(
        "Low N:P",
        "Moderate N:P",
        "High N:P"
      ),
      include.lowest = TRUE
    )
  )


#Add molar N:P nutrient regime classifications to trajectory dataset
traj_metrics <- traj_metrics %>%
  left_join(
    NP_Regimes %>%
      select(Site, Nutrient_group_molL_ratio),
    by = "Site"
  )


###4A.IndicSpecies##############################################################

##Create separate biomass datasets per round/season & aggregate dup genera (sum biomass)
#start with main dataset
zoop_biomass <- taxa.biomass.dataNEW %>%
  select(Round, Site, Date, Zone, unique_id, Taxon, TotalBiomassPerSample)

zoop_biomass <- zoop_biomass %>% #remove round 0
  filter(Round !=0)

#filter for rare taxa
taxon_occurence <- zoop_biomass %>% 
  group_by(Taxon) %>%
  summarise(
    n_samples = sum(TotalBiomassPerSample > 0, na.rm = TRUE),
    total_samples = n_distinct(unique_id),
    occurrence = n_samples / total_samples
  )
keep_taxa <- taxon_occurence %>%
  filter(occurrence >= 0.10) %>%
  pull(Taxon)
zoop_biomass_filtered <- zoop_biomass %>%
  filter(Taxon %in% keep_taxa)

#separate into seasons (avoids needing to account for multiple sampling rounds)
#june
zoop_biomass_june <- zoop_biomass_filtered %>%
  filter(Round == 1)
zoop_biomass_june <- zoop_biomass_june %>%
  group_by(Site, Round, Taxon) %>%
  summarise(
    TotalBiomass = sum(TotalBiomassPerSample, na.rm = TRUE),
    .groups = "drop"
  )
zoop_biomass_june_wide <- zoop_biomass_june %>% #switch from long dataset to wide
  select(Site, Taxon, TotalBiomass) %>%
  pivot_wider(
    names_from = Taxon,
    values_from = TotalBiomass,
    values_fill = 0
  )
zoop_biomass_june_wide <- zoop_biomass_june_wide %>% #sample id to rownames
  column_to_rownames("Site")

#july
zoop_biomass_july <- zoop_biomass_filtered %>%
  filter(Round == 2)
zoop_biomass_july <- zoop_biomass_july %>%
  group_by(Site, Round, Taxon) %>%
  summarise(
    TotalBiomass = sum(TotalBiomassPerSample, na.rm = TRUE),
    .groups = "drop"
  )
zoop_biomass_july_wide <- zoop_biomass_july %>% #switch from long dataset to wide
  select(Site, Taxon, TotalBiomass) %>%
  pivot_wider(
    names_from = Taxon,
    values_from = TotalBiomass,
    values_fill = 0
  )
zoop_biomass_july_wide <- zoop_biomass_july_wide %>% #sample id to rownames
  column_to_rownames("Site")

#august
zoop_biomass_august <- zoop_biomass_filtered %>%
  filter(Round == 3)
zoop_biomass_august <- zoop_biomass_august %>%
  group_by(Site, Round, Taxon) %>%
  summarise(
    TotalBiomass = sum(TotalBiomassPerSample, na.rm = TRUE),
    .groups = "drop"
  )
zoop_biomass_august_wide <- zoop_biomass_august %>% #switch from long dataset to wide
  select(Site, Taxon, TotalBiomass) %>%
  pivot_wider(
    names_from = Taxon,
    values_from = TotalBiomass,
    values_fill = 0
  )
zoop_biomass_august_wide <- zoop_biomass_august_wide %>% #sample id to rownames
  column_to_rownames("Site")

#NP Regimes averaged across rounds per site
#add sites to rownames
NP_Regimes <- NP_Regimes %>% #site to rownames
  column_to_rownames("Site")

##IndicSpecies Analysis
#Round 1/ June
June_indval <- multipatt(
  zoop_biomass_june_wide,
  NP_Regimes$Nutrient_group_molL_ratio,
  func = "IndVal.g",
  control = how(nperm = 999)
)
June_indval

#Round 2/ July
July_indval <- multipatt(
  zoop_biomass_july_wide,
  NP_Regimes$Nutrient_group_molL_ratio,
  func = "IndVal.g",
  control = how(nperm = 999)
)
July_indval

#Round 3/ August
August_indval <- multipatt(
  zoop_biomass_august_wide,
  NP_Regimes$Nutrient_group_molL_ratio,
  func = "IndVal.g",
  control = how(nperm = 999)
)
August_indval



###4B.TITAN######################################################################

####Continuous gradients#########################################################
#install.packages("TITAN2")
library(TITAN2)

#add n:p mg/l
NP_Regimes <- NP_Regimes %>%
  mutate(
    NP.Ratio_mgL = Total.Nitrogen / Total.Phosphorus
  )

NP_Regimes <- NP_Regimes %>%
  column_to_rownames("Site")

#aggregate biomass to site level
comm_site<- zoop_biomass %>%
  group_by(Site, Taxon) %>%
  summarise(
    TotalBiomass = mean(TotalBiomassPerSample, na.rm = TRUE), #mean biomass rather than total
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Taxon,
    values_from = TotalBiomass,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(comm_site) <- comm_site$Site #Sites to rownames
comm_site$Site <- NULL

#filter taxa for titan (taxa > or = 3 samples)
comm_TITAN <- comm_site[, colSums(comm_site > 0) >= 3, drop = FALSE]

sort(colSums(comm_TITAN > 0)) #check

#align nutrient data to community matrix
nutrient_TITAN <- NP_Regimes[
  match(rownames(comm_TITAN), rownames(NP_Regimes)),
  ,
  drop = FALSE
]

identical( #check
  rownames(comm_TITAN),
  rownames(nutrient_TITAN)
)
nrow(comm_TITAN)
nrow(nutrient_TITAN)

#create vectors
TN_vec <- nutrient_TITAN$Total.Nitrogen
TP_vec <- nutrient_TITAN$Total.Phosphorus
NP_vec <- nutrient_TITAN$NP.Ratio_mgL

length(TN_vec) #check
length(TP_vec)
length(NP_vec)

nrow(comm_TITAN)

sum(is.na(TN_vec)) #check
sum(is.na(TP_vec))
sum(is.na(NP_vec))

#Run TITAN for TN
#####TN##########################################################################
set.seed(123)
Titan_TN <- titan(
  env = TN_vec,
  txa = comm_TITAN,
  minSplt = 5,
  numPerm = 999,
  boot = TRUE,
  nBoot = 999,
  imax = FALSE
)

summary(Titan_TN)
plot_sumz(Titan_TN)
plot_taxa(Titan_TN)

Titan_TN$sppmax

##print top 5 taxa candidates
Titan_TN_topTaxa <- as.data.frame(Titan_TN$sppmax) %>%
  rownames_to_column(var = "Taxon") %>%
  mutate(
    across(
      c(ienv.cp, zenv.cp, maxgrp, obsiv.prob,
        zscore, purity, reliability, filter),
      as.numeric
    )
  )
#screening thresholds for top 5
min_zscore <- 1.50
max_perm_p <- 0.12
min_purity <- 0.70
min_reliability <- 0.30

#select & rank taxa
candidate_taxa_TN <- Titan_TN_topTaxa %>%
  dplyr::filter(
    zscore >= min_zscore,
    obsiv.prob <= max_perm_p,
    purity >= min_purity,
    reliability >= min_reliability
  ) %>%
  mutate(
    Direction = case_when(
      maxgrp == 1 ~ "Decreaser",
      maxgrp == 2 ~ "Increaser",
      TRUE ~ "Unknown"
    )
  ) %>%
  arrange(
    desc(zscore),
    obsiv.prob,
    desc(purity),
    desc(reliability)
  ) %>%
  slice_head(n = 5) %>%
  transmute(
    Taxon,
    Direction,
    `TN change point` = zenv.cp,
    `Z-score` = zscore,
    `p-value` = obsiv.prob,
    Purity = purity,
    `Reliability` = reliability
  ) %>%
  mutate(
    `TN change point` = round(`TN change point`, 2),
    `Z-score` = round(`Z-score`, 2),
    `p-value` = signif(`p-value`, 3),
    Purity = round(Purity, 3),
    `Reliability` = round(`Reliability`, 3)
  )

print(candidate_taxa_TN)

#Run TITAN for TP
#####TP##########################################################################
set.seed(123)
Titan_TP <- titan(
  env = TP_vec,
  txa = comm_TITAN,
  minSplt = 5,
  numPerm = 999,
  boot = TRUE,
  nBoot = 999,
  imax = FALSE
)

summary(Titan_TP)
plot_sumz(Titan_TP)
plot_taxa(Titan_TP)

##print top 5 taxa candidates
Titan_TP_topTaxa <- as.data.frame(Titan_TP$sppmax) %>%
  rownames_to_column(var = "Taxon") %>%
  mutate(
    across(
      c(ienv.cp, zenv.cp, maxgrp, obsiv.prob,
        zscore, purity, reliability, filter),
      as.numeric
    )
  )
#screening thresholds for top 5
min_zscore <- 1.50
max_perm_p <- 0.12
min_purity <- 0.70
min_reliability <- 0.30

#select & rank taxa
candidate_taxa_TP <- Titan_TP_topTaxa %>%
  dplyr::filter(
    zscore >= min_zscore,
    obsiv.prob <= max_perm_p,
    purity >= min_purity,
    reliability >= min_reliability
  ) %>%
  mutate(
    Direction = case_when(
      maxgrp == 1 ~ "Decreaser",
      maxgrp == 2 ~ "Increaser",
      TRUE ~ "Unknown"
    )
  ) %>%
  arrange(
    desc(zscore),
    obsiv.prob,
    desc(purity),
    desc(reliability)
  ) %>%
  slice_head(n = 5) %>%
  transmute(
    Taxon,
    Direction,
    `TP change point` = zenv.cp,
    `Z-score` = zscore,
    `p-value` = obsiv.prob,
    Purity = purity,
    `Reliability` = reliability
  ) %>%
  mutate(
    `TP change point` = round(`TP change point`, 2),
    `Z-score` = round(`Z-score`, 2),
    `p-value` = signif(`p-value`, 3),
    Purity = round(Purity, 3),
    `Reliability` = round(`Reliability`, 3)
  )

print(candidate_taxa_TP)


#Run TITAN for NP
#####NP##########################################################################
set.seed(123)
Titan_NP <- titan(
  env = NP_vec,
  txa = comm_TITAN,
  minSplt = 5,
  numPerm = 999,
  boot = TRUE,
  nBoot = 999,
  imax = FALSE
)

summary(Titan_NP)
plot_sumz(Titan_NP)
plot_taxa(Titan_NP)

##print top 5 taxa candidates
Titan_NP_topTaxa <- as.data.frame(Titan_NP$sppmax) %>%
  rownames_to_column(var = "Taxon") %>%
  mutate(
    across(
      c(ienv.cp, zenv.cp, maxgrp, obsiv.prob,
        zscore, purity, reliability, filter),
      as.numeric
    )
  )
#screening thresholds for top 5
min_zscore <- 1.50
max_perm_p <- 0.12
min_purity <- 0.70
min_reliability <- 0.30

#select & rank taxa
candidate_taxa_NP <- Titan_NP_topTaxa %>%
  dplyr::filter(
    zscore >= min_zscore,
    obsiv.prob <= max_perm_p,
    purity >= min_purity,
    reliability >= min_reliability
  ) %>%
  mutate(
    Direction = case_when(
      maxgrp == 1 ~ "Decreaser",
      maxgrp == 2 ~ "Increaser",
      TRUE ~ "Unknown"
    )
  ) %>%
  arrange(
    desc(zscore),
    obsiv.prob,
    desc(purity),
    desc(reliability)
  ) %>%
  slice_head(n = 5) %>%
  transmute(
    Taxon,
    Direction,
    `NP change point` = zenv.cp,
    `Z-score` = zscore,
    `p-value` = obsiv.prob,
    Purity = purity,
    `Reliability` = reliability
  ) %>%
  mutate(
    `NP change point` = round(`NP change point`, 2),
    `Z-score` = round(`Z-score`, 2),
    `p-value` = signif(`p-value`, 3),
    Purity = round(Purity, 3),
    `Reliability` = round(`Reliability`, 3)
  )

print(candidate_taxa_NP)

####Continuous Figure#####################################################
make_continuous_plots <- function(data, nutrient_name) {
  
  ggplot(
    filter(data, Nutrient == nutrient_name),
    aes(
      x = Change_point,
      y = Taxon,
      color = Direction
    )
  ) +
    geom_segment(
      aes(
        x = Lower,
        xend = Upper,
        yend = Taxon
      ),
      linewidth = 0.7,
      alpha = 0.7
    ) +
    geom_point(size = 3) +
    scale_color_manual(values = c(
      "Decreaser" = "#D66B35",
      "Increaser" = "#2878A5"
    )) +
    labs(
      title = nutrient_name,
      x = "Estimated nutrient change point (mg/L)",
      y = "Taxon",
      color = "Response direction"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 12,
        hjust = 0.5
      ),
      axis.text.y = element_text(
        face = "italic",
        size = 10
      ),
      axis.text.x = element_text(size = 9),
      axis.title = element_text(size = 10),
      axis.line = element_line(
        colour = "black",
        linewidth = 0.3
      ),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.3
      ),
      legend.position = "bottom"
    )
}

#create plots and remove axis titles to only have 1 x and 1 y axis title
Figure1_TN <- make_continuous_plots(titan_top5, "TN")
Figure1_TN <- Figure1_TN +
  labs(x = NULL, y = NULL)

Figure1_TP <- make_continuous_plots(titan_top5, "TP")
Figure1_TP <- Figure1_TP +
  labs(x = NULL)

Figure1_NP <- make_continuous_plots(titan_top5, "N:P")
Figure1_NP <- Figure1_NP +
  labs(y = NULL)

#patchwork plot
ContinuousTitan <- Figure1_TN /
  Figure1_TP /
  Figure1_NP +
  plot_layout(
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )
ContinuousTitan

#increase plot box line width
ContinuousTitan <- ContinuousTitan &
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.65,
      linetype = "solid"
    )
  )
ContinuousTitan

ggsave("ContinuousTitanPanel.png", width = 6, height = 6, units = "in", dpi = 300)

####Complementary temporal TITAN#################################################
#create a biomass matrix for each sampling round

temp_comm_matrix <- function(zoop_biomass, round_number) {
  
  comm <- zoop_biomass %>%
    filter(Round == round_number) %>%
    group_by(Site, Taxon) %>%
    summarise(
      mean_biomass = mean(TotalBiomassPerSample, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = Taxon,
      values_from = mean_biomass,
      values_fill = 0
    ) %>%
    as.data.frame()
  
  rownames(comm) <- comm$Site
  comm$Site <- NULL
  
  as.matrix(comm)
}

#matrices for each sampling round
R_1 <- temp_comm_matrix(zoop_biomass, 1)
R_2 <- temp_comm_matrix(zoop_biomass, 2)
R_3 <- temp_comm_matrix(zoop_biomass, 3)

#species occurrence filter (retain taxa occurring in at least three sites)
R_1_filt <- R_1[, colSums(R_1 > 0) >=3, drop = FALSE]
R_2_filt <- R_2[, colSums(R_2 > 0) >=3, drop = FALSE]
R_3_filt <- R_3[, colSums(R_3 > 0) >=3, drop = FALSE]

#direct comparison of same taxa across sampling rounds, identify taxa that passed
#the filter in all three rounds
common_taxa <- Reduce(
  intersect,
  list(
    colnames(R_1_filt),
    colnames(R_2_filt),
    colnames(R_3_filt)
  )
)

R_1_filt <- R_1_filt[, common_taxa, drop = FALSE]
R_2_filt <- R_2_filt[, common_taxa, drop = FALSE]
R_3_filt <- R_3_filt[, common_taxa, drop = FALSE]

#check
dim(R_1_filt)
dim(R_2_filt)
dim(R_3_filt)

identical(colnames(R_1_filt), colnames(R_2_filt))
identical(colnames(R_2_filt), colnames(R_3_filt))

##create nutrient datasets for round x nutrient comparisons
#TN
TN_R1 <- Nutrients %>%
  filter(Round == 1) %>%
  select(Site, Total.Nitrogen)

TN_R2 <- Nutrients %>%
  filter(Round == 2) %>%
  select(Site, Total.Nitrogen)

TN_R3 <- Nutrients %>%
  filter(Round == 3) %>%
  select(Site, Total.Nitrogen)

# #check
# nrow(TN_R1)
# nrow(TN_R2)
# nrow(TN_R3)
# 
# head(TN_R1)
# head(TN_R2)
# head(TN_R3)

#TP
TP_R1 <- Nutrients %>%
  filter(Round == 1) %>%
  select(Site, Total.Phosphorus)

TP_R2 <- Nutrients %>%
  filter(Round == 2) %>%
  select(Site, Total.Phosphorus)

TP_R3 <- Nutrients %>%
  filter(Round == 3) %>%
  select(Site, Total.Phosphorus)

#NP
NP_R1 <- Nutrients %>%
  filter(Round == 1) %>%
  select(Site, NP.ratio)

NP_R2 <- Nutrients %>%
  filter(Round == 2) %>%
  select(Site, NP.ratio)

NP_R3 <- Nutrients %>%
  filter(Round == 3) %>%
  select(Site, NP.ratio)


##align nutrient values with community matrices
##TN
align_TN <- function(comm, Nutrients) {
  site_match <- match(rownames(comm), Nutrients$Site)
  
  if (anyNA(site_match)) {
    stop("Some community matrix sites are missing from the nutrient dataset.")
  }
  
  TN <- Nutrients$Total.Nitrogen[site_match]
  names(TN) <- rownames(comm)
  
  
  return(TN)
}

#apply
TN_R1 <- align_TN(R_1_filt, TN_R1)
TN_R2 <- align_TN(R_2_filt, TN_R2)
TN_R3 <- align_TN(R_3_filt, TN_R3)

#verify
identical(rownames(R_1_filt), names(TN_R1))
identical(rownames(R_2_filt), names(TN_R2))
identical(rownames(R_3_filt), names(TN_R3))

sum(is.na(TN_R1))
sum(is.na(TN_R2))
sum(is.na(TN_R3))


##TP
align_TP <- function(comm, Nutrients) {
  site_match <- match(rownames(comm), Nutrients$Site)
  
  if (anyNA(site_match)) {
    stop("Some community matrix sites are missing from the nutrient dataset.")
  }
  
  TP <- Nutrients$Total.Phosphorus[site_match]
  names(TP) <- rownames(comm)
  
  
  return(TP)
}

#apply
TP_R1 <- align_TP(R_1_filt, TP_R1)
TP_R2 <- align_TP(R_2_filt, TP_R2)
TP_R3 <- align_TP(R_3_filt, TP_R3)

#verify
identical(rownames(R_1_filt), names(TP_R1))
identical(rownames(R_2_filt), names(TP_R2))
identical(rownames(R_3_filt), names(TP_R3))

sum(is.na(TP_R1))
sum(is.na(TP_R2))
sum(is.na(TP_R3))


##NP
align_NP <- function(comm, Nutrients) {
  site_match <- match(rownames(comm), Nutrients$Site)
  
  if (anyNA(site_match)) {
    stop("Some community matrix sites are missing from the nutrient dataset.")
  }
  
  NP <- Nutrients$NP.ratio[site_match]
  names(NP) <- rownames(comm)
  
  
  return(NP)
}

#apply
NP_R1 <- align_NP(R_1_filt, NP_R1)
NP_R2 <- align_NP(R_2_filt, NP_R2)
NP_R3 <- align_NP(R_3_filt, NP_R3)

#verify
identical(rownames(R_1_filt), names(NP_R1))
identical(rownames(R_2_filt), names(NP_R2))
identical(rownames(R_3_filt), names(NP_R3))

sum(is.na(NP_R1))
sum(is.na(NP_R2))
sum(is.na(NP_R3))



###Run TITAN2
set.seed(123)

#####TN##########################################################################
######R1#########################################################################
Titan_R1_TN <- titan(
  env = TN_R1,
  txa = R_1_filt,
  minSplt = 3,
  numPerm = 999,
  boot = TRUE,
  nBoot = 999,
  imax = FALSE
)
summary(Titan_R1_TN)
Titan_R1_TN$sppmax

##print top 5 taxa candidates
Titan_R1_TN_topTaxa <- as.data.frame(Titan_R1_TN$sppmax) %>%
  rownames_to_column(var = "Taxon") %>%
  mutate(
    across(
      c(ienv.cp, zenv.cp, maxgrp, obsiv.prob,
        zscore, purity, reliability, filter),
      as.numeric
    )
  )
#screening thresholds for top 5
min_zscore <- 1.50
max_perm_p <- 0.12
min_purity <- 0.70
min_reliability <- 0.30

#select & rank taxa
candidate_taxa_R1_TN <- Titan_R1_TN_topTaxa %>%
  dplyr::filter(
    zscore >= min_zscore,
    obsiv.prob <= max_perm_p,
    purity >= min_purity,
    reliability >= min_reliability
  ) %>%
  mutate(
    Direction = case_when(
      maxgrp == 1 ~ "Decreaser",
      maxgrp == 2 ~ "Increaser",
      TRUE ~ "Unknown"
    )
  ) %>%
  arrange(
    desc(zscore),
    obsiv.prob,
    desc(purity),
    desc(reliability)
  ) %>%
  slice_head(n = 5) %>%
  transmute(
    Taxon,
    Direction,
    `TN change point` = zenv.cp,
    `Z-score` = zscore,
    `p-value` = obsiv.prob,
    Purity = purity,
    `Reliability` = reliability
  ) %>%
  mutate(
    `TN change point` = round(`TN change point`, 2),
    `Z-score` = round(`Z-score`, 2),
    `p-value` = signif(`p-value`, 3),
    Purity = round(Purity, 3),
    `Reliability` = round(`Reliability`, 3)
  )

print(candidate_taxa_R1_TN)

######R2#########################################################################
Titan_R2_TN <- titan(
  env = TN_R2,
  txa = R_2_filt,
  minSplt = 3,
  numPerm = 999,
  boot = TRUE,
  nBoot = 999,
  imax = FALSE
)
summary(Titan_R2_TN)
Titan_R2_TN$sppmax

##print top 5 taxa candidates
Titan_R2_TN_topTaxa <- as.data.frame(Titan_R2_TN$sppmax) %>%
  rownames_to_column(var = "Taxon") %>%
  mutate(
    across(
      c(ienv.cp, zenv.cp, maxgrp, obsiv.prob,
        zscore, purity, reliability, filter),
      as.numeric
    )
  )

#screening thresholds for top 5
min_zscore <- 1.50
max_perm_p <- 0.12
min_purity <- 0.70
min_reliability <- 0.30

#select & rank taxa
candidate_taxa_R2_TN <- Titan_R2_TN_topTaxa %>%
  dplyr::filter(
    zscore >= min_zscore,
    obsiv.prob <= max_perm_p,
    purity >= min_purity,
    reliability >= min_reliability
  ) %>%
  mutate(
    Direction = case_when(
      maxgrp == 1 ~ "Decreaser",
      maxgrp == 2 ~ "Increaser",
      TRUE ~ "Unknown"
    )
  ) %>%
  arrange(
    desc(zscore),
    obsiv.prob,
    desc(purity),
    desc(reliability)
  ) %>%
  slice_head(n = 5) %>%
  transmute(
    Taxon,
    Direction,
    `TN change point` = zenv.cp,
    `Z-score` = zscore,
    `p-value` = obsiv.prob,
    Purity = purity,
    `Reliability` = reliability
  ) %>%
  mutate(
    `TN change point` = round(`TN change point`, 2),
    `Z-score` = round(`Z-score`, 2),
    `p-value` = signif(`p-value`, 3),
    Purity = round(Purity, 3),
    `Reliability` = round(`Reliability`, 3)
  )

print(candidate_taxa_R2_TN)

######R3#########################################################################
Titan_R3_TN <- titan(
  env = TN_R3,
  txa = R_3_filt,
  minSplt = 3,
  numPerm = 999,
  boot = TRUE,
  nBoot = 999,
  imax = FALSE
)
summary(Titan_R3_TN)
Titan_R3_TN$sppmax

##print top 5 taxa candidates
Titan_R3_TN_topTaxa <- as.data.frame(Titan_R3_TN$sppmax) %>%
  rownames_to_column(var = "Taxon") %>%
  mutate(
    across(
      c(ienv.cp, zenv.cp, maxgrp, obsiv.prob,
        zscore, purity, reliability, filter),
      as.numeric
    )
  )

#screening thresholds for top 5
min_zscore <- 1.50
max_perm_p <- 0.12
min_purity <- 0.70
min_reliability <- 0.30

#select & rank taxa
candidate_taxa_R3_TN <- Titan_R3_TN_topTaxa %>%
  dplyr::filter(
    zscore >= min_zscore,
    obsiv.prob <= max_perm_p,
    purity >= min_purity,
    reliability >= min_reliability
  ) %>%
  mutate(
    Direction = case_when(
      maxgrp == 1 ~ "Decreaser",
      maxgrp == 2 ~ "Increaser",
      TRUE ~ "Unknown"
    )
  ) %>%
  arrange(
    desc(zscore),
    obsiv.prob,
    desc(purity),
    desc(reliability)
  ) %>%
  slice_head(n = 5) %>%
  transmute(
    Taxon,
    Direction,
    `TN change point` = zenv.cp,
    `Z-score` = zscore,
    `p-value` = obsiv.prob,
    Purity = purity,
    `Reliability` = reliability
  ) %>%
  mutate(
    `TN change point` = round(`TN change point`, 2),
    `Z-score` = round(`Z-score`, 2),
    `p-value` = signif(`p-value`, 3),
    Purity = round(Purity, 3),
    `Reliability` = round(`Reliability`, 3)
  )

print(candidate_taxa_R3_TN)


#####TP##########################################################################
######R1#########################################################################
Titan_R1_TP <- titan(
  env = TP_R1,
  txa = R_1_filt,
  minSplt = 3,
  numPerm = 999,
  boot = TRUE,
  nBoot = 999,
  imax = FALSE
)
summary(Titan_R1_TP)
Titan_R1_TP$sppmax

##print top 5 taxa candidates
Titan_R1_TP_topTaxa <- as.data.frame(Titan_R1_TP$sppmax) %>%
  rownames_to_column(var = "Taxon") %>%
  mutate(
    across(
      c(ienv.cp, zenv.cp, maxgrp, obsiv.prob,
        zscore, purity, reliability, filter),
      as.numeric
    )
  )

#screening thresholds for top 5
min_zscore <- 1.50
max_perm_p <- 0.12
min_purity <- 0.70
min_reliability <- 0.30

#select & rank taxa
candidate_taxa_R1_TP <- Titan_R1_TP_topTaxa %>%
  dplyr::filter(
    zscore >= min_zscore,
    obsiv.prob <= max_perm_p,
    purity >= min_purity,
    reliability >= min_reliability
  ) %>%
  mutate(
    Direction = case_when(
      maxgrp == 1 ~ "Decreaser",
      maxgrp == 2 ~ "Increaser",
      TRUE ~ "Unknown"
    )
  ) %>%
  arrange(
    desc(zscore),
    obsiv.prob,
    desc(purity),
    desc(reliability)
  ) %>%
  slice_head(n = 5) %>%
  transmute(
    Taxon,
    Direction,
    `TP change point` = zenv.cp,
    `Z-score` = zscore,
    `p-value` = obsiv.prob,
    Purity = purity,
    `Reliability` = reliability
  ) %>%
  mutate(
    `TP change point` = round(`TP change point`, 2),
    `Z-score` = round(`Z-score`, 2),
    `p-value` = signif(`p-value`, 3),
    Purity = round(Purity, 3),
    `Reliability` = round(`Reliability`, 3)
  )

print(candidate_taxa_R1_TP)

######R2#########################################################################
Titan_R2_TP <- titan(
  env = TP_R2,
  txa = R_2_filt,
  minSplt = 3,
  numPerm = 999,
  boot = TRUE,
  nBoot = 999,
  imax = FALSE
)
summary(Titan_R2_TP)
Titan_R2_TP$sppmax

##print top 5 taxa candidates
Titan_R2_TP_topTaxa <- as.data.frame(Titan_R2_TP$sppmax) %>%
  rownames_to_column(var = "Taxon") %>%
  mutate(
    across(
      c(ienv.cp, zenv.cp, maxgrp, obsiv.prob,
        zscore, purity, reliability, filter),
      as.numeric
    )
  )

#screening thresholds for top 5
min_zscore <- 1.50
max_perm_p <- 0.12
min_purity <- 0.70
min_reliability <- 0.30

#select & rank taxa
candidate_taxa_R2_TP <- Titan_R2_TP_topTaxa %>%
  dplyr::filter(
    zscore >= min_zscore,
    obsiv.prob <= max_perm_p,
    purity >= min_purity,
    reliability >= min_reliability
  ) %>%
  mutate(
    Direction = case_when(
      maxgrp == 1 ~ "Decreaser",
      maxgrp == 2 ~ "Increaser",
      TRUE ~ "Unknown"
    )
  ) %>%
  arrange(
    desc(zscore),
    obsiv.prob,
    desc(purity),
    desc(reliability)
  ) %>%
  slice_head(n = 5) %>%
  transmute(
    Taxon,
    Direction,
    `TP change point` = zenv.cp,
    `Z-score` = zscore,
    `p-value` = obsiv.prob,
    Purity = purity,
    `Reliability` = reliability
  ) %>%
  mutate(
    `TP change point` = round(`TP change point`, 2),
    `Z-score` = round(`Z-score`, 2),
    `p-value` = signif(`p-value`, 3),
    Purity = round(Purity, 3),
    `Reliability` = round(`Reliability`, 3)
  )

print(candidate_taxa_R2_TP)

######R3#########################################################################
Titan_R3_TP <- titan(
  env = TP_R3,
  txa = R_3_filt,
  minSplt = 3,
  numPerm = 999,
  boot = TRUE,
  nBoot = 999,
  imax = FALSE
)
summary(Titan_R3_TP)
Titan_R3_TP$sppmax

##print top 5 taxa candidates
Titan_R3_TP_topTaxa <- as.data.frame(Titan_R3_TP$sppmax) %>%
  rownames_to_column(var = "Taxon") %>%
  mutate(
    across(
      c(ienv.cp, zenv.cp, maxgrp, obsiv.prob,
        zscore, purity, reliability, filter),
      as.numeric
    )
  )

#screening thresholds for top 5
min_zscore <- 1.50
max_perm_p <- 0.12
min_purity <- 0.70
min_reliability <- 0.30

#select & rank taxa
candidate_taxa_R3_TP <- Titan_R3_TP_topTaxa %>%
  dplyr::filter(
    zscore >= min_zscore,
    obsiv.prob <= max_perm_p,
    purity >= min_purity,
    reliability >= min_reliability
  ) %>%
  mutate(
    Direction = case_when(
      maxgrp == 1 ~ "Decreaser",
      maxgrp == 2 ~ "Increaser",
      TRUE ~ "Unknown"
    )
  ) %>%
  arrange(
    desc(zscore),
    obsiv.prob,
    desc(purity),
    desc(reliability)
  ) %>%
  slice_head(n = 5) %>%
  transmute(
    Taxon,
    Direction,
    `TP change point` = zenv.cp,
    `Z-score` = zscore,
    `p-value` = obsiv.prob,
    Purity = purity,
    `Reliability` = reliability
  ) %>%
  mutate(
    `TP change point` = round(`TP change point`, 2),
    `Z-score` = round(`Z-score`, 2),
    `p-value` = signif(`p-value`, 3),
    Purity = round(Purity, 3),
    `Reliability` = round(`Reliability`, 3)
  )

print(candidate_taxa_R3_TP)


#####TN:TP#######################################################################
######R1#########################################################################
Titan_R1_NP <- titan(
  env = NP_R1,
  txa = R_1_filt,
  minSplt = 3,
  numPerm = 999,
  boot = TRUE,
  nBoot = 999,
  imax = FALSE
)
summary(Titan_R1_NP)
Titan_R1_NP$sppmax

##print top 5 taxa candidates
Titan_R1_NP_topTaxa <- as.data.frame(Titan_R1_NP$sppmax) %>%
  rownames_to_column(var = "Taxon") %>%
  mutate(
    across(
      c(ienv.cp, zenv.cp, maxgrp, obsiv.prob,
        zscore, purity, reliability, filter),
      as.numeric
    )
  )

#screening thresholds for top 5
min_zscore <- 1.50
max_perm_p <- 0.12
min_purity <- 0.70
min_reliability <- 0.30

#select & rank taxa
candidate_taxa_R1_NP <- Titan_R1_NP_topTaxa %>%
  dplyr::filter(
    zscore >= min_zscore,
    obsiv.prob <= max_perm_p,
    purity >= min_purity,
    reliability >= min_reliability
  ) %>%
  mutate(
    Direction = case_when(
      maxgrp == 1 ~ "Decreaser",
      maxgrp == 2 ~ "Increaser",
      TRUE ~ "Unknown"
    )
  ) %>%
  arrange(
    desc(zscore),
    obsiv.prob,
    desc(purity),
    desc(reliability)
  ) %>%
  slice_head(n = 5) %>%
  transmute(
    Taxon,
    Direction,
    `NP change point` = zenv.cp,
    `Z-score` = zscore,
    `p-value` = obsiv.prob,
    Purity = purity,
    `Reliability` = reliability
  ) %>%
  mutate(
    `NP change point` = round(`NP change point`, 2),
    `Z-score` = round(`Z-score`, 2),
    `p-value` = signif(`p-value`, 3),
    Purity = round(Purity, 3),
    `Reliability` = round(`Reliability`, 3)
  )

print(candidate_taxa_R1_NP)

######R2#########################################################################
Titan_R2_NP <- titan(
  env = NP_R2,
  txa = R_2_filt,
  minSplt = 3,
  numPerm = 999,
  boot = TRUE,
  nBoot = 999,
  imax = FALSE
)
summary(Titan_R2_NP)
Titan_R2_NP$sppmax

##print top 5 taxa candidates
Titan_R2_NP_topTaxa <- as.data.frame(Titan_R2_NP$sppmax) %>%
  rownames_to_column(var = "Taxon") %>%
  mutate(
    across(
      c(ienv.cp, zenv.cp, maxgrp, obsiv.prob,
        zscore, purity, reliability, filter),
      as.numeric
    )
  )

#screening thresholds for top 5
min_zscore <- 1.50
max_perm_p <- 0.12
min_purity <- 0.70
min_reliability <- 0.30

#select & rank taxa
candidate_taxa_R2_NP <- Titan_R2_NP_topTaxa %>%
  dplyr::filter(
    zscore >= min_zscore,
    obsiv.prob <= max_perm_p,
    purity >= min_purity,
    reliability >= min_reliability
  ) %>%
  mutate(
    Direction = case_when(
      maxgrp == 1 ~ "Decreaser",
      maxgrp == 2 ~ "Increaser",
      TRUE ~ "Unknown"
    )
  ) %>%
  arrange(
    desc(zscore),
    obsiv.prob,
    desc(purity),
    desc(reliability)
  ) %>%
  slice_head(n = 5) %>%
  transmute(
    Taxon,
    Direction,
    `NP change point` = zenv.cp,
    `Z-score` = zscore,
    `p-value` = obsiv.prob,
    Purity = purity,
    `Reliability` = reliability
  ) %>%
  mutate(
    `NP change point` = round(`NP change point`, 2),
    `Z-score` = round(`Z-score`, 2),
    `p-value` = signif(`p-value`, 3),
    Purity = round(Purity, 3),
    `Reliability` = round(`Reliability`, 3)
  )

print(candidate_taxa_R2_NP)

######R3#########################################################################
Titan_R3_NP <- titan(
  env = NP_R3,
  txa = R_3_filt,
  minSplt = 3,
  numPerm = 999,
  boot = TRUE,
  nBoot = 999,
  imax = FALSE
)
summary(Titan_R3_NP)
Titan_R3_NP$sppmax

##print top 5 taxa candidates
Titan_R3_NP_topTaxa <- as.data.frame(Titan_R3_NP$sppmax) %>%
  rownames_to_column(var = "Taxon") %>%
  mutate(
    across(
      c(ienv.cp, zenv.cp, maxgrp, obsiv.prob,
        zscore, purity, reliability, filter),
      as.numeric
    )
  )

#screening thresholds for top 5
min_zscore <- 1.50
max_perm_p <- 0.12
min_purity <- 0.70
min_reliability <- 0.30

#select & rank taxa
candidate_taxa_R3_NP <- Titan_R3_NP_topTaxa %>%
  dplyr::filter(
    zscore >= min_zscore,
    obsiv.prob <= max_perm_p,
    purity >= min_purity,
    reliability >= min_reliability
  ) %>%
  mutate(
    Direction = case_when(
      maxgrp == 1 ~ "Decreaser",
      maxgrp == 2 ~ "Increaser",
      TRUE ~ "Unknown"
    )
  ) %>%
  arrange(
    desc(zscore),
    obsiv.prob,
    desc(purity),
    desc(reliability)
  ) %>%
  slice_head(n = 5) %>%
  transmute(
    Taxon,
    Direction,
    `NP change point` = zenv.cp,
    `Z-score` = zscore,
    `p-value` = obsiv.prob,
    Purity = purity,
    `Reliability` = reliability
  ) %>%
  mutate(
    `NP change point` = round(`NP change point`, 2),
    `Z-score` = round(`Z-score`, 2),
    `p-value` = signif(`p-value`, 3),
    Purity = round(Purity, 3),
    `Reliability` = round(`Reliability`, 3)
  )

print(candidate_taxa_R3_NP)


####Temporal Figures############################################################
make_temporal_plot <- function(data, nutrient_name) {
  
  ggplot(
    filter(data, Nutrient == nutrient_name),
    aes(
      x = Month,
      y = Change_point,
      group = Taxon,
      color = Direction
    )
  ) +
    geom_errorbar(
      aes(ymin = Lower, ymax = Upper),
      width = 0.15,
      linewidth = 0.5,
      alpha = 0.6
    ) +
    geom_line(
      color = "grey50",
      linewidth = 0.5
    ) +
    geom_point(size = 2.5) +
    facet_wrap(
      ~ Taxon,
      scales = "free_y",
      ncol = 5
    ) +
    scale_color_manual(values = c(
      "Decreaser" = "#D66B35",
      "Increaser" = "#2878A5"
    )) +
    labs(
      title = nutrient_name,
      x = "Sampling month",
      y = "Estimated change point (mg/L)",
      color = "Response direction"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 13,
        hjust = 0.5
      ),
      strip.text = element_text(
        face = "italic",
        size = 8
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 8
      ),
      axis.text.y = element_text(size = 8),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.4
      ),
      panel.spacing = unit(0.5, "lines"),
      legend.position = "bottom"
    )
}

Figure2_TN <- make_temporal_plot(
  temporal_plot_data, "TN"
)
Figure2_TN
ggsave("TemporalTITAN_TN.png", width = 7, height = 7, units = "in", dpi = 300)


Figure2_TP <- make_temporal_plot(
  temporal_plot_data, "TP"
)
Figure2_TP
ggsave("TemporalTITAN_TP.png", width = 7, height = 7, units = "in", dpi = 300)


Figure2_NP <- make_temporal_plot(
  temporal_plot_data, "N:P"
)
Figure2_NP
ggsave("TemporalTITAN_NP.png", width = 7, height = 7, units = "in", dpi = 300)


###Wetland site map#############################################################
#Map libraries
library(sf)
library(ggspatial)
library(prettymapr)
library(raster) #masks select function from dplyr
library(tigris)
library(patchwork)
library(rosm)
library(tmap)

####All Wetland Sites####
AllwetlandLocs = read.csv("Wetland_sites.csv", header = T, sep = ",")
wetlands_sf <- st_as_sf(
  AllwetlandLocs,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)
options(tigris_use_cache = TRUE)

#classifying wetland type
wetlands_sf <- wetlands_sf %>%
  mutate(
    WetlandType = case_when(
      ISU_Name %in% c("DD15N", "SS", "DD8", "DD65", "LP", "WW", "FCC2",
                      "KS", "RS") ~ "CREP",
      ISU_Name %in% c("UH1", "UH2") ~ "Restored for habitat",
      ISU_Name %in% c("Bjork") ~ "Remnant"
    )
  )
# Iowa state boundary
iowa <- states(cb = TRUE) |>
  dplyr::filter(NAME == "Iowa") |>
  st_transform(4326)
# Iowa counties
iowa_counties <- counties(state = "IA", cb = TRUE) |>
  st_transform(4326)

Sampling_Map8.5 <- ggplot() +
  geom_sf(data = iowa, fill = NA, color = "transparent") +
  # Basemap
  annotation_map_tile(type = "cartolight", zoom = 9, alpha = 0.5) +
  coord_sf(
    xlim = c(-96.8, -90.1),
    ylim = c(40.3, 43.6),
    expand = FALSE
  ) +
  # Counties
  geom_sf(data = iowa_counties, fill = NA, color = "grey50") +
  # State boundary
  geom_sf(data = iowa, fill = NA, color = "black") +
  # Wetland sites
  geom_sf(
    data = wetlands_sf,
    aes(shape = WetlandType),
    size = 2.5
  ) +
  # #scale_color_manual(
  #   values = c(
  #     "CREP" = "red",
  #     "Remnant" = "green",
  #     "Restored for habitat" = "blue"
  #   )
  # ) +
  # Labels
  geom_text_repel(
    data = wetlands_sf,
    aes(x = lon, y = lat, label = ISU_Name),
    size = 2.5
  ) +
  annotation_scale(
    location = "bl",
    height = unit(0.20, "cm"),
    width_hint = 0.20,
    text_cex = 0.65
  ) +
  annotation_north_arrow(
    location = "tr",
    height = unit(0.75, "cm"),
    width = unit (0.75, "cm"),
    style = north_arrow_fancy_orienteering(
      text_size = 6
    )
  ) +
  coord_sf(
    xlim = c(-96, -92),
    ylim = c(41.8, 43.5),
    expand = FALSE
  ) +
  labs(
    x = "Longitude",
    y = "Latitude",
    shape = "Wetland Type"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom"
  )
Sampling_Map8.5
#ggsave("Sampling_Map9.png", width = 7, height = 4, units = "in", dpi = 300)

####Extent Inset Box Map####
#same limits as all wetland map
xmin <- -96
xmax <- -92
ymin <- 41.8
ymax <- 43.5

#turn into polygon
study_box <- st_as_sfc(
  st_bbox(
    c(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    ),
    crs = 4326
  )
)
options(tigris_use_cache = TRUE)

# Iowa state boundary
iowa <- states(cb = TRUE) |>
  filter(NAME == "Iowa") |>
  st_transform(4326)

IowaStudyPlot <- ggplot() +
  # Iowa outline
  geom_sf(
    data = iowa,
    fill = "white",
    color = "black",
    linewidth = 0.8
  ) +  
  # Study area box
  geom_sf(
    data = study_box,
    fill = NA,
    color = "red",
    linewidth = 2
  ) +
  theme_minimal() +
  labs(
    title = "Study Area Extent in Iowa",
    x = "Longitude",
    y = "Latitude"
  )

IowaStudyPlot
#ggsave("IowaStudyPlot.png", width = 7, height = 4, dpi = 300)
