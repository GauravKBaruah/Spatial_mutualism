
rm(list=ls())
require(deSolve) ## for integrating ordinary differential equations
require(tidyverse) ## for efficient data manipulation & plotting
require(cowplot) ## for arranging plots in a grid
library(dplyr)
library(readr)
library(ggplot2)
library(grid)
library(gridExtra)

load("Mutualism_data.RData")



fact %>% 
  ggplot(aes(x=factor(Strength_mutualism), y= (Biomass.animal+Biomass.plant), color=interaction))+
  geom_boxplot()+xlab("Mutualism strength")+ylab("Community abundance")+
  theme_classic()+facet_wrap(.~individual.variation,ncol=2)+ylim(c(0,30))

fact %>% 
  ggplot(aes(x=factor(Strength_mutualism), y= (trait.matching), color=interaction))+
  geom_boxplot()+xlab("Mutualism strength")+ylab("Trait matching")+
  theme_classic()+facet_wrap(.~individual.variation,ncol=2)

# 
# fact %>% 
#   ggplot(aes(x=factor(Strength_mutualism), y= (CV.trait), color=interaction))+
#   geom_boxplot()+xlab("Mutualism strength")+ylab("Trait patterning")+
#   theme_classic()+facet_wrap(.~individual.variation,ncol=2)

fact %>% 
  ggplot(aes(x=factor(Strength_mutualism), y= (richness), color=interaction))+
  geom_boxplot()+xlab("Mutualism strength")+ylab("Species richness")+
  theme_classic()+facet_wrap(.~individual.variation,ncol=2)


