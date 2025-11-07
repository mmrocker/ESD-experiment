
title: "_ESD_plots"
author: "MMRocker"
date: "13 August 2025"
output: html_document


#set up working file & R

getwd() #check working directory


library(ggplot2)
library(readxl)
library(lubridate)
require(grid)
library(gridExtra)
library(reshape2)
library(plyr)
library(ggpattern)
library(dplyr)


df<-read_excel("ESDExpt_PopPrev.xlsx",sheet="ESD_long")



df$Sampling_Date <- mdy(df$Sampling_Date)
df$Sampling_Day <- as.numeric(df$Sampling_Date)
df$Sampling_Day <- as.numeric(df$Sampling_Day) - 19905
df$Container <- as.factor(df$Container)
df$Treatment <- as.factor(df$Treatment)
df$Table <- as.factor(df$Table)
df$Source_Pop <- as.factor(df$Source_Pop)
df$ESD_Stage_photo <- as.factor(df$ESD_Stage_photo)
df$Molt_Stage <- as.factor(df$Molt_Stage)
df$ESDxMolt_Stage <- as.factor(df$ESDxMolt_Stage)
df$Source_Pop <- factor(df$Source_Pop, levels = c("ME", "MA", "CT"))


df$TrtPop <- paste(df$Treatment, df$Source_Pop)
df$TrtPop <- factor(df$TrtPop, levels = c("cold ME",  "cold MA", "cold CT",
                                          "warm ME","warm MA", "warm CT" ))
df$ESD_Stage_photo <- factor(df$ESD_Stage_photo, levels = c( "Dead","Healthy", "Mild", "Moderate", 
                                                "Severe"))
df$Molt_Stage <- factor(df$Molt_Stage, levels = c("Dead", "Molt", "Healthy",
                                                  "Mild", "Moderate", "Severe"))
df$ESDxMolt_Stage <- factor(df$ESDxMolt_Stage, levels = c("Dead", "Molt-Healthy", "Healthy",
                                                      "Molt-Mild",
                                                      "Mild", "Moderate", "Severe"))



#rename & reorder states
df$Source_Pop <- revalue(df$Source_Pop, c("CT" = "Niantic Bay, CT",
                                                        "MA" = "Buzzards Bay, MA",
                                                        "ME" = "Casco Bay, ME"))
df$Source_Pop <- factor(df$Source_Pop, levels = c("Casco Bay, ME",
                                                  "Buzzards Bay, MA",
                                                  "Niantic Bay, CT"))

#rename temp treatments
df$Treatment <- revalue(df$Treatment, c("cold" = "cold (13°C)",
                                        "warm" = "warm (20°C)"))
                        
#rename & reorder Sampling TimePoints
df$Sampling_TP <- df$Sampling_Day
df$Sampling_Day <- as.factor(df$Sampling_Day)
df$Sampling_Day <- revalue(df$Sampling_Day, c("140" = "141"))
df$Sampling_TP <- as.factor(df$Sampling_TP)
df$Sampling_TP <- revalue(df$Sampling_TP, c("1" = "1",
                                          "29" = "2",
                                          "57" = "3",
                                          "87" = "4",
                                          "115" = "5",
                                          "140" = "6",
                                          "141" = "6"))
df$Sampling_TP <- as.factor(df$Sampling_TP)


df$Molt <- df$Molt_Stage
df$Molt <- revalue(df$Molt, c("Dead" = "Dead", "Healthy" = "none", 
                                    "Molt" = "Molted", "Mild" = "none", 
                                    "Moderate" = "none", "Severe" = "none"))

df$Molt <- factor(df$Molt, levels = c("Molted", "none", "Dead"))

df$ESD_Stage_photo2 <- df$ESD_Stage_photo
df$ESD_Stage_photo2 <- revalue(df$ESD_Stage_photo2, c("Dead" = "Dead", "Healthy" = "Healthy", 
                              "Mild" = "Healthy", "Moderate" = "Moderate", "Severe" = "Severe"))


df2 <- df %>%
  filter(!TAG %in% c("3781", "3746")) #remove lobsters that molted within first 3 weeks of expt (n=2 -- 3781 & 3746)


#disease plot
PLOT7 <- ggplot(df2, aes(x = Sampling_Day, fill = ESD_Stage_photo, pattern = Molt)) +
  geom_bar_pattern(position = "fill", pattern_fill= "black") +
  scale_fill_manual(
    values = c("black", viridis::magma(4, direction = 1, begin = 0.97, end = 0.4)),
    labels = c("Healthy" = "No disease") ) +
  scale_pattern_manual(
    values = c("circle", "none", "none"),
    labels = c("none" = " ", "Dead" = " ") ) +
  facet_grid(rows = vars(Treatment), cols = vars(Source_Pop)) +
  labs(
    y = "Proportion of lobster",
    x = "Day of experiment",
    fill = "ESD Index" ) +
  theme_bw(base_size = 12) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(),
    legend.position = "right",
    legend.box = "vertical",
    legend.box.just = "left",         # Aligns stacked legends left
    legend.justification = "left",    # Aligns entire legend box left
    legend.text.align = 0,            # Left-aligns text within legend
    legend.key = element_rect(fill = NA, color = NA) ) +
  guides(
    fill = guide_legend(order = 1),
    pattern = guide_legend(
      order = 2,
      title = NULL,
      override.aes = list(fill = NA)
    ) )


PLOT7

ggsave("PopulationPrevalence_Fig3.png",
       device = "png",
       height = 4, # in inches, but you can set the units
       width = 8,
       dpi = 300)


#active disease plot
PLOT8 <- ggplot(df2, aes(x = Sampling_Day, fill = ESD_Stage_photo2, pattern = Molt)) +
  geom_bar_pattern(position = "fill") +
  scale_fill_manual(
    values = c("black", viridis::magma(3, direction = 1, begin = 0.97, end = 0.4)),
    labels = c("Healthy" = "No ACTIVE disease") ) +
  scale_pattern_manual(
    values = c("pch", "none", "none"),
    labels = c("none" = " ", "Dead" = " ") ) +
  facet_grid(rows = vars(Treatment), cols = vars(Source_Pop)) +
  labs(
    y = "Proportion of lobster",
    x = "Day of experiment",
    fill = "ESD Index" ) +
  theme_bw(base_size = 12) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(),
    legend.position = "right",
    legend.box = "vertical",
    legend.box.just = "left",         # Aligns stacked legends left
    legend.justification = "left",    # Aligns entire legend box left
    legend.text.align = 0,            # Left-aligns text within legend
    legend.key = element_rect(fill = NA, color = NA) ) +
  guides(
    fill = guide_legend(order = 1),
    pattern = guide_legend(
      order = 2,
      title = NULL,
      override.aes = list(fill = NA)
    ) )


PLOT8
