---
title: "ESD survival & molt_plots"
author: "MMRocker"
date: "11 August 2025"
output: html_document
---


#set up working file & R

getwd() #check working directory
#setwd("/Users/jc229544/Desktop/")

library(ggplot2)
library(readxl)
library(lubridate)
require(grid)
library(gridExtra)
library(reshape2)
library(stringr)
library(tidyr)
library(viridis)
library(dplyr)
library(plyr)
library(survival)
library(coxme)
library(MuMIn)
library(cowplot)
library(patchwork)



st<-read_excel("ESDExpt_MoltEvents.xlsx",sheet="Survival")

st<-as.data.frame(st)

st <- subset(st, TAG!="3741") #remove early acclimation death from analysis
st2 <- subset(st, Premolt_Days>21) #remove lobsters that molted within first 3 weeks of expt (n=2) & death w/o molt (n=1) 

st$Treatment <- factor(st$Treatment, levels = c("CT cold", "MA cold", "ME cold", 
                                          "CT warm", "MA warm", "ME warm"))

st2$Source_Pop <- factor(st2$Source_Pop, levels = c("ME", "MA", "CT"))

st2$MoltStage1 <- factor(st2$MoltStage1, levels = c("C4", "D0"))


st2 %>%  
  ggplot(aes(y = factor(Source_Pop), x = Premolt_Days, fill = Temperature)) +
  geom_boxplot() +     # add colour to boxplots
  geom_point(alpha = 0.5) +                # alpha = transparency
  facet_wrap(~ MoltStage1, ncol = 1) +       # spread by Moltstage
  scale_fill_manual(values=c("blue", "red")) +
  theme(legend.position = "none") +         # remove legend
  xlab("Days pre-molt") +                            # label x-axis
  ylab("Source Pop") +         # label y-axis
  ggtitle("") # add title

###viewing count table
ftable(st2$Temperature, st2$Source_Pop, st2$MoltStage1, st2$Molt)


####too many variables, model DID NOT converge: AICc - 727.9556
coxph(Surv(Premolt_Days,Molt) ~ MoltStage1 * Temperature * Source_Pop * maxcoverage,
      data=st2)
AICc(coxph(Surv(Premolt_Days,Molt) ~ MoltStage1 * Temperature * Source_Pop * maxcoverage,
      data=st2))

#dredging full model
BIGGESTBOY <- coxph(Surv(Premolt_Days,Molt) ~ Source_Pop * Temperature * MoltStage1 + maxcoverage,
                    data=st2, na.action = "na.fail")
dredge(BIGGESTBOY)

#dredging with limited interactions
BIGBOY <- coxph(Surv(Premolt_Days,Molt) ~ Source_Pop + Temperature + MoltStage1 + maxcoverage
                + Source_Pop:Temperature + Source_Pop:MoltStage1 + MoltStage1:Temperature,
                data=st2, na.action = "na.fail")

dredge(BIGBOY)

####SAME MODEL RESULTS

#model factors
#****BEST dredge1 model: AIC - 166.5713
#model factors: pop, temp, molt stage, coverage, molt x pop
AICc(coxph(Surv(Premolt_Days,Molt) ~  Source_Pop + Temperature + MoltStage1 + maxcoverage
           + MoltStage1:Source_Pop,
           data=st2))
coxph(Surv(Premolt_Days,Molt) ~ Source_Pop + Temperature + MoltStage1 + maxcoverage
      + MoltStage1:Source_Pop ,
      data=st2)
#significant model factors: max cov, molt stage, temp

#****BEST dredge2 model: AIC - 166.7346
#model factors: pop, temp, molt stage, coverage
AICc(coxph(Surv(Premolt_Days,Molt) ~  Source_Pop + Temperature + MoltStage1 + maxcoverage,
           data=st2))
coxph(Surv(Premolt_Days,Molt) ~ Source_Pop + Temperature + MoltStage1 + maxcoverage,
      data=st2)
#significant model factors: max cov, molt stage, temp, popCT, popMA


#****BEST dredge3 model: AIC - 168.5148; pop x molt interaction; max cov, moltD0, temp sig
#model factors: pop, temp, molt stage
AICc(coxph(Surv(Premolt_Days,Molt) ~  Source_Pop + Temperature + MoltStage1 ,
           data=st2))
coxph(Surv(Premolt_Days,Molt) ~ Source_Pop + Temperature + MoltStage1 ,
      data=st2)
#significant model factors: molt stage, temp, popCT, popMA



####BEST MODEL!!!! w/o interaction for visual/figures
BESTmodel <- coxph(Surv(Premolt_Days,Molt) ~ Source_Pop + Temperature + MoltStage1 + maxcoverage,
                   data=st2)

summary(BESTmodel)

library(survminer)
ggcoxdiagnostics(BESTmodel, type = "martingale") #for non-linearity and outliers
ggcoxdiagnostics(BESTmodel, type = "schoenfeld") #test proportional hazard assumption
ggcoxdiagnostics(BESTmodel, type = "deviance") #assess goodness of fit & ID poorly fit observations
ggcoxdiagnostics(BESTmodel, type = "score") #evaluate influence of observations on parameter estimates

ftest <- cox.zph(BESTmodel)
ggcoxzph(ftest)

#strata() allows different MoltStage levels to have different survival functions/non proportionality with respect to molt stage 
nonPHmodel <- coxph(Surv(Premolt_Days,Molt) ~ Source_Pop + Temperature + strata(MoltStage1) + maxcoverage,
                  data=st2)
summary(nonPHmodel) 
ggcoxdiagnostics(nonPHmodel, type = "martingale") #for non-linearity and outliers
ggcoxdiagnostics(nonPHmodel, type = "schoenfeld") #test proportional hazard assumption
ggcoxdiagnostics(nonPHmodel, type = "deviance") #assess goodness of fit & ID poorly fit observations
ggcoxdiagnostics(nonPHmodel, type = "score") #evaluate influence of observations on parameter estimates

ftest <- cox.zph(nonPHmodel)
ggcoxzph(ftest)

noMOLTmodel <- coxph(Surv(Premolt_Days,Molt) ~ Source_Pop + Temperature +  maxcoverage,
                    data=st2)

AICc(nonPHmodel) #better AIC!!!
AICc(BESTmodel)
AICc(noMOLTmodel)

anova(noMOLTmodel, BESTmodel, nonPHmodel, test = "LRT") #adding strata(molt) improve model

summary(nonPHmodel) 
bh <- basehaz(nonPHmodel, centered = FALSE)
ggplot(bh, aes(x = time, y = hazard, color = strata)) +
  geom_step() +
  labs(y = "Cumulative baseline hazard")


####predicting model to plot curves based on set parameters
###predict: 1. Pop; 2. Temp; 3. Molt stage; 4. Disease coverage (0, 10, 20, 30, 40, 50%)


#create prediction dataframe
pred_dat50_10 <- expand.grid(Molt = 1,
                             Premolt_Days = 1:140,
                             Source_Pop = c("ME","CT","MA"),
                             Temperature = c("cold","warm"),
                             MoltStage1 = c("C4", "D0"),
                             maxcoverage = c(0,.10,.20,.30,.40,.50))



preds50_10 <- predict(nonPHmodel, newdata = pred_dat50_10,
                      type = "survival", se.fit = TRUE)
pred_dat50_10$prob <- preds50_10$fit
pred_dat50_10$se.fit <- preds50_10$se.fit
pred_dat50_10$lcl <- preds50_10$fit - 1.96*preds50_10$se.fit
pred_dat50_10$ucl <- preds50_10$fit + 1.96*preds50_10$se.fit
pred_dat50_10$prob_1 <- 1 - pred_dat50_10$prob
pred_dat50_10$lcl_1 <- 1 - pred_dat50_10$lcl
pred_dat50_10$ucl_1 <- 1 - pred_dat50_10$ucl
pred_dat50_10$prob_1_100 <- 100 * pred_dat50_10$prob_1
pred_dat50_10$maxcoverage100 <- 100* pred_dat50_10$maxcoverage

summary(pred_dat50_10)

#visual predictions
fitnew <- survfit(nonPHmodel, newdata = pred_dat50_10)
plot(fitnew, conf.int = TRUE)


#rename & reorder states
pred_dat50_10$Source_Pop <- revalue(pred_dat50_10$Source_Pop, c("CT" = "Connecticut",
                                                          "MA" = "Massachusetts",
                                                          "ME" = "Maine"))
pred_dat50_10$Source_Pop <- factor(pred_dat50_10$Source_Pop, levels = c("Maine", "Massachusetts", "Connecticut" ))

pred_dat50_10$Temperature <- as.factor(pred_dat50_10$Temperature)
pred_dat50_10$Temperature <- factor(pred_dat50_10$Temperature, levels = c("cold",  "warm"))
pred_dat50_10$Temperature <- revalue(pred_dat50_10$Temperature, c("cold" = "cold (13°C)",
                                                                  "warm" = "warm (20°C)"))

st2$Temperature <- revalue(st2$Temperature, c("cold" = "cold (13°C)",
                                                                  "warm" = "warm (20°C)"))
st2$Source_Pop <- revalue(st2$Source_Pop, c("CT" = "Connecticut",
                                                                "MA" = "Massachusetts",
                                                                "ME" = "Maine"))
st2$Source_Pop <- factor(st2$Source_Pop, levels = c("Maine", "Massachusetts", "Connecticut" ))


pred_dat50_10$MoltStage1 <- as.factor(pred_dat50_10$MoltStage1)
pred_dat50_10$MoltStage1 <- factor(pred_dat50_10$MoltStage1, levels = c("C4",  "D0"))




#results visualisation
###3graphs
MoltFig <- ggplot(pred_dat50_10 
       |> filter(Source_Pop == c("Maine"),
                 Temperature == "warm (20°C)", 
                 # MoltStage1 == "C4",
                 maxcoverage100 == 30), aes(x = Premolt_Days, y = prob_1)) + 
  geom_vline(xintercept=(c(1, 29, 57, 87, 115, 140,141)),
             color = "gray", linewidth=0.5)+
  geom_ribbon(aes(ymin = lcl_1, ymax =  ucl_1, fill = MoltStage1), alpha = 0.2) +
  geom_line(aes(col = MoltStage1), linewidth=1,) +
  scale_color_viridis_d(option="cividis",begin = .20, end = 0.95, name = NULL) +
  scale_fill_viridis_d(option="cividis",begin = .20, end = 0.95, name = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 10) +
  labs(
    #title='Molt stage',
    x='Day of experiment',
    y='Probabality of molting',
  ) +
  theme(legend.position = c(0.01, 0.99),
        legend.justification = c(0,1),
        legend.text = element_text(size = 8),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_line())

MoltFig

###2
TempFig <- ggplot(pred_dat50_10 
                  |> filter(Source_Pop == "Maine",
                            #Temperature == "warm (20°C)", 
                            MoltStage1 == "C4",
                            maxcoverage100 == 30), aes(x = Premolt_Days, y = prob_1)) + 
  geom_vline(xintercept=(c(1, 29, 57, 87, 115, 140,141)),
             color = "gray", linewidth=0.5)+
  geom_ribbon(aes(ymin = lcl_1, ymax =  ucl_1, fill = Temperature), alpha = 0.2) +
  geom_line(aes(col = Temperature), linewidth=1) +
  scale_color_manual(values=c("blue", "red"), name = NULL) +
  scale_fill_manual(values=c("blue", "red"), name = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 10) +
  labs(
    #title='Temperature',
    x=NULL,
    y=NULL
    ) +
  theme(legend.position = c(0.01, 0.99),
        legend.justification = c(0,1),
        legend.text = element_text(size = 8),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_line())

TempFig


# Define colors
star_colors <- viridis::magma(3, direction = 1, begin = 0.4, end = 0.97)
names(star_colors) <- c("Maine", "Massachusetts", "Connecticut")

###1
SourceFig <- ggplot(pred_dat50_10 
                  |> filter(#Source_Pop == "Maine",
                            Temperature == "warm (20°C)", 
                            MoltStage1 == "C4",
                            maxcoverage100 == 30), aes(x = Premolt_Days, y = prob_1)) + 
  geom_vline(xintercept=(c(1, 29, 57, 87, 115, 140,141)),
             color = "gray", linewidth=0.5)+
  geom_ribbon(aes(ymin = lcl_1, ymax =  ucl_1, fill = Source_Pop), alpha = 0.2) +
  geom_line(aes(col = Source_Pop), linewidth=1,) +
  scale_color_manual(values = star_colors, labels = c("Casco Bay, ME", "Buzzards Bay, MA", "Niantic Bay, CT"), name = NULL) +
  scale_fill_manual(values = star_colors, labels = c("Casco Bay, ME", "Buzzards Bay, MA", "Niantic Bay, CT"), name = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 10) +
  labs(
    #title='Source population',
    x=NULL,
    y='Probabality of molting') +
  theme(legend.position = c(0.01, 0.99),
        legend.justification = c(0,1),
        legend.text = element_text(size = 8),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_line())

SourceFig

###4
DiseaseFig <- ggplot(pred_dat50_10 
                  |> filter(Source_Pop == "Maine",
                            Temperature == "warm (20°C)", 
                            MoltStage1 == "C4",
                            #maxcoverage100 == 30
                            ), aes(x = Premolt_Days, y = prob_1)) + 
  geom_vline(xintercept=(c(1, 29, 57, 87, 115, 140,141)),
             color = "gray", linewidth=0.5)+
  geom_ribbon(aes(ymin = lcl_1, ymax =  ucl_1, fill = as.factor(maxcoverage100)), alpha = 0.2) +
  geom_line(aes(col = as.factor(maxcoverage100)), linewidth=1,) +
  scale_color_viridis_d(option="magma",begin = .97, end = 0.6, name = NULL) +
  scale_fill_viridis_d(option="magma",begin = .97, end = 0.6, name = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 10) +
  labs(
    #title='Maximum disease severity',
    x='Day of experiment',
    y=NULL
    ) +
  guides(color=guide_legend(reverse = TRUE),
         fill=guide_legend(reverse = TRUE)) +
  theme(legend.position = c(0.01, 0.99),
        legend.justification = c(0,1),
        legend.text = element_text(size = 8),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line())

DiseaseFig


# Create a separate "title" plots with a gray background
tMolt <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Molt stage", size = 3.5) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "gray90", color = NA),
    plot.margin = margin(2, 0, 2, 0)
  )

tTemp <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Temperature", size = 3.5) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "gray90", color = NA),
    plot.margin = margin(2, 0, 2, 0)
  )

tSource <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Source population", size = 3.5) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "gray90", color = NA),
    plot.margin = margin(2, 0, 2, 0)
  )

tDisease <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Maximum disease severity", size = 3.5) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "gray90", color = NA),
    plot.margin = margin(2, 0, 2, 0)
  )



# Combine using patchwork
pMolt <- tMolt / MoltFig +
  plot_layout(heights = c(0.08, 1))  # Adjust height ratio of title to plot
pMolt

pTemp <- tTemp / TempFig +
  plot_layout(heights = c(0.08, 1))  # Adjust height ratio of title to plot
pTemp

pSource <- tSource / SourceFig +
  plot_layout(heights = c(0.08, 1))  # Adjust height ratio of title to plot
pSource

pDisease <- tDisease / DiseaseFig +
  plot_layout(heights = c(0.08, 1))  # Adjust height ratio of title to plot
pDisease

gridMOLT <- plot_grid(pSource, pTemp, pDisease, pMolt,
                      ncol = 2#,  # 2 x 2 layout
                      #labels = c("B", "C", "D", "E")
)

gridMOLT


ggsave("MoltProbability_Fig5.png",
       device = "png",
       height = 5, # in inches, but you can set the units
       width = 7,
       dpi = 300)


