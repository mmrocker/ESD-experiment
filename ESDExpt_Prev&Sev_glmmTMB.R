---
title: "zero-inflated general linear mixed model for ESD"
author: "MMRocker"
date: "Aug 2025"
output: html_document
---


#set up working file & R

getwd() #check working directory



#current citation for glmmTMB is:
#  Brooks ME, Kristensen K, van Benthem KJ, Magnusson A, Berg
#CW, Nielsen A, Skaug HJ, Maechler M, Bolker BM(2017). “glmmTMB
#Balances Speed and Flexibility Among Packages for Zero-inflated
#Generalized Linear Mixed Modeling.” The R Journal, 9(2), 378–
#400. doi:10.32614/RJ-2017-066

#glmmTMBpackage_BenBolker.pdf for basics

#fit a zero-inflated Beta response by specifying ziformula
#https://stackoverflow.com/questions/73347804/glmm-with-beta-distribution-and-lots-of-zeros-in-y-variable


#Load required packages:
library(openxlsx)
library(glmmTMB)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(plyr)
library(emmeans)
library(patchwork)
library(MuMIn)
library(patchwork)
library(cowplot)
library(DHARMa)




#response variable: ESD % coverage 
#factors: Treatment (temp), Source population, TAG/Container (individuals), Table
#random effect: individuals nested in tables (measured repeatedly)
#distribution: beta
#model type: zero-inflated generalized linear mixed model (ZIGLMM)

#data set & manipulation
#upload data
df<-read.xlsx("ESDExpt_PrevalenceAndSeverity.xlsx", sheet="COVERAGE_long")


#define factors and numbers
df <- as.data.frame(df)
df$TAG<-as.factor(df$TAG)
df$Table=as.factor(df$Table)
df$Table=factor(df$Table, levels = c("A", "B", "C", "D", "E", "F"))
df$Treatment=as.factor(df$Treatment)
df$Source_Pop<-as.factor(df$Source_Pop)
df$Source_Pop=factor(df$Source_Pop, levels = c("ME","MA","CT"))
df$Container<-as.factor(df$Container)
df$MoltStage<-as.factor(df$MoltStage)

df$disease <- ifelse(df$percentCoverage == 0, 0, 1) #create presence/absence of disease

df2 <- subset(df, !MoltStage=="D2") #remove lobsters that molted within first 3 weeks of expt (n=2) 
df2$MoltStage<-factor(df2$MoltStage, levels = c("C4","D0"))

summary(df2)
dim(df2)

df_NO_NAs <- df2[!rowSums(is.na(df2[11])), ] #remove entries with NAs (D08 lobster; 2x post-death timepoints)
dim(df_NO_NAs)
summary(df_NO_NAs)


#basic glmmTMB fit — a zero-inflated with a single zero-inflation 
  #parameter applying to all observations (ziformula~1). 
  #(Excluding zero-inflation is glmmTMB’s default: to exclude it explicitly, use ziformula ~0.)

#percent coverage ~ temperature treatment * source population * DAY with individual nested in the tables as a random effect

#create full model
BIGGEST <- glmmTMB(
  percentCoverage ~ Treatment * Source_Pop * MoltStage * poly(Day,2) + (1 | Table/TAG),
  data = df_NO_NAs,
  ziformula = ~1,         # intercept-only zero-inflation
  family = beta_family(),
  na.action = "na.fail"
)

#dredge and get top model
dd <- dredge(BIGGEST)


#model without any 3- & 4-way interactions for over parameterization & model convergence issues
BIG <- glmmTMB(
  percentCoverage ~ Treatment + Source_Pop + MoltStage + poly(Day,2) + 
  Treatment:Source_Pop + Treatment:MoltStage + Treatment:poly(Day,2) +
  Source_Pop:MoltStage + Source_Pop:poly(Day,2) +
  MoltStage:poly(Day,2) + (1 | Table/TAG),
  data = df_NO_NAs,
  ziformula = ~1,         # intercept-only zero-inflation
  family = beta_family(),
  na.action = "na.fail"
)

#dredge and get top model
dd2 <- dredge(BIG)


#top model from dredge (simplest model within delta 2 AICc)
BESTmodel <- glmmTMB(
  percentCoverage ~ Source_Pop + Treatment + MoltStage + poly(Day,2) + 
    Source_Pop:poly(Day,2) + Treatment:poly(Day,2) + MoltStage:poly(Day,2) + 
    (1 | Table/TAG),
  data = df_NO_NAs, 
  ziformula = ~1,         # intercept-only zero-inflation
  family = beta_family()
)

summary(BESTmodel)
BESTmodel

#allowing for zero-inflation
BESTmodel2 <- glmmTMB(
  percentCoverage ~ Source_Pop + Treatment + MoltStage + poly(Day,2) + 
    Source_Pop:poly(Day,2) + Treatment:poly(Day,2) + MoltStage:poly(Day,2) + 
    (1 | Table/TAG),
  data = df_NO_NAs, 
  ziformula = ~.,         # zero-inflation
  family = beta_family()
)

summary(BESTmodel2)
BESTmodel2


AICc(BESTmodel) 
AICc(BESTmodel2)

#testing assumptions & model fit
sim_res <- simulateResiduals(BESTmodel2, n = 1000)
plot(sim_res)
plotResiduals(sim_res, df_NO_NAs$Day)
plotResiduals(sim_res, df_NO_NAs$Source_Pop)
plotResiduals(sim_res, df_NO_NAs$Treatment)
plotResiduals(sim_res, df_NO_NAs$MoltStage)
testDispersion(sim_res)
testZeroInflation(sim_res)


#Generate predictions
df_NO_NAs$predicted <- predict(BESTmodel2, type = "response")  # Predicted values on the response scale
df_NO_NAs$predictedC <- predict(BESTmodel2, type = "conditional")  # Predicted values (Severity (%)) on the response scale
df_NO_NAs$predictedZI <- predict(BESTmodel2, type = "zprob")  # Predicted values (Prevalence (probability)) on the response scale

#plot predictions vs observed - full model
ggplot(df_NO_NAs, aes(x = percentCoverage, y = predicted)) +
  geom_point(color = "blue", alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Observed Values", y = "Predicted Values", 
       title = "Observed vs. Predicted Values") +
  theme_minimal()

#plot predictions vs observed - conditional model - severity only
ggplot(df_NO_NAs, aes(x = percentCoverage, y = predictedC)) +
  geom_point(color = "blue", alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Observed Values", y = "Predicted Values - conditional model", 
       title = "Observed vs. Predicted Values - conditional model") +
  theme_minimal()


#plot predictions vs observed - zero-inflated model - prevalence only
ggplot(df_NO_NAs, aes(x = as.factor(disease), y = (1-predictedZI))) +
  geom_boxplot() +
  geom_point(color = "blue", alpha = 0.7) +
  #geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Observed Values", y = "Predicted Values - zero-inflated model model", 
       title = "Observed vs. Predicted Values - zero-inflated model model") +
  theme_minimal()



ColdHot <- c("#0000FF", "#FF0000") #blue,  red
LightColdHot <- c("#CCCCFF",  "#FFCCCC") #light blue, light red




#create one visual of model for each of the Source*Temp*Molt; this is the model fit for each
#MA cold C4, MA warm C4, CT cold C4, CT warm C4, ME cold C4, ME warm C4,
#CT cold D0, CT warm D0, ME cold D0, ME warm D0 (***no D0 in MA)
#overlay this on plot of points

#MA COLD C4
#create data frame for modeling
df_MAcold = data.frame(Day = 1:141,
                       Table = "A",
                       TAG = 3777,
                       Source_Pop = "MA",
                       MoltStage = "C4",
                       Treatment = "cold")
#View(df_MAcold)

#predict response of model and se and add to data frame
df_MAcold$fit = predict(BESTmodel2, 
                 newdata = df_MAcold,
                 re.form = ~0,
                 type = "response")

df_MAcold$se = predict(BESTmodel2, 
                        newdata = df_MAcold,
                        re.form = ~0,
                        type = "response",
                       se.fit = T)$se.fit

df_MAcold$fitC = predict(BESTmodel2, 
                        newdata = df_MAcold,
                        re.form = ~0,
                        type = "conditional")

df_MAcold$seC = predict(BESTmodel2, 
                       newdata = df_MAcold,
                       re.form = ~0,
                       type = "conditional",
                       se.fit = T)$se.fit

df_MAcold$fitZI = predict(BESTmodel2, 
                        newdata = df_MAcold,
                        re.form = ~0,
                        type = "zprob")

df_MAcold$seZI = predict(BESTmodel2, 
                       newdata = df_MAcold,
                       re.form = ~0,
                       type = "zprob",
                       se.fit = T)$se.fit
#Full model 
ggplot(data = df_MAcold,
       aes(x = Day,
           y = fit)) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = fit - se, ymax = fit + se),fill="lightblue" ) +
  geom_line(col="blue")

#Conditional component - severity
ggplot(data = df_MAcold,
       aes(x = Day,
           y = fitC)) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = fitC - seC, ymax = fitC + seC),fill="lightblue" ) +
  geom_line(col="blue")

#ZI component - prevalence
ggplot(data = df_MAcold,
       aes(x = Day,
           y = (1-fitZI))) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = (1-fitZI) - seZI, ymax = (1-fitZI) + seZI),fill="lightblue" ) +
  geom_line(col="blue")



#MA WARM C4
#create data frame for modeling
df_MAwarm = data.frame(Day = 1:141,
                       Table = "B",
                       TAG = 3765,
                       Source_Pop = "MA",
                       MoltStage = "C4",
                       Treatment = "warm")
#View(df_MAwarm)
#predict response of model and se and add to data frame
df_MAwarm$fit = predict(BESTmodel2, 
                        newdata = df_MAwarm,
                        re.form = ~0,
                        type = "response")

df_MAwarm$se = predict(BESTmodel2, 
                       newdata = df_MAwarm,
                       re.form = ~0,
                       type = "response",
                       se.fit = T)$se.fit

df_MAwarm$fitC = predict(BESTmodel2, 
                         newdata = df_MAwarm,
                         re.form = ~0,
                         type = "conditional")

df_MAwarm$seC = predict(BESTmodel2, 
                        newdata = df_MAwarm,
                        re.form = ~0,
                        type = "conditional",
                        se.fit = T)$se.fit

df_MAwarm$fitZI = predict(BESTmodel2, 
                          newdata = df_MAwarm,
                          re.form = ~0,
                          type = "zprob")

df_MAwarm$seZI = predict(BESTmodel2, 
                         newdata = df_MAwarm,
                         re.form = ~0,
                         type = "zprob",
                         se.fit = T)$se.fit

#Full model
ggplot(data = df_MAwarm,
       aes(x = Day,
           y = fit)) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = fit - se, ymax = fit + se),fill="pink" ) +
  geom_line(col="red")

#Conditional component - severity
ggplot(data = df_MAwarm,
       aes(x = Day,
           y = fitC)) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = fitC - seC, ymax = fitC + seC),fill="pink" ) +
  geom_line(col="red")

#ZI component - prevalence
ggplot(data = df_MAwarm,
       aes(x = Day,
           y = (1-fitZI))) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = (1-fitZI) - seZI, ymax = (1-fitZI) + seZI),fill="pink" ) +
  geom_line(col="red")


#ME COLD C4
#create data frame for modeling
df_MEcold = data.frame(Day = 1:141,
                       Table = "C",
                       TAG = 3738,
                       Source_Pop = "ME",
                       MoltStage = "C4",
                       Treatment = "cold")
#View(df_MEcold)
#predict response of model and se and add to data frame
df_MEcold$fit = predict(BESTmodel2, 
                        newdata = df_MEcold,
                        re.form = ~0,
                        type = "response")

df_MEcold$se = predict(BESTmodel2, 
                       newdata = df_MEcold,
                       re.form = ~0,
                       type = "response",
                       se.fit = T)$se.fit

df_MEcold$fitC = predict(BESTmodel2, 
                         newdata = df_MEcold,
                         re.form = ~0,
                         type = "conditional")

df_MEcold$seC = predict(BESTmodel2, 
                        newdata = df_MEcold,
                        re.form = ~0,
                        type = "conditional",
                        se.fit = T)$se.fit

df_MEcold$fitZI = predict(BESTmodel2, 
                          newdata = df_MEcold,
                          re.form = ~0,
                          type = "zprob")

df_MEcold$seZI = predict(BESTmodel2, 
                         newdata = df_MEcold,
                         re.form = ~0,
                         type = "zprob",
                         se.fit = T)$se.fit
#Full model 
ggplot(data = df_MEcold,
       aes(x = Day,
           y = fit)) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = fit - se, ymax = fit + se),fill="lightblue" ) +
  geom_line(col="blue")

#Conditional component - severity
ggplot(data = df_MEcold,
       aes(x = Day,
           y = fitC)) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = fitC - seC, ymax = fitC + seC),fill="lightblue" ) +
  geom_line(col="blue")

#ZI component - prevalence
ggplot(data = df_MEcold,
       aes(x = Day,
           y = (1-fitZI))) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = (1-fitZI) - seZI, ymax = (1-fitZI) + seZI),fill="lightblue" ) +
  geom_line(col="blue")




#ME WARM C4
#create data frame for modeling
df_MEwarm = data.frame(Day = 1:141,
                       Table = "E",
                       TAG = 3763,
                       Source_Pop = "ME",
                       MoltStage = "C4",
                       Treatment = "warm")
#View(df_MEwarm)
#predict response of model and se and add to data frame
df_MEwarm$fit = predict(BESTmodel2, 
                        newdata = df_MEwarm,
                        re.form = ~0,
                        type = "response")

df_MEwarm$se = predict(BESTmodel2, 
                       newdata = df_MEwarm,
                       re.form = ~0,
                       type = "response",
                       se.fit = T)$se.fit

df_MEwarm$fitC = predict(BESTmodel2, 
                         newdata = df_MEwarm,
                         re.form = ~0,
                         type = "conditional")

df_MEwarm$seC = predict(BESTmodel2, 
                        newdata = df_MEwarm,
                        re.form = ~0,
                        type = "conditional",
                        se.fit = T)$se.fit

df_MEwarm$fitZI = predict(BESTmodel2, 
                          newdata = df_MEwarm,
                          re.form = ~0,
                          type = "zprob")

df_MEwarm$seZI = predict(BESTmodel2, 
                         newdata = df_MEwarm,
                         re.form = ~0,
                         type = "zprob",
                         se.fit = T)$se.fit

#Full model
ggplot(data = df_MEwarm,
       aes(x = Day,
           y = fit)) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = fit - se, ymax = fit + se),fill="pink" ) +
  geom_line(col="red")

#Conditional component - severity
ggplot(data = df_MEwarm,
       aes(x = Day,
           y = fitC)) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = fitC - seC, ymax = fitC + seC),fill="pink" ) +
  geom_line(col="red")

#ZI component - prevalence
ggplot(data = df_MEwarm,
       aes(x = Day,
           y = (1-fitZI))) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = (1-fitZI) - seZI, ymax = (1-fitZI) + seZI),fill="pink" ) +
  geom_line(col="red")



#CT COLD C4
#create data frame for modeling
df_CTcold = data.frame(Day = 1:141,
                       Table = "F",
                       TAG = 3782,
                       Source_Pop = "CT",
                       MoltStage = "C4",
                       Treatment = "cold")
#View(df_CTcold)
#predict response of model and se and add to data frame
df_CTcold$fit = predict(BESTmodel2, 
                        newdata = df_CTcold,
                        re.form = ~0,
                        type = "response")

df_CTcold$se = predict(BESTmodel2, 
                       newdata = df_CTcold,
                       re.form = ~0,
                       type = "response",
                       se.fit = T)$se.fit

df_CTcold$fitC = predict(BESTmodel2, 
                         newdata = df_CTcold,
                         re.form = ~0,
                         type = "conditional")

df_CTcold$seC = predict(BESTmodel2, 
                        newdata = df_CTcold,
                        re.form = ~0,
                        type = "conditional",
                        se.fit = T)$se.fit

df_CTcold$fitZI = predict(BESTmodel2, 
                          newdata = df_CTcold,
                          re.form = ~0,
                          type = "zprob")

df_CTcold$seZI = predict(BESTmodel2, 
                         newdata = df_CTcold,
                         re.form = ~0,
                         type = "zprob",
                         se.fit = T)$se.fit
#Full model 
ggplot(data = df_CTcold,
       aes(x = Day,
           y = fit)) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = fit - se, ymax = fit + se),fill="lightblue" ) +
  geom_line(col="blue")

#Conditional component - severity
ggplot(data = df_CTcold,
       aes(x = Day,
           y = fitC)) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = fitC - seC, ymax = fitC + seC),fill="lightblue" ) +
  geom_line(col="blue")

#ZI component - prevalence
ggplot(data = df_CTcold,
       aes(x = Day,
           y = (1-fitZI))) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = (1-fitZI) - seZI, ymax = (1-fitZI) + seZI),fill="lightblue" ) +
  geom_line(col="blue")



#CT WARM C4
#create data frame for modeling
df_CTwarm = data.frame(Day = 1:141,
                       Table = "D",
                       TAG = 3759,
                       Source_Pop = "CT",
                       MoltStage = "C4",
                       Treatment = "warm")
#View(df_CTwarm)
#predict response of model and se and add to data frame
df_CTwarm$fit = predict(BESTmodel2, 
                        newdata = df_CTwarm,
                        re.form = ~0,
                        type = "response")

df_CTwarm$se = predict(BESTmodel2, 
                       newdata = df_CTwarm,
                       re.form = ~0,
                       type = "response",
                       se.fit = T)$se.fit

df_CTwarm$fitC = predict(BESTmodel2, 
                         newdata = df_CTwarm,
                         re.form = ~0,
                         type = "conditional")

df_CTwarm$seC = predict(BESTmodel2, 
                        newdata = df_CTwarm,
                        re.form = ~0,
                        type = "conditional",
                        se.fit = T)$se.fit

df_CTwarm$fitZI = predict(BESTmodel2, 
                          newdata = df_CTwarm,
                          re.form = ~0,
                          type = "zprob")

df_CTwarm$seZI = predict(BESTmodel2, 
                         newdata = df_CTwarm,
                         re.form = ~0,
                         type = "zprob",
                         se.fit = T)$se.fit

#Full model
ggplot(data = df_CTwarm,
       aes(x = Day,
           y = fit)) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = fit - se, ymax = fit + se),fill="pink" ) +
  geom_line(col="red")

#Conditional component - severity
ggplot(data = df_CTwarm,
       aes(x = Day,
           y = fitC)) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = fitC - seC, ymax = fitC + seC),fill="pink" ) +
  geom_line(col="red")

#ZI component - prevalence
ggplot(data = df_CTwarm,
       aes(x = Day,
           y = (1-fitZI))) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = (1-fitZI) - seZI, ymax = (1-fitZI) + seZI),fill="pink" ) +
  geom_line(col="red")






#ME COLD D0
#create data frame for modeling
df_MEcoldD0 = data.frame(Day = 1:141,
                       Table = "C",
                       TAG = 3738,
                       Source_Pop = "ME",
                       MoltStage = "D0",
                       Treatment = "cold")
#View(df_MEcold)
#predict response of model and se and add to data frame
df_MEcoldD0$fit = predict(BESTmodel2, 
                        newdata = df_MEcoldD0,
                        re.form = ~0,
                        type = "response")

df_MEcoldD0$se = predict(BESTmodel2, 
                       newdata = df_MEcoldD0,
                       re.form = ~0,
                       type = "response",
                       se.fit = T)$se.fit

df_MEcoldD0$fitC = predict(BESTmodel2, 
                         newdata = df_MEcoldD0,
                         re.form = ~0,
                         type = "conditional")

df_MEcoldD0$seC = predict(BESTmodel2, 
                        newdata = df_MEcoldD0,
                        re.form = ~0,
                        type = "conditional",
                        se.fit = T)$se.fit

df_MEcoldD0$fitZI = predict(BESTmodel2, 
                          newdata = df_MEcoldD0,
                          re.form = ~0,
                          type = "zprob")

df_MEcoldD0$seZI = predict(BESTmodel2, 
                         newdata = df_MEcoldD0,
                         re.form = ~0,
                         type = "zprob",
                         se.fit = T)$se.fit
#Full model 
ggplot(data = df_MEcoldD0,
       aes(x = Day,
           y = fit)) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = fit - se, ymax = fit + se),fill="lightblue" ) +
  geom_line(col="blue")

#Conditional component - severity
ggplot(data = df_MEcoldD0,
       aes(x = Day,
           y = fitC)) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = fitC - seC, ymax = fitC + seC),fill="lightblue" ) +
  geom_line(col="blue")

#ZI component - prevalence
ggplot(data = df_MEcoldD0,
       aes(x = Day,
           y = (1-fitZI))) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = (1-fitZI) - seZI, ymax = (1-fitZI) + seZI),fill="lightblue" ) +
  geom_line(col="blue")




#ME WARM D0
#create data frame for modeling
df_MEwarmD0 = data.frame(Day = 1:141,
                       Table = "E",
                       TAG = 3763,
                       Source_Pop = "ME",
                       MoltStage = "D0",
                       Treatment = "warm")
#View(df_MEwarm)
#predict response of model and se and add to data frame
df_MEwarmD0$fit = predict(BESTmodel2, 
                        newdata = df_MEwarmD0,
                        re.form = ~0,
                        type = "response")

df_MEwarmD0$se = predict(BESTmodel2, 
                       newdata = df_MEwarmD0,
                       re.form = ~0,
                       type = "response",
                       se.fit = T)$se.fit

df_MEwarmD0$fitC = predict(BESTmodel2, 
                         newdata = df_MEwarmD0,
                         re.form = ~0,
                         type = "conditional")

df_MEwarmD0$seC = predict(BESTmodel2, 
                        newdata = df_MEwarmD0,
                        re.form = ~0,
                        type = "conditional",
                        se.fit = T)$se.fit

df_MEwarmD0$fitZI = predict(BESTmodel2, 
                          newdata = df_MEwarmD0,
                          re.form = ~0,
                          type = "zprob")

df_MEwarmD0$seZI = predict(BESTmodel2, 
                         newdata = df_MEwarmD0,
                         re.form = ~0,
                         type = "zprob",
                         se.fit = T)$se.fit

#Full model
ggplot(data = df_MEwarmD0,
       aes(x = Day,
           y = fit)) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = fit - se, ymax = fit + se),fill="pink" ) +
  geom_line(col="red")

#Conditional component - severity
ggplot(data = df_MEwarmD0,
       aes(x = Day,
           y = fitC)) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = fitC - seC, ymax = fitC + seC),fill="pink" ) +
  geom_line(col="red")

#ZI component - prevalence
ggplot(data = df_MEwarmD0,
       aes(x = Day,
           y = (1-fitZI))) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = (1-fitZI) - seZI, ymax = (1-fitZI) + seZI),fill="pink" ) +
  geom_line(col="red")



#CT COLD D0
#create data frame for modeling
df_CTcoldD0 = data.frame(Day = 1:141,
                       Table = "F",
                       TAG = 3782,
                       Source_Pop = "CT",
                       MoltStage = "D0",
                       Treatment = "cold")
#View(df_CTcold)
#predict response of model and se and add to data frame
df_CTcoldD0$fit = predict(BESTmodel2, 
                        newdata = df_CTcoldD0,
                        re.form = ~0,
                        type = "response")

df_CTcoldD0$se = predict(BESTmodel2, 
                       newdata = df_CTcoldD0,
                       re.form = ~0,
                       type = "response",
                       se.fit = T)$se.fit

df_CTcoldD0$fitC = predict(BESTmodel2, 
                         newdata = df_CTcoldD0,
                         re.form = ~0,
                         type = "conditional")

df_CTcoldD0$seC = predict(BESTmodel2, 
                        newdata = df_CTcoldD0,
                        re.form = ~0,
                        type = "conditional",
                        se.fit = T)$se.fit

df_CTcoldD0$fitZI = predict(BESTmodel2, 
                          newdata = df_CTcoldD0,
                          re.form = ~0,
                          type = "zprob")

df_CTcoldD0$seZI = predict(BESTmodel2, 
                         newdata = df_CTcoldD0,
                         re.form = ~0,
                         type = "zprob",
                         se.fit = T)$se.fit
#Full model 
ggplot(data = df_CTcoldD0,
       aes(x = Day,
           y = fit)) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = fit - se, ymax = fit + se),fill="lightblue" ) +
  geom_line(col="blue")

#Conditional component - severity
ggplot(data = df_CTcoldD0,
       aes(x = Day,
           y = fitC)) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = fitC - seC, ymax = fitC + seC),fill="lightblue" ) +
  geom_line(col="blue")

#ZI component - prevalence
ggplot(data = df_CTcoldD0,
       aes(x = Day,
           y = (1-fitZI))) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = (1-fitZI) - seZI, ymax = (1-fitZI) + seZI),fill="lightblue" ) +
  geom_line(col="blue")



#CT WARM D0
#create data frame for modeling
df_CTwarmD0 = data.frame(Day = 1:141,
                       Table = "D",
                       TAG = 3759,
                       Source_Pop = "CT",
                       MoltStage = "D0",
                       Treatment = "warm")
#View(df_CTwarm)
#predict response of model and se and add to data frame
df_CTwarmD0$fit = predict(BESTmodel2, 
                        newdata = df_CTwarmD0,
                        re.form = ~0,
                        type = "response")

df_CTwarmD0$se = predict(BESTmodel2, 
                       newdata = df_CTwarmD0,
                       re.form = ~0,
                       type = "response",
                       se.fit = T)$se.fit

df_CTwarmD0$fitC = predict(BESTmodel2, 
                         newdata = df_CTwarmD0,
                         re.form = ~0,
                         type = "conditional")

df_CTwarmD0$seC = predict(BESTmodel2, 
                        newdata = df_CTwarmD0,
                        re.form = ~0,
                        type = "conditional",
                        se.fit = T)$se.fit

df_CTwarmD0$fitZI = predict(BESTmodel2, 
                          newdata = df_CTwarmD0,
                          re.form = ~0,
                          type = "zprob")

df_CTwarmD0$seZI = predict(BESTmodel2, 
                         newdata = df_CTwarmD0,
                         re.form = ~0,
                         type = "zprob",
                         se.fit = T)$se.fit

#Full model
ggplot(data = df_CTwarmD0,
       aes(x = Day,
           y = fit)) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = fit - se, ymax = fit + se),fill="pink" ) +
  geom_line(col="red")

#Conditional component - severity
ggplot(data = df_CTwarmD0,
       aes(x = Day,
           y = fitC)) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = fitC - seC, ymax = fitC + seC),fill="pink" ) +
  geom_line(col="red")

#ZI component - prevalence
ggplot(data = df_CTwarmD0,
       aes(x = Day,
           y = (1-fitZI))) +
  facet_wrap(~ MoltStage) +
  geom_ribbon(aes(ymin = (1-fitZI) - seZI, ymax = (1-fitZI) + seZI),fill="pink" ) +
  geom_line(col="red")







#merge fitted data into one dataset
df_MAcold2 <- subset(df_MAcold[,c(1,4:12)])
df_MAwarm2 <- subset(df_MAwarm[,c(1,4:12)])
df_CTcold2 <- subset(df_CTcold[,c(1,4:12)])
df_CTwarm2 <- subset(df_CTwarm[,c(1,4:12)])
df_MEcold2 <- subset(df_MEcold[,c(1,4:12)])
df_MEwarm2 <- subset(df_MEwarm[,c(1,4:12)])

df_CTcoldD02 <- subset(df_CTcoldD0[,c(1,4:12)])
df_CTwarmD02 <- subset(df_CTwarmD0[,c(1,4:12)])
df_MEcoldD02 <- subset(df_MEcoldD0[,c(1,4:12)])
df_MEwarmD02 <- subset(df_MEwarmD0[,c(1,4:12)])

df_merge = merge(df_MAcold2, df_MAwarm2, by=c("Day", "Source_Pop", "MoltStage", "Treatment", "fit", "se",
                                               "fitC", "seC",
                                               "fitZI", "seZI"),all=T)
df_merge = merge(df_merge, df_CTcold2, by=c("Day", "Source_Pop", "MoltStage", "Treatment", "fit", "se",
                                            "fitC", "seC",
                                            "fitZI", "seZI"),all=T)
df_merge = merge(df_merge, df_CTwarm2, by=c("Day", "Source_Pop", "MoltStage", "Treatment", "fit", "se",
                                            "fitC", "seC",
                                            "fitZI", "seZI"),all=T)
df_merge = merge(df_merge, df_MEcold2, by=c("Day", "Source_Pop", "MoltStage", "Treatment", "fit", "se",
                                            "fitC", "seC",
                                            "fitZI", "seZI"),all=T)
df_merge = merge(df_merge, df_MEwarm2, by=c("Day", "Source_Pop", "MoltStage", "Treatment", "fit", "se",
                                            "fitC", "seC",
                                            "fitZI", "seZI"),all=T)
df_merge = merge(df_merge, df_CTcoldD02, by=c("Day", "Source_Pop", "MoltStage", "Treatment", "fit", "se",
                                            "fitC", "seC",
                                            "fitZI", "seZI"),all=T)
df_merge = merge(df_merge, df_CTwarmD02, by=c("Day", "Source_Pop", "MoltStage", "Treatment", "fit", "se",
                                            "fitC", "seC",
                                            "fitZI", "seZI"),all=T)
df_merge = merge(df_merge, df_MEcoldD02, by=c("Day", "Source_Pop", "MoltStage", "Treatment", "fit", "se",
                                            "fitC", "seC",
                                            "fitZI", "seZI"),all=T)
df_merge = merge(df_merge, df_MEwarmD02, by=c("Day", "Source_Pop", "MoltStage", "Treatment", "fit", "se",
                                            "fitC", "seC",
                                            "fitZI", "seZI"),all=T)



df_merge$fitZI_1 <- (1-df_merge$fitZI) #reverse the probability for graphing

#rename & reorder states
names(df_merge)[names(df_merge) == "Treatment"] <- "Temperature"
df_merge$Source_Pop <- revalue(df_merge$Source_Pop, c("CT" = "Connecticut",
                                                      "MA" = "Massachusetts",
                                                      "ME" = "Maine"))
df_merge$Source_Pop <- factor(df_merge$Source_Pop, levels = c("Maine", "Massachusetts", "Connecticut"))
df_merge$MoltStage <- factor(df_merge$MoltStage, levels = c("C4", "D0"))


#rename & reorder states
names(df_NO_NAs)[names(df_NO_NAs) == "Treatment"] <- "Temperature"
df_NO_NAs$Source_Pop <- revalue(df_NO_NAs$Source_Pop, c("CT" = "Connecticut",
                                                        "MA" = "Massachusetts",
                                                        "ME" = "Maine"))
df_NO_NAs$Source_Pop <- factor(df_NO_NAs$Source_Pop, levels = c("Maine", "Massachusetts", "Connecticut"))

#rename & reorder states
df_merge$Source_Pop <- revalue(df_merge$Source_Pop, c("Connecticut" = "Niantic Bay, CT",
                                                      "Massachusetts" = "Buzzards Bay, MA",
                                                      "Maine" = "Casco Bay, ME"))
df_merge$Source_Pop <- factor(df_merge$Source_Pop, levels = c("Casco Bay, ME",
                                                              "Buzzards Bay, MA",
                                                              "Niantic Bay, CT"))

df_NO_NAs$Source_Pop <- revalue(df_NO_NAs$Source_Pop, c("Connecticut" = "Niantic Bay, CT",
                                                              "Massachusetts" = "Buzzards Bay, MA",
                                                              "Maine" = "Casco Bay, ME"))
df_NO_NAs$Source_Pop <- factor(df_NO_NAs$Source_Pop, levels = c("Casco Bay, ME",
                                                                      "Buzzards Bay, MA",
                                                                      "Niantic Bay, CT"))


# Define colors
star_colors <- viridis::magma(3, direction = 1, begin = 0.4, end = 0.97)
names(star_colors) <- c("Casco Bay, ME", "Buzzards Bay, MA", "Niantic Bay, CT")

###df_merge ==> model predictions












ZI3 <- ggplot(df_merge, aes(x = Day, y = fitZI_1)) + #predicted
  geom_vline(xintercept=(c(1, 29, 57, 87, 115, 140,141)),
             color = "gray", linewidth=0.5)+
  #geom_line(data=df_NO_NAs, aes(x=Day, y= disease, group = TAG), col="lightgray") + #lines connecting actual
  geom_point(data=df_NO_NAs, aes(x = Day, y = disease, group = Source_Pop, col = Source_Pop), shape=16, size = 1, 
             position = position_jitter(width = 2, height = 0.05)) + #actual prevalence
  geom_ribbon(aes(ymin = fitZI_1 - seZI, ymax = fitZI_1 + seZI, fill = Source_Pop), alpha = 0.5) + #predicted
  geom_line(aes(col = Source_Pop), linewidth=1) +
  facet_grid(rows = vars(MoltStage), cols = vars(Temperature),
             labeller = labeller(
               MoltStage = c(
                 "C4" = "C4 (intermolt)",
                 "D0" = "D0 (plateau phase)"
               ),
               Temperature = c(
                 "cold" = "cold (13°C)",
                 "warm" = "warm (20°C)"
               )
             )) +
  scale_color_manual(values = star_colors, labels = c("Casco Bay, ME", "Buzzards Bay, MA", "Niantic Bay, CT"), name = NULL) +
  scale_fill_manual(values = star_colors, labels = c("Casco Bay, ME", "Buzzards Bay, MA", "Niantic Bay, CT"), name = NULL) +
  scale_x_continuous( breaks = c(0, 50, 100, 140)) +
  theme_bw(base_size = 11)+
  labs(
    #title='conditional model fit + se with observed (color + line) values',
    x = NULL,
    y='ESD prevalence') + 
  theme(legend.position = "none",
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_line(),
        axis.text.x = element_blank())
ZI3

COND3 <-ggplot(df_merge, aes(x = Day, y = 100*(fitC))) + #predicted
  geom_vline(xintercept=(c(1, 29, 57, 87, 115, 140,141)),
             color = "gray", linewidth=0.5)+
  geom_line(data = df_NO_NAs, aes(x=Day, y= 100*(percentCoverage), group = TAG), col="lightgray") + #lines connecting actual
  geom_point(data = df_NO_NAs, aes(x=Day, y= 100*(percentCoverage), group = Source_Pop, col = Source_Pop), shape=16, size =1) + #actual
  geom_ribbon(aes(ymin = 100*(fitC - seC), ymax = 100*(fitC + seC), fill = Source_Pop), alpha = 0.5) + #predicted
  geom_line(aes(col = Source_Pop), linewidth=1) +
  facet_grid(rows = vars(MoltStage), cols = vars(Temperature),
             labeller = labeller(
               MoltStage = c(
                 "C4" = "C4 (intermolt)",
                 "D0" = "D0 (plateau phase)"
               ),
               Temperature = c(
                 "cold" = "cold (13°C)",
                 "warm" = "warm (20°C)"
               )
             )) +
  scale_color_manual(values = star_colors, labels = c("Casco Bay, ME", "Buzzards Bay, MA", "Niantic Bay, CT"), name = NULL) +
  scale_fill_manual(values = star_colors, labels = c("Casco Bay, ME", "Buzzards Bay, MA", "Niantic Bay, CT"), name = NULL) +
  scale_x_continuous(breaks = c(0, 50, 100, 140)) +
  labs(
    #title='conditional model fit + se with observed (color + line) values',
    x='Day of experiment',
    y='ESD severity (%)') +
  coord_cartesian(ylim = c(0, 50)) +
  theme_bw(base_size = 11) +
  theme(legend.position = c(0.15, 0.3),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_line())
COND3

( (ZI3) / (COND3)) 

ggsave("PrevSev_FigS2.png",
       device = "png",
       height = 8, # in inches, but you can set the units
       width = 8,
       dpi = 300)



####SEVERITY
###1
PopFig <- ggplot(df_merge 
                  |> filter(#Source_Pop == "Casco Bay, ME",
                            Temperature == "warm", 
                            MoltStage == "C4"
                  ), aes(x = Day, y = 100*(fitC))) + 
  geom_vline(xintercept=(c(1, 29, 57, 87, 115, 140,141)),
             color = "gray", linewidth=0.5)+
  geom_line(data = df_NO_NAs
            |> filter(#Source_Pop == "Casco Bay, ME",
              Temperature == "warm", 
              MoltStage == "C4"
            ), aes(x=Day, y= 100*(percentCoverage), group = TAG), col="lightgray", alpha = 0.5) + #lines connecting actual
  geom_point(data = df_NO_NAs
             |> filter(#Source_Pop == "Casco Bay, ME",
               Temperature == "warm", 
               MoltStage == "C4"
             ), aes(x = Day, y = 100*(percentCoverage), group = Source_Pop, col = Source_Pop), alpha = 0.7, size = .8) + #actual severity
  geom_ribbon(aes(ymin = 100*(fitC - seC), ymax = 100*(fitC + seC), fill = Source_Pop), alpha = 0.3) +
  geom_line(aes(col = Source_Pop), linewidth=1, alpha = 1.5) +
  scale_color_manual(values = star_colors, labels = c("Casco Bay, ME", "Buzzards Bay, MA", "Niantic Bay, CT"), name = NULL) +
  scale_fill_manual(values = star_colors, labels = c("Casco Bay, ME", "Buzzards Bay, MA", "Niantic Bay, CT"), name = NULL) +
  theme_bw(base_size = 10) +
  coord_cartesian(ylim = c(0, 40)) +
  labs(
    #title='Source population',
    x='Day of experiment',
    y='Disease severity (%)'
  ) +
  theme(legend.position = c(0.01, 0.99),
        legend.justification = c(0,1),
        legend.text = element_text(size = 8),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line())

PopFig


###2
TempFig <- ggplot(df_merge 
                  |> filter(Source_Pop == "Casco Bay, ME",
                            #Temperature == "warm", 
                            MoltStage == "C4"
                  ), aes(x = Day, y = 100*(fitC))) + 
  geom_vline(xintercept=(c(1, 29, 57, 87, 115, 140,141)),
             color = "gray", linewidth=0.5)+
  geom_line(data = df_NO_NAs
            |> filter(Source_Pop == "Casco Bay, ME",
              #Temperature == "warm", 
              MoltStage == "C4"
            ), aes(x=Day, y= 100*(percentCoverage), group = TAG), col="lightgray", alpha = 0.5) + #lines connecting actual
  geom_point(data = df_NO_NAs
             |> filter(Source_Pop == "Casco Bay, ME",
               #Temperature == "warm", 
               MoltStage == "C4"
             ), aes(x = Day, y = 100*(percentCoverage), group = Temperature, col = Temperature), alpha = 0.7, size = .8) + #actual severity
  geom_ribbon(aes(ymin = 100*(fitC - seC), ymax = 100*(fitC + seC), fill = Temperature), alpha = 0.2) +
  geom_line(aes(col = Temperature), linewidth=1,) +
  scale_color_manual(values=c("blue", "red"), labels = c("cold" = "cold (13°C)","warm" = "warm (20°C)"), name = NULL) +
  scale_fill_manual(values=c("blue", "red"), labels = c("cold" = "cold (13°C)","warm" = "warm (20°C)"), name = NULL) +
  theme_bw(base_size = 10) +
  coord_cartesian(ylim = c(0, 40)) +
  labs(
    #title='Temperature',
    x='Day of experiment',
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

###3
MoltFig <- ggplot(df_merge 
                     |> filter(Source_Pop == "Casco Bay, ME",
                               Temperature == "warm", 
                               #MoltStage1 == "C4"
                     ), aes(x = Day, y = 100*(fitC))) + 
  geom_vline(xintercept=(c(1, 29, 57, 87, 115, 140,141)),
             color = "gray", linewidth=0.5)+
  geom_line(data = df_NO_NAs
            |> filter(Source_Pop == "Casco Bay, ME",
              Temperature == "warm", 
              #MoltStage == "C4"
            ), aes(x=Day, y= 100*(percentCoverage), group = TAG), col="lightgray", alpha = 0.5) + #lines connecting actual
  geom_point(data = df_NO_NAs
             |> filter(Source_Pop == "Casco Bay, ME",
               Temperature == "warm", 
               #MoltStage == "C4"
             ), aes(x = Day, y = 100*(percentCoverage), group = MoltStage, col = MoltStage), alpha = 0.7, size = .8) + #actual severity
  geom_ribbon(aes(ymin = 100*(fitC - seC), ymax = 100*(fitC + seC), fill = MoltStage), alpha = 0.2) +
  geom_line(aes(col = MoltStage), linewidth=1,) +
  scale_color_viridis_d(option="cividis",begin = .20, end = 0.95, labels = c("C4" = "C4 (intermolt)", "D0" = "D0 (plateau phase)"), name = NULL) +
  scale_fill_viridis_d(option="cividis",begin = .20, end = 0.95, labels = c("C4" = "C4 (intermolt)", "D0" = "D0 (plateau phase)"), name = NULL) +
  theme_bw(base_size = 10) +
  coord_cartesian(ylim = c(0, 40)) +
  labs(
    #title='Molt Stage',
    x='Day of experiment',
    y=NULL #'Disease severity (%)'
  ) +
  theme(legend.position = c(0.01, 0.99),
        legend.justification = c(0,1),
        legend.text = element_text(size = 8),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line())

MoltFig

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


# Combine using patchwork
pMolt <- tMolt / MoltFig +
  plot_layout(heights = c(0.08, 1))  # Adjust height ratio of title to plot
pMolt

pTemp <- tTemp / TempFig +
  plot_layout(heights = c(0.08, 1))  # Adjust height ratio of title to plot
pTemp

pSource <- tSource / PopFig +
  plot_layout(heights = c(0.08, 1))  # Adjust height ratio of title to plot
pSource


gridSEVERITY <- plot_grid(pSource, pTemp, pMolt, 
                      ncol = 2#,  # 2 x 2 layout
                      #labels = c("B", "C", "D", "E")
)

gridSEVERITY



####PREVALENCE
###1
PopFig2 <- ggplot(df_merge 
                 |> filter(#Source_Pop == "Casco Bay, ME",
                   Temperature == "warm", 
                   MoltStage == "C4"
                 ), aes(x = Day, y = fitZI_1)) + 
  geom_vline(xintercept=(c(1, 29, 57, 87, 115, 140,141)),
             color = "gray", linewidth=0.5)+
  geom_point(data = df_NO_NAs
             |> filter(#Source_Pop == "Casco Bay, ME",
                       Temperature == "warm", 
                       MoltStage == "C4"
             ), aes(x = Day, y = disease, group = Source_Pop, col = Source_Pop), shape=16, size = 1, 
             position = position_jitter(width = 2, height = 0.05)) + #actual prevalence
  geom_ribbon(aes(ymin = fitZI_1 - seZI, ymax = fitZI_1 + seZI, fill = Source_Pop), alpha = 0.2) +
  geom_line(aes(col = Source_Pop), linewidth=1,) +
  scale_color_manual(values = star_colors, labels = c("Casco Bay, ME", "Buzzards Bay, MA", "Niantic Bay, CT"), name = NULL) +
  scale_fill_manual(values = star_colors, labels = c("Casco Bay, ME", "Buzzards Bay, MA", "Niantic Bay, CT"), name = NULL) +
  theme_bw(base_size = 10) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    #title='Source population',
    x= NULL,'Day of experiment',
    y='Disease  prevalence'
  ) +
  theme(legend.position = 'none', # c(0.01, 0.99),
        legend.justification = c(0,1),
        legend.text = element_text(size = 8),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line())

PopFig2


###2
TempFig2 <- ggplot(df_merge 
                  |> filter(Source_Pop == "Casco Bay, ME",
                            #Temperature == "warm", 
                            MoltStage == "C4"
                  ), aes(x = Day, y = fitZI_1)) + 
  geom_vline(xintercept=(c(1, 29, 57, 87, 115, 140,141)),
             color = "gray", linewidth=0.5)+
  geom_point(data = df_NO_NAs
             |> filter(Source_Pop == "Casco Bay, ME",
               #Temperature == "warm", 
               MoltStage == "C4"
             ), aes(x = Day, y = disease, group = Temperature, col = Temperature), shape=16, size = 1, 
             position = position_jitter(width = 2, height = 0.05)) + #actual prevalence
  geom_ribbon(aes(ymin = fitZI_1 - seZI, ymax = fitZI_1 + seZI, fill = Temperature), alpha = 0.2) +
  geom_line(aes(col = Temperature), linewidth=1,) +
  scale_color_manual(values=c("blue", "red"), name = NULL) +
  scale_fill_manual(values=c("blue", "red"), name = NULL) +
  theme_bw(base_size = 10) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    #title='Temperature',
    x= NULL,'Day of experiment',
    y=NULL
  ) +
  theme(legend.position = 'none', # c(0.01, 0.99),
        legend.justification = c(0,1),
        legend.text = element_text(size = 8),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line())

TempFig2

###3
MoltFig2 <- ggplot(df_merge 
                  |> filter(Source_Pop == "Casco Bay, ME",
                            Temperature == "warm", 
                            #MoltStage1 == "C4"
                  ), aes(x = Day, y = fitZI_1)) + 
  geom_vline(xintercept=(c(1, 29, 57, 87, 115, 140,141)),
             color = "gray", linewidth=0.5)+
  geom_point(data = df_NO_NAs
             |> filter(Source_Pop == "Casco Bay, ME",
               Temperature == "warm", 
               #MoltStage == "C4"
             ), aes(x = Day, y = disease, group = MoltStage, col = MoltStage), shape=16, size = 1, 
             position = position_jitter(width = 2, height = 0.05)) + #actual prevalence
  geom_ribbon(aes(ymin = fitZI_1 - seZI, ymax = fitZI_1 + seZI, fill = MoltStage), alpha = 0.2) +
  geom_line(aes(col = MoltStage), linewidth=1,) +
  scale_color_viridis_d(option="cividis",begin = .20, end = 0.95, name = NULL) +
  scale_fill_viridis_d(option="cividis",begin = .20, end = 0.95, name = NULL) +
  theme_bw(base_size = 10) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    #title='Molt Stage',
    x= NULL, #'Day of experiment',
    y= NULL #'Disease prevalence'
  ) +
  theme(legend.position = 'none', #c(0.01, 0.99),
        legend.justification = c(0,1),
        legend.text = element_text(size = 8),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line())

MoltFig2


# Combine using patchwork
pMolt2 <- tMolt / MoltFig2 +
  plot_layout(heights = c(0.08, 1))  # Adjust height ratio of title to plot
pMolt2

pTemp2 <- tTemp / TempFig2 +
  plot_layout(heights = c(0.08, 1))  # Adjust height ratio of title to plot
pTemp2

pSource2 <- tSource / PopFig2 +
  plot_layout(heights = c(0.08, 1))  # Adjust height ratio of title to plot
pSource2


gridPREVALENCE <- plot_grid(pSource2, pTemp2, pMolt2, 
                          ncol = 2#,  # 2 x 2 layout
                          #labels = c("B", "C", "D", "E")
)

gridPREVALENCE



gridPreSev <- plot_grid(pSource2, pTemp2, pMolt2,
                        pSource, pTemp, pMolt,
                        nrow = 2#,  # 2 x 3 layout
                        #labels = c("B", "C", "D", "E")
)

gridPreSev


ggsave("PrevSev_Fig4.png",
       device = "png",
       height = 6, # in inches, but you can set the units
       width = 9,
       dpi = 300)

