
# Clean up the R environment
rm(list=ls())

# Load libraries
library(lmerTest)
library(corrplot)

# Read in the data
data<-read.csv("Long_data.csv")
data$Diet<-as.factor(data$Diet)
head(data)
data$n_X.pup<-data$n_total - data$n_pup
data$n_X.eclo<-data$n_pup - data$n_eclo

# Make the HP diet the reference group
contrasts(data$Diet)<-contr.treatment(levels(data$Diet), base=2)

# Fit the models
model_G<-glmer(cbind(n_pup, n_X.pup) ~ Diet + (1|DGRP_number), family="binomial", data=data, control=glmerControl(optimizer = "bobyqa"))
summary(model_G)

model_GxE<-glmer(cbind(n_pup, n_X.pup) ~ Diet + (Diet|DGRP_number), family="binomial", data=data, control=glmerControl(optimizer = "bobyqa"))
summary(model_GxE)
anova(model_G, model_GxE)

# What about for ecolsion
model_G<-glmer(cbind(n_eclo, n_X.eclo) ~ Diet + (1|DGRP_number), family="binomial", data=data, control=glmerControl(optimizer = "bobyqa"))
summary(model_G)

model_GxE<-glmer(cbind(n_eclo, n_X.eclo) ~ Diet + (Diet|DGRP_number), family="binomial", data=data, control=glmerControl(optimizer = "bobyqa"))
summary(model_GxE)
anova(model_G, model_GxE)
