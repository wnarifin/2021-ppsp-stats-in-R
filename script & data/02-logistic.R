# About ====
#' title: "Statistical Analysis Using RStudio"
#' subtitle: "Logistic Regression"
#' author: "Wan Nor Arifin"

# Library ====
# official CRAN
library(foreign)
library(tidyverse)
library(gtsummary)
library(ggplot2)
library(ggpubr)
library(GGally)
library(psych)
library(rsq)
library(broom)
library(ResourceSelection)
# custom function
# desc_cat()
source("https://raw.githubusercontent.com/wnarifin/medicalstats-in-R/master/functions/desc_cat_fun.R")

# Data ====
coronary = read.dta("coronary.dta")
str(coronary)

## descriptive ----
coronary %>% select(-id) %>%
  tbl_summary(by = cad,
              statistic = all_continuous() ~ "{mean} ({sd})", 
              digits = all_continuous() ~ 1)

# SLogR ====
# example: cad ~ gender

## explore ----
# descriptive
coronary %>% select(gender, cad) %>% desc_cat()
with(coronary, table(gender, cad)) %>% print() %>% 
  prop.table(margin = 2) * 100

## fit slogr ----
slg_cad = glm(cad ~ gender, data = coronary, family = binomial)
summary(slg_cad)
tidy(slg_cad, conf.int = TRUE, exponentiate = TRUE)  # OR

# MLogR ====
# cad ~ many factors
names(coronary)  # recall what we have

## explore ----
# descriptive
coronary %>% select(-id, -cad, -race, -gender) %>% describeBy(., coronary$cad)  # numerical
coronary %>% select(race, gender) %>% by(., coronary$cad, desc_cat)  # categorical
# or just use tbl_summary()
# plots
coronary %>% select(-id, -cad) %>% by(., coronary$cad, ggpairs)  # or using built-in graphics

## fit univariable / slogr ----
# again, you can do one-by-one univariable on your own
slg_cad0 = glm(cad ~ 1, data = coronary, family = binomial)
summary(slg_cad0)
add1(slg_cad0, scope = ~ sbp + dbp + chol + age + bmi + race + gender, test = "LRT")
# all sig. except bmi & race
# choose only vars p-value < 0.25 to proceed in MLogR

## fit multivariable / mlogr ----
# race not included, p-value > 0.25
# sbp not included, redundant
mlg_cad = glm(cad ~ dbp + chol + age + bmi + gender, 
              data = coronary, family = binomial)
summary(mlg_cad)

# stepwise
# both
mlg_cad_stepboth = step(mlg_cad, direction = "both")
summary(mlg_cad_stepboth)  # cad ~ dbp + gender
# forward
mlg_cad_stepforward = step(slg_cad0, 
                           scope = ~ dbp + chol + age + bmi + gender, 
                           direction = "forward")
summary(mlg_cad_stepforward)  # cad ~ dbp + gender
# backward
mlg_cad_stepback = step(mlg_cad, direction = "backward")
summary(mlg_cad_stepback)  # cad ~ dbp + gender

# mlg_cad1: cad ~ dbp + gender
mlg_cad_sel = glm(cad ~ dbp + gender, data = coronary, family = binomial)
summary(mlg_cad_sel)
tidy(mlg_cad_sel, conf.int = TRUE)
tbl_regression(mlg_cad_sel)
rsq(mlg_cad_sel)

# multicollinearity, MC
# by looking at the estimates and standard errors, SEs
# when SE > Estimate -- MC problem
# How large? Relatively large... not specific in book.
# Sometimes, the estimates are unusually large, i.e. large ORs
# illogical -- also indicates MC problem
summary(mlg_cad_sel)  # all SEs < Estimates/Coefficients

# interaction, *
summary(glm(cad ~ dbp*gender, data = coronary, family = binomial))
# insig. *

# model fit
# 1. Hosmer-Lemeshow test
hl_cad = hoslem.test(mlg_cad_sel$y, mlg_cad_sel$fitted.values)
hl_cad  # does not fit, slightly... ideally > 0.05
# usually because small number of variables in the model
cbind(hl_cad$observed, hl_cad$expected)
# 2. classification table
coronary$cad_prob = mlg_cad_sel$fitted.values  # save predicted probabilities
coronary$cad_pred = ifelse(coronary$cad_prob > 0.5, "cad", "no cad") %>% as.factor()
table(coronary$cad, coronary$cad_pred)
# correctly classified %
100 * sum(coronary$cad == coronary$cad_pred) / length(coronary$cad)  # = 80%
# 3. Receiver Operating Characteristic (ROC) curve
roc_cad = epiDisplay::lroc(mlg_cad_sel)
roc_cad$auc
# we assume the model fit based on 2/3 criteria

# rename the selected model to final model
mlg_cad_final = mlg_cad_sel

# Present the results and interpret ====
summary(mlg_cad_final)
tib_mlg = tidy(mlg_cad_final, conf.int = TRUE, exponentiate = TRUE)  # in Odds ratio
rsq(mlg_cad_final, adj = T)
# R-squared is usually reported for linear regression.
# But R-squared is also available for GLM, in our case logistic regression.
# This is usually known as pseudo-R-squared. In GLM, it is made possible by the
# work of Zhang (2016), the author of "rsq" package.

# use kable to come up with nice table
knitr::kable(tib_mlg, format = "simple")
tbl_regression(mlg_cad_final)
# export it to a csv file for use later
write.csv(tib_mlg, "mlg_final.csv")

# Prediction ====

# predict probability of getting CAD ----
coronary$cad_prob = predict(mlg_cad_final, type = "response")  # in probability
# converted from logit, by adding type = "response"
head(coronary$cad_prob)
# you can also use mlg_cad_final$fitted.values
# but as we will see below, we need predict() for new data

## one subject ----
# dbp = 110, gender = man
predict(mlg_cad_final, list(dbp = 110, gender = "man"), type = "response")
# probability > 0.5 = cad

## many subjects ----
# more data points
new_data = data.frame(dbp = c(100, 110, 120, 100, 110, 120),
                      gender = c("man", "man", "man", "woman", "woman", "woman"))
new_data
predict(mlg_cad_final, new_data, type = "response")
new_data$prob_cad = predict(mlg_cad_final, new_data, type = "response")
new_data
new_data$pred_cad = cut(new_data$prob_cad, breaks = c(-Inf, 0.5, Inf),
                        labels = c("no cad", "cad"))
# alternative to ifelse
new_data

# Exercise ====
# - obtain OR for 5mmHg increase in DBP
