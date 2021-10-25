# About ====
#' title: "Statistical Analysis Using RStudio"
#' subtitle: "Linear Regression"
#' author: "Wan Nor Arifin"

# Library ====
# official CRAN
library(foreign)
library(tidyverse)
library(gtsummary)
library(ggplot2)
library(ggpubr)
library(GGally)
library(psych)  # for describe
library(rsq)
library(broom)
library(car)  # for vif
# custom function
# desc_cat()
source("https://raw.githubusercontent.com/wnarifin/medicalstats-in-R/master/functions/desc_cat_fun.R")

# Data ====
coronary = read.dta("coronary.dta")
str(coronary)

## descriptive ----
tbl_summary(coronary)  # median IQR
coronary %>% select(-id) %>%
  tbl_summary(statistic = all_continuous() ~ "{mean} ({sd})", 
              digits = all_continuous() ~ 1)
# customization: http://www.danieldsjoberg.com/gtsummary/index.html
# or a different lecture...

# SLR ====
# example: chol ~ dbp

## explore ----
# descriptive
coronary %>% select(chol, dbp) %>% describe()
# plots
hist_chol = ggplot(coronary, aes(chol)) + geom_histogram(color = "blue", fill = "white")
hist_dbp = ggplot(coronary, aes(dbp)) + geom_histogram(color = "red", fill = "white")
bplot_chol = ggplot(coronary, aes(chol)) + geom_boxplot(color = "blue", )
bplot_dbp = ggplot(coronary, aes(dbp)) + geom_boxplot(color = "red")
ggarrange(hist_chol, bplot_chol, hist_dbp, bplot_dbp)

## fit slr ----
# fit glm, chol ~ dbp
slr_chol = glm(chol ~ dbp, data = coronary)
summary(slr_chol)
tidy(slr_chol, conf.int = TRUE)  # broom package
rsq(slr_chol, adj = T)  # adjusted R2 - penalized for number of predictor p
# scatter plot
plot_slr = ggplot(coronary, aes(x = dbp, y = chol)) + geom_point() + geom_smooth(method = lm)
plot_slr
# plot(chol ~ dbp, data = coronary)
# abline(slr_chol)  # or using built-in graphics

# MLR ====
# chol ~ many factors
names(coronary)  # recall what we have

## explore ----
# descriptive
coronary %>% select(-id, -cad, -race, -gender) %>% describe()  # numerical
coronary %>% select(race, gender) %>% desc_cat()  # categorical
# or just use
# coronary %>% select(race, gender) %>% tbl_summary()
# plots
coronary %>% select(-id, -cad) %>% ggpairs()
# coronary %>% select(-id, -cad) %>% plot()  # or using built-in graphics

## fit univariable / slr ----
# can do univariable one-by-one on your own
slr_chol0 = glm(chol ~ 1, data = coronary)
summary(slr_chol0)
add1(slr_chol0, scope = ~ sbp + dbp + age + bmi + race + gender, test = "LRT")
# all sig. except gender
# may choose only vars p-value < 0.25 to proceed in MLR / expert judgement

## fit multivariable / mlr----
# modeling considerations:
# - select either sbp / dbp! redundant based on plot before, highly correlated
# - gender not sig., may exclude
# - for exercise reason, exclude age
mlr_chol = glm(chol ~ dbp + bmi + race, data = coronary)
summary(mlr_chol)
rsq(mlr_chol, adj = T)

# stepwise
# important to know stepwise/automatic selection is meant for exploratory analysis
# for confirmatory analysis, expert opinion in variable selection is preferable
# both
mlr_chol_stepboth = step(mlr_chol, direction = "both")
summary(mlr_chol_stepboth)
# forward
mlr_chol_stepforward = step(slr_chol0, scope = ~ dbp + bmi + race, 
                            direction = "forward")
summary(mlr_chol_stepforward)  # same with both
# backward
mlr_chol_stepback = step(mlr_chol, direction = "backward")
summary(mlr_chol_stepback)  # same with both

# select: chol ~ dbp + race
mlr_chol_sel = glm(chol ~ dbp + race, data = coronary)
summary(mlr_chol_sel)
tidy(mlr_chol_sel, conf.int = TRUE)
tbl_regression(mlr_chol_sel)
rsq(mlr_chol_sel)

# Multicollinearity, MC
# by Variance Inflation Factor (VIF)
vif(mlr_chol_sel)  # all < 10

# Interaction, *
summary(glm(chol ~ dbp * race, data = coronary))  # dbp * race not sig.
# in R, it is easy to fit interaction by *
# dbp * race will automatically include all vars involved i.e. equal to
# glm(chol ~ dbp + race + dbp:race, data = coronary)
# use : to just include interaction

# model fit: residuals
rraw_chol = resid(mlr_chol_sel)  # unstandardized
hist(rraw_chol)
boxplot(rraw_chol)  # normally distributed

rstd_chol = rstandard(mlr_chol_sel)  # standardized residuals
pstd_chol = scale(predict(mlr_chol_sel))  # standardized predicted values
plot(rstd_chol ~ pstd_chol, xlab = "Std predicted", ylab = "Std residuals")
abline(0, 0)  # normal, linear, equal variance

plot(rraw_chol ~ coronary$dbp, xlab = "DBP", ylab = "Raw Residuals")
abline(0, 0)

# rename the selected model
mlr_chol_final = mlr_chol_sel

# Present the results and interpret ====
tib_mlr = tidy(mlr_chol_final, conf.int = TRUE)
rsq(mlr_chol_final, adj = T)
# use kable to come up with nice table
knitr::kable(tib_mlr, format = "simple")
tbl_regression(mlr_chol_final)
# export it to a csv file for use later
write.csv(tib_mlr, "mlr_final.csv")

# Prediction ====

## predict cholesterol level ----
coronary$pred_chol = predict(mlr_chol_final)
# coronary = coronary %>% mutate(pred_chol = predict(mlr_chol_final))
head(coronary)

## one subject ----
# dbp = 90, race = indian
predict(mlr_chol_final, list(dbp = 90, race = "indian"))

## many subjects ----
new_data = data.frame(dbp = c(90, 90, 90), race = c("malay", "chinese", "indian"))
new_data
predict(mlr_chol_final, new_data)
new_data$pred_chol = predict(mlr_chol_final, new_data)
new_data

# Exercise ====
# - replace dbp with sbp
# - replace race with age
# - obtain coefficient for 5mmHg increase in DBP
# Answer:
# coronary$dbp5 = coronary$dbp/5
# mlr_chol_sel1 = glm(chol ~ dbp5 + race, data = coronary)
# summary(mlr_chol_sel1)
# also can just multiply by 5
