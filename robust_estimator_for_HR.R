###########################################################
############ Robust Variance Estimator #########################
###########################################################
# R codes for variance estimator for log-hazard ratio when using matching with replacement 
# Reference: Austin PC, Cafri G. Variance estimation when using propensity-score matching with replacement with survival or time-to-event outcomes. Statistics in Medicine 2020;39:1623–40. doi:10.1002/sim.8502

# load package
library(survival)

# load data, the dataset of nearest neighbor matching with replacement within specified caliper widths in a 1:5 variable ratio.
data <- readRDS("data.RDS")
# 25339 rows, including 4301 denosumab users and 21038 controls   

# pair.id, matched cluster id;
# subject.id, subject's original id;
# unique.id, unique id for each row.
# treat, exposure (0 for control and 1 for denosumab)

# Cox model with variance estimator that accounts for clustering
# within matched pairs.
cox.nnm.pair <- coxph(Surv(time=start, time2 = follow_up_5, event=event_5) ~ factor(treatment) + cluster(pair.id), data = data)
summary(cox.nnm.pair)

# Cox model with variance estimator that accounts for clustering
# within subjects.
cox.nnm.subject <- coxph(Surv(time=start, time2 = follow_up_5, event=event_5) ~ factor(treatment) + cluster(subjid), data = data)
summary(cox.nnm.subject)

# Cox model with variance estimator that treats all subjects as
# independent.
cox.nnm.cross <- coxph(Surv(time=start, time2 = follow_up_5, event=event_5) ~ factor(treatment) + cluster(unique.id), data=data)
summary(cox.nnm.cross)

# calculations
beta.nnm <- cox.nnm.pair$coef
var.pair.nnm <- cox.nnm.pair$var
var.subject.nnm <- cox.nnm.subject$var
var.cross.nnm <- cox.nnm.cross$var
K.pair.nnm <- nrow(sub[sub$group=="dmab",])
K.subject.nnm <- length(unique(sub$subjid))
K.cross.nnm <- nrow(sub)

# Variance estimator for the log-hazard ratio.
var.nnm <- (K.pair.nnm/(K.pair.nnm-1))*var.pair.nnm +
  (K.subject.nnm/(K.subject.nnm-1))*var.subject.nnm -
  (K.cross.nnm/(K.cross.nnm-1))*var.cross.nnm

# Standard error of estimated log-hazard ratio.
se.nnm <- sqrt(var.nnm)
exp(beta.nnm)
# Upper and lower endpoints of a 95% confidence interval for the hazard ratio
ci.lower.nnm <- exp(beta.nnm - 1.96*se.nnm)
ci.upper.nnm <- exp(beta.nnm + 1.96*se.nnm)
# save results
estimates_cis <- c(exp(beta.nnm),ci.lower.nnm, ci.upper.nnm)
estimates_cis
###########################################################
