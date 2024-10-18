################################################################################
# Purpose: THESIS
# Steps  : Step 1: Load data and libraries
#             - The dataset to be analyzed, contains data from study IMPower 110.
#               It compares the effect of Atezolizumab VS platinum based chemotherapy, 
#               with 572 patients for 38 months. 
#         
# Author : Kalatzi Marilena
# Date   : August 2024
################################################################################

################################ LOAD PACKAGES #################################
library(readODS)
library(survival)
library(flexsurv)
library(ggplot2)
library(survminer)
library(RColorBrewer)
library(viridis)
library(dplyr)
library(muhaz)
library(tidyverse)
library(kableExtra)
library(magrittr)
library(knitr)
################################################################################

################################ LOAD DATA #####################################
# Define a function to read specific rows
LoadData = function(file, sheet, start_row, end_row) {
  # Read the entire sheet
  data = read_ods(file, sheet = sheet)
  # Extract the specific rows
  specific_rows = data[start_row:end_row, ]
  return(specific_rows)
}

FilePath = "C:/Users/kalat/OneDrive/Desktop/Thesis/Publications on 1L NSCLC/Dummy survival ML-NMR study data.ods"
# FilePath = "Dummy survival ML-NMR study data.ods"

AtezoDf = LoadData(FilePath, sheet = 3, start_row = 3771, end_row = 4324)

AtezoDf = subset(AtezoDf, select = - Data_type)
head(AtezoDf)

AtezoDf = data.frame(AtezoDf)
AtezoDf$Arm = factor(AtezoDf$Arm, labels = c("Chemotherapy", "Atezolizumab"))
################################################################################

# ---- Comparator 6 study - Atezolizumab VS Chemotherapy ----

# ---- Kaplan Meier Estimator

ChemoArm = AtezoDf[AtezoDf$Arm == "Chemotherapy",]
AtezoArm = AtezoDf[AtezoDf$Arm == "Atezolizumab",]

OverallKM = survfit(Surv(Time, Event)~1, data = AtezoDf)
ggsurvplot(OverallKM, 
           data             = AtezoDf,
           censor           = T,
           risk.table       = T,
           ylim             = c(0,1),
           xlab             = "Time (in months)",
           title            = "Kaplan-Meier Survival Curve for Overall Survival",
           break.time.by    = 5,
)
OverallKM
# Median survival 16.2 months

TreatKM = survfit(Surv(Time, Event)~Arm, data = AtezoDf)
ggsurvplot(TreatKM, 
           data             = AtezoDf,
           censor           = T,
           risk.table       = T,
           ylim             = c(0,1),
           title            = "Kaplan-Meier with treatment as predictor",
           xlab             = "Time (in months)",
           break.time.by    = 5,
           surv.scale       = 'percent',
           legend.labs      = c('Chemotherapy', 'Atezolizumab'),
           legend.title     = 'Treatment',
           cumevents        = F,
           pval             = T,
           surv.median.line = "hv",
           pval.method      = T)
TreatKM
# Median survival for ChemoArmgroup: 14.2
# Median survival for AtezoArmgroup: 17.8

# ---- Fit PH Cox Model

CoxModel = coxph(Surv(Time, Event)~Arm, 
                 data = AtezoDf,
                 ties = "breslow",robust = F)

cox.zph(CoxModel)

# ---- Shoenfeld residuals
CoxModelSchoenfeld = ggcoxzph(cox.zph(CoxModel))
CoxModelSchoenfeld

# do not reject proportionality test

# ---- Log Cumulative Hazard Plots

Atezolog = survfit(Surv(log(Time),Event)~Arm, data = AtezoDf)

plot(Atezolog, 
     fun  = "cloglog", 
     main = "log-log Survival Curves for PH Assumption Verification",
     xlab = "log(Time)",
     ylab = "log(-log(S(t)))",
     lwd  = 2,
     lty  = c(2,1),
     col  = c("coral","cyan3")
     )
legend(x      = "topleft", 
       xpd    = TRUE,  
       cex    = 1,  
       bty    = "n",
       lty    = c(2,1), 
       legend = c("Chemotherapy", "Atezolizumab"), 
       col    = c("coral","cyan3")
       )

# not showing  proportionality

# ---- Plot of Smoothed hazard curves for the 2 arms 

# Chemotherapy arm curve

ChemoHazard = muhaz(times    = ChemoArm$Time,
                    delta    = ChemoArm$Event,
                    min.time = 0, 
                    max.time = max(ChemoArm$Time)
)

plot(ChemoHazard$est.grid, 
     ChemoHazard$haz.est, 
     type = 'l', 
     xlab = 'Time (months)', 
     ylab = 'Smoothed hazard',
     lwd  = 2, 
     lty  = 1, 
     col  = 'coral', 
     main = 'Chemotherapy Hazard',
     )

# Atezolizumab arm curve

AtezoHazard = muhaz(times    = AtezoArm$Time,
                    delta    = AtezoArm$Event,
                    min.time = 0, 
                    max.time = max(AtezoArm$Time)
)

plot(AtezoHazard$est.grid, 
     AtezoHazard$haz.est, 
     type = 'l', 
     xlab = 'Time (months)', 
     ylab = 'Smoothed hazard',
     lwd  = 2, 
     lty  = 1, 
     col  = 'cyan3', 
     main = 'Atezolizumab Hazard',
)

# PH assumption is not verified - we should model each arm separately

# ---- FSEA ----

# ---- Chemotherapy Arm ----

# ---- Kaplan Meier Estimator

ChemoKM = survfit(Surv(Time, Event)~1, data = ChemoArm)
ggsurvplot(ChemoKM, 
           data             = ChemoArm,
           censor           = T,
           risk.table       = T,
           ylim             = c(0,1),
           title            = "Kaplan-Meier Survival Curve for Chemotherapy Arm",
           xlab             = "Time (in months)",
           break.time.by    = 5,
           surv.median.line = "hv"
           
)
ChemoKM
# Median survival for Chemotherapy group 14.2 months

# ---- Fit Parametric Models
DistributionNames = c('Exponential' ,
                      'Weibull'     , 
                      'Gompertz'    ,
                      'Log-Logistic', 
                      'Log-Normal'  , 
                      'Generalized Gamma'
                      )
 
# ---- Chemotherapy Parametric fit

ChemoExp = flexsurvreg(Surv(Time, Event) ~ 1, data = ChemoArm, dist = "exp")

ChemoWeib = flexsurvreg(Surv(Time, Event) ~ 1, data = ChemoArm, dist = "weibull")

ChemoGomp = flexsurvreg(Surv(Time, Event) ~ 1, data = ChemoArm, dist = "gompertz")

ChemoLogis = flexsurvreg(Surv(Time, Event) ~ 1, data = ChemoArm, dist = "llogis")

ChemoNorm = flexsurvreg(Surv(Time, Event) ~ 1, data = ChemoArm, dist = "lnorm")

ChemoGen = flexsurvreg(Surv(Time, Event) ~ 1, data = ChemoArm, dist = "gengamma")

# Set up a colors palette for the survival curves of the parametric models
colors = brewer.pal(n = 8, name = "Dark2")

plot(ChemoKM$time,
     ChemoKM$surv,
     col = c(0,0), ylim = c(0,1),
     xlab=("Time (in months)"), main="Chemotherapy Arm Parametric Fit", ylab="Survival Probability")  
summary(ChemoExp)[[1]] %$% lines(time, est, col= colors[1], lwd = 2)
summary(ChemoWeib)[[1]] %$% lines(time, est, col=colors[2], lwd = 2)
summary(ChemoGomp)[[1]] %$% lines(time, est, col=colors[3], lwd = 2)
summary(ChemoLogis)[[1]] %$% lines(time, est, col=colors[4], lwd = 2)
summary(ChemoNorm)[[1]] %$% lines(time, est, col=colors[5], lwd = 2)
summary(ChemoGen)[[1]] %$% lines(time, est, col=colors[6], lwd = 2)
lines(ChemoKM$time,
      ChemoKM$surv,col = c(1,0),lwd=2,t="s")  
legend("bottomleft",
       c("KM",DistributionNames),
       col=c(1,colors[1:6]), 
       lty=1, 
       lwd=2,
       xpd    = TRUE,  
       cex    = 1,  
       bty    = "n")

ChemoAic = AIC(ChemoExp,
               ChemoWeib,
               ChemoGomp,
               ChemoLogis,
               ChemoNorm,
               ChemoGen
               )

models = list(ChemoExp,
              ChemoWeib,
              ChemoGomp,
              ChemoLogis,
              ChemoNorm,
              ChemoGen
)
model_names = c("ChemoExp",
                "ChemoWeib",
                "ChemoGomp",
                "ChemoLogis",
                "ChemoNorm",
                "ChemoGen"
)

all_bic = NULL
for (i in 1:length(models)){
  bic = BIC(models[[i]])
  all_bic = c(all_bic,bic)
}
ChemoBic = data.frame(model_names,all_bic)

ChemoAic
ChemoBic

ChemoIC = data.frame(Model = DistributionNames, AIC = ChemoAic$AIC, BIC = ChemoBic$all_bic)
ChemoIC

kable(ChemoIC, 
      col.names = c("Model","AIC", "BIC"), 
      align     = "c",
      caption   = "Parametric Fit for Chemo Arm: Information Criteria"
      ) %>% 
  kable_styling(bootstrap_options = "striped")

#  Log-logistic model has the lowest AIC and BIC values 

# ---- Hazard plot

plot(muhaz(times=ChemoArm$Time, delta = ChemoArm$Event, min.time=0, max.time=max(ChemoArm$Time)),xlab='Months', ylab='Hazard', main='Chemotherapy Arm Hazards',lwd = 2)
plot(ChemoExp, type='hazard', add=T, ci=F, col=colors[1])
plot(ChemoWeib, type='hazard', add=T, ci=F, col=colors[2])
plot(ChemoGomp, type='hazard', add=T, ci=F, col=colors[3])
plot(ChemoLogis, type='hazard', add=T, ci=F, col=colors[4])
plot(ChemoNorm, type='hazard', add=T, ci=F, col=colors[5])
plot(ChemoGen, type='hazard', add=T, ci=F, col=colors[6])

legend('topleft', lty=1, lwd=2, col=c('black',colors[1:6]),
                                              xpd    = TRUE,  
                                              cex    = 1,  
                                              bty    = "n",
       legend=c('Kernel estimate', DistributionNames))

# ---- Cumulative hazard plot 
plot(ChemoKM, lwd=1, fun='cumhaz', xlab='Months', ylab='Cumulative hazard', main='Chemotherapy Arm Cumulative Hazard')

t=seq(from=0, to=max(ChemoKM$time), length=1000)

lines(x=t, Hexp(t, rate=ChemoExp$res[1], log=F), lwd=2, col=colors[1])
lines(x=t, Hweibull(t, shape=ChemoWeib$res[1, 1], scale=ChemoWeib$res[2, 1], log=F), lwd=2, col=colors[2])
lines(x=t, Hgompertz(t, shape=ChemoGomp$res[1, 1], rate=ChemoGomp$res[2, 1], log=F), lwd=2, col=colors[3])
lines(x=t, Hllogis(t, shape=ChemoLogis$res[1, 1], scale=ChemoLogis$res[2, 1], log=F), lwd=2, col=colors[4])
lines(x=t, Hlnorm(t, meanlog=ChemoNorm$res[1, 1], sdlog=ChemoNorm$res[2, 1], log=F), lwd=2, col=colors[5])
lines(x=t, Hgengamma(t, mu=ChemoGen$res[1, 1], sigma=ChemoGen$res[2, 1], Q=ChemoGen$res[3, 1]), lwd=2, col=colors[6])

legend('topleft', lty=1, lwd=2, col=c('black',colors[1:6]),
       xpd    = TRUE,  
       cex    = 1,  
       bty    = "n",
       legend=c('Kaplan Meier', DistributionNames))

# ---- Log-cumulative hazard plot 
plot(ChemoKM, lwd=1, fun='cloglog', xlab='Log-Months', ylab='Cumulative hazard', main='Chemotherapy Arm Log-Cumulative Hazard')

t=seq(from=0, to=max(ChemoKM$time), length=1000)

lines(x=t, Hexp(t, rate=ChemoExp$res[1], log=T), lwd=2, col=colors[1])
lines(x=t, Hweibull(t, shape=ChemoWeib$res[1, 1], scale=ChemoWeib$res[2, 1], log=T), lwd=2, col=colors[2])
lines(x=t, Hgompertz(t, shape=ChemoGomp$res[1, 1], rate=ChemoGomp$res[2, 1], log=T), lwd=2, col=colors[3])
lines(x=t, Hllogis(t, shape=ChemoLogis$res[1, 1], scale=ChemoLogis$res[2, 1], log=T), lwd=2, col=colors[4])
lines(x=t, Hlnorm(t, meanlog=ChemoNorm$res[1, 1], sdlog=ChemoNorm$res[2, 1], log=T), lwd=2, col=colors[5])
lines(x=t, log(Hgengamma(t, mu=ChemoGen$res[1, 1], sigma=ChemoGen$res[2, 1], Q=ChemoGen$res[3, 1])), lwd=2, col=colors[6])


legend('topleft', lty=1, lwd=2, col=c('black',colors[1:6]),
       xpd    = TRUE,  
       cex    = 1,  
       bty    = "n",
       legend=c('Kaplan Meier', DistributionNames))


# ---- Extrapolation Parametric Survival Curves

Years = 20
ExtrapHorizon = 12*Years 
times = seq(0,ExtrapHorizon, by = 1)

ExtrapChemoExp = summary(ChemoExp, t = times)[[1]]
ExtrapChemoWeib = summary(ChemoWeib, t = times)[[1]]
ExtrapChemoGomp = summary(ChemoGomp, t = times)[[1]]
ExtrapChemoLogis = summary(ChemoLogis, t = times)[[1]]
ExtrapChemoNorm = summary(ChemoNorm, t = times)[[1]]
ExtrapChemoGen = summary(ChemoGen, t = times)[[1]]

# Plot the extrapolated survival curves
plot(NULL, xlim = c(0,12*Years), ylim = c(0, 1), 
     xlab = "Time (Years)", ylab = "Survival Probability", 
     main = "Extrapolated Survival Curves for Chemotherapy Group",
     xaxt = "n")
axis(side=1, at=12*seq(from=0, to=ExtrapHorizon, by=1), labels=as.character(0:ExtrapHorizon))
lines(ChemoKM, lwd = 2)
lines(est~time, data = ExtrapChemoExp, type = "l", col = colors[1],lwd = 2)
lines(est~time, data = ExtrapChemoWeib, type = "l", col = colors[2],lwd = 2)
lines(est~time, data = ExtrapChemoGomp, type = "l", col = colors[3],lwd = 2)
lines(est~time, data = ExtrapChemoLogis, type = "l", col = colors[4],lwd = 2)
lines(est~time, data = ExtrapChemoNorm, type = "l", col = colors[5],lwd = 2)
lines(est~time, data = ExtrapChemoGen, type = "l", col = colors[6],lwd = 2)

legend("topright", 
       legend = c("KM",DistributionNames),
       col = c(1,colors),
       lty = 1,
       lwd = 2,
       bty = "n")

# ---- Parametric models Mean Survival Times

MeanExp = round(mean_exp(rate = ChemoExp$res[1]), 1)    
MeanWeib = round(mean_weibull(shape=ChemoWeib$res[1, 1], scale=ChemoWeib$res[2, 1]), 1)                            # Weibull                        # Log-logistic
MeanGomp = round(mean_gompertz(shape=ChemoGomp$res[1, 1], rate=ChemoGomp$res[2, 1]), 1)     
MeanLogis = round(mean_llogis(shape=ChemoLogis$res[1, 1], scale=ChemoLogis$res[2, 1]), 1)   
MeanNorm = round(mean_lnorm(meanlog=ChemoNorm$res[1, 1], sdlog=ChemoNorm$res[2, 1]), 1)                        # Log-Normal
MeanGen = round(mean_gengamma(mu = ChemoGen$res[1,1], sigma = ChemoGen$res[2,1],Q = ChemoGen$res[3,1]),1)

ChemoParamMean = c(MeanExp, MeanWeib, MeanGomp, MeanLogis, MeanNorm, MeanGen)
ChemoParamMeanY = round(c(MeanExp, MeanWeib, MeanGomp, MeanLogis, MeanNorm, MeanGen)/12,2)
ChemoParamMeandf = data.frame(Distribution                  = DistributionNames, 
                            `Mean Survival Time (Months)` = ChemoParamMean,
                            `Mean Survival Time (Years)`  = ChemoParamMeanY)

kable(ChemoParamMeandf, 
      col.names = c("Distribution","Mean Survival Time (Months)", "Mean Survival Time (Years)"),
      caption = "Mean Survival Times for 20 year horizon of parametric distributions",
      align = "c") %>% 
  kable_styling(bootstrap_options = "striped")

# Survival curve for each parametric model with shaded AUC

# Define a function to plot survival curves with AUC
plot_auc_curve = function(times, survival_probs, col, model_name) {
  plot(times, survival_probs, type = "l", col = col, lwd = 2, ylim = c(0, 1), xlim = c(0, max(times)),
       xlab = "Time (Years)", ylab = "Survival Probability", main = paste("AUC for", model_name, "of Chemo"),xaxt = "n")
  axis(side=1, at=12*seq(from=0, to=ExtrapHorizon, by=1), labels=as.character(0:ExtrapHorizon))
  polygon(c(times, rev(times)), c(survival_probs, rep(0, length(survival_probs))), col = adjustcolor(col, alpha.f = 0.3), border = NA)
}

# Plot each survival curve with AUC
par(mfrow = c(2, 3)) 
plot_auc_curve(ExtrapChemoExp$time, ExtrapChemoExp$est, colors[1], "Exponential")
plot_auc_curve(ExtrapChemoWeib$time, ExtrapChemoWeib$est, colors[2], "Weibull")
plot_auc_curve(ExtrapChemoGomp$time, ExtrapChemoGomp$est, colors[3], "Gompertz")
plot_auc_curve(ExtrapChemoLogis$time, ExtrapChemoLogis$est, colors[4], "Log-logistic")
plot_auc_curve(ExtrapChemoNorm$time, ExtrapChemoNorm$est, colors[5], "Log-normal")
plot_auc_curve(ExtrapChemoGen$time, ExtrapChemoGen$est, colors[6], "Generalized Gamma")
dev.off()

# ---- Spline Modeling

ChemoSpline1 = flexsurvspline(Surv(Time, Event) ~ 1, data = ChemoArm, k = 1, scale = "hazard")
ChemoSpline2 = flexsurvspline(Surv(Time, Event) ~ 1, data = ChemoArm, k = 2, scale = "hazard")
ChemoSpline3 = flexsurvspline(Surv(Time, Event) ~ 1, data = ChemoArm, k = 3, scale = "hazard")
ChemoSplineDef = flexsurvspline(Surv(Time, Event) ~ 1, data = ChemoArm, knots = log(c(19)), scale = "hazard")
ChemoSplineDef2 = flexsurvspline(Surv(Time, Event) ~ 1, data = ChemoArm, knots = log(c(19,29)), scale = "hazard")
ChemoSplineDef3 = flexsurvspline(Surv(Time, Event) ~ 1, data = ChemoArm, knots = log(c(11,19,29)), scale = "hazard")

SplineAIC = AIC(ChemoSpline1, ChemoSpline2, ChemoSpline3,ChemoSplineDef,ChemoSplineDef2,ChemoSplineDef3)

SplineModels = list(ChemoSpline1,
                 ChemoSpline2,
                 ChemoSpline3,
                 ChemoSplineDef,
                 ChemoSplineDef2,
                 ChemoSplineDef3)

all_spline_bic = NULL
for (i in 1:length(SplineModels)){
  bic = BIC(SplineModels[[i]])
  all_spline_bic = c(all_spline_bic,bic)
}

SplineNames = c("1 knot", 
                "2 knots", 
                "3 knots", 
                "Knot in month 19",
                "Knots in months 19, 29",
                "Knots in months 11, 19, 29")

ChemoSplineBic = data.frame(SplineNames,all_spline_bic)


plot(ChemoKM$time,
     ChemoKM$surv,
     col = c(0,0), ylim = c(0,1),
     xlab=("Time (in months)"), main="Chemotherapy Arm Spline Fit", ylab="Survival Probability")  
summary(ChemoSpline1)[[1]] %$% lines(time, est, col= colors[1], lwd = 2)
summary(ChemoSpline2)[[1]] %$% lines(time, est, col=colors[2], lwd = 2)
summary(ChemoSpline3)[[1]] %$% lines(time, est, col=colors[3], lwd = 2)
summary(ChemoSplineDef)[[1]] %$% lines(time, est, col=colors[4], lwd = 2)
summary(ChemoSplineDef2)[[1]] %$% lines(time, est, col=colors[5], lwd = 2)
summary(ChemoSplineDef3)[[1]] %$% lines(time, est, col=colors[6], lwd = 2)
lines(ChemoKM$time,
      ChemoKM$surv,col = c(1,0),lwd=2,t="s")  
legend("bottomleft",
       c("KM",SplineNames),
       col=c(1,colors[1:6]), 
       lty=1, 
       lwd=2,
       xpd    = TRUE,  
       cex    = 1,  
       bty    = "n")

# Spline model with specified knots in months 11, 19 and 29 is the best model
# among all 12 models fitted (standard parametric and spline based)

ChemoSplineIC = data.frame(Model = SplineNames, AIC = SplineAIC$AIC, BIC = ChemoSplineBic$all_spline_bic)
ChemoSplineIC

kable(ChemoSplineIC, 
      col.names = c("Model","AIC", "BIC"), 
      align = "c",
      caption = "Spline Fit for Chemo Arm: Information Criteria") %>% 
  kable_styling(bootstrap_options = "striped")

All_IC = data.frame(Model = c(model_names,SplineNames), 
                    AIC = c(ChemoAic$AIC, SplineAIC$AIC),
                    BIC = c(ChemoBic$all_bic, ChemoSplineBic$all_spline_bic))

All_IC

All_IC_sorted = All_IC[order(All_IC$AIC), ]
All_IC_sorted

BIC_IC_sorted = All_IC[order(All_IC$BIC),]

BIC_IC_sorted


plot(ChemoKM$time,
     ChemoKM$surv,
     col = c(0,0), ylim = c(0,1),
     xlab=("Time (in months)"), main="Chemotherapy Arm Top Fits", ylab="Survival Probability")
summary(ChemoLogis)[[1]] %$% lines(time, est, col=colors[4], lwd = 2)
summary(ChemoNorm)[[1]] %$% lines(time, est, col=colors[5], lwd = 2)
summary(ChemoSplineDef)[[1]] %$% lines(time, est, col="skyblue3", lwd = 2)
summary(ChemoSpline2)[[1]] %$% lines(time, est, col="firebrick", lwd = 2)
lines(ChemoKM$time,
      ChemoKM$surv,col = c(1,0),lwd=2,t="s")  
legend("bottomleft",
       c("KM",DistributionNames[c(4,5)],SplineNames[c(4,2)]),
       col=c(1,colors[4:5],c("skyblue3","firebrick")), 
       lty=1, 
       lwd=2,
       xpd    = TRUE,  
       cex    = 1,  
       bty    = "n")


# ---- Extrapolation for Splines

Years = 20
ExtrapHorizon = 12*Years 
times = seq(0,ExtrapHorizon, by = 1)

ExtrapSpline1 = summary(ChemoSpline1, t = times)[[1]]
ExtrapSpline2 = summary(ChemoSpline2, t = times)[[1]]
ExtrapSpline3 = summary(ChemoSpline3, t = times)[[1]]
ExtrapSpline1_Custom = summary(ChemoSplineDef, t = times)[[1]]
ExtrapSpline2_Custom = summary(ChemoSplineDef2, t = times)[[1]]
ExtrapSpline3_Custom = summary(ChemoSplineDef3, t = times)[[1]]

# Plot the extrapolated survival curves
plot(NULL, xlim = c(0,12*Years), ylim = c(0, 1), 
     xlab = "Time (Years)", ylab = "Survival Probability", 
     main = "Extrapolated Survival Curves for Chemotherapy Group",
     xaxt = "n")
axis(side=1, at=12*seq(from=0, to=ExtrapHorizon, by=1), labels=as.character(0:ExtrapHorizon))
lines(ChemoKM, lwd = 2)
lines(est~time, data = ExtrapSpline1, type = "l", col = colors[1],lwd = 2)
lines(est~time, data = ExtrapSpline2, type = "l", col = colors[2],lwd = 2)
lines(est~time, data = ExtrapSpline3, type = "l", col = colors[3],lwd = 2)
lines(est~time, data = ExtrapSpline1_Custom, type = "l", col = colors[4],lwd = 2)
lines(est~time, data = ExtrapSpline2_Custom, type = "l", col = colors[5],lwd = 2)
lines(est~time, data = ExtrapSpline3_Custom, type = "l", col = colors[6],lwd = 2)

legend("topright", 
       legend = c("KM",SplineNames),
       col = c(1,colors),
       lty = 1,
       lwd = 2,
       bty = "n")

# Mean survival times
ChemoSpline1Mean = mean_survspline1(gamma0=ChemoSpline1$res[1, 1],
                 gamma1=ChemoSpline1$res[2, 1],
                 gamma2=ChemoSpline1$res[3, 1],
                 knots=ChemoSpline1$knots)        

ChemoSpline2Mean = mean_survspline2(gamma0=ChemoSpline2$res[1, 1],
                                    gamma1=ChemoSpline2$res[2, 1],
                                    gamma2=ChemoSpline2$res[3, 1],
                                    gamma3=ChemoSpline2$res[4, 1],
                                    knots=ChemoSpline2$knots)     

ChemoSpline3Mean = mean_survspline3(gamma0=ChemoSpline3$res[1, 1],
                                    gamma1=ChemoSpline3$res[2, 1],
                                    gamma2=ChemoSpline3$res[3, 1],
                                    gamma3=ChemoSpline3$res[4, 1],
                                    gamma4=ChemoSpline3$res[5, 1],
                                    knots=ChemoSpline3$knots)

ChemoSplineDefMean = mean_survspline1(gamma0=ChemoSplineDef$res[1, 1],
                                      gamma1=ChemoSplineDef$res[2, 1],
                                      gamma2=ChemoSplineDef$res[3, 1],
                                      knots=ChemoSplineDef$knots)

ChemoSplineDef2Mean = mean_survspline2(gamma0=ChemoSplineDef2$res[1, 1],
                                       gamma1=ChemoSplineDef2$res[2, 1],
                                       gamma2=ChemoSplineDef2$res[3, 1],
                                       gamma3=ChemoSplineDef2$res[4, 1],
                                       knots=ChemoSplineDef2$knots) 

ChemoSplineDef3Mean = mean_survspline3(gamma0=ChemoSplineDef3$res[1, 1],
                                       gamma1=ChemoSplineDef3$res[2, 1],
                                       gamma2=ChemoSplineDef3$res[3, 1],
                                       gamma3=ChemoSplineDef3$res[4, 1],
                                       gamma4=ChemoSplineDef3$res[5, 1],
                                       knots=ChemoSplineDef3$knots)

# ---- Extrapolation Parametric Survival Curves

Years = 20
ExtrapHorizon = 12*Years 
times = seq(0,ExtrapHorizon, by = 1)

ExtrapSplineDef = summary(ChemoSplineDef, t = times)[[1]]
ExtrapSpline2   = summary(ChemoSpline2, t = times)[[1]]

# Plot the extrapolated survival curves
plot(NULL, xlim = c(0,12*Years), ylim = c(0, 1), 
     xlab = "Time (Years)", ylab = "Survival Probability", 
     main = "Extrapolated Survival Curves for Chemotherapy Group Best Fits",
     xaxt = "n")
axis(side=1, at=12*seq(from=0, to=ExtrapHorizon, by=1), labels=as.character(0:ExtrapHorizon))
lines(ChemoKM, lwd = 2)
lines(est~time, data = ExtrapChemoLogis, type = "l", col = colors[4],lwd = 2)
lines(est~time, data = ExtrapChemoNorm, type = "l", col = colors[5],lwd = 2)
lines(est~time, data = ExtrapSplineDef, type = "l", col = "skyblue3", lwd = 2)
lines(est~time, data = ExtrapSpline2, type = "l", col = "firebrick", lwd = 2)
legend("topright",
       c("KM",DistributionNames[c(4,5)],SplineNames[c(4,2)]),
       col=c(1,colors[4:5],c("skyblue3","firebrick")), 
       lty=1, 
       lwd=2,
       xpd    = TRUE,  
       cex    = 1,  
       bty    = "n")


# ---- Survival probabilities for 1, 5 and 10 years
# Time points of interest in months 

time_points = c(12, 60, 120) 

SurvExp = summary(ChemoExp, t = time_points)[[1]]$est * 100
SurvWeib = summary(ChemoWeib, t = time_points)[[1]]$est * 100
SurvGomp = summary(ChemoGomp, t = time_points)[[1]]$est * 100
SurvLogis = summary(ChemoLogis, t = time_points)[[1]]$est * 100
SurvNorm = summary(ChemoNorm, t = time_points)[[1]]$est * 100
SurvGen = summary(ChemoGen, t = time_points)[[1]]$est * 100

# Combine results into a data frame
ChemoParamSurv = data.frame(
  Distribution = c("Exponential", "Weibull", "Gompertz", "Log-Logistic", "Log-Normal", "Generalized Gamma"),
  `1-Year Survival Probability (%)` = c(SurvExp[1], SurvWeib[1], SurvGomp[1], SurvLogis[1], SurvNorm[1], SurvGen[1]),
  `5-Year Survival Probability (%)` = c(SurvExp[2], SurvWeib[2], SurvGomp[2], SurvLogis[2], SurvNorm[2], SurvGen[2]),
  `10-Year Survival Probability (%)` = c(SurvExp[3], SurvWeib[3], SurvGomp[3], SurvLogis[3], SurvNorm[3], SurvGen[3])
)

kable(ChemoParamSurv, 
      col.names = c("Distribution", "1-Year Survival Probability (%)", "5-Year Survival Probability (%)", "10-Year Survival Probability (%)"),
      caption = "Predicted Survival Probabilities for Different Parametric Distributions",
      align = "c") %>% 
  kable_styling(bootstrap_options = "striped")

# ---- Spline Survival Probabilities

SurvChemoSpline1 = summary(ChemoSpline1, t = time_points)[[1]]$est * 100
SurvChemoSpline2 = summary(ChemoSpline2, t = time_points)[[1]]$est * 100
SurvChemoSpline3 = summary(ChemoSpline3, t = time_points)[[1]]$est * 100
SurvChemoSplineDef = summary(ChemoSplineDef, t = time_points)[[1]]$est * 100
SurvChemoSplineDef2 = summary(ChemoSplineDef2, t = time_points)[[1]]$est * 100
SurvChemoSplineDef3 = summary(ChemoSplineDef3, t = time_points)[[1]]$est * 100

# ---- Parametric and Spline models Survival Probabilities

# Spline Model Results
ChemoSplineSurv = data.frame(
  Distribution = c("Spline 1", "Spline 2", "Spline 3", "Spline Def (log(19))", "Spline Def (log(19,29))", "Spline Def (log(11,19,29))"),
  `1-Year Survival Probability (%)` = c(SurvChemoSpline1[1], SurvChemoSpline2[1], SurvChemoSpline3[1], SurvChemoSplineDef[1], SurvChemoSplineDef2[1], SurvChemoSplineDef3[1]),
  `5-Year Survival Probability (%)` = c(SurvChemoSpline1[2], SurvChemoSpline2[2], SurvChemoSpline3[2], SurvChemoSplineDef[2], SurvChemoSplineDef2[2], SurvChemoSplineDef3[2]),
  `10-Year Survival Probability (%)` = c(SurvChemoSpline1[3], SurvChemoSpline2[3], SurvChemoSpline3[3], SurvChemoSplineDef[3], SurvChemoSplineDef2[3], SurvChemoSplineDef3[3])
)

# Combine Parametric and Spline Results
ChemoCombSurv = rbind(ChemoParamSurv, ChemoSplineSurv)

kable(ChemoCombSurv, 
      col.names = c("Distribution", "1-Year Survival Probability (%)", "5-Year Survival Probability (%)", "10-Year Survival Probability (%)"),
      caption = "Predicted Survival Probabilities for Chemotherapy Arm (Parametric and Spline Models)",
      align = "c") %>% 
  kable_styling(bootstrap_options = "striped")

# ---- Atezolizumab Arm ----

# ---- Kaplan Meier Estimator

AtezoKM = survfit(Surv(Time, Event)~1, data = AtezoArm)
ggsurvplot(AtezoKM, 
           data             = AtezoArm,
           censor           = T,
           risk.table       = T,
           ylim             = c(0,1),
           xlab             = "Time (in months)",
           title            = "Kaplan-Meier Survival Curve for Atezolizumab Arm",
           break.time.by    = 5,
           palette = "cyan3"
)
AtezoKM
# Median survival for Atezolizumab group 17.8 months

# ---- Fit Parametric Models

# ---- Atezolizumab Parametric fit

AtezoExp = flexsurvreg(Surv(Time, Event) ~ 1, data = AtezoArm, dist = "exp")

AtezoWeib = flexsurvreg(Surv(Time, Event) ~ 1, data = AtezoArm, dist = "weibull")

AtezoGomp = flexsurvreg(Surv(Time, Event) ~ 1, data = AtezoArm, dist = "gompertz")

AtezoLogis = flexsurvreg(Surv(Time, Event) ~ 1, data = AtezoArm, dist = "llogis")

AtezoNorm = flexsurvreg(Surv(Time, Event) ~ 1, data = AtezoArm, dist = "lnorm")

AtezoGen = flexsurvreg(Surv(Time, Event) ~ 1, data = AtezoArm, dist = "gengamma")

plot(AtezoKM$time,
     AtezoKM$surv,
     col = c(0,0), ylim = c(0,1),
     xlab=("Time (in months)"), main="Atezolizumab Arm", ylab="Survival Probability")  
summary(AtezoExp)[[1]] %$% lines(time, est, col= colors[1], lwd = 2)
summary(AtezoWeib)[[1]] %$% lines(time, est, col=colors[2], lwd = 2)
summary(AtezoGomp)[[1]] %$% lines(time, est, col=colors[3], lwd = 2)
summary(AtezoLogis)[[1]] %$% lines(time, est, col=colors[4], lwd = 2)
summary(AtezoNorm)[[1]] %$% lines(time, est, col=colors[5], lwd = 2)
summary(AtezoGen)[[1]] %$% lines(time, est, col=colors[6], lwd = 2)
lines(AtezoKM$time,
      AtezoKM$surv,col = c(1,0),lwd=2,t="s")  
legend("bottomleft",
       c("KM",DistributionNames),
       col=c(1,colors[1:6]), 
       lty=1, 
       lwd=2,
       xpd    = TRUE,  
       cex    = 1,  
       bty    = "n")

AtezoAic = AIC(AtezoExp,
               AtezoWeib,
               AtezoGomp,
               AtezoLogis,
               AtezoNorm,
               AtezoGen
)

models = list(AtezoExp,
              AtezoWeib,
              AtezoGomp,
              AtezoLogis,
              AtezoNorm,
              AtezoGen
)

model_names = c("AtezoExp",
                "AtezoWeib",
                "AtezoGomp",
                "AtezoLogis",
                "AtezoNorm",
                "AtezoGen"
)

all_bic = NULL
for (i in 1:length(models)){
  bic = BIC(models[[i]])
  all_bic = c(all_bic,bic)
}
AtezoBic = data.frame(model_names,all_bic)

AtezoAic
AtezoBic

AtezoIC = data.frame(Model = DistributionNames, AIC = AtezoAic$AIC, BIC = AtezoBic$all_bic)
AtezoIC

# Find the row with the minimum AIC and BIC
min_aic_row = which.min(AtezoIC$AIC)
min_bic_row = which.min(AtezoIC$BIC)

# Create the table with bold formatting for the rows with minimum AIC and BIC
kable(AtezoIC, col.names = c("Model","AIC", "BIC"), 
      align = "c", 
      caption = "Atezolizumab Arm Parametric Fit") %>% 
  kable_styling(bootstrap_options = "striped") %>%
  row_spec(min_aic_row, bold = TRUE, color = "red") %>%
  row_spec(min_bic_row, bold = TRUE, color = "red")


BICs_spl=data.table::setDT(as.data.frame(all_bic), keep.rownames='BIC')
BICs_spl_long=BICs_spl %>% gather(Model, BIC)
G2=ggplot(AtezoBic, aes(x=model_names, y=all_bic, col=model_names)) +
  geom_segment(aes(x=reorder(model_names, all_bic),
                   xend=model_names,
                   y=min(all_bic) - 10,
                   yend=all_bic),
               size=10) +
  theme_bw() + coord_flip() + ggtitle('Pembrolizumab OS') +
  guides(col='none')


# Log-normal model has the lowest AIC and BIC values 

# ---- Hazard plots
dev.off()
plot(muhaz(times=AtezoArm$Time, delta = AtezoArm$Event, min.time=0, max.time=max(AtezoArm$Time)),xlab='Months', ylab='Hazard', main='Atezolizumab Arm Hazard', lwd = 2,ylim = c(0,0.06))
plot(AtezoExp, type='hazard', add=T, ci=F, col=colors[1])
plot(AtezoWeib, type='hazard', add=T, ci=F, col=colors[2])
plot(AtezoGomp, type='hazard', add=T, ci=F, col=colors[3])
plot(AtezoLogis, type='hazard', add=T, ci=F, col=colors[4])
plot(AtezoNorm, type='hazard', add=T, ci=F, col=colors[5])
plot(AtezoGen, type='hazard', add=T, ci=F, col=colors[6])

legend('bottomleft', 
       lty=1, 
       lwd=2, 
       col=c('black', colors[1:6]),
       legend=c('Kernel estimate', DistributionNames),
       bty = "n")

# ---- Cumulative hazard plot 
plot(AtezoKM, lwd=1, fun='cumhaz', xlab='Months', ylab='Cumulative hazard', main='Atezolizumab Arm Cumulative Hazard')

t=seq(from=0, to=max(AtezoKM$time), length=1000)

lines(x=t, Hexp(t, rate=AtezoExp$res[1], log=F), lwd=2, col=colors[1])
lines(x=t, Hweibull(t, shape=AtezoWeib$res[1, 1], scale=AtezoWeib$res[2, 1], log=F), lwd=2, col=colors[2])
lines(x=t, Hgompertz(t, shape=AtezoGomp$res[1, 1], rate=AtezoGomp$res[2, 1], log=F), lwd=2, col=colors[3])
lines(x=t, Hllogis(t, shape=AtezoLogis$res[1, 1], scale=AtezoLogis$res[2, 1], log=F), lwd=2, col=colors[4])
lines(x=t, Hlnorm(t, meanlog=AtezoNorm$res[1, 1], sdlog=AtezoNorm$res[2, 1], log=F), lwd=2, col=colors[5])
lines(x=t, Hgengamma(t, mu=AtezoGen$res[1, 1], sigma=AtezoGen$res[2, 1], Q=AtezoGen$res[3, 1]), lwd=2, col=colors[6])

legend('topleft', lty=1, lwd=2, col=c('black',colors[1:6]),
       xpd    = TRUE,  
       cex    = 1,  
       bty    = "n",
       legend=c('Kaplan Meier', DistributionNames))

# ---- Log-cumulative hazard plot 
plot(AtezoKM, lwd=1, fun='cloglog', xlab='Log-Months', ylab='Cumulative hazard', main='Atezolizumab Arm Log-Cumulative Hazard')

t=seq(from=0, to=max(AtezoKM$time), length=1000)

lines(x=t, Hexp(t, rate=AtezoExp$res[1], log=T), lwd=2, col=colors[1])
lines(x=t, Hweibull(t, shape=AtezoWeib$res[1, 1], scale=AtezoWeib$res[2, 1], log=T), lwd=2, col=colors[2])
lines(x=t, Hgompertz(t, shape=AtezoGomp$res[1, 1], rate=AtezoGomp$res[2, 1], log=T), lwd=2, col=colors[3])
lines(x=t, Hllogis(t, shape=AtezoLogis$res[1, 1], scale=AtezoLogis$res[2, 1], log=T), lwd=2, col=colors[4])
lines(x=t, Hlnorm(t, meanlog=AtezoNorm$res[1, 1], sdlog=AtezoNorm$res[2, 1], log=T), lwd=2, col=colors[5])
lines(x=t, log(Hgengamma(t, mu=AtezoGen$res[1, 1], sigma=AtezoGen$res[2, 1], Q=AtezoGen$res[3, 1])), lwd=2, col=colors[6])


legend('topleft', lty=1, lwd=2, col=c('black',colors[1:6]),
       xpd    = TRUE,  
       cex    = 1,  
       bty    = "n",
       legend=c('Kaplan Meier', DistributionNames))


# ---- Extrapolation Parametric Survival Curves

ExtrapHorizon = 12*Years
times = seq(0,ExtrapHorizon, by = 1)

ExtrapAtezoExp = summary(AtezoExp, t = times)[[1]]
ExtrapAtezoWeib = summary(AtezoWeib, t = times)[[1]]
ExtrapAtezoGomp = summary(AtezoGomp, t = times)[[1]]
ExtrapAtezoLogis = summary(AtezoLogis, t = times)[[1]]
ExtrapAtezoNorm = summary(AtezoNorm, t = times)[[1]]
ExtrapAtezoGen = summary(AtezoGen, t = times)[[1]]

# Plot the extrapolated survival curves
plot(NULL, xlim = c(0,12*Years), ylim = c(0, 1), 
     xlab = "Time (Years)", ylab = "Survival Probability", 
     main = "Extrapolated Survival Curves for Atezolizumab Group",
     xaxt = "n")
axis(side=1, at=12*seq(from=0, to=ExtrapHorizon, by=1), labels=as.character(0:ExtrapHorizon))
lines(AtezoKM, lwd = 2)
lines(est~time, data = ExtrapAtezoExp, type = "l", col = colors[1],lwd = 2)
lines(est~time, data = ExtrapAtezoWeib, type = "l", col = colors[2],lwd = 2)
lines(est~time, data = ExtrapAtezoGomp, type = "l", col = colors[3],lwd = 2)
lines(est~time, data = ExtrapAtezoLogis, type = "l", col = colors[4],lwd = 2)
lines(est~time, data = ExtrapAtezoNorm, type = "l", col = colors[5],lwd = 2)
lines(est~time, data = ExtrapAtezoGen, type = "l", col = colors[6],lwd = 2)
legend("topright", 
       legend = c("KM",DistributionNames),
       col = c(1,colors[1:6]),
       lty = 1,
       lwd = 2,
       bty = "n")

# ---- Parametric models Mean Survival Times

AtezoMeanExp = round(mean_exp(rate = AtezoExp$res[1]), 1)    
AtezoMeanWeib = round(mean_weibull(shape=AtezoWeib$res[1, 1], scale=AtezoWeib$res[2, 1]), 1)                            # Weibull                        # Log-logistic
AtezoMeanGomp = round(mean_gompertz(shape=AtezoGomp$res[1, 1], rate=AtezoGomp$res[2, 1]), 1)   
# integral divergent

AtezoMeanLogis = round(mean_llogis(shape=AtezoLogis$res[1, 1], scale=AtezoLogis$res[2, 1]), 1) 
AtezoMeanNorm = round(mean_lnorm(meanlog=AtezoNorm$res[1, 1], sdlog=AtezoNorm$res[2, 1]), 1)   
AtezoMeanGen = round(mean_gengamma(mu = AtezoGen$res[1,1], sigma = AtezoGen$res[2,1],Q = AtezoGen$res[3,1]),1)

AtezoParamMean = c(AtezoMeanExp, AtezoMeanWeib, AtezoMeanLogis, AtezoMeanNorm, AtezoMeanGen)
AtezoParamMeanY = c(AtezoMeanExp, AtezoMeanWeib, AtezoMeanLogis, AtezoMeanNorm, AtezoMeanGen)/12
AtezoParamMeanDF = data.frame(Distribution                  = c(DistributionNames[1:2],DistributionNames[4:6]),
                            `Mean Survival Time (Months)` = AtezoParamMean,
                            `Mean Survival Time (Years)`  = round(AtezoParamMeanY,2))

kable(AtezoParamMeanDF, 
      col.names = c("Distribution","Mean Survival Time (Months)", "Mean Survival Time (Years)"), 
      align = "c",
      caption = "Atezolizumab Arm Mean Survival Times of Parametric Fits for 20 year horizon") %>% 
  kable_styling(bootstrap_options = "striped")

# Survival curve for each parametric model with shaded AUC

# Define a function to plot survival curves with AUC
plot_auc_curve = function(times, survival_probs, col, model_name) {
  plot(times, survival_probs, type = "l", col = col, lwd = 2, ylim = c(0, 1), xlim = c(0, max(times)),
       xlab = "Time (Years)", ylab = "Survival Probability", main = paste("AUC for", model_name, "of Atezolizumab"),xaxt = "n")
  axis(side=1, at=12*seq(from=0, to=ExtrapHorizon, by=1), labels=as.character(0:ExtrapHorizon))
  polygon(c(times, rev(times)), c(survival_probs, rep(0, length(survival_probs))), col = adjustcolor(col, alpha.f = 0.3), border = NA)
}

# Plot each survival curve with AUC
par(mfrow = c(2, 3))
plot_auc_curve(ExtrapAtezoExp$time, ExtrapAtezoExp$est, colors[1], "Exponential")
plot_auc_curve(ExtrapAtezoWeib$time, ExtrapAtezoWeib$est, colors[2], "Weibull")
plot_auc_curve(ExtrapAtezoGomp$time, ExtrapAtezoGomp$est, colors[3], "Gompertz")
plot_auc_curve(ExtrapAtezoLogis$time, ExtrapAtezoLogis$est, colors[4], "Log-logistic")
plot_auc_curve(ExtrapAtezoNorm$time, ExtrapAtezoNorm$est, colors[5], "Log-normal")
plot_auc_curve(ExtrapAtezoGen$time, ExtrapAtezoGen$est, colors[6], "Generalized Gamma")
dev.off()

# ---- Spline Modeling

AtezoSpline1 = flexsurvspline(Surv(Time, Event) ~ 1, data = AtezoArm, k = 1, scale = "hazard")
AtezoSpline2 = flexsurvspline(Surv(Time, Event) ~ 1, data = AtezoArm, k = 2, scale = "hazard")
AtezoSpline3 = flexsurvspline(Surv(Time, Event) ~ 1, data = AtezoArm, k = 3, scale = "hazard")
AtezoSplineDef = flexsurvspline(Surv(Time, Event) ~ 1, data = AtezoArm, knots = log(c(31)), scale = "hazard")
AtezoSplineDef2 = flexsurvspline(Surv(Time, Event) ~ 1, data = AtezoArm, knots = log(c(15,31)), scale = "hazard")
AtezoSplineDef3 = flexsurvspline(Surv(Time, Event) ~ 1, data = AtezoArm, knots = log(c(1,15,31)), scale = "hazard")

AtezoSplAIC = AIC(AtezoSpline1, AtezoSpline2, AtezoSpline3,AtezoSplineDef, AtezoSplineDef2,AtezoSplineDef3)

AtezoSplineModels = list(AtezoSpline1,
                    AtezoSpline2,
                    AtezoSpline3,
                    AtezoSplineDef,
                    AtezoSplineDef2,
                    AtezoSplineDef3)

atezo_spline_bic = NULL
for (i in 1:length(AtezoSplineModels)){
  bic = BIC(AtezoSplineModels[[i]])
  atezo_spline_bic = c(atezo_spline_bic,bic)
}


AtezoSplineNames = c("1 knot", 
                     "2 knots", 
                     "3 knots", 
                     "Knot in month 31",
                     "Knots in months 15, 31",
                     "Knots in months 1, 15, 31")

AtezoSplineBic = data.frame(AtezoSplineNames,atezo_spline_bic)
AtezoAllIC = data.frame(Model = c(DistributionNames, AtezoSplineNames), 
                           AIC = c(AtezoAic$AIC, AtezoSplAIC$AIC),
                           BIC = c(AtezoBic$all_bic, AtezoSplineBic$atezo_spline_bic))
plot(AtezoKM$time,
     AtezoKM$surv,
     col = c(0,0), ylim = c(0,1),
     xlab=("Time (in months)"), main="Atezolizumab Arm Spline Fit", ylab="Survival Probability")  
summary(AtezoSpline1)[[1]] %$% lines(time, est, col= colors[1], lwd = 2)
summary(AtezoSpline2)[[1]] %$% lines(time, est, col=colors[2], lwd = 2)
summary(AtezoSpline3)[[1]] %$% lines(time, est, col=colors[3], lwd = 2)
summary(AtezoSplineDef)[[1]] %$% lines(time, est, col=colors[4], lwd = 2)
summary(AtezoSplineDef2)[[1]] %$% lines(time, est, col=colors[5], lwd = 2)
summary(AtezoSplineDef3)[[1]] %$% lines(time, est, col=colors[6], lwd = 2)
lines(AtezoKM$time,
      AtezoKM$surv,col = c(1,0),lwd=2,t="s")  
legend("bottomleft",
       c("KM",AtezoSplineNames),
       col=c(1,colors[1:6]), 
       lty=1, 
       lwd=2,
       xpd    = TRUE,  
       cex    = 1,  
       bty    = "n")

# Spline model with specified knots in months 11, 19 and 29 is the best model
# among all 12 models fitted (standard parametric and spline based)

AtezoMinAIC = which.min(AtezoAllIC$AIC)
AtezoMinBIC = which.min(AtezoAllIC$BIC)

kable(AtezoAllIC, 
      col.names = c("Model","AIC", "BIC"), 
      align = "c",
      caption = "All Fits for Atezolizumab Arm: Information Criteria") %>% 
  kable_styling(bootstrap_options = "striped")%>%
  row_spec(AtezoMinAIC, bold = TRUE, col = "blue") %>%
  row_spec(AtezoMinBIC, bold = TRUE, col = "red")  %>%
add_footnote("Rows highlighted in blue: model with smallest AIC.", 
                                                                notation = "none")%>%
add_footnote("Rows highlighted in red: model with smallest BIC", 
               notation = "none")

# Remove the specified rows and reorder by AIC in ascending order
AtezoAllIC_filtered = AtezoAllIC %>%
  filter(!Model %in% c("Knots in months 15, 31", "Knots in months 1, 15, 31")) %>%
  arrange(AIC)

AtezoMinAIC = which.min(AtezoAllIC_filtered$AIC)
AtezoMinBIC = which.min(AtezoAllIC_filtered$BIC)
# Display the table
kable(AtezoAllIC_filtered, 
      col.names = c("Model","AIC", "BIC"), 
      align = "c",
      caption = "Filtered Fits for Atezolizumab Arm: Information Criteria") %>% 
  kable_styling(bootstrap_options = "striped") %>%
  row_spec(AtezoMinAIC, bold = TRUE, col = "blue")


plot(AtezoKM$time,
     AtezoKM$surv,
     col = c(0,0), ylim = c(0,1),
     xlab=("Time (in months)"), main="Atezolizumab Arm Top Fits", ylab="Survival Probability")
summary(AtezoGomp)[[1]] %$% lines(time, est, col=colors[3], lwd = 2)
summary(AtezoLogis)[[1]] %$% lines(time, est, col=colors[4], lwd = 2)
summary(AtezoNorm)[[1]] %$% lines(time, est, col=colors[5], lwd = 2)
summary(AtezoSplineDef)[[1]] %$% lines(time, est, col="firebrick", lwd = 2)
lines(AtezoKM$time,
      AtezoKM$surv,col = c(1,0),lwd=2,t="s")  
legend("bottomleft",
       c("KM",DistributionNames[c(5,4,3)],AtezoSplineNames[4]),
       col=c(1,colors[c(5,4,3)],c("firebrick")), 
       lty=1, 
       lwd=2,
       xpd    = TRUE,  
       cex    = 1,  
       bty    = "n")


# ---- Extrapolation Parametric Survival Curves

Years = 20
ExtrapHorizon = 12*Years 
times = seq(0,ExtrapHorizon, by = 1)

ExtrapAtSpline = summary(AtezoSplineDef, t = times)[[1]]

# Plot the extrapolated survival curves
plot(NULL, xlim = c(0,12*Years), ylim = c(0, 1), 
     xlab = "Time (Years)", ylab = "Survival Probability", 
     main = "Extrapolated Survival Curves for Atezotherapy Group",
     xaxt = "n")
axis(side=1, at=12*seq(from=0, to=ExtrapHorizon, by=1), labels=as.character(0:ExtrapHorizon))
lines(AtezoKM, lwd = 2)
lines(est~time, data = ExtrapAtezoGomp, type = "l", col = colors[3],lwd = 2)
lines(est~time, data = ExtrapAtezoLogis, type = "l", col = colors[4],lwd = 2)
lines(est~time, data = ExtrapAtezoNorm, type = "l", col = colors[5],lwd = 2)
lines(est~time, data = ExtrapAtSpline, type = "l", col = "firebrick", lwd = 2)
legend("topright", 
       legend = c("KM",DistributionNames[c(5,4,3)], AtezoSplineNames[4]),
       col = c(1,colors[c(5,4,3)],"firebrick"),
       lty = 1,
       lwd = 2,
       bty = "n")


# ---- Extrapolation Splines

ExtrapHorizon = 12*Years
times = seq(0,ExtrapHorizon, by = 1)

ExtrapAtSpline1 = summary(AtezoSpline1, t = times)[[1]]
ExtrapAtSpline2 = summary(AtezoSpline2, t = times)[[1]]
ExtrapAtSpline3 = summary(AtezoSpline3, t = times)[[1]]
ExtrapAtSplineDef = summary(AtezoSplineDef, t = times)[[1]]
ExtrapAtSplineDef2 = summary(AtezoSplineDef2, t = times)[[1]]
ExtrapAtSplineDef3 = summary(AtezoSplineDef3, t = times)[[1]]

# Plot the extrapolated survival curves
plot(NULL, xlim = c(0,12*Years), ylim = c(0, 1), 
     xlab = "Time (Years)", ylab = "Survival Probability", 
     main = "Extrapolated Survival Curves for Atezolizumab Group",
     xaxt = "n")
axis(side=1, at=12*seq(from=0, to=ExtrapHorizon, by=1), labels=as.character(0:ExtrapHorizon))
lines(AtezoKM, lwd = 2)
lines(est~time, data = ExtrapAtSpline1, type = "l", col = colors[1],lwd = 2)
lines(est~time, data = ExtrapAtSpline2, type = "l", col = colors[2],lwd = 2)
lines(est~time, data = ExtrapAtSpline3, type = "l", col = colors[3],lwd = 2)
lines(est~time, data = ExtrapAtSplineDef, type = "l", col = colors[4],lwd = 2)
lines(est~time, data = ExtrapAtSplineDef2, type = "l", col = colors[5],lwd = 2)
lines(est~time, data = ExtrapAtSplineDef3, type = "l", col = colors[6],lwd = 2)
legend("topright", 
       legend = c("KM",DistributionNames),
       col = c(1,colors[1:6]),
       lty = 1,
       lwd = 2,
       bty = "n")

# Mean survival times
AtezoSpline1Mean = mean_survspline1(gamma0=AtezoSpline1$res[1, 1],
                                    gamma1=AtezoSpline1$res[2, 1],
                                    gamma2=AtezoSpline1$res[3, 1],
                                    knots=AtezoSpline1$knots)        

AtezoSpline2Mean = mean_survspline2(gamma0=AtezoSpline2$res[1, 1],
                                    gamma1=AtezoSpline2$res[2, 1],
                                    gamma2=AtezoSpline2$res[3, 1],
                                    gamma3=AtezoSpline2$res[4, 1],
                                    knots=AtezoSpline2$knots)     

AtezoSpline3Mean = mean_survspline3(gamma0=AtezoSpline3$res[1, 1],
                                    gamma1=AtezoSpline3$res[2, 1],
                                    gamma2=AtezoSpline3$res[3, 1],
                                    gamma3=AtezoSpline3$res[4, 1],
                                    gamma4=AtezoSpline3$res[5, 1],
                                    knots=AtezoSpline3$knots)

AtezoSplineDefMean = mean_survspline1(gamma0=AtezoSplineDef$res[1, 1],
                                      gamma1=AtezoSplineDef$res[2, 1],
                                      gamma2=AtezoSplineDef$res[3, 1],
                                      knots=AtezoSplineDef$knots)

AtezoSplineDef2Mean = mean_survspline2(gamma0=AtezoSplineDef2$res[1, 1],
                                       gamma1=AtezoSplineDef2$res[2, 1],
                                       gamma2=AtezoSplineDef2$res[3, 1],
                                       gamma3=AtezoSplineDef2$res[4, 1],
                                       knots=AtezoSplineDef2$knots) 

AtezoSplineDef3Mean = mean_survspline3(gamma0=AtezoSplineDef3$res[1, 1],
                                       gamma1=AtezoSplineDef3$res[2, 1],
                                       gamma2=AtezoSplineDef3$res[3, 1],
                                       gamma3=AtezoSplineDef3$res[4, 1],
                                       gamma4=AtezoSplineDef3$res[5, 1],
                                       knots=AtezoSplineDef3$knots)


# ---- Summaries

# Summaries
# Survival probabilities
years         = c(5,10,15,20)
SummaryTimes  = 12*years            
MeansAtExp    = round(summary(AtezoExp, type='survival', t=SummaryTimes)[[1]][, 2], 3)
MeansAtWeib   = round(summary(AtezoWeib, type='survival', t=SummaryTimes)[[1]][, 2], 3)
MeansAtGomp   = round(summary(AtezoGomp, type='survival', t=SummaryTimes)[[1]][, 2], 3)
MeansAtLogis  = round(summary(AtezoLogis, type='survival', t=SummaryTimes)[[1]][, 2], 3)
MeansAtNorm   = round(summary(AtezoNorm, type='survival', t=SummaryTimes)[[1]][, 2], 3)
MeansAtGen    = round(summary(AtezoGen, type='survival', t=SummaryTimes)[[1]][, 2], 3)
MeansAtSpline = round(summary(AtezoSplineDef, type='survival', t=SummaryTimes)[[1]][, 2], 3)

AtezoMeansDF = data.frame(Year              = years, 
                          Exponential       = 100*MeansAtExp,
                          Weibull           = 100*MeansAtWeib,
                          Gompertz          = 100*MeansAtGomp,
                          Logistic          = 100*MeansAtLogis,
                          Normal            = 100*MeansAtNorm,
                          GenGamma          = 100*MeansAtGen,
                          `Spline month 31` = 100*MeansAtSpline
                          )

kable(AtezoMeansDF, col.names = c("Year",DistributionNames, "Spline w/ knot in 31"), align = "c") %>% 
  kable_styling(bootstrap_options = "striped")

# Median survival times
round(summary(AtezoExp, type='median')[[1]]['est'], 1)
round(summary(AtezoWeib, type='median')[[1]]['est'], 1)
round(summary(AtezoGomp, type='median')[[1]]['est'], 1)
round(summary(AtezoLogis, type='median')[[1]]['est'], 1)
round(summary(AtezoNorm, type='median')[[1]]['est'], 1)
round(summary(AtezoGen, type='median')[[1]]['est'], 1)
round(summary(AtezoSplineDef, type='median')[[1]]['est'], 1)

# ---- Survival Probabilities

SurvAtezoExp = summary(AtezoExp, t = time_points)[[1]]$est * 100
SurvAtezoWeib = summary(AtezoWeib, t = time_points)[[1]]$est * 100
SurvAtezoGomp = summary(AtezoGomp, t = time_points)[[1]]$est * 100
SurvAtezoLogis = summary(AtezoLogis, t = time_points)[[1]]$est * 100
SurvAtezoNorm = summary(AtezoNorm, t = time_points)[[1]]$est * 100
SurvAtezoGen = summary(AtezoGen, t = time_points)[[1]]$est * 100

# Combine results into a data frame
AtezoParamSurv = data.frame(
  Distribution = c("Exponential", "Weibull", "Gompertz", "Log-Logistic", "Log-Normal", "Generalized Gamma"),
  `1-Year Survival Probability (%)` = c(SurvAtezoExp[1], SurvAtezoWeib[1], SurvAtezoGomp[1], SurvAtezoLogis[1], SurvAtezoNorm[1], SurvAtezoGen[1]),
  `5-Year Survival Probability (%)` = c(SurvAtezoExp[2], SurvAtezoWeib[2], SurvAtezoGomp[2], SurvAtezoLogis[2], SurvAtezoNorm[2], SurvAtezoGen[2]),
  `10-Year Survival Probability (%)` = c(SurvAtezoExp[3], SurvAtezoWeib[3], SurvAtezoGomp[3], SurvAtezoLogis[3], SurvAtezoNorm[3], SurvAtezoGen[3])
)

kable(AtezoParamSurv, 
      col.names = c("Distribution", "1-Year Survival Probability (%)", "5-Year Survival Probability (%)", "10-Year Survival Probability (%)"),
      caption = "Predicted Survival Probabilities for Atezo Arm for Different Parametric Distributions",
      align = "c") %>% 
  kable_styling(bootstrap_options = "striped")

# Spline models

SurvAtezoSpline1 = summary(AtezoSpline1, t = time_points)[[1]]$est * 100
SurvAtezoSpline2 = summary(AtezoSpline2, t = time_points)[[1]]$est * 100
SurvAtezoSpline3 = summary(AtezoSpline3, t = time_points)[[1]]$est * 100
SurvAtezoSplineDef = summary(AtezoSplineDef, t = time_points)[[1]]$est * 100
SurvAtezoSplineDef2 = summary(AtezoSplineDef2, t = time_points)[[1]]$est * 100
SurvAtezoSplineDef3 = summary(AtezoSplineDef3, t = time_points)[[1]]$est * 100

# Combine results into a data frame for splines
AtezoSplineSurv = data.frame(
  Distribution = c("Spline 1", "Spline 2", "Spline 3", "Spline Def (log(31))", "Spline Def (log(15,31))", "Spline Def (log(1,15,31))"),
  `1-Year Survival Probability (%)` = c(SurvAtezoSpline1[1], SurvAtezoSpline2[1], SurvAtezoSpline3[1], SurvAtezoSplineDef[1], SurvAtezoSplineDef2[1], SurvAtezoSplineDef3[1]),
  `5-Year Survival Probability (%)` = c(SurvAtezoSpline1[2], SurvAtezoSpline2[2], SurvAtezoSpline3[2], SurvAtezoSplineDef[2], SurvAtezoSplineDef2[2], SurvAtezoSplineDef3[2]),
  `10-Year Survival Probability (%)` = c(SurvAtezoSpline1[3], SurvAtezoSpline2[3], SurvAtezoSpline3[3], SurvAtezoSplineDef[3], SurvAtezoSplineDef2[3], SurvAtezoSplineDef3[3])
)

kable(AtezoSplineSurv, 
      col.names = c("Spline Model", "1-Year Survival Probability (%)", "5-Year Survival Probability (%)", "10-Year Survival Probability (%)"),
      caption = "Predicted Survival Probabilities for Atezo Arm Spline Models",
      align = "c") %>% 
  kable_styling(bootstrap_options = "striped")

# Combine both parametric and spline results into one data frame
AtezoCombSurv = rbind(AtezoParamSurv, AtezoSplineSurv)

kable(AtezoCombSurv, 
      col.names = c("Model/Distribution", "1-Year Survival Probability (%)", "5-Year Survival Probability (%)", "10-Year Survival Probability (%)"),
      caption = "Predicted Survival Probabilities for Atezo Arm (Parametric and Spline Models)",
      align = "c") %>% 
  kable_styling(bootstrap_options = "striped")

