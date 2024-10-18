################################################################################
# Purpose: THESIS
#         
# Author : Kalatzi Marilena
# Date   : August 2024
################################################################################

rm(list = ls())

################################ LOAD PACKAGES #################################
library(survextrap)
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
library(readxl)
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

AtezoDf$Time = AtezoDf$Time /12

ChemoArm = AtezoDf[AtezoDf$Arm == "Chemotherapy",]
AtezoArm = AtezoDf[AtezoDf$Arm == "Atezolizumab",]

# ---- DEFINE M-SPLINE AND PRIORS ----

# ---- DEFAULT M-SPLINE SPECIFICATION
# Set up the default spline to allow changes in the hazard for up to 20 years.
# After that, we don't expect to see other changes, or we are not interested.

mspline = mspline_spec(Surv(Time, Event) ~ 1, data=AtezoDf, df=10, add_knots=20)

# ---- HAZARD SCALE PARAMETER (η) PRIOR
# for prior η to be specified, we think about the median age (64) of the trial
# "and we want to imply a prior mean survival of 25 years after diagnosis
# but with variance chosen so that mean survival could be as high as 100 years 
# after diagnosis"

prior_hscale = p_meansurv(median=25, upper=100, mspline=mspline)

# function prior_haz_const() translates a normal prior for η to the corresponding
# beliefs about survival. Using this function we can check that the lower limit 
# set is sensible for our data.

prior_haz_const(mspline, prior_hscale = prior_hscale)

# ---- HAZARD VARIABILITY PARAMETER (σ) PRIOR
# "prior for σ is chosen so that the highest hazard values over the 20 year horizon
# (i.e 90% quantile for some reason?) are expected to be about ρ = 2 times the 
# lowest values (10% quantile), with a vague 95% credible interval between 1 and 16"
# "prior_haz_sd() uses simulation to estimate the beliefs implied by a particular
# Gamma prior for σ, jointly with the prior specified for hazard scale. 
# Gamma(2,5) is chosen throug trial and error to achieve a value around 2 in the 
# quantity returned (hr represents the wanted result for ρ)"

set.seed(178)
prior_hsd = p_gamma(2, 4)
prior_haz_sd(mspline = mspline,
             prior_hsd = prior_hsd,
             prior_hscale = prior_hscale)

# ---- HAZARD RATIO FOR TREATMENT EFFECT

prior_loghr = p_hr(median=1, upper=50)
prior_hr(prior_loghr)

# ---- HAZARD RATIO VARIABILITY PARAMETER τ FOR NON-PROPORTIONAL HAZARDS MODEL

# governs the size of departure from proportional hazards, 
# i.e. the variability in the hazard ratio over time

set.seed(1)
prior_hrsd = p_gamma(2, 3)
prior_hr_sd(mspline = mspline,                              
            prior_hsd = prior_hsd,
            prior_loghr = prior_loghr,
            prior_hscale = prior_hscale,
            prior_hrsd = prior_hrsd,
            formula = ~ Arm,
            nonprop = ~ Arm,
            newdata = data.frame(Arm=1), 
            newdata0 = data.frame(Arm=0))

# ---- MODELING -----
# ---- Single Treatment Group ----

# options(mc.cores = 2)
# mod_con = survextrap(Surv(Time/12, Event) ~ 1, data=ChemoArm, mspline=mspline,
#                       prior_hscale=prior_hscale, prior_hsd = prior_hsd)
# 
# # mod_con5 = survextrap(Surv(Time/12, Event) ~ 1, data=ChemoArm, 
# #                        prior_hscale=prior_hscale, prior_hsd = prior_hsd)
# 
# surv_const = survival(mod_con5, tmax=20) %>% mutate(model="Constant hazard")
# surv_ipd = survival(mod_con, tmax=20) %>% mutate(model="Varying hazard")
# surv_single = rbind(surv_const, surv_ipd)
# haz_const = hazard(mod_con5, tmax=20) %>% mutate(model="Constant hazard")
# haz_ipd = hazard(mod_con, tmax=20) %>% mutate(model="Varying hazard")
# haz_single = rbind(haz_const, haz_ipd)
# 
# ps = ggplot(surv_single, aes(x=t, y=median, 
#                               group=model, col=model, fill=model)) + 
#   geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, colour=NA) +
#   geom_line(lwd=1.3) + 
#   geom_step(data=mod_con$km, aes(x=time, y=surv), lwd=1.3,
#             inherit.aes = FALSE) +
#   theme_minimal() + 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   xlab("Years after diagnosis") + 
#   ylab("Survival probability") +
#   scale_color_viridis(discrete=TRUE) + 
#   scale_fill_viridis(discrete=TRUE) + 
#   theme(legend.position = c(0.6, 0.8)) + labs(col=NULL, fill=NULL) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   geom_vline(xintercept = max(ChemoArm$Time/12)) +
#   geom_vline(xintercept = mod_con5$mspline$iknots, col="gray80", lty=2)
# 
# ph = ggplot(haz_single, aes(x=t, y=median, 
#                              group=model, col=model, fill=model)) + 
#   geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2, colour=NA) +
#   geom_line(lwd=1.3) + 
#   theme_minimal() + xlab("Years after diagnosis") + ylab("Hazard") +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   scale_color_viridis(discrete=TRUE) + 
#   scale_fill_viridis(discrete=TRUE) + 
#   geom_vline(xintercept = max(ChemoArm$Time/12)) +
#   geom_vline(xintercept = mod_con5$mspline$iknots, col="gray80", lty=2) + 
#   theme(legend.position = "none")
# 
# grid::grid.newpage()
# grid::grid.draw(cbind(ggplotGrob(ps), 
#                       ggplotGrob(ph)))

# ---- DETERMINING THE OPTIMAL NUMBER OF KNOTS
options(mc.cores = 2)
chains = 4; iter = 2000
dfs = 5:12
rescomp = as.data.frame(matrix(nrow=length(dfs), ncol=7))
for (i in seq_along(dfs)){
  msp = mspline_spec(Surv(Time, Event) ~ 1, data=ChemoArm, df=dfs[i])
  mod = survextrap(Surv(Time/12, Event) ~ 1, data=ChemoArm, mspline=msp,
                    prior_hscale=prior_hscale, prior_hsd = prior_hsd)
  rescompi = cbind(list(df = dfs[i]), 
                    list(looic = mod$loo$estimates["looic","Estimate"]),
                    rmst(mod,t=5))
  rescomp[i,] = rescompi
}
names(rescomp) = names(rescompi)
rescomp %>% 
  mutate(looic = round(looic, 1), 
         rmf = sprintf("%s (%s, %s)", round(median,2), 
                       round(`lower`, 2), round(`upper`, 2))) %>%
  select(df, looic, rmf) %>%
  knitr::kable(col.names = c("df", "LOOIC", "Restricted mean survival (5 years)"))


# ---- Compare different models

# ---- Proportional Hazards model
chains = 2
iter = 2000
TreatPh = survextrap(Surv(Time, Event) ~ Arm, data=AtezoDf, 
                     mspline=mspline, 
                     chains=chains, iter=iter,
                     prior_hscale=prior_hscale, prior_hsd = prior_hsd,
                     prior_loghr=prior_loghr)

# ---- Non-proportional hazards model

TreatNonPh = survextrap(Surv(Time, Event) ~ Arm, data=AtezoDf, 
                     mspline=mspline, 
                     chains=chains, iter=iter, 
                     nonprop = T, #non-proportionality
                     prior_hscale=prior_hscale, prior_hsd = prior_hsd,
                     prior_loghr=prior_loghr)

# ---- Fit arms separately 

ChemoMod = survextrap(Surv(Time, Event) ~ 1, data=ChemoArm, mspline=mspline,
                                   prior_hscale=prior_hscale, chains = chains,
                      iter = iter, prior_hsd = prior_hsd)

AtezoMod = survextrap(Surv(Time, Event) ~ 1, data=AtezoArm, mspline=mspline, 
                      chains=chains, iter=iter,
                      prior_hscale=prior_hscale, prior_hsd = prior_hsd)

# ---- VISUALIZATION
years = 20
times = years
SurvPh = survival(TreatPh, tmax=times) %>% mutate(Model="Proportional hazards")

SurvNonPh = survival(TreatNonPh, tmax=times) %>% mutate(Model="Non-proportional hazards")

SurvChemo = survival(ChemoMod, tmax=times) %>% 
  mutate(Model="Separate fits", Arm="Chemotherapy") %>% relocate(Arm)

SurvAtezo = survival(AtezoMod, tmax=times) %>% 
  mutate(Model="Separate fits", Arm="Atezolizumab") %>% relocate(Arm)

SurvModels = rbind(SurvPh, SurvNonPh, SurvChemo, SurvAtezo)

pp = ggplot(SurvModels, aes(x=t, y=median, col=Model, lty=Arm)) + 
  geom_step(data=TreatPh$km, aes(x=time, y=surv, lty=Arm), lwd=1.3,
            inherit.aes = FALSE) +
  geom_line(lwd=1.1) + 
  labs(lty="Treatment", col=NULL) +
  ggtitle("Model Comparison without external data")+
  theme_minimal() +
  theme(legend.position = c(0.6, 0.8), 
        legend.box = "vertical",  
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size = 16),
        axis.title = element_text(size=14),  # Increase the axis title text size
        axis.text = element_text(size=12),  # Increase the axis text size
        plot.title = element_text(size=16, face="bold"),
        legend.box.just = "left") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Years after diagnosis") + 
  ylab("Survival probability") +
  scale_color_viridis(discrete=TRUE) + 
  geom_vline(xintercept = max(ChemoArm$Time))

# ---- four different graphs

# Proportional Hazards Model
ggplot(SurvPh, aes(x=t, y=median, col=Arm)) + 
  geom_step(data=TreatPh$km, aes(x=time, y=surv, lty = Arm), lwd=1.3, inherit.aes = FALSE) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=Arm), alpha=0.2, linetype=0) +  # Shaded CI for Chemotherapy
  geom_line(lwd=1.1) + 
  labs(lty="Treatment", col="Arm") + 
  ggtitle("Proportional Hazards Model") +
  theme_minimal() +
  theme(legend.position = c(0.6, 0.8), 
        legend.box = "vertical",  
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(size=14),  # Increase the axis title text size
        axis.text = element_text(size=12),  # Increase the axis text size
        plot.title = element_text(size=16, face="bold"),
        legend.box.just = "left") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Years after diagnosis") + 
  ylab("Survival probability") +
  scale_color_manual(values=c("Chemotherapy" = "#FDE725", "Atezolizumab" = "#440154")) +  # Yellow for Chemo, Purple for Atezo
  scale_fill_manual(values=c("Chemotherapy" = "#FDE725", "Atezolizumab" = "#440154")) +  # Match CI fill colors
  geom_vline(xintercept = max(ChemoArm$Time))

# Non-Proportional Hazards Model


# Proportional Hazards Model
ggplot(SurvNonPh, aes(x=t, y=median, col=Arm)) + 
  geom_step(data=TreatNonPh$km, aes(x=time, y=surv, lty = Arm), lwd=1.3, inherit.aes = FALSE) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=Arm), alpha=0.2, linetype=0) +  # Shaded CI for Chemotherapy
  geom_line(lwd=1.1) + 
  labs(lty="Treatment", col="Arm") + 
  ggtitle("Non Proportional Hazards Model") +
  theme_minimal() +
  theme(legend.position = c(0.6, 0.8), 
        legend.box = "vertical",  
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(size=14),  # Increase the axis title text size
        axis.text = element_text(size=12),  # Increase the axis text size
        plot.title = element_text(size=16, face="bold"),
        legend.box.just = "left") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Years after diagnosis") + 
  ylab("Survival probability") +
  scale_color_manual(values=c("Chemotherapy" = "#FDE725", "Atezolizumab" = "#440154")) +  # Yellow for Chemo, Purple for Atezo
  scale_fill_manual(values=c("Chemotherapy" = "#FDE725", "Atezolizumab" = "#440154")) +  # Match CI fill colors
  geom_vline(xintercept = max(ChemoArm$Time))


# Proportional and Non-proportional

# Combine the datasets from Proportional and Non-Proportional Hazards models
SurvCombined <- rbind(SurvPh %>% mutate(Model = "Proportional Hazards"),
                      SurvNonPh %>% mutate(Model = "Non-Proportional Hazards"))

# Combined Plot: Proportional vs Non-Proportional Hazards
ggplot(SurvCombined, aes(x=t, y=median, col=Model, lty=Arm)) + 
  geom_step(data=TreatPh$km, aes(x=time, y=surv, lty = Arm), lwd=1.3, inherit.aes = FALSE) +
  geom_line(lwd=1.1) + 
  labs(lty="Treatment", col="Model") + 
  ggtitle("Proportional vs Non-Proportional Hazards Model") +
  theme_minimal() +
  theme(legend.position = c(0.6, 0.8), 
        legend.box = "vertical",  
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(size=14),  # Increase the axis title text size
        axis.text = element_text(size=12),  # Increase the axis text size
        plot.title = element_text(size=16, face="bold"),
        legend.box.just = "left") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Years after diagnosis") + 
  ylab("Survival probability") +
  scale_color_viridis(discrete=TRUE) + 
  geom_vline(xintercept = max(ChemoArm$Time))


# Chemotherapy Plot with Confidence Intervals and Yellow Color
ggplot(SurvChemo, aes(x=t, y=median, col=Arm)) + 
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=Arm), alpha=0.2, linetype=0) +  # Shaded CI for Chemotherapy
  geom_step(data=ChemoMod$km, aes(x=time, y=surv), lwd=1.3, inherit.aes = FALSE, lty=2) +  # Step curve for KM estimates
  geom_line(lwd=1.1) + 
  labs(lty="Treatment", col="Arm", fill="Arm") +  # Add fill legend for CI
  ggtitle("Chemotherapy Arm") +
  theme_minimal() +
  theme(legend.position = c(0.6, 0.8), 
        legend.box = "vertical",  
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12),
        axis.title = element_text(size=14),  # Increase axis title size
        axis.text = element_text(size=12),  # Increase axis text size
        plot.title = element_text(size=16, face="bold"),
        legend.box.just = "left") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Years after diagnosis") + 
  ylab("Survival probability") +
  scale_color_manual(values=c("Chemotherapy" = "#FDE725")) +  # Assign yellow color to Chemotherapy
  scale_fill_manual(values=c("Chemotherapy" = "#FDE725")) +  # Assign yellow color to CI
  geom_vline(xintercept = max(ChemoArm$Time))

# Atezolizumab Plot with Confidence Intervals and Purple Color
ggplot(SurvAtezo, aes(x=t, y=median, col=Arm)) + 
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=Arm), alpha=0.2, linetype=0) +  # Shaded CI for Atezolizumab
  geom_step(data=AtezoMod$km, aes(x=time, y=surv), lwd=1.3, inherit.aes = FALSE, lty = 2) +  # Step curve for KM estimates
  geom_line(lwd=1.1) + 
  labs(lty="Treatment", col="Arm", fill="Arm") +  # Add fill legend for CI
  ggtitle("Atezolizumab Arm") +
  theme_minimal() +
  theme(legend.position = c(0.6, 0.8), 
        legend.box = "vertical",  
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12),
        axis.title = element_text(size=14),  # Increase axis title size
        axis.text = element_text(size=12),  # Increase axis text size
        plot.title = element_text(size=16, face="bold"),
        legend.box.just = "left") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Years after diagnosis") + 
  ylab("Survival probability") +
  scale_color_manual(values=c("Atezolizumab" = "#440154")) +  # Assign purple color to Atezolizumab
  scale_fill_manual(values=c("Atezolizumab" = "#440154")) +  # Assign purple color to CI
  geom_vline(xintercept = max(ChemoArm$Time))


# ---- Model Summaries

ModelComp <- function(mod){
  rm <- rmst(mod, niter=1000, newdata=data.frame(Arm="Chemotherapy"), t=c(3,10,20)) %>% 
    rename(rm_med="median", rm_lower="lower", rm_upper="upper") %>% 
    select(t, rm_med, rm_lower, rm_upper)
  rm2 <- rmst(mod, niter=1000, newdata=data.frame(Arm="Atezolizumab"), t=c(3,10,20)) %>% 
    rename(rm_med2="median", rm_lower2="lower", rm_upper2="upper") %>% 
    select(rm_med2, rm_lower2, rm_upper2)
  cbind(name=deparse(substitute(mod)), rm, rm2, 
        list(looic = mod$loo$estimates["looic","Estimate"]))
}

RmChemo <- rmst(ChemoMod,t=c(3,10,20)) %>%
  rename(rm_med ="median", rm_lower ="lower", rm_upper ="upper") %>% 
  select(t, rm_med, rm_lower, rm_upper)

RmAtezo <- rmst(AtezoMod, t=c(3,10,20)) %>%
  rename(rm_med2="median", rm_lower2="lower", rm_upper2="upper") %>% 
  select(rm_med2, rm_lower2, rm_upper2)

# irmst_sep <- rmst(AtezoMod, t=c(36,120), sample=TRUE) - rmst(ChemoMod, t=c(36,120), sample=TRUE) 
# 
# irmst_sep <- irmst_sep %>%
#   posterior::summarise_draws(median,
#                              ~quantile(.x, probs=c(0.025, 0.975))) %>%
#   rename(t=variable, ir_med="median",ir_lower="2.5%", ir_upper="97.5%") %>% 
#   select(ir_med, ir_lower, ir_upper)

SepArmLooic <- ChemoMod$loo$estimates["looic","Estimate"] +
  AtezoMod$loo$estimates["looic","Estimate"]
SepArmComp <- cbind(list(name="mod_sep"), RmChemo, RmAtezo, list(looic=SepArmLooic))

TreatComp <- rbind(ModelComp(TreatPh), ModelComp(TreatNonPh), SepArmComp)
TreatCompf <- TreatComp %>%
  arrange(t) %>%
  mutate(rm = sprintf("%s (%s,%s)", round(rm_med,2),
                      round(rm_lower,2), round(rm_upper,2)),
         rm2 = sprintf("%s (%s,%s)", round(rm_med2,2),
                      round(rm_lower2,2), round(rm_upper2,2)),
         looic = round(looic),
         model = forcats::fct_recode(name, 
                                     "(a) Proportional hazards"="TreatPh",
                                     "(b) Non-proportional hazards"="TreatNonPh",
                                     "(c) Separate arms"="mod_sep"
         )) %>%
  select(model, t, rm, rm2, looic) %>%
  tibble::remove_rownames()

knitr::kable(TreatCompf, col.names=c("Model","Time horizon (years)",
                                   "RMST Chemotherapy (CI)",
                                   "RMST Atezolizumab (CI)",
                                   "LOOIC"))%>% 
  kable_styling(bootstrap_options = "striped")


kable(TreatCompf, 
      align = "c", 
      col.names = c("Model", "Time horizon (years)", "RMS Chemotherapy (CI)", "RMR Atezolizumab (CI)", "LOOIC"),
      caption = "RMST for the 3 fitted models, in 3 different timepoints") %>%
  kable_styling(bootstrap_options = "striped") %>%
  row_spec(3, extra_css = "border-bottom: 1px solid blue;") %>%
  row_spec(6, extra_css = "border-bottom: 1px solid blue;")


# ---- Survival Probabilities
# Define time points for 1-year, 5-year, and 10-year survival probabilities
# comparison_times <- c(12, 60, 120)  # Use comparison_times instead of times
comparison_times = c(1,5,10)
# Filter for these time points and include CIs
SurvPh_1_5_10 <- SurvPh %>%
  filter(t %in% comparison_times) %>%
  select(t, median, lower, upper) %>%
  mutate(Model="Proportional hazards")

SurvNonPh_1_5_10 <- SurvNonPh %>%
  filter(t %in% comparison_times) %>%
  select(t, median, lower, upper) %>%
  mutate(Model="Non-proportional hazards")

SurvChemo_1_5_10 <- SurvChemo %>%
  filter(t %in% comparison_times) %>%
  select(t, median, lower, upper) %>%
  mutate(Model="Chemotherapy")

SurvAtezo_1_5_10 <- SurvAtezo %>%
  filter(t %in% comparison_times) %>%
  select(t, median, lower, upper) %>%
  mutate(Model="Atezolizumab")

# Combine all models into one dataframe
SurvModels_1_5_10 <- rbind(SurvPh_1_5_10, SurvNonPh_1_5_10, SurvChemo_1_5_10, SurvAtezo_1_5_10)

# Summarize the data to handle duplicates by averaging the values for median, lower, and upper
SurvModels_1_5_10_summarized <- SurvModels_1_5_10 %>%
  group_by(Model, t) %>%
  summarise(
    median = mean(median, na.rm = TRUE),
    lower = mean(lower, na.rm = TRUE),
    upper = mean(upper, na.rm = TRUE)
  ) %>%
  ungroup()


##################
# Convert the tibble to a data frame to use reshape()
SurvModels_1_5_10_summarized_df <- as.data.frame(SurvModels_1_5_10_summarized)

# Reshape from long to wide format using reshape()
SurvTable_CI <- reshape(SurvModels_1_5_10_summarized_df,
                        timevar = "t",  # Use the 't' column as the time variable
                        idvar = "Model",  # Use 'Model' as the identifier variable
                        direction = "wide")  # Wide format

SurvTable_CI[, -1] <- SurvTable_CI[, -1] * 100

colnames(SurvTable_CI) <- c("Model", 
                            "1-Year Median (%)", "1-Year Lower CI (%)", "1-Year Upper CI (%)",
                            "5-Year Median (%)", "5-Year Lower CI (%)", "5-Year Upper CI (%)",
                            "10-Year Median (%)", "10-Year Lower CI (%)", "10-Year Upper CI (%)")


SurvTable_CI <- SurvTable_CI %>%
  mutate(`1-Year Median (CI)` = sprintf("%.2f (%.2f, %.2f)", `1-Year Median (%)`, `1-Year Lower CI (%)`, `1-Year Upper CI (%)`),
         `5-Year Median (CI)` = sprintf("%.2f (%.2f, %.2f)", `5-Year Median (%)`, `5-Year Lower CI (%)`, `5-Year Upper CI (%)`),
         `10-Year Median (CI)` = sprintf("%.2f (%.2f, %.2f)", `10-Year Median (%)`, `10-Year Lower CI (%)`, `10-Year Upper CI (%)`))

SurvTable_CI_cleaned <- SurvTable_CI %>%
  select(Model, `1-Year Median (CI)`, `5-Year Median (CI)`, `10-Year Median (CI)`)


kable(SurvTable_CI_cleaned, 
      col.names = c("Model", 
                    "1-Year Median (CI)", "5-Year Median (CI)", "10-Year Median (CI)"),
      caption = "Predicted Survival Probabilities with Confidence Intervals (in %)") %>%
  kable_styling(bootstrap_options = "striped", full_width = F) %>%
  row_spec(0, bold = TRUE)  # To emphasize the header (optional)

# ---- External Registry Data for Chemotherapy group

ExternalPath = "C:/Users/kalat/OneDrive/Desktop/Thesis/Lung Cancer Statistics/DeathsByAgeUK.xlsx"
External = read_excel(ExternalPath, sheet = 4)

ExternalDf = as.data.frame(External)
colnames(ExternalDf) = c("start","stop","r","n","treat")

ExternalDf$treat = factor(ExternalDf$treat)
ExternalDf$n = round(as.integer(ExternalDf$n)/10000)
ExternalDf$r = round(as.integer(ExternalDf$r)/10000)
# ExternalDf$start = ExternalDf$start*12
# ExternalDf$stop = ExternalDf$stop*12


ExternalDf <- ExternalDf %>%
  mutate(
    haz = -log(r / n)/5,  # Hazard estimate
    haz_lower = -log(qbeta(0.975, r, n - r))/5,  # Lower bound of hazard
    haz_upper = -log(qbeta(0.025, r, n - r))/5   # Upper bound of hazard
  )

# View the updated ExternalDf with hazard and CIs
ExternalDf
mspline_registry = mspline_spec(Surv(Time, Event) ~ 1, data=AtezoDf, df=10,  add_knots=c(8,13, 25))
# mspline = list(add_knots=c(4, 10, 25))
# cetux_months = cetux_seer
# cetux_months$start = cetux_months$start*12
# cetux_months$stop  = cetux_months$stop*12
tail(AtezoDf)

ChemoModRegistry = survextrap(Surv(Time, Event) ~ 1, 
                              data=ChemoArm, 
                              external = ExternalDf,
                              mspline=mspline_registry,
                              prior_hscale=prior_hscale, 
                              prior_hsd = prior_hsd,
                              chains = 2)
ChemoModRegistry

# test external data
ChemoModRegistry = survextrap(Surv(Time, Event) ~ 1, 
                              data=ChemoArm, 
                              external = cetux_seer,
                              mspline=mspline_registry,
                              prior_hscale=prior_hscale, 
                              prior_hsd = prior_hsd,
                              chains = 2)
# 

plot(ChemoModRegistry)

times = 25

# Extract survival data for Kaplan-Meier curve
km_data <- ChemoMod$km
# Extract survival data from both models
SurvChemo <- survival(ChemoMod) %>%
  mutate(Model = "Without external data")
SurvChemoRegistry <- survival(ChemoModRegistry, tmax = times) %>%
  mutate(Model = "With external data")

# Extract hazard data from both models
HazChemo <- hazard(ChemoMod) %>%
  mutate(Model = "Without external data")
tail(HazChemo)
HazChemoRegistry <- hazard(ChemoModRegistry, tmax = times) %>%
  mutate(Model = "With external data")

# Combine survival and hazard data into respective data frames
SurvData <- rbind(SurvChemo, SurvChemoRegistry)
HazData <- rbind(HazChemo, HazChemoRegistry)

# Plot survival comparison with CIs and KM line
ggplot(SurvData, aes(x = t, y = median, col = Model, fill = Model)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) + 
  geom_step(data = km_data, aes(x = time, y = surv), color = "black", linetype = "dashed", size = 1.2, inherit.aes = FALSE) +
  labs(title = "Survival Curves Comparison with Kaplan-Meier",
       x = "Time (Years)", 
       y = "Survival Probability") +
  theme_minimal() + 
  theme(legend.position = c(0.6, 0.8), 
        legend.box = "vertical",  
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(size=14),  # Increase the axis title text size
        axis.text = element_text(size=12),  # Increase the axis text size
        plot.title = element_text(size=16, face="bold")) +
  scale_x_continuous(labels = function(x) round(x, 0)) +  # Convert months to years
  scale_color_viridis(discrete = TRUE) + 
  scale_fill_viridis(discrete = TRUE) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Plot hazard comparison with CIs
ggplot(HazData, aes(x = t, y = median, col = Model, fill = Model)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) + 
  labs(title = "Hazard Curves Comparison",
       x = "Time (Years)", 
       y = "Hazard Rate") +
  theme_minimal() + 
  theme(legend.position = c(0.6, 0.8), 
        legend.box = "vertical",  
        legend.title = element_text(size=14, face="bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(size=14),  # Increase the axis title text size
        axis.text = element_text(size=12),  # Increase the axis text size
        plot.title = element_text(size=16, face="bold")) +
  scale_color_viridis(discrete = TRUE) + 
  scale_fill_viridis(discrete = TRUE) + 
  geom_vline(xintercept = max(ChemoArm$Time)) +  # KM cutoff for ChemoMod
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(HazData, aes(x = t, y = median, col = Model, fill = Model)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
  geom_line(size = 1.2) +
  geom_step(data = ExternalDf, aes(x = start, y = haz), col = "gray30", inherit.aes = FALSE) +
  geom_step(data = ExternalDf, aes(x = start, y = haz_lower), col = "gray70", linetype = 2, inherit.aes = FALSE) +
  geom_step(data = ExternalDf, aes(x = start, y = haz_upper), col = "gray70", linetype = 2, inherit.aes = FALSE) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Years after diagnosis") + 
  ylab("Hazard") +
  scale_color_viridis(discrete = TRUE) + 
  scale_fill_viridis(discrete = TRUE) + 
  geom_vline(xintercept = max(ChemoArm$Time)) +  # KM cutoff for ChemoMod
  theme(legend.position = c(0.4, 0.9), legend.background = element_blank()) +
  labs(col = NULL, fill = NULL)

# test with cetux_seer

ggplot(HazData, aes(x = t, y = median, col = Model, fill = Model)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
  geom_line(size = 1.2) +
  geom_step(data = cetux_seer, aes(x = start, y = haz), col = "gray30", inherit.aes = FALSE) +
  geom_step(data = cetux_seer, aes(x = start, y = haz_lower), col = "gray70", linetype = 2, inherit.aes = FALSE) +
  geom_step(data = cetux_seer, aes(x = start, y = haz_upper), col = "gray70", linetype = 2, inherit.aes = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Years after diagnosis") + 
  ylab("Hazard") +
  scale_color_viridis(discrete = TRUE) + 
  scale_fill_viridis(discrete = TRUE) + 
  geom_vline(xintercept = max(ChemoArm$Time)) +  # KM cutoff for ChemoMod
  theme(legend.position = c(0.4, 0.9), legend.background = element_blank()) +
  labs(col = NULL, fill = NULL)

cetux_seer

#################################
