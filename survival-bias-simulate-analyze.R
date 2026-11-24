#==============================================================================
# R Code for: 
# "Aging by the numbers: Truncation by death is not just missing data"
# =============================================================================
# Simulates data for baseline plus 3 follow-up visits
# with truncation by death
# X = baseline continuous covariate (standardized age)
# Y0 = pre-exposure cognition (3MS)
# A = exposure (= 1 if CKD, = 0 if no CKD)
# W = post-exposure continuous covariate (change in CVD risk)
# Y1, Y2, Y3 = cognition outcome (3MS) at follow-up 1, 2, 3
# S1, S2, S3 = survival status at follow-up 1, 2, 3
# Y1_OBS, Y2_OBS, Y3_OBS = NA if truncated by death;  = Y1, Y2, Y3 otherwise
#
# Data analyzed using:
# Independence Estimating Equations (incomplete and full cohort) -
#     [Goal 1: Generalize to survivors]
# Principal Stratification (Survivor Average Causal Effects) and 
#     Weighted Estimating Equations (incomplete and full cohort) -
#     [Goal 2: Infer causality]
# Mixed-effects models (incomplete and full cohort) -
#     [Goal 3: Estimate “immortal cohort” association]
#
# Deliverables:
# Simulated dataset (Truncation_By_Death.csv)
# Longitudinal estimates (Figure2.pdf)
#=============================================================================

# Load Required Libraries
# install.packages(c("survey", "plyr", "dplyr", "nlme", "multcomp", "ggplot2", "gridExtra"))

library(survey)    # WEE to handle truncation by death using inverse probability weighting
library(plyr) 	   # Data manipulation and visualization
library(dplyr) 	   # Data manipulation and visualization
library(nlme)      # Mixed-effects models (Immortal Cohort approach)
library(multcomp)  # Linear combinations of coefficients
library(ggplot2)   # Figure 2
library(gridExtra) # Figure 2 arrangement

set.seed(06012026)
n_samples <- 5000

# ==============================================================================
# Generate Simulated Dataset (CKD Exposure & Cognition [3MS] Outcome)
# ==============================================================================

SD_Y0 <- 1
SD_Y1 <- 1
SD_Y2 <- 3
SD_Y3 <- 3
MU_Y0 <- 80

sim_data <- tibble(
  id = 1:n_samples,

  # Baseline covariate (standardized age):  X ~ normal(0, sd=1)
  X = rnorm(n_samples, mean = 0, sd = 1),

  # Pre-exposure outcome (3MS):  Y0 ~ normal(MU_Y0, sd=SD_Y0)
  Y0_temp = MU_Y0 + rnorm(n_samples, mean = 0, sd = SD_Y0),  
  Y0 = pmax(0, pmin(Y0_temp, 100)),

) %>%
  rowwise() %>%
  mutate(
    # Exposure (1 = yes CKD, 0 = no CKD): logit(E[A | X, Y0]) =  a0 + a*1X + a2*scale(Y0) + a3*X*scale(Y0)
    P_A = plogis(0 + 0.5*X - 0.25*(Y0-80)/SD_Y0 + 0.75*X*(Y0-80)/SD_Y0),
    A = rbinom(1, 1, P_A), 

    # Post-exposure outcome-death confounder, W 
    #E[W | A, X, Y0] =  h0 + h1X + h2*scale(Y0) + h3*A
    W = 0 + 0.75*X - 0.25*(Y0-80)/SD_Y0 + 1*A +  rnorm(1, mean = 0, sd = 2),

    # Mortality status (1 = alive for year 1, 0 = dead)
    # logit(E[S1 | X, Y0, A, W]) =  g01 + g11*X + g21*Y0 + g31*W + g41*A 
    P_S1 = plogis(1.2 - 0.2*X + 0.2*(Y0-80)/SD_Y0 - 1*W - 1*A),
    S1 = rbinom(1, 1, P_S1), 

    # Outcome at year 1 (also, outcome-death confounder for year 2)
    # E[Y1 | Y0, X, A, W] = b01 + b11*Y0 + b21*X + b31*A + b41*W + b51*A*X + b61*W*A
    Y1_temp = Y0 - 3 - 1*X - 4*A - 4*W - 1*A*X + 3*W*A  +  rnorm(1, mean = 0, sd = SD_Y1),
    Y1 = pmax(0, pmin(Y1_temp, 100)),

    # Mortality status (1 = alive for year 2, 0 = dead)
    # logit(E[S2 | X, Y1, A, W]) =  g02 + g12*X + g22*Y1 + g32*W + g42*A 
    P_S2 = plogis(1.3 - 0.2*X + 0.1*(Y1-75)/SD_Y1 - 1*W - 1*A),
    S2 = rbinom(1, 1, P_S2)*S1, 

    # Outcome at year 2 (also, outcome-death confounder for year 3)
    # E[Y2 | Y1, X, A, W] = b02 + b12*Y1 + b22*X + b32*A + b42*W + b52*A*X + b62*W*A
    Y2_temp = Y1 - 4 - 1*X - 4*A - 3*W - 1*A*X + 2*W*A  +  rnorm(1, mean = 0, sd = SD_Y2),
    Y2 = pmax(0, pmin(Y2_temp, 100)),

    # Mortality status (1 = alive for visit 3, 0 = dead)
    # logit(E[S3 | X, Y2, A, W]) =  g02 + g12*X + g22*Y2 + g32*W + g42*A 
    P_S3 = plogis(1.1 - 0.2*X + 0.1*(Y2-68)/SD_Y2 - 1*W - 1*A),
    S3 = rbinom(1, 1, P_S3)*S2*S1,

    # Outcome at visit 3 (final visit)
    # E[Y3 | Y2, X, A, W] = b03 + b13*Y1 + b23*X + b33*A + b43*W + b53*A*X + b63*W*A
    Y3_temp = Y2 - 5 - 1*X - 4*A - 2*W - 1*A*X + 1*W*A  +  rnorm(1, mean = 0, sd = SD_Y3),
    Y3 = pmax(0, pmin(Y3_temp, 100)),

    # Observed outcome: Undefined (NA) if the individual dies before visit
    Y0_OBS = Y0,
    Y1_OBS = ifelse(S1 == 1, Y1, NA_real_),
    Y2_OBS = ifelse(S2 == 1, Y2, NA_real_),
    Y3_OBS = ifelse(S3 == 1, Y3, NA_real_)
  ) %>%
  ungroup()  %>%
  dplyr::select(-c(Y0_temp, Y1_temp, Y2_temp, Y3_temp, P_A, P_S1, P_S2, P_S3))


# save simulated data
write.csv(sim_data, "Truncation_By_Death.csv", row.names=F)

######################################################
# uncomment to read and analye dataset used in paper #
# sim_data <- read.csv("Truncation_By_Death.csv")    #
######################################################

# for use in models that use IPW to control for confounding 
sim_data$PSCORE <- predict(glm(A ~ X*scale(Y0), family = "binomial", data = sim_data), type = "response")
sim_data$IPW_A <- mean(1-sim_data$A)*(1 - sim_data$A)/(1 - sim_data$PSCORE) +  mean(sim_data$A)*sim_data$A/sim_data$PSCORE 

# Reshape observed data into a long-format dataset for standard longitudinal (mixed) models
sim_obsdata_long <- sim_data %>%
  tidyr::pivot_longer(
    cols = c(Y0_OBS, Y1_OBS, Y2_OBS, Y3_OBS), 
    names_to = "TIME", 
    values_to = "Y"
  ) %>%
  mutate(YEAR = 1*(TIME == "Y1_OBS") + 2*(TIME == "Y2_OBS") + 3*(TIME == "Y3_OBS"),
         S = ifelse(is.na(Y), 0, 1)) 
 
# Reshape complete (immortal) data into a long-format dataset for standard longitudinal (mixed) models
sim_compdata_long <- sim_data %>%
  tidyr::pivot_longer(
    cols = c(Y0_OBS, Y1, Y2, Y3), 
    names_to = "TIME", 
    values_to = "Y"
  ) %>%
  mutate(YEAR = 1*(TIME == "Y1") + 2*(TIME == "Y2") + 3*(TIME == "Y3"),
         S = ifelse(is.na(Y), 0, 1)) 

# mean (SD) baseline cognition score
mean_Y0 <- mean(sim_data$Y0)
sd_Y0 <- sd(sim_data$Y0)

# number with CKD
table(sim_data$A)

# proportion with CKD
prop.table(table(sim_data$A))

# number of CKD status and year 1 survival status
table(sim_data$A, sim_data$S1)

# proportion of year 1 survival status by CKD status
prop.table(table(sim_data$A, sim_data$S1),margin=1)

# number of CKD status and year 2 survival status
table(sim_data$A, sim_data$S2)

# proportion of year 2 survival status by CKD status
prop.table(table(sim_data$A, sim_data$S2),margin=1)

# number of CKD status and year 3 survival status
table(sim_data$A, sim_data$S3)

# proportion of year 3 survival status by CKD status
prop.table(table(sim_data$A, sim_data$S3),margin=1)

# change in CVD risk by CKD
aggregate(W ~ A, data = sim_data, FUN = mean, na.rm = TRUE)
aggregate(W ~ A, data = sim_data, FUN = sd, na.rm = TRUE)

# ==============================================================================
# Goal 1: Conditioning on Survival (Independence Estimating Equations)
# ==============================================================================
# An analysis conditioned on being alive.

# IEE model
sim_obsdata_long$YEAR.f <- factor(sim_obsdata_long$YEAR)
IEE_des <- svydesign(id=~id, weights=~1, data=sim_obsdata_long) 
IEE_model <- svyglm(Y ~ A*YEAR.f + Y0 + X, design=subset(IEE_des,YEAR>0))
cat("\n--- Goal 1:Conditioning on Survival IEE Results ---\n")
print(summary(IEE_model)$coefficients)

# adjusted mean difference at Y1
IEE_lincom_MD_Y1 <- glht(IEE_model, linfct = c("A  = 0"))
confint(IEE_lincom_MD_Y1)$confint

# adjusted mean difference at Y2
IEE_lincom_MD_Y2 <- glht(IEE_model, linfct = c("A + A:YEAR.f2 = 0"))
confint(IEE_lincom_MD_Y2)$confint

# adjusted mean difference at Y3
IEE_lincom_MD_Y3 <- glht(IEE_model, linfct = c("A + A:YEAR.f3 = 0"))
confint(IEE_lincom_MD_Y3)$confint

# dataset to derive adjusted means
grid <- expand.grid(
  A = c(0,1),
  YEAR.f = levels(factor(sim_obsdata_long$YEAR[sim_obsdata_long$YEAR>0]))
)

grid$X <- 0
grid$Y0 <- mean_Y0
Y0.ests <- t.test(sim_data$Y0)

K_matrix <- diag(nrow(grid)) # Identity matrix to evaluate each fitted value
rownames(K_matrix) <- paste(grid$A, grid$YEAR.f, sep = "_")

X_matrix_adjust <- model.matrix(~  A*YEAR.f + Y0 + X, grid)
IEE_glht <- glht(IEE_model, linfct = K_matrix %*% X_matrix_adjust)
IEE_ci_results <- confint(IEE_glht, calpha = univariate_calpha())$confint

# dataset of adjusted means and 95% CIs
IEE_results_df <- data.frame(
  A = grid$A,
  YEAR.f = as.numeric(grid$YEAR.f),
  Estimate = IEE_ci_results[,"Estimate"],
  Lower_CI = IEE_ci_results[,"lwr"],
  Upper_CI = IEE_ci_results[,"upr"]
)


IEE_results_df <- rbind( c(0,0, mean_Y0, Y0.ests$conf.int[1], Y0.ests$conf.int[2]),
                         c(1,0, mean_Y0, Y0.ests$conf.int[1], Y0.ests$conf.int[2]),
                          IEE_results_df
                       )

IEE_results_df$A <- factor(IEE_results_df$A,labels=c("No","Yes"))


# ==============================================================================
# Goal 1 Complete data comparison:  (Independence Estimating Equations)
# ==============================================================================

# IEE model of complete data
sim_compdata_long$YEAR.f <- factor(sim_compdata_long$YEAR)
IEE_complete_des <- svydesign(id=~id, weights=~1, data=sim_compdata_long) 
IEE_complete_model <- svyglm(Y ~ A*YEAR.f + Y0 + X, design=subset(IEE_complete_des,YEAR>0))
cat("\n--- IEE with complete data ---\n")
print(summary(IEE_complete_model)$coefficients)

# adjusted mean difference at Y1
IEE_complete_lincom_MD_Y1 <- glht(IEE_complete_model, linfct = c("A  = 0"))
confint(IEE_complete_lincom_MD_Y1)$confint

# adjusted mean difference at Y2
IEE_complete_lincom_MD_Y2 <- glht(IEE_complete_model, linfct = c("A + A:YEAR.f2 = 0"))
confint(IEE_complete_lincom_MD_Y2)$confint

# adjusted mean difference at Y3
IEE_complete_lincom_MD_Y3 <- glht(IEE_complete_model, linfct = c("A + A:YEAR.f3 = 0"))
confint(IEE_complete_lincom_MD_Y3)$confint

# adjusted means from complete data
IEE_complete_glht <- glht(IEE_complete_model, linfct = K_matrix %*% X_matrix_adjust)
IEE_complete_ci_results <- confint(IEE_complete_glht, calpha = univariate_calpha())$confint

# dataset of adjusted means and 95% CIs from complete data 
IEE_complete_results_df <- data.frame(
  A = grid$A,
  YEAR.f = as.numeric(grid$YEAR.f),
  Estimate = IEE_complete_ci_results[,"Estimate"],
  Lower_CI = IEE_complete_ci_results[,"lwr"],
  Upper_CI = IEE_complete_ci_results[,"upr"]
)

IEE_complete_results_df <- rbind( c(0,0, mean_Y0, Y0.ests$conf.int[1], Y0.ests$conf.int[2]),
                         c(1,0, mean_Y0, Y0.ests$conf.int[1], Y0.ests$conf.int[2]),
                          IEE_complete_results_df
                       )

IEE_complete_results_df$A <- factor(IEE_complete_results_df$A,labels=c("No","Yes"))


# ==============================================================================
# Goal 2a: Inferring Causality (Weighted Independence Estimating Equations)
# ==============================================================================
# A weighted analysis that accounts for confounding and survival bias.
# This estimates an exposure causal effect that is a function of SACE as in
# Shardell et al Biostatistics 2015; Tchetgen Tchetgen et al Am J Epidemiol 2012 

# Create lagged variables in long dataset
sim_obsdata_long$S <- as.numeric(!is.na(sim_obsdata_long$Y)) # longitudinal survival status
sim_obsdata_long <- sim_obsdata_long %>%
  arrange(id, YEAR) %>%      
  group_by(id) %>%          
  mutate(
    S_lag = lag(S, n = 1),   # Create the lag variable for S
    Y_lag = lag(Y, n = 1)    # Create the lag variable for Y
  ) %>%
  ungroup()                    

# subset data based on survivors at previous visit
# model death given survival at previous visit and compute estimated hazard probability
S_mod_sim_obsdata_long <- subset(sim_obsdata_long, YEAR>0 & S_lag==1)
S_mod <- glm(I(1-S)~YEAR.f*(X + scale(Y_lag) + W + A), S_mod_sim_obsdata_long, family="binomial")
S_mod_sim_obsdata_long$h_pred  <-  predict(S_mod, type = "response", na.action = na.pass)

# compute predicted survival probability using estimated hazards 
S_mod_sim_obsdata_long <- S_mod_sim_obsdata_long %>%
arrange(id, YEAR) %>%   
  group_by(id) %>%            
  mutate(cumSurv = cumprod(1-h_pred)) %>%
  ungroup()  

# merge data with predicted survival probility
sim_obsdata_long2 <- merge(sim_obsdata_long, S_mod_sim_obsdata_long[,c("id","YEAR","h_pred","cumSurv")], by=c("id","YEAR"),all=T)
sim_obsdata_long2 <- sim_obsdata_long2[with(sim_obsdata_long2, order(id, YEAR)), ]

# product of IPW for survival and confounding
sim_obsdata_long2$cumSurv[sim_obsdata_long2$YEAR == 0] <- 1
sim_obsdata_long2$IPW_A[sim_obsdata_long2$YEAR == 0] <- 1
sim_obsdata_long2$IPW_S <- mean(sim_obsdata_long2$S) * sim_obsdata_long2$S/sim_obsdata_long2$cumSurv
sim_obsdata_long2$IPW_WEE <- sim_obsdata_long2$IPW_A * sim_obsdata_long2$IPW_S
sim_obsdata_long2$IPW_WEE[sim_obsdata_long2$S == 0] <- 0

# WEE model
WEE_des <- svydesign(id=~id, weights=~IPW_WEE, data=sim_obsdata_long2) 
WEE_model <- svyglm(Y ~ A*YEAR.f, design=subset(WEE_des,YEAR>0))
cat("\n--- Goal 2: Inferring Causality WEE Results ---\n")
print(summary(WEE_model)$coefficients)

# standardized mean difference at Y1
WEE_lincom_MD_Y1 <- glht(WEE_model, linfct = c("A  = 0"))
confint(WEE_lincom_MD_Y1)$confint

# standardized mean difference at Y2
WEE_lincom_MD_Y2 <- glht(WEE_model, linfct = c("A + A:YEAR.f2 = 0"))
confint(WEE_lincom_MD_Y2)$confint

# standardized mean difference at Y3
WEE_lincom_MD_Y3 <- glht(WEE_model, linfct = c("A + A:YEAR.f3 = 0"))
confint(WEE_lincom_MD_Y3)$confint

# standardized means
X_matrix <- model.matrix(~  A*YEAR.f, grid)
WEE_glht <- glht(WEE_model, linfct = K_matrix %*% X_matrix)
WEE_ci_results <- confint(WEE_glht, calpha = univariate_calpha())$confint

# dataset of standardized means and 95% CIs
WEE_results_df <- data.frame(
  A = grid$A,
  YEAR.f = as.numeric(grid$YEAR.f),
  Estimate = WEE_ci_results[,"Estimate"],
  Lower_CI = WEE_ci_results[,"lwr"],
  Upper_CI = WEE_ci_results[,"upr"]
)

WEE_results_df <- rbind( c(0,0, mean_Y0, Y0.ests$conf.int[1], Y0.ests$conf.int[2]),
                         c(1,0, mean_Y0, Y0.ests$conf.int[1], Y0.ests$conf.int[2]),
                          WEE_results_df
                       )

WEE_results_df$A <- factor(WEE_results_df$A,labels=c("No","Yes"))


# ==============================================================================
# Goal 2b: Inferring Causality via SACE using Principal Stratification
# ==============================================================================
# Estimates the Survivor Average Causal Effect (SACE) for the "always survivors"
# principal stratum group using an alternative inverse probability weighting approach.
# A weighted analysis that accounts for confounding and survival bias as in 
# Tchetgen Tchetgen Stat Med 2014 

# creating a dataset of all exposed participant
S_mod_sim_expdata_long <- S_mod_sim_obsdata_long
S_mod_sim_expdata_long$A <- 1 

# estimated hazard if exposed
S_mod_sim_obsdata_long$h_pred_exp  <-  predict(S_mod, type = "response", newdata=S_mod_sim_expdata_long, na.action = na.pass)

# compute predicted survival probability using estimated hazards if exposed
S_mod_sim_obsdata_long <- S_mod_sim_obsdata_long %>%
arrange(id, YEAR) %>%   
  group_by(id) %>%           
  mutate(cumSurv_exp = cumprod(1-h_pred_exp)) %>%
  ungroup()  

# merge data with predicted survival probility if exposed
sim_obsdata_long3 <- merge(sim_obsdata_long2, S_mod_sim_obsdata_long[,c("id","YEAR","h_pred_exp","cumSurv_exp")], by=c("id","YEAR"),all=T)
sim_obsdata_long3 <- sim_obsdata_long3[with(sim_obsdata_long3, order(id, YEAR)), ]

# product of IPW for survival ratio (if exposed/actual exposure) and confounding
sim_obsdata_long3$cumSurv_exp[sim_obsdata_long3$YEAR == 0] <- 1
sim_obsdata_long3$IPW_S_exp <- sim_obsdata_long3$S * sim_obsdata_long3$cumSurv_exp/sim_obsdata_long3$cumSurv
sim_obsdata_long3$IPW_SACE <- sim_obsdata_long3$IPW_A * sim_obsdata_long3$IPW_S_exp
sim_obsdata_long3$IPW_SACE[sim_obsdata_long3$S == 0] <- 0

# SACE model
SACE_des <- svydesign(id=~id, weights=~IPW_SACE, data=sim_obsdata_long3) 
SACE_model <- svyglm(Y ~ A*YEAR.f, design=subset(SACE_des,YEAR>0))
cat("\n--- Goal 2: Inferring Causality SACE Results ---\n")
print(summary(SACE_model)$coefficients)

# standardized mean difference at Y1
SACE_lincom_MD_Y1 <- glht(SACE_model, linfct = c("A  = 0"))
confint(SACE_lincom_MD_Y1)$confint

# standardized mean difference at Y2
SACE_lincom_MD_Y2 <- glht(SACE_model, linfct = c("A + A:YEAR.f2 = 0"))
confint(SACE_lincom_MD_Y2)$confint

# standardized mean difference at Y3
SACE_lincom_MD_Y3 <- glht(SACE_model, linfct = c("A + A:YEAR.f3 = 0"))
confint(SACE_lincom_MD_Y3)$confint

# standardized means
SACE_glht <- glht(SACE_model, linfct = K_matrix %*% X_matrix)
SACE_ci_results <- confint(SACE_glht, calpha = univariate_calpha())$confint

# dataset of standardized means and 95% CIs
SACE_results_df <- data.frame(
  A = grid$A,
  YEAR.f = as.numeric(grid$YEAR.f),
  Estimate = SACE_ci_results[,"Estimate"],
  Lower_CI = SACE_ci_results[,"lwr"],
  Upper_CI = SACE_ci_results[,"upr"]
)

SACE_results_df <- rbind( c(0,0, mean_Y0, Y0.ests$conf.int[1], Y0.ests$conf.int[2]),
                         c(1,0, mean_Y0, Y0.ests$conf.int[1], Y0.ests$conf.int[2]),
                          SACE_results_df
                       )

SACE_results_df$A <- factor(SACE_results_df$A,labels=c("No","Yes"))


# ==============================================================================
# Goal 2 Complete data comparison: Weighted Estimating Equations (for confounding)
# ==============================================================================

# WEE model of complete data (weighting for confounding)
WEE_complete_des <- svydesign(id=~id, weights=~IPW_A, data=sim_compdata_long) 
WEE_complete_model <- svyglm(Y ~ A*YEAR.f, design=subset(WEE_complete_des,YEAR>0))
cat("\n--- WEE with complete data ---\n")
print(summary(WEE_complete_model)$coefficients)

# mean difference at Y1
WEE_complete_lincom_MD_Y1 <- glht(WEE_complete_model, linfct = c("A  = 0"))
confint(WEE_complete_lincom_MD_Y1)$confint

# mean difference at Y2
WEE_complete_lincom_MD_Y2 <- glht(WEE_complete_model, linfct = c("A + A:YEAR.f2 = 0"))
confint(WEE_complete_lincom_MD_Y2)$confint

# mean difference at Y3
WEE_complete_lincom_MD_Y3 <- glht(WEE_complete_model, linfct = c("A + A:YEAR.f3 = 0"))
confint(WEE_complete_lincom_MD_Y3)$confint

# standardized means for complete data
WEE_complete_glht <- glht(WEE_complete_model, linfct = K_matrix %*% X_matrix)
WEE_complete_ci_results <- confint(WEE_complete_glht, calpha = univariate_calpha())$confint

# dataset of standardized means and 95% CIs from complete data
WEE_complete_results_df <- data.frame(
  A = grid$A,
  YEAR.f = as.numeric(grid$YEAR.f),
  Estimate = WEE_complete_ci_results[,"Estimate"],
  Lower_CI = WEE_complete_ci_results[,"lwr"],
  Upper_CI = WEE_complete_ci_results[,"upr"]
)


WEE_complete_results_df <- rbind( c(0,0, mean_Y0, Y0.ests$conf.int[1], Y0.ests$conf.int[2]),
                         c(1,0, mean_Y0, Y0.ests$conf.int[1], Y0.ests$conf.int[2]),
                          WEE_complete_results_df
                       )

WEE_complete_results_df$A <- factor(WEE_complete_results_df$A,labels=c("No","Yes"))


# ==============================================================================
# Goal 3: The Immortal Cohort Framework (Linear Mixed Model)
# ==============================================================================
# Conventional mixed effects model that implicitly handles all sources of attrition 
# as ignorable dropout assuming missing at random instead of structural death.

# Mixed-effects model
mixed_model <- lme(
  fixed =Y ~ A*YEAR.f + Y0 + X, 
  random = ~ 1  | id, 
  data = subset(sim_obsdata_long, YEAR>0), 
  na.action = na.omit
)
cat("\n--- Goal 3: Linear Mixed Model (Immortal Cohort) Results ---\n")
print(summary(mixed_model)$tTable)

# mean difference at Y1
mixed_lincom_MD_Y1 <- glht(mixed_model, linfct = c("A  = 0"))
confint(mixed_lincom_MD_Y1)$confint

# mean difference at Y2
mixed_lincom_MD_Y2 <- glht(mixed_model, linfct = c("A + A:YEAR.f2 = 0"))
confint(mixed_lincom_MD_Y2)$confint

# mean difference at Y3
mixed_lincom_MD_Y3 <- glht(mixed_model, linfct = c("A + A:YEAR.f3 = 0"))
confint(mixed_lincom_MD_Y3)$confint

# adjusted means
mixed_glht <- glht(mixed_model, linfct = K_matrix %*% X_matrix_adjust)
mixed_ci_results <- confint(mixed_glht, calpha = univariate_calpha())$confint

# dataset of adjusted means and 95% CIs
mixed_results_df <- data.frame(
  A = grid$A,
  YEAR.f = as.numeric(grid$YEAR.f),
  Estimate = mixed_ci_results[,"Estimate"],
  Lower_CI = mixed_ci_results[,"lwr"],
  Upper_CI = mixed_ci_results[,"upr"]
)

mixed_results_df <- rbind( c(0,0, mean_Y0, Y0.ests$conf.int[1], Y0.ests$conf.int[2]),
                         c(1,0, mean_Y0, Y0.ests$conf.int[1], Y0.ests$conf.int[2]),
                          mixed_results_df
                       )

mixed_results_df$A <- factor(mixed_results_df$A,labels=c("No","Yes"))


# ==============================================================================
# Goal 3 Complete data comparison: Linear Mixed Model
# ==============================================================================

# Mixed-effects model of complete data
mixed_complete_model <- lme(
  fixed =Y ~ A*YEAR.f + Y0 + X, 
  random = ~ 1  | id, 
  data = subset(sim_compdata_long, YEAR>0), 
  na.action = na.omit
)
cat("\n--- Linear Mixed Model with complete data ---\n")
print(summary(mixed_complete_model)$tTable)

# mean difference at Y1
mixed_complete_lincom_MD_Y1 <- glht(mixed_complete_model, linfct = c("A  = 0"))
confint(mixed_complete_lincom_MD_Y1)$confint

# mean difference at Y2
mixed_complete_lincom_MD_Y2 <- glht(mixed_complete_model, linfct = c("A + A:YEAR.f2 = 0"))
confint(mixed_complete_lincom_MD_Y2)$confint

# mean difference at Y3
mixed_complete_lincom_MD_Y3 <- glht(mixed_complete_model, linfct = c("A + A:YEAR.f3 = 0"))
confint(mixed_complete_lincom_MD_Y3)$confint

# adjusted means
mixed_complete_glht <- glht(mixed_complete_model, linfct = K_matrix %*% X_matrix_adjust)
mixed_complete_ci_results <- confint(mixed_complete_glht, calpha = univariate_calpha())$confint

# dataset of adjusted means and 95% CIs from complete data
mixed_complete_results_df <- data.frame(
  A = grid$A,
  YEAR.f = as.numeric(grid$YEAR.f),
  Estimate = mixed_complete_ci_results[,"Estimate"],
  Lower_CI = mixed_complete_ci_results[,"lwr"],
  Upper_CI = mixed_complete_ci_results[,"upr"]
)

mixed_complete_results_df <- rbind( c(0,0, mean_Y0, Y0.ests$conf.int[1], Y0.ests$conf.int[2]),
                         c(1,0, mean_Y0, Y0.ests$conf.int[1], Y0.ests$conf.int[2]),
                          mixed_complete_results_df
                       )

mixed_complete_results_df$A <- factor(mixed_complete_results_df$A,labels=c("No","Yes"))


# ===================
#    Figure 2 
# ===================

## IEE plot data
IEE_complete_results_df$Model <- " 2A. IEE - Complete Data"
IEE_results_df$Model <- "2B. IEE - Survivors"
IEE_all_df <- rbind(IEE_complete_results_df, IEE_results_df)
IEE_all_df$Model <- factor(IEE_all_df$Model, levels=unique(IEE_all_df$Model))


## WEE and SACE plot data
WEE_complete_results_df$Model <- "2C. WEE - Complete Data"
SACE_results_df$Model <- "2D. SACE"
WEE_results_df$Model <- "2E. WEE - Survivors"
WEE_all_df <- rbind(WEE_complete_results_df, SACE_results_df, WEE_results_df)
WEE_all_df$Model <- factor(WEE_all_df$Model, levels=unique(WEE_all_df$Model))


## Mixed Model plot data
mixed_complete_results_df$Model <- "2F. Mixed Model - Complete Data"
mixed_results_df$Model <- "2G. Mixed Model - Survivors"
mixed_all_df <- rbind(mixed_complete_results_df, mixed_results_df)
mixed_all_df$Model <- factor(mixed_all_df$Model, levels=unique(mixed_all_df$Model))


## Plot y-axis limits
y_min <- 10*floor(min(IEE_all_df$Lower_CI,WEE_all_df$Lower_CI, mixed_all_df$Lower_CI)/10) 
y_max <- 10*ceiling(max(IEE_all_df$Upper_CI,WEE_all_df$Upper_CI, mixed_all_df$Upper_CI)/10)

## Create labels for plotting
add_labels <- function(df){

  df$Estimate_lbl <- sprintf("%.1f", df$Estimate)

  df$CI_lbl <- paste0(
    sprintf("%.1f", df$Lower_CI),
    ", ",
    sprintf("%.1f", df$Upper_CI)
  )

  df$CI_width <- df$Upper_CI - df$Lower_CI

  ## vertical spacing
  df$offset <- pmax(1.7, 0.15 * df$CI_width)

  ## common horizontal location for estimate and CI
  df$x_lbl <- ifelse(
    df$A == "No",
    as.numeric(df$YEAR.f) - 0.05,
    as.numeric(df$YEAR.f) - 0.10
  )

  ## estimate location
  df$Estimate_y <- ifelse(
    df$A == "No",
    df$Upper_CI + df$offset,
    df$Lower_CI - df$offset
  )

  ## CI above estimate for CKD=No
  ## CI below estimate for CKD=Yes
  df$CI_y <- ifelse(
    df$A == "No",
    df$Estimate_y + 1.1 * df$offset,
    df$Estimate_y - 1.1 * df$offset
  )

  df
}

## update datasets with label coordinates
IEE_all_df   <- add_labels(IEE_all_df)
WEE_all_df   <- add_labels(WEE_all_df)
mixed_all_df <- add_labels(mixed_all_df)


## IEE plots
 p_IEE_all <- ggplot(IEE_all_df, aes(x = YEAR.f, color = A, fill = A)) +
       geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.10) +
       geom_line(aes(y=Estimate), linewidth = 1, , alpha=0.5) +
       geom_point(aes(y = Estimate), size = 1, alpha=0.2) +
       geom_text(aes(x = x_lbl, y = Estimate_y, label = Estimate_lbl),
             hjust = 0.5, size = 3.1, fontface = "bold", show.legend = FALSE) +
       geom_text(aes(x = x_lbl, y = CI_y, label = CI_lbl),
             hjust = 0.5, size = 2.5, show.legend = FALSE) +
       labs(x = "Year", y = "Mean 3MS", color = "CKD", fill = "CKD") +
       scale_color_manual(values=c("No" = "blue","Yes"="red"))+
       scale_fill_manual(values=c("No" = "blue","Yes"="red"))+
       coord_cartesian(ylim = c(y_min, y_max),xlim = c(-0.25, 3.10)) + 
       facet_wrap(~ Model, ncol = 2) +
       theme(strip.text = element_text(hjust = 0)) +
       theme_minimal() +
       theme(legend.position = "bottom", legend.direction = "horizontal", strip.text = element_text(hjust = 0))

## WEE and SACE plots
 p_WEE_all <- ggplot(WEE_all_df, aes(x = YEAR.f, color = A, fill = A)) +
       geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.10) +
       geom_line(aes(y=Estimate), linewidth = 1, , alpha=0.5) +
       geom_point(aes(y = Estimate), size = 1, alpha=0.2) +
       geom_text(aes(x = x_lbl, y = Estimate_y, label = Estimate_lbl),
             hjust = 0.5, size = 3.1, fontface = "bold", show.legend = FALSE) +
       geom_text(aes(x = x_lbl, y = CI_y, label = CI_lbl),
             hjust = 0.5, size = 2.5, show.legend = FALSE) +
       labs(x = "Year", y = "Mean 3MS", color = "CKD", fill = "CKD") +
       scale_color_manual(values=c("No" = "blue","Yes"="red"))+
       scale_fill_manual(values=c("No" = "blue","Yes"="red"))+
       coord_cartesian(ylim = c(y_min, y_max),xlim = c(-0.25, 3.10)) + 
       facet_wrap(~ Model, ncol = 3) +
       theme(strip.text = element_text(hjust = 0)) + 
       theme_minimal() +
       theme(legend.position = "bottom", legend.direction = "horizontal", strip.text = element_text(hjust = 0))

## Mixed Model plots
 p_mixed_all <- ggplot(mixed_all_df, aes(x = YEAR.f, color = A, fill = A)) +
       geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.10) +
       geom_line(aes(y=Estimate), linewidth = 1, , alpha=0.5) +
       geom_point(aes(y = Estimate), size = 1, alpha=0.2) +
       geom_text(aes(x = x_lbl, y = Estimate_y, label = Estimate_lbl),
             hjust = 0.5, size = 3.1, fontface = "bold", show.legend = FALSE) +
       geom_text(aes(x = x_lbl, y = CI_y, label = CI_lbl),
             hjust = 0.5, size = 2.5, show.legend = FALSE) +
       labs(x = "Year", y = "Mean 3MS", color = "CKD", fill = "CKD") +
       scale_color_manual(values=c("No" = "blue","Yes"="red"))+
       scale_fill_manual(values=c("No" = "blue","Yes"="red"))+
       coord_cartesian(ylim = c(y_min, y_max),xlim = c(-0.25, 3.10)) +  
       facet_wrap(~ Model, ncol = 2) +
       theme_minimal() +
       theme(legend.position = "bottom", legend.direction = "horizontal", strip.text = element_text(hjust = 0)) 

get_legend <- function(a.gplot){
       tmp <- ggplot_gtable(ggplot_build(a.gplot))
       leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
       legend <- tmp$grobs[[leg]]
       return(legend)}

shared_legend <- get_legend(p_IEE_all)


## Save figure

pdf("Figure2.pdf", width = 8, height = 10) 

Fig2_alt <- grid.arrange(p_IEE_all + theme(legend.position = "none"), 
                         p_WEE_all + theme(legend.position = "none"),
                         p_mixed_all + theme(legend.position = "none"), nrow=3, bottom = shared_legend)
dev.off()



