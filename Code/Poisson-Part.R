library(readxl)
library(tidyverse)
library(car)
library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)
library(MASS)
data <- read_excel("Data/carotene.xlsx")
summary(data)

#Betaplasma är count variabel, units av betaplasma är poisson fördelat och 
# och vi modeller de med våra kovariater

#Factors 
data <- data %>%
  mutate(smokstat = factor(smokstat,
                           levels = c(1, 2, 3),
                           labels = c("Never", "Former", "Current")))

data <- data %>%
  mutate(sex = factor(sex,
                      levels = c(1, 2),
                      labels = c("Male", "Female")))
# Set baseline
data$sex <- relevel(data$sex, ref = "Female")

data <- data %>%
  mutate(vituse = factor(vituse,
                         levels = c(1, 2, 3),
                         labels = c("Often", "Rarely", "Never")))

# Set baseline
data$vituse <- relevel(data$vituse, ref = "Never")



################## Start of Poisson ##########################

#Frequency table over betaplasma levels
data |> count(betaplasma) -> betaplasma_table
betaplasma_table


#Add proportions of observed levels
betaplasma_table |> mutate(observed = n/sum(n)) -> betaplasma_table
betaplasma_table


# Average
mu <- mean(data$betaplasma)
mu 
# => 189.9238

# Variance
var(data$betaplasma)
# => 33477.52 

# NOTICE! Poisson is not ideal since mu != var

# Expected frequencies
betaplasma_table |> 
  mutate(expected = dpois(betaplasma, mu)) -> betaplasma_table
betaplasma_table

#Plot expected and observed
ggplot(betaplasma_table) +
  geom_col(aes(x = betaplasma + 0.2, y = observed, fill = "Observed"), 
           width = 0.4) +
  geom_col(aes(x = betaplasma - 0.2, y = expected, fill = "Expected"), 
           width = 0.4) +
  xlim(0, 500) +
  labs(y = "Probability", x = "units (levels) of betaplasma", fill = "",
       title = "Distribution of betaplasma levels vs Poisson Distribution") 
#OBS! Ta bort raden xlim(0, 500) om du vill se den faktiska fördelningen - vi har en lång svans upp till över 1000 pga höga betaplasma värden 

##############################################################################
#Alternative way?
#In the example they get the same plot - but I get a different result. 
#Don't know why...
pivot_longer(betaplasma_table, cols = c("observed", "expected"), 
             values_to = "Proportion", 
             names_to = "variable") -> betaplasma_long
betaplasma_long

ggplot(betaplasma_long, aes(x = betaplasma, y = Proportion, fill = variable)) +
  geom_col(position = "dodge") +
  labs(y = "Probability", x = "units (levels) of betaplasma", fill = "",
       title = "Distribution of betaplasma levels: Poisson?") +
  theme(text = element_text(size = 14))
###########################################################################

#Mean and variance by program
data |> summarise(mu = mean(betaplasma),
                   s2 = var(betaplasma))
  # RESULT: mu = 190, s^2 = 33478
  # NOTICE: Mu is far from the variance, which indicates Poisson is a poor fit. 
  # This is called OVERDISPERSION - that the variance of our data is much greater
  # than the variance of the distribution - THIS is what causes the poor fit. 
  # Since Poisson requires E(Y) = mu, Var(Y) = Mu, 

#Model testing
  #-> Only main terms, without calories - refer to project 1 and 2. 
model <- glm(betaplasma ~ bmi + age + fat + cholesterol + fiber + 
               alcohol + betadiet + smokstat + sex + vituse, family = "poisson", data = data)
summary(model)

#Rate ratios and confidence intervals
beta <- model$coefficients
RR <- exp(beta)
ci_beta = confint(model)
ci_RR = exp(ci_beta)

cbind(beta = beta, ci_beta, RR = RR, ci_RR) |> round(digits = 2)
  #NOTICE: The betas have almost no relevance


#Estimated mean with c.i.
data_pred <- cbind(
  data,
  xb = predict(model, se.fit = TRUE))
data_pred |> mutate(
  xb.residual.scale = NULL,
  xb.lwr = xb.fit - 1.96*xb.se.fit,
  xb.upr = xb.fit + 1.96*xb.se.fit,
  muhat = exp(xb.fit),
  mu.lwr = exp(xb.lwr),
  mu.upr = exp(xb.upr)) -> data_pred

glimpse(data_pred)

#Plot estimated mean with its c.i. - don't know if this serves a purpose here -
#since there are many varibles to plot against. 

#For example against bmi
ggplot(data_pred, aes(bmi, betaplasma)) +
  geom_point() +
  geom_line(aes(y = muhat), linewidth = 1) +
  geom_ribbon(aes(ymin = mu.lwr, ymax = mu.upr), alpha = 0.1) +
  labs(title = "Expected number of awards",
       caption = "95% confidence interval",
       color = "program")

#Or against bmi but for different genders
ggplot(data_pred, aes(bmi, betaplasma, color = sex)) +
  geom_point() +
  geom_line(aes(y = muhat), linewidth = 1) +
  geom_ribbon(aes(ymin = mu.lwr, ymax = mu.upr), alpha = 0.1) +
  labs(title = "Expected number of awards",
       caption = "95% confidence interval",
       color = "program") +
  facet_wrap(~ sex)





#MODELS (OBS all excl calories)
###############################################################################

#Model 1. Null Model
Model_1 <- glm(betaplasma ~ 1, family = "poisson", data = data)

#Model 2. Full model = Main terms (10) + interaction (66), total 76 variables (+ 1 if coeffiecients bcz of intercept)
Model_2 <- glm(betaplasma ~ (bmi + age + fat + cholesterol + fiber + 
                                 alcohol + betadiet + smokstat + sex + vituse)^2,
               family = "poisson", data = data)

#Model 3. Backward Elimination, start from the full model 
Model_3 <- step(Model_2,
                      direction = "backward",
                      scope = list(lower = formula(Model_1),
                      upper = formula(Model_2)),
                      k = log(nobs(Model_1)))

log(nrow(data)) == log(nobs(Model_1))

#Model 4. Forward Selection, start from the null model 
Model_4 <- step(Model_1,
                direction = "forward",
                scope = list(lower = formula(Model_1),
                             upper = formula(Model_2)),
                k = log(nobs(Model_1)))

#Model 5. Stepwise regression from model 3
Model_5 <- step(Model_3,
                direction = "both",
                scope = list(lower = formula(Model_1),
                upper = formula(Model_2)),
                k = log(nobs(Model_3)))

#Model 6. Stepwise regression from model 4
Model_6 <- step(Model_4,
                direction = "both",
                scope = list(lower = formula(Model_1),
                upper = formula(Model_2)),
                k = log(nobs(Model_4)))

summary(Model_1)
summary(Model_2)
summary(Model_3)
summary(Model_4)
summary(Model_5)
summary(Model_6)

############ OBS KODTEST #######################################################
#Här testade jag så jag använder selektionsmetoderna rätt genom att switcha till
#negativ bionomial fördelning och jag får då små modeller. 
model_test <- glm.nb(betaplasma ~ (bmi + age + fat + cholesterol + fiber + 
                                     alcohol + betadiet + smokstat + sex + vituse)^2,
                     data = data)
model_test_null <- glm.nb(betaplasma ~ 1,
                          data = data)
model_test_backward <- step(model_test,
                direction = "backward",
                scope = list(lower = formula(model_test_null),
                             upper = formula(model_test)),
                k = log(nobs(model_test_null)),
                trace = 1)

#RESULT:
#Step:  AIC=3818.4
#betaplasma ~ bmi + cholesterol + betadiet + vituse

# Df Deviance    AIC
#<none>             335.50 3818.4
#- cholesterol  1   346.78 3823.9
#- betadiet     1   354.98 3832.1
#- vituse       2   360.95 3832.3
#- bmi          1   363.32 3840.5
################################################################################



#See if we can reduce vituse - the categorical variable, into two categories:
#------------------------------------------------------------------------------

# Transform three levels to two, but assigning vituse category "Never" into "Never" 
# and if the category is not "Never" it's assigned into "Used". 
new_data <- data
new_data$vituseReduced <- ifelse(new_data$vituse == "Never", "Never", "Used")
new_data$vituseReduced <- factor(new_data$vituseReduced)
#Then we create a new model, based on this two-factor collapsed version (reduced):
Model_5_Reduced <- glm(betaplasma ~ bmi + age + fat + cholesterol + fiber + alcohol + 
                                       betadiet + smokstat + sex + vituseReduced + bmi:fat + bmi:cholesterol + 
                                       bmi:fiber + bmi:alcohol + bmi:betadiet + bmi:smokstat + bmi:sex + 
                                       bmi:vituseReduced + age:fat + age:cholesterol + age:fiber + age:betadiet + 
                                       age:smokstat + age:vituseReduced + fat:cholesterol + fat:alcohol + 
                                       fat:betadiet + fat:smokstat + fat:sex + fat:vituseReduced + cholesterol:fiber + 
                                       cholesterol:alcohol + cholesterol:betadiet + cholesterol:smokstat + 
                                       fiber:alcohol + fiber:betadiet + fiber:smokstat + fiber:sex + 
                                       fiber:vituseReduced + alcohol:betadiet + alcohol:smokstat + alcohol:sex + 
                                       alcohol:vituseReduced + betadiet:vituseReduced + smokstat:sex + smokstat:vituseReduced + 
                                       sex:vituseReduced,
                        family = "poisson", data = new_data)
#Save variables
sum_1 <- summary(Model_5)
sum_2 <- summary(Model_5_Reduced)
deviance_1 <- sum_1$deviance
deviance_2 <- sum_2$deviance
df_1 <- sum_1$df.residual
df_2 <- sum_2$df.residual

#Calculate differences
D_diff <- deviance_2 - deviance_1
df_diff <- df_2 - df_1
chi2_alpha <- qchisq(p = 1 - 0.05, df = df_diff)
Pvalue <- pchisq(q = D_diff, df = df_diff, lower.tail = FALSE)

#I create a table and remove the vertical title using rownames, just to clean it up! 
table_1e <- cbind(D_diff, df_diff, chi2_alpha, Pvalue)
rownames(table_1e) <- c(" ")
table_1e

# 1. Name of test: Partial likelihood ratio test 
# 2. Null Hypothesis, H_0: Beta_often = Beta_Rarely. In words, that the model with 
#    one dummy variable Beta_Used is enough. 
# 3. The value of the test statistic is D = 2543.28
# 4. The distribution is Chi-squared
# 5. The P-value is P = 0
# 6. The conclusion is to with great certainty reject H_0 since D_diff >> Chi2_alpha, 
# and since Pvalue < 0.05.

#CONCLUSION: Keep the original model_5 and throw away model_5_reduced. 



#BIC AND AIC 
df1 <- data.frame(variable = names(Model_1$coefficients),
                  b_model1 = Model_1$coefficients, row.names = NULL)
df2 <- data.frame(variable = names(Model_2$coefficients),
                  b_model2 = Model_2$coefficients, row.names = NULL)
df3 <- data.frame(variable = names(Model_3$coefficients),
                  b_model3 = Model_3$coefficients, row.names = NULL)
df4 <- data.frame(variable = names(Model_4$coefficients),
                  b_model4 = Model_4$coefficients, row.names = NULL)
df5 <- data.frame(variable = names(Model_5$coefficients),
                  b_model5 = Model_5$coefficients, row.names = NULL)
df6 <- data.frame(variable = names(Model_6$coefficients),
                  b_model6 = Model_6$coefficients, row.names = NULL)
All_Models <- full_join(df1, df2) |> full_join(df3) |> full_join(df4) |> full_join(df5) |> full_join(df6)

#Writes a table:
table_3b <- kable(All_Models, caption = "Table 3(b): Estimated β-parameters in each of the six models")
table_3b
#Exported to excel:
#write.csv(All_Models, "Table_3b_ny.csv", row.names = FALSE)

# Calculate McFadden’s adjusted pseudo R2, AIC and BIC for all models from Table.3(b), 
# and indicate which model is best, according to each of these criteria.

aic <- AIC(Model_1, Model_2, Model_3, Model_4, Model_5, Model_6)
bic <- BIC(Model_1, Model_2, Model_3, Model_4, Model_5, Model_6)

#Create dataframe for the AIC- and BIC
collect.AICetc <- data.frame(aic, bic, D = c(Model_1$deviance, Model_2$deviance, Model_3$deviance, Model_4$deviance, Model_5$deviance, Model_6$deviance)) 

#Remove unnecessary df.1 column
collect.AICetc |> mutate(df.1 = NULL, D0 = Model_1$deviance,
                         p = df - 1) -> collect.AICetc

#Calculate Psuedo R & McFadden's adjusted psuedo R2
collect.AICetc |> mutate(
  pseudoR2 = 1 - D/D0,
  pseudoR2adj = 1 - (D + p)/D0) -> collect.AICetc

#Show result
collect.AICetc

#Model with best (lowest) AIC: Model_4 (Forward Selection)
#Model with best (lowest) BIC: Either Model_3, 5 or 6 - they appear the same
#( Model with best (highest) R^2_adj: Model_2 )


# Leverage 
data_pred |> mutate(v_3 = hatvalues(Model_3)) -> data_pred
data_pred |> mutate(v_4 = hatvalues(Model_4)) -> data_pred

#Standardized deviance residuals
data_infl_3 <- influence(Model_3)
data_infl_4 <- influence(Model_4)
#These gives deviance and pearson residuals!
#Extract deviance residuals and standardize 
data_pred |> mutate(devres_3 = data_infl_3$dev.res,
                     std.devres_3 = devres_3/sqrt(1 - v_3),
                    devres_4 = data_infl_4$dev.res,
                    std.devres_4 = devres_4/sqrt(1 - v_4)) -> data_pred

#Cooks D
data_pred |> mutate(D_3 = cooks.distance(Model_3)) -> data_pred
data_pred |> mutate(D_4 = cooks.distance(Model_4)) -> data_pred




# Standardiserade residualer mot förväntat värde, dvs µ^hat

ggplot(data_pred, aes(x = fitted(Model_3), y = std.devres_3)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  geom_hline(yintercept = -2, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 3, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -3, linetype = "dashed", color = "blue") +
  labs(x = "Expected values", y = "Standardized deviance-residuals",
       title = "Residuals plot – Model 3")

ggplot(data_pred, aes(x = fitted(Model_4), y = std.devres_4)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  geom_hline(yintercept = -2, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 3, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -3, linetype = "dashed", color = "blue") +
  labs(x = "Expected values", y = "Standardized deviance-residuals",
       title = "Residuals plot – Model 4")

# Leverage mot standardiserade residualer, färg enligt cooks D

ggplot(data_pred, aes(x = v_3, y = std.devres_3, color = D_3)) +
  geom_point(alpha = 0.6) +
  scale_size_continuous(name = "Cook’s D") +
  labs(x = "Leverage", y = "Standardized deviance residuals",
       title = "Influence plot – Model 3") +
  theme_minimal()

ggplot(data_pred, aes(x = v_4, y = std.devres_4, color = D_4)) +
  geom_point(alpha = 0.6) +
  scale_size_continuous(name = "Cook’s D") +
  labs(x = "Leverage", y = "Standardized deviance residuals",
       title = "Influence plot – Model 4") +
  theme_minimal()


#Cooks D

ggplot(data_pred, aes(x = 1:nrow(data_pred), y = D_3)) +
  geom_bar(stat = "identity") +
  labs(x = "Observation", y = "Cook’s D", title = "Cook’s D – Model 3") +
  geom_hline(yintercept = 4/nrow(data_pred), linetype = "dotted", color = "red")

#Cooks D

ggplot(data_pred, aes(x = 1:nrow(data_pred), y = D_4)) +
  geom_bar(stat = "identity") +
  labs(x = "Observation", y = "Cook’s D", title = "Cook’s D – Model 4") +
  geom_hline(yintercept = 4/nrow(data_pred), linetype = "dotted", color = "red")


##################################

# HYPOTHESIS TESTING

##################################

D_diff_3 <- Model_3$null.deviance - Model_3$deviance
D_diff_4 <- Model_4$null.deviance - Model_4$deviance

df_diff_3 <- Model_3$df.null - Model_3$df.residual
df_diff_4 <- Model_4$df.null - Model_4$df.residual

chi2_alpha_3 <- qchisq(p = 1 - 0.05, df = df_diff_3)
chi2_alpha_4 <- qchisq(p = 1 - 0.05, df = df_diff_4)
Pvalue_3 <- pchisq(q = D_diff_3, df = df_diff_3, lower.tail = FALSE)
Pvalue_4 <- pchisq(q = D_diff_4, df = df_diff_4, lower.tail = FALSE)

#I create a table and remove the vertical title using rownames, just to clean it up! 
table_3 <- cbind(D_diff_3, df_diff_3, chi2_alpha_3, Pvalue_3)
table_4 <- cbind(D_diff_4, df_diff_4, chi2_alpha_4, Pvalue_4)
rownames(table_3) <- c(" ")
rownames(table_4) <- c(" ")
table_3
table_4

#RESULT:
#> table_3
#D_diff_3 df_diff_3 chi2_alpha_3 Pvalue_3
#18505.63        66     85.96491        0
#> table_4
#D_diff_4 df_diff_4 chi2_alpha_4 Pvalue_4
#18530.88        73     93.94534        0





















#------------------------------------------------------------------------------
# Calculate the standardized deviance residuals for the model with the best AIC 
# and the model with the best BIC, from 3(c). 





# Standardized residuals for respective model 
res_Model_3 <- rstandard(Model_3, type = "deviance")
res_Model_4 <- rstandard(Model_4, type = "deviance")

#Linear predictor
lp_Model_3 <- predict(Model_3, type = "link")
lp_Model_4 <- predict(Model_4, type = "link")




# Task:
# Plot the standardized deviance residuals, with suitable reference lines, 
# against the linear predictors for each of the two models. Also make QQ-plots 
# for the residuals. Discuss which of the models has the best behaved residuals


#UPDATED: Second try for plot 1 in order to get color-coded 

#Plot 1:
Model_5_infl <- influence(Stepwise_From_Backward)
Model_6_infl <- influence(Stepwise_From_Forward_Reduced)

data_pred <- cbind(data,
                   xbeta5 = predict(Stepwise_From_Backward),
                   xbeta6 = predict(Stepwise_From_Forward_Reduced),
                   v5 = Model_5_infl$hat,
                   v6 = Model_6_infl$hat)


data_pred |> mutate(devresid5 = Model_5_infl$dev.res,
                    stddevresid5 = devresid5/sqrt(1 - v5)) -> data_pred
data_pred |> mutate(devresid6 = Model_6_infl$dev.res,
                    stddevresid6 = devresid6/sqrt(1 - v6)) -> data_pred

ggplot(data_pred, aes(x = xbeta5, 
                      y = stddevresid5, 
                      color = as.factor(lowplasma_01))) +
  geom_point() +
  geom_hline(yintercept = c(-3, -2, 0, 2, 3), 
             linetype = c("dotted", "dashed", "solid", "dashed", "dotted"),
             linewidth = 1) +
  labs(title = "Standardized deviance residuals vs linear predictor",
       x = "Linear predictor, xb", y = "Standardized deviance residuals, devstd",
       color = "Y")

ggplot(data_pred, aes(x = xbeta6, 
                      y = stddevresid6, 
                      color = as.factor(lowplasma_01))) +
  geom_point() +
  geom_hline(yintercept = c(-3, -2, 0, 2, 3), 
             linetype = c("dotted", "dashed", "solid", "dashed", "dotted"),
             linewidth = 1) +
  labs(title = "Standardized deviance residuals vs linear predictor",
       x = "Linear predictor, xb", y = "Standardized deviance residuals, devstd",
       color = "Y")



