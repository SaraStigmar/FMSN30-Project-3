library(readxl)
library(tidyverse)
library(car)
library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)
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
  labs(y = "Probability", x = "units (levels) of betaplasma", fill = "",
       title = "Distribution of betaplasma levels: Poisson?") 

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

#MODELS

#Only main terms, without calories - refer to project 1 and 2. 
model <- glm(betaplasma ~ bmi + age + fat + cholesterol + fiber + 
               alcohol + betadiet + smokstat + sex + vituse, family = "poisson", data = data)
summary(model)

#Rate ratios and confidence intervals
beta <- model$coefficients
RR <- exp(beta)
ci_beta = confint(model)
ci_RR = exp(ci_beta)

cbind(beta = beta, ci_beta, RR = RR, ci_RR) |> round(digits = 2)

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
