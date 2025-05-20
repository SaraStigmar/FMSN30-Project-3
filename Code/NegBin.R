source("Code/Base.R")
library(ggeffects)

# Select only numeric columns
numeric_df <- data %>% select_if(is.numeric)

# Calculate correlations between betaplasma and other variables
correlations <- sapply(numeric_df, function(x) cor(numeric_df$betaplasma, x, use = "complete.obs"))
print(correlations)
# fiber and bmi has strongest correlation

#betaplasma_glmnb_fiber <- glm.nb(betaplasma ~ fiber, data=data)
#summary(betaplasma_glmnb_fiber)

betaplasma_glmnb_bmi <- glm.nb(betaplasma ~ bmi, data=data)
summary(betaplasma_glmnb_bmi)

# Confidence intervals for beta
beta <- betaplasma_glmnb_bmi$coefficients
ci_beta <- confint(betaplasma_glmnb_bmi)
cbind(beta = beta, ci_beta) |> round(digits=2)

cbind(RR = exp(beta), exp(ci_beta)) |> round(digits=2)


#estimate means
betaplasma_bmi_pred <- cbind(
  data,
  xb = predict(betaplasma_glmnb_bmi, se.fit = TRUE))
betaplasma_bmi_pred |> mutate(
  xb.residual.scale = NULL,
  xb.lwr = xb.fit - 1.96*xb.se.fit,
  xb.upr = xb.fit + 1.96*xb.se.fit,
  muhat = exp(xb.fit),
  mu.lwr = exp(xb.lwr),
  mu.upr = exp(xb.upr)
) -> betaplasma_bmi_pred

ggplot(betaplasma_bmi_pred, aes(bmi, betaplasma)) +
  geom_point() +
  geom_line(aes(y = muhat), linewidth=1) +
  geom_ribbon(aes(ymin = mu.lwr, ymax = mu.upr), alpha=0.1)

#Likelihood ratio test
betaplasma_null <- glm.nb(betaplasma ~1, data = data)
betaplasma_null$null.deviance
betaplasma_glmnb_bmi$null.deviance

anova(betaplasma_null,betaplasma_glmnb_bmi) #better than null

#Full model without calories
fullmodel <- glm.nb(betaplasma ~ bmi + age + fat + cholesterol + fiber + 
               alcohol + betadiet + smokstat + sex + vituse, data = data)
summary(fullmodel) #Many insignificant terms

#Confidence intervals
beta <- fullmodel$coefficients
ci_beta = confint(fullmodel)
cbind(beta = beta, ci_beta) |> round(digits = 2)

#Estimate means with confidence intervals
fullmodel_pred <- cbind(
  data,
  xb = predict(fullmodel, se.fit = TRUE))
fullmodel_pred |> 
  mutate(
    xb.residual.scale = NULL,
    xb.lwr = xb.fit - 1.96*xb.se.fit,
    xb.upr = xb.fit + 1.96*xb.se.fit,
    muhat = exp(xb.fit),
    mu.lwr = exp(xb.lwr),
    mu.upr = exp(xb.upr)) ->
  fullmodel_pred

ggplot(fullmodel_pred, aes(bmi,betaplasma, color = vituse)) +
  geom_point() +
  geom_line(aes(y = muhat), linewidth = 1) +
  geom_ribbon(aes(ymin = mu.lwr, ymax = mu.upr), alpha = 0.1) +
  facet_wrap(~ vituse)


pred_bmi <- ggeffects::ggpredict(fullmodel, terms = "bmi [all]")
plot(pred_bmi)

ggplot(data, aes(x = bmi, y = betaplasma, color = vituse)) +
  geom_point(alpha = 0.4) +  # Observations
  # Add model predictions
  geom_line(data = pred_bmi, aes(x = x, y = predicted, color = group), size = 1) +
  geom_ribbon(data = pred_bmi, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), 
              alpha = 0.2, inherit.aes = FALSE) +
  labs(y = "Betaplasma", x = "BMI", title = "Observed data and model-predicted means") +
  theme_minimal()


#Selection
Fullscope_Model <- glm.nb(betaplasma ~ (bmi + age + fat + cholesterol + fiber + 
                                        alcohol + betadiet + smokstat + sex + vituse)^2,data = data)

vars <- c("bmi", "age", "fat", "cholesterol", "fiber", 
          "alcohol", "betadiet", "smokstat", "sex", "vituse")

# Create formula string with main effects and all 2-way interactions
scope_formula <- as.formula(
  paste("~ (", paste(vars, collapse = " + "), ")^2")
)

nullmodel <- glm.nb(betaplasma ~ 1, data = data)

n <- nrow(data)  # sample size
forward_model <- step(nullmodel,
                  scope = list(lower = ~1, upper = formula(Fullscope_Model)),
                  direction = "forward",
                  trace = TRUE,
                  k = log(n))  # BIC instead of AIC

backward_model <- step(Fullscope_Model,
                       scope = list(lower = ~1, upper = formula(Fullscope_Model)),
                       direction = "backward",
                       trace = TRUE,
                       k = log(n))  # BIC

stepwise_model_n <- step(nullmodel,
                       scope = list(lower = ~1, upper = formula(Fullscope_Model)),
                       direction = "both",
                       trace = TRUE,
                       k = log(n))  # BIC

stepwise_model_f <- step(Fullscope_Model,
                         scope = list(lower = ~1, upper = formula(Fullscope_Model)),
                         direction = "both",
                         trace = TRUE,
                         k = log(n))  # BIC



summary(forward_model)
summary(backward_model)
summary(stepwise_model_n)
summary(stepwise_model_f)

#All the same
new_data <- data

# Transform three levels to two, but assigning vituse category "Never" into "Never" 
# and if the category is not "Never" it's assigned into "Used". 
new_data$vituse_reduced <- ifelse(new_data$vituse == "Never", "Never", "Used")
new_data$vituse_reduced <- factor(new_data$vituse_reduced)

model_reduced <- glm.nb(betaplasma ~ bmi + cholesterol + betadiet + vituse_reduced, 
                                     data = new_data)
model_reduced$coefficients

#Save variables
sum_1 <- summary(forward_model)
sum_2 <- summary(model_reduced)
deviance_1 <- sum_1$deviance
deviance_2 <- sum_2$deviance
df_1 <- sum_1$df.residual
df_2 <- sum_2$df.residual

#Calculate differences
D_diff <- deviance_2 - deviance_1
df_diff <- df_2 - df_1
chi2_alpha <- qchisq(p = 1 - 0.05, df = df_diff)
Pvalue <- pchisq(q = D_diff, df = df_diff, lower.tail = FALSE) #Cant conlude that full model is better than reduced

#Is this model better than the simple
anova(betaplasma_glmnb_bmi,model_reduced)

#AIC, BIC, pseudo R2

data.frame(AIC(nullmodel, model_reduced, forward_model), 
           BIC(nullmodel, model_reduced, forward_model),
           D = c(nullmodel$deviance, model_reduced$deviance, forward_model$deviance),
           D0 = c(nullmodel$null.deviance, model_reduced$null.deviance, forward_model$null.deviance),
           p = c(nullmodel$df.null - nullmodel$df.residual,
                 model_reduced$df.null - model_reduced$df.residual,
                 forward_model$df.null - forward_model$df.residual)) |> 
  mutate(df.1 = NULL) -> collect.AICetc
collect.AICetc


collect.AICetc <- data.frame(
  AIC(nullmodel, model_reduced, forward_model, betaplasma_glmnb_bmi), 
  BIC(nullmodel, model_reduced, forward_model, betaplasma_glmnb_bmi),
  D = c(nullmodel$deviance, model_reduced$deviance, forward_model$deviance, betaplasma_glmnb_bmi$deviance),
  D0 = c(nullmodel$null.deviance, model_reduced$null.deviance, forward_model$null.deviance, betaplasma_glmnb_bmi$null.deviance),
  p = c(nullmodel$df.null - nullmodel$df.residual,
        model_reduced$df.null - model_reduced$df.residual,
        forward_model$df.null - forward_model$df.residual,
        betaplasma_glmnb_bmi$df.null - betaplasma_glmnb_bmi$df.residual)
)

collect.AICetc

collect.AICetc |> mutate(
  R2 = 1 - D/D0,
  R2.adj = 1 - (D + p)/D0) |> round(digits = 3)

betaplasma_pred <- cbind(
  new_data,
  xb = predict(model_reduced, se.fit = TRUE))
betaplasma_pred |> 
  mutate(
    xb.residual.scale = NULL,
    xb.lwr = xb.fit - 1.96*xb.se.fit,
    xb.upr = xb.fit + 1.96*xb.se.fit,
    muhat = exp(xb.fit),
    mu.lwr = exp(xb.lwr),
    mu.upr = exp(xb.upr)) ->
  betaplasma_pred


betaplasma_infl <- influence(model_reduced)
betaplasma_pred |> mutate(
  v = hatvalues(model_reduced),
  devres = betaplasma_infl$dev.res,
  std.devres = devres/sqrt(1 - v)) ->
  betaplasma_pred

betaplasma_pois_glm <- glm(betaplasma ~ bmi + cholesterol + betadiet + vituse_reduced, family = "poisson", data = new_data)
betaplasma_pois_infl <- influence(betaplasma_pois_glm)

betaplasma_pred |> mutate(
  xb_pois = predict(betaplasma_pois_glm),
  v_pois = hatvalues(betaplasma_pois_glm),
  devres_pois = betaplasma_pois_infl$dev.res,
  std.devres_pois = devres_pois/sqrt(1 - v_pois)) -> betaplasma_pred

betaplasma_pred |> summarise(
  min_nb = min(std.devres),
  max_nb = max(std.devres),
  min_pois = min(std.devres_pois),
  max_pois = max(std.devres_pois))

ggplot(betaplasma_pred, aes(x = xb_pois)) +
  geom_point(aes(y = std.devres_pois), size = 2) +
  geom_hline(yintercept = c(-3, -2, 0, 2, 3), 
             linetype = c("dotted", "dashed", "solid", "dashed", "dotted"), 
             linewidth = 1) +
  expand_limits(y = c(-22.5, 49)) +
  labs(y = "std dev.res", x = "xb",
       title = "Betaplasma: Poisson model")

ggplot(betaplasma_pred, aes(x = xb.fit)) +
  geom_point(aes(y = std.devres), size = 2) +
  geom_hline(yintercept = c(-3, -2, 0, 2, 3), 
             linetype = c("dotted", "dashed", "solid", "dashed", "dotted"), 
             linewidth = 1) +
  expand_limits(y = c(-22.5, 49)) +
  labs(y = "std dev.res", x = "xb",
       title = "Betaplasma: Negbin model")
#For this model negative binomial clearly better

D_pois <- -2*logLik(betaplasma_pois_glm)[1]
D_pois

D_nb <- -2*logLik(model_reduced)[1]
D_nb

D_diff <- D_pois - D_nb
D_diff

qchisq(1 - 0.05, 1)
pchisq(D_diff, 1, lower.tail = FALSE)
