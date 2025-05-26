source("Code/Base.R")
source("Code/Functions/boot_negbin.R")
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
# Add predicted values with standard errors
betaplasma_bmi_pred <- cbind(
  data,
  xb = predict(betaplasma_glmnb_bmi, se.fit = TRUE)
)

# Add confidence intervals and fitted means
betaplasma_bmi_pred <- betaplasma_bmi_pred |> mutate(
  xb.residual.scale = NULL,
  xb.lwr = xb.fit - 1.96 * xb.se.fit,
  xb.upr = xb.fit + 1.96 * xb.se.fit,
  muhat = exp(xb.fit),
  mu.lwr = exp(xb.lwr),
  mu.upr = exp(xb.upr)
)

# Extract influence measures and hat values
bmi_infl <- influence(betaplasma_glmnb_bmi)
v <- hatvalues(betaplasma_glmnb_bmi)

# Calculate deviance and Pearson residuals
devres <- bmi_infl$dev.res
pearson_res <- residuals(betaplasma_glmnb_bmi, type = "pearson")

# Add residuals and standardized residuals to dataframe
betaplasma_bmi_pred <- betaplasma_bmi_pred |> mutate(
  v = v,
  devres = devres,
  std.devres = devres / sqrt(1 - v),
  pearsonres = pearson_res,
  std.pearsonres = pearson_res / sqrt(1 - v)
)

ggplot(betaplasma_bmi_pred, aes(bmi, betaplasma)) +
  geom_point() +
  geom_line(aes(y = muhat), linewidth=1) +
  geom_ribbon(aes(ymin = mu.lwr, ymax = mu.upr), alpha=0.1)

#Likelihood ratio test
betaplasma_null <- glm.nb(betaplasma ~1, data = data)
betaplasma_null$null.deviance
betaplasma_glmnb_bmi$null.deviance

anova(betaplasma_null,betaplasma_glmnb_bmi) #better than null


#Residuals:
#Standard deviance
ggplot(betaplasma_bmi_pred, aes(x = xb.fit)) +
  geom_point(aes(y = std.devres), size = 2) +
  geom_hline(yintercept = c(-3, -2, 0, 2, 3), 
             linetype = c("dotted", "dashed", "solid", "dashed", "dotted"), 
             linewidth = 1) +
  expand_limits(y = c(-4.5, 7.5)) +
  labs(y = "std dev.res", x = "xb",
       title = "Betaplasma - BMI: Negbin model (Std Dev Residuals)")

# Plot using standardized Pearson residuals
ggplot(betaplasma_bmi_pred, aes(x = xb.fit)) +
  geom_point(aes(y = std.pearsonres), size = 2) +
  geom_hline(yintercept = c(-3, -2, 0, 2, 3), 
             linetype = c("dotted", "dashed", "solid", "dashed", "dotted"), 
             linewidth = 1) +
  expand_limits(y = c(-4.5, 7.5)) +
  labs(y = "std Pearson res.", x = "xb",
       title = "Betaplasma - BMI: Negbin model (Std Pearson Residuals)")

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


# Extract influence measures and hat values
full_infl <- influence(fullmodel)
v <- hatvalues(fullmodel)

# Calculate deviance and Pearson residuals
devres <- full_infl$dev.res
pearson_res <- residuals(fullmodel, type = "pearson")

# Add residuals and standardized residuals to dataframe
fullmodel_pred <- fullmodel_pred |> mutate(
  v = v,
  devres = devres,
  std.devres = devres / sqrt(1 - v),
  pearsonres = pearson_res,
  std.pearsonres = pearson_res / sqrt(1 - v)
)


#Residuals:
#Standard deviance
ggplot(fullmodel_pred, aes(x = xb.fit)) +
  geom_point(aes(y = std.devres), size = 2) +
  geom_hline(yintercept = c(-3, -2, 0, 2, 3), 
             linetype = c("dotted", "dashed", "solid", "dashed", "dotted"), 
             linewidth = 1) +
  expand_limits(y = c(-4.5, 7.5)) +
  labs(y = "std dev.res", x = "xb",
       title = "Betaplasma - Full Model: Negbin model (Std Dev Residuals)")

# Plot using standardized Pearson residuals
ggplot(fullmodel_pred, aes(x = xb.fit)) +
  geom_point(aes(y = std.pearsonres), size = 2) +
  geom_hline(yintercept = c(-3, -2, 0, 2, 3), 
             linetype = c("dotted", "dashed", "solid", "dashed", "dotted"), 
             linewidth = 1) +
  expand_limits(y = c(-4.5, 7.5)) +
  labs(y = "std Pearson res.", x = "xb",
       title = "Betaplasma - Full model: Negbin model (Std Pearson Residuals)")


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
Fullscope_Model1 <- glm.nb(betaplasma ~ (bmi + age + fat + cholesterol + fiber + 
                                        alcohol + betadiet + smokstat + sex + vituse)^2,data = data, control = glm.control(maxit=1000))
Fullscope_Model <- glm.nb(betaplasma ~ bmi + age + fat + cholesterol + fiber + 
                                          alcohol + betadiet + smokstat + sex + vituse,data = data, control = glm.control(maxit=1000))
summary(Fullscope_Model)
vars <- c("bmi", "age", "fat", "cholesterol", "fiber", 
          "alcohol", "betadiet", "smokstat", "sex", "vituse")

# Create formula string with main effects and all 2-way interactions
scope_formula <- as.formula(
  paste("~ (", paste(vars, collapse = " + "), ")^2")
)

nullmodel <- glm.nb(betaplasma ~ 1, data = data)

n <- nrow(data)  # sample size
forward_model <- step(nullmodel,
                  scope = list(lower = ~1, upper = formula(Fullscope_Model1)),
                  direction = "forward",
                  trace = TRUE,
                  k = log(n))  # BIC instead of AIC

backward_model <- step(Fullscope_Model,
                       scope = list(lower = ~1, upper = formula(Fullscope_Model)),
                       direction = "backward",
                       trace = TRUE,
                       k = log(n))  # BIC

stepwise_model_n <- step(nullmodel,
                       scope = list(lower = ~1, upper = formula(Fullscope_Model1)),
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

step_model <- glm.nb(betaplasma ~(bmi + cholesterol + betadiet + vituse)^2, data = data)

forward_model <- step(nullmodel,
                      scope = list(lower = ~1, upper = formula(step_model)),
                      direction = "forward",
                      trace = TRUE,
                      k = log(n))  # BIC instead of AIC

backward_model <- step(step_model,
                       scope = list(lower = ~1, upper = formula(step_model)),
                       direction = "backward",
                       trace = TRUE,
                       k = log(n))  # BIC

stepwise_model_n <- step(nullmodel,
                         scope = list(lower = ~1, upper = formula(step_model)),
                         direction = "both",
                         trace = TRUE,
                         k = log(n))  # BIC

stepwise_model_f <- step(step_model,
                         scope = list(lower = ~1, upper = formula(step_model)),
                         direction = "both",
                         trace = TRUE,
                         k = log(n)) 



#All of them result in the same

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


data.frame(AIC(nullmodel, model_reduced, forward_model,fullmodel), 
           BIC(nullmodel, model_reduced, forward_model,fullmodel),
           D = c(nullmodel$deviance, model_reduced$deviance, forward_model$deviance,fullmodel$deviance),
           D0 = c(nullmodel$null.deviance, model_reduced$null.deviance, forward_model$null.deviance,fullmodel$null.deviance),
           p = c(nullmodel$df.null - nullmodel$df.residual,
                 model_reduced$df.null - model_reduced$df.residual,
                 forward_model$df.null - forward_model$df.residual,
                 fullmodel$df.null - fullmodel$df.residual)) |> 
  mutate(df.1 = NULL) -> collect.AICetc
collect.AICetc


collect.AICetc <- data.frame(
  AIC(nullmodel, model_reduced, forward_model, betaplasma_glmnb_bmi,fullmodel), 
  BIC(nullmodel, model_reduced, forward_model, betaplasma_glmnb_bmi,fullmodel),
  D = c(nullmodel$deviance, model_reduced$deviance, forward_model$deviance, betaplasma_glmnb_bmi$deviance,fullmodel$deviance),
  D0 = c(nullmodel$null.deviance, model_reduced$null.deviance, forward_model$null.deviance, betaplasma_glmnb_bmi$null.deviance,fullmodel$null.deviance),
  p = c(nullmodel$df.null - nullmodel$df.residual,
        model_reduced$df.null - model_reduced$df.residual,
        forward_model$df.null - forward_model$df.residual,
        betaplasma_glmnb_bmi$df.null - betaplasma_glmnb_bmi$df.residual,
        fullmodel$df.null - fullmodel$df.residual)
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

pearson_res <- residuals(model_reduced, type = "pearson")

# Add residuals and standardized residuals to dataframe
betaplasma_pred <- betaplasma_pred |> mutate(
  v = v,
  devres = devres,
  std.devres = devres / sqrt(1 - v),
  pearsonres = pearson_res,
  std.pearsonres = pearson_res / sqrt(1 - v)
)

ggplot(betaplasma_pred, aes(x = xb.fit)) +
  geom_point(aes(y = std.devres), size = 2) +
  geom_hline(yintercept = c(-3, -2, 0, 2, 3), 
             linetype = c("dotted", "dashed", "solid", "dashed", "dotted"), 
             linewidth = 1) +
  expand_limits(y = c(-4.5, 7.5)) +
  labs(y = "std dev.res", x = "xb",
       title = "Betaplasma - Final: Negbin model (Std Dev Residuals)")

# Plot using standardized Pearson residuals
ggplot(betaplasma_pred, aes(x = xb.fit)) +
  geom_point(aes(y = std.pearsonres), size = 2) +
  geom_hline(yintercept = c(-3, -2, 0, 2, 3), 
             linetype = c("dotted", "dashed", "solid", "dashed", "dotted"), 
             linewidth = 1) +
  expand_limits(y = c(-4.5, 7.5)) +
  labs(y = "std Pearson res.", x = "xb",
       title = "Betaplasma - Final: Negbin model (Std Pearson Residuals)")



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
       title = "Betaplasma: Poisson model (Std. Dev residuals")

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
pchisq(D_diff, 1, lower.tail = FALSE) #P-value 0


#Prediction intervals for just BMI



bmi_seq <- seq(min(new_data$bmi), max(new_data$bmi), length.out = 100)

pred_grid <- expand.grid(
  bmi = bmi_seq,
  cholesterol = median(new_data$cholesterol, na.rm = TRUE),
  betadiet = median(new_data$betadiet, na.rm = TRUE),
  vituse_reduced = levels(new_data$vituse_reduced)
)

boot_results <- boot.nb(
  model = model_reduced,
  odata = new_data,
  newdata = pred_grid,
  N = 1000,           # Adjust for speed vs accuracy
  p = 0.95
)

pred_grid$pred <- boot_results$pred
pred_grid$lower <- boot_results$lower
pred_grid$upper <- boot_results$upper


ggplot(pred_grid, aes(x = bmi, color = vituse_reduced)) +
  geom_line(aes(y = pred)) +
  geom_line(aes(y = lower), linetype = "dashed") +
  geom_line(aes(y = upper), linetype = "dashed") +
  facet_wrap(~ vituse_reduced) +
  labs(y = "Predicted betaplasma", title = "Bootstrap Prediction Intervals")


# Combine original data with model predictions (mean fit)
new_data$fit <- predict(model_reduced, type = "response")
# Calculate confidence intervals for the mean prediction
link <- predict(model_reduced, type = "link", se.fit = TRUE)
critval <- qnorm(0.975)
new_data$mu.lwr <- model_reduced$family$linkinv(link$fit - critval * link$se.fit)
new_data$mu.upr <- model_reduced$family$linkinv(link$fit + critval * link$se.fit)



# Compute confidence intervals for the mean prediction
link <- predict(model_reduced, type = "link", se.fit = TRUE)
critval <- qnorm(0.975)
new_data$mu.lwr <- model_reduced$family$linkinv(link$fit - critval * link$se.fit)
new_data$mu.upr <- model_reduced$family$linkinv(link$fit + critval * link$se.fit)


cooks_d <- cooks.distance(model_reduced)

# Plot Cook's distance
plot(cooks_d, 
     type = "h", 
     main = "Cook's Distance", 
     ylab = "Cook's Distance", 
     xlab = "Observation", 
     col = "blue", 
     lwd = 2)

# Add a reference line for identifying influential observations
abline(h = 4 / length(cooks_d), col = "red", lty = 2)

# Optionally label points above the threshold
threshold <- 4 / length(cooks_d)
influential <- which(cooks_d > threshold)
text(influential, cooks_d[influential], labels = influential, pos = 3, cex = 0.8)


library(dplyr)
library(ggplot2)

# Calculate dynamic threshold
n <- nobs(model_reduced)            # Number of observations
p <- length(coef(model_reduced))   # Number of parameters (incl. intercept)
leverage_threshold <- 2 * p / n    # You can also try 3*p/n for a stricter cutoff

# Prepare data for plotting
diagnostics_df <- data %>%
  mutate(
    fitted = fitted(model_reduced),
    Dcook = cooks.distance(model_reduced),
    leverage = hatvalues(model_reduced),
    std_resid = rstandard(model_reduced),
    high_resid = abs(std_resid) > 2,
    high_leverage = leverage > leverage_threshold
  )

# Plot: Cook's Distance vs Fitted Values
ggplot(diagnostics_df, aes(x = fitted, y = Dcook)) +
  geom_point(aes(color = factor(high_leverage))) +
  geom_point(data = filter(diagnostics_df, high_leverage),
             aes(x = fitted, y = Dcook, color = "High leverage"),
             size = 3) +
  geom_point(data = filter(diagnostics_df, high_resid),
             aes(x = fitted, y = Dcook, color = "High residuals"),
             size = 3, shape = 21, fill = "orange", stroke = 1.2) +
  geom_hline(yintercept = 4 / n, linewidth = 1, linetype = "dashed") +
  labs(
    title = "Cook's Distance vs Fitted Values",
    subtitle = paste("High leverage threshold:", round(leverage_threshold, 4), 
                     " | Dashed line: Cook's D > 4/n"),
    x = "Fitted values", y = "Cook's Distance", color = "Highlight"
  ) +
  scale_color_manual(values = c("TRUE" = "green", 
                                "High leverage" = "green", 
                                "High residuals" = "orange")) +
  theme(legend.position = "top")







library(dplyr)
library(ggplot2)

# Calculate dynamic threshold
n <- nobs(model_reduced)            # Number of observations
p <- length(coef(model_reduced))   # Number of parameters (incl. intercept)
leverage_threshold <- 2 * p / n    # You can also try 3*p/n for a stricter cutoff

# Prepare data for plotting
diagnostics_df <- data %>%
  mutate(
    fitted = fitted(model_reduced),
    Dcook = cooks.distance(model_reduced),
    leverage = hatvalues(model_reduced),
    std_resid = rstandard(model_reduced),
    high_resid = abs(std_resid) > 2,
    high_leverage = leverage > leverage_threshold
  )

# Plot: Cook's Distance vs Fitted Values
ggplot(diagnostics_df, aes(x = fitted, y = Dcook)) +
  geom_point(aes(color = factor(high_leverage))) +
  geom_point(data = filter(diagnostics_df, high_leverage),
             aes(x = fitted, y = Dcook, color = "High leverage"),
             size = 3) +
  geom_point(data = filter(diagnostics_df, high_resid),
             aes(x = fitted, y = Dcook, color = "High residuals"),
             size = 3, shape = 21, fill = "orange", stroke = 1.2) +
  geom_hline(yintercept = 4 / n, linewidth = 1, linetype = "dashed") +
  labs(
    title = "Cook's Distance vs Fitted Values",
    subtitle = paste("High leverage threshold:", round(leverage_threshold, 4), 
                     " | Dashed line: Cook's D > 4/n"),
    x = "Fitted values", y = "Cook's Distance", color = "Highlight"
  ) +
  scale_color_manual(values = c("TRUE" = "green", 
                                "High leverage" = "green", 
                                "High residuals" = "orange")) +
  theme(legend.position = "top")



library(ggplot2)
library(MASS)
library(dplyr)

# Assuming model_reduced is your fitted glm.nb model
n <- nobs(model_reduced)

# Create a data frame with diagnostics
influence_data <- data.frame(
  leverage = hatvalues(model_reduced),
  fitted = fitted(model_reduced),
  cooksD = cooks.distance(model_reduced)
)

# Optional: add row numbers for labeling influential points
influence_data$obs <- seq_len(n)

# Define a dynamic threshold for high leverage
p <- length(coef(model_reduced))  # number of parameters
high_leverage_threshold <- 2 * p / n

# Plot: Fitted Values vs Leverage with bubble size for Cook's distance
ggplot(influence_data, aes(x = fitted, y = leverage)) +
  geom_point(aes(size = cooksD), alpha = 0.6) +
  geom_hline(yintercept = high_leverage_threshold, linetype = "dashed", color = "darkgrey") +
  scale_size_continuous(name = "Cook's Distance", range = c(1, 6)) +
  labs(title = "Fitted Values vs Leverage",
       subtitle = paste0("Dashed line: high leverage threshold (", round(high_leverage_threshold, 3), ")"),
       x = "Fitted Values",
       y = "Leverage") +
  theme_minimal()
