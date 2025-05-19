source("Code/Base.R")

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

anova(betaplasma_null,betaplasma_glmnb_bmi)
