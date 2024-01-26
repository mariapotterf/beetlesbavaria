

library(ggplot2)
library(ggeffects)
library(MASS)



set.seed(123) # For reproducibility
n <- 100 # Number of observations
df <- data.frame(
  predictor1 = runif(n, 1, 10),
  predictor2 = runif(n, 1, 5)
)
df$response_var <- rnbinom(n, size = 1, mu = 2 + df$predictor1 * 0.3 - df$predictor2 * 0.1)

# Fit a Negative Binomial model
model <- glm.nb(response_var ~ predictor1 + predictor2, data = df)

summary(model)

# Assuming 'model' is your glm.nb model
p1 <- ggpredict(model, terms = "predictor1")
p2 <- ggpredict(model, terms = "predictor2")

# Plot

ggplot(p1, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high)) +
  geom_line() +
  geom_ribbon(alpha = 0.1) +
  theme_minimal()

ggplot(p2, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high)) +
  geom_line() +
  geom_ribbon(alpha = 0.1) +
  theme_minimal()


