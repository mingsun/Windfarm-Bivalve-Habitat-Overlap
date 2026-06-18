# parametric Gaussian-shaped curve
# SI ~ a * exp(-0.5 * ((mid - opt) / sigma)^2)
# where a is the maximum SI value, opt is the optimal value of the predictor variable, and sigma controls the width of the curve

# For environmental variables where the empirical SI response was expected to be unimodal,
# we developed a Gaussian response curve fitted to the binned SI values.
# The model estimated the maximum suitability, the environmental optimum, and the response width,
# thereby constraining the final SI relationship to a single bell-shaped curve.
# Fitted values were then rescaled to 0–1 and used as the final SI curve.

library(dplyr)
library(minpack.lm)
library(rlang)

fit_bell_SI <- function(SI_df, SI_col) {
  SI_col <- enquo(SI_col)

  fit_df <- SI_df |>
    transmute(
      mid = mid,
      SI = !!SI_col
    ) |>
    filter(
      !is.na(mid),
      !is.na(SI),
      is.finite(mid),
      is.finite(SI)
    )

  start_a <- max(fit_df$SI, na.rm = TRUE)
  start_opt <- fit_df$mid[which.max(fit_df$SI)]
  start_sigma <- sd(fit_df$mid, na.rm = TRUE)

  bell_fit <- nlsLM(
    SI ~ a * exp(-0.5 * ((mid - opt) / sigma)^2),
    data = fit_df,
    start = list(
      a = start_a,
      opt = start_opt,
      sigma = start_sigma
    ),
    lower = c(
      a = 0,
      opt = min(fit_df$mid, na.rm = TRUE),
      sigma = 0.001
    ),
    upper = c(
      a = 2,
      opt = max(fit_df$mid, na.rm = TRUE),
      sigma = Inf
    )
  )

  return(bell_fit)
}
