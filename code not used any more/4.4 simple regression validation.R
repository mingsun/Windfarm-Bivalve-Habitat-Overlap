# surfclam
data <- read.csv("../results/HSI/surfclam/SC_tows_with_composite_SI.csv") # slope = 0.575, R2 = 0.167, p < 0.001

# quahog
data <- read.csv("../results/HSI/quahog/OQ_tows_with_composite_SI.csv") # slope = 0.158, R2 = 0.108, p < 0.001

# scallop
data <- read.csv("../results/HSI/scallop/SS_tows_with_composite_SI.csv") # slope = 0.239, R2 = 0.108, p < 0.001


data <- data |>
  mutate(log.N = log1p(NperMsq)) |>
  dplyr::select(NperMsq, log.N)


# fit linear model
lm_fit <- lm(HSI_weighted ~ log.N, data = data)

# extract stats
slope <- coef(lm_fit)[2]
r2 <- summary(lm_fit)$r.squared
pval <- summary(lm_fit)$coefficients[2, 4]

# format p-value nicely
p_text <- ifelse(pval < 0.001, "p < 0.001", paste0("p = ", round(pval, 3)))

# plot
ggplot(data, aes(x = NperMsq, y = HSI_weighted)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = paste0(
      "Slope = ",
      round(slope, 3),
      "\nR² = ",
      round(r2, 3),
      "\n",
      p_text
    ),
    hjust = 1.1,
    vjust = 1.5,
    size = 5
  ) +
  theme_bw()
