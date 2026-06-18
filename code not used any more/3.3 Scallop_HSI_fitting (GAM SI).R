library(tidyverse)
library(classInt)
library(mgcv)

data <- read.csv("../data/scallop data/SS_tows_MA_final.csv") |>
  filter(YEAR >= 2013) |> # this is more recent and in August
  filter(!is.na(NperMsq) & !is.na(TEMP) & !is.na(DEPTH) & !is.na(SED_PHI_MEAN))

# check for correlation among predictors
cor_df <- data %>% dplyr::select(TEMP, DEPTH, SED_PHI_MEAN)
cor_mat <- cor(cor_df, use = "complete.obs")
cor_mat # all < 0.3, so no strong correlation

# check for multicollinearity using VIF
library(usdm)
vif(cor_df) # all < 2, so no strong multicollinearity

remove(cor_df, cor_mat)

k = 6 # number of knots for GAM smoothing

# 1. Build Suitability Indices (SIs) ----

## 1.1 SI_sediment ----

# use Fisher's natural breaks to determine sediment bins, which minimizes within-bin variance and maximizes between-bin variance
sediment <- as.numeric(data$SED_PHI_MEAN)
sed_brks <- classIntervals(sediment, 20, style = "fisher")$brks
sed_brks[1] <- sed_brks[1] - .Machine$double.eps # ensure lowest value falls in first bin

# create bins
data <- data |>
  mutate(
    sediment_bins = cut(SED_PHI_MEAN, breaks = sed_brks, include.lowest = TRUE)
  )

sediment_SI <- data |>
  group_by(sediment_bins) |>
  summarize(abundance = mean(NperMsq, na.rm = TRUE), .groups = "drop") |>
  mutate(
    SI_sediment = {
      mn <- min(abundance, na.rm = TRUE)
      mx <- max(abundance, na.rm = TRUE)
      if (mx == mn) {
        rep(0, length(abundance))
      } else {
        (abundance - mn) / (mx - mn)
      }
    }
  ) |>
  mutate(
    mid = (sed_brks[-1] + sed_brks[-length(sed_brks)])[as.integer(
      sediment_bins
    )] /
      2
  )

# smooth SI across bin midpoints using GAM
sed_gam <- gam(SI_sediment ~ s(mid, k = k), data = sediment_SI)

# predict SI_sediment across bins and re-scale to 0-1
sediment_SI <- sediment_SI |>
  mutate(
    SI_sediment = predict(sed_gam, newdata = sediment_SI)
  ) |>
  mutate(
    SI_sediment = {
      mn <- min(SI_sediment, na.rm = TRUE)
      mx <- max(SI_sediment, na.rm = TRUE)
      if (mx == mn) {
        rep(0, length(SI_sediment))
      } else {
        (SI_sediment - mn) / (mx - mn)
      }
    }
  )

data <- left_join(
  data,
  dplyr::select(sediment_SI, sediment_bins, SI_sediment),
  by = "sediment_bins"
)

saveRDS(sediment_SI, "../results/HSI/scallop/scallop_sediment_SI.rds")

remove(sediment, sed_brks, sediment_SI, sed_gam)

## ------------------------------------------------------------- ##

## 1.2 depth ----

# use Fisher's natural breaks to determine sediment bins, which minimizes within-bin variance and maximizes between-bin variance
depth <- as.numeric(data$DEPTH)
depth_brks <- classIntervals(depth, 20, style = "fisher")$brks
depth_brks[1] <- depth_brks[1] - .Machine$double.eps # ensure lowest value falls in first bin

# create bins
data <- data |>
  mutate(
    depth_bins = cut(DEPTH, breaks = depth_brks, include.lowest = TRUE)
  )

depth_SI <- data |>
  group_by(depth_bins) |>
  summarize(abundance = mean(NperMsq, na.rm = TRUE), .groups = "drop") |>
  mutate(
    SI_depth = {
      mn <- min(abundance, na.rm = TRUE)
      mx <- max(abundance, na.rm = TRUE)
      if (mx == mn) {
        rep(0, length(abundance))
      } else {
        (abundance - mn) / (mx - mn)
      }
    }
  ) |>
  mutate(
    mid = (depth_brks[-1] + depth_brks[-length(depth_brks)])[as.integer(
      depth_bins
    )] /
      2
  )

# smooth SI across bin midpoints using GAM
depth_gam <- gam(SI_depth ~ s(mid, k = k), data = depth_SI)

# predict SI_depth across bins and re-scale to 0-1
depth_SI <- depth_SI |>
  mutate(
    SI_depth = predict(depth_gam, newdata = depth_SI)
  ) |>
  mutate(
    SI_depth = {
      mn <- min(SI_depth, na.rm = TRUE)
      mx <- max(SI_depth, na.rm = TRUE)
      if (mx == mn) {
        rep(0, length(SI_depth))
      } else {
        (SI_depth - mn) / (mx - mn)
      }
    }
  )


data <- left_join(
  data,
  dplyr::select(depth_SI, depth_bins, SI_depth),
  by = "depth_bins"
)

saveRDS(depth_SI, "../results/HSI/scallop/scallop_depth_SI.rds")

remove(depth, depth_brks, depth_SI, depth_gam)

## ------------------------------------------------------------- ##

## 1.3 temperature ----

# use Fisher's natural breaks to determine sediment bins, which minimizes within-bin variance and maximizes between-bin variance
temp <- as.numeric(data$TEMP)
temp_brks <- classIntervals(temp, 20, style = "fisher")$brks
temp_brks[1] <- temp_brks[1] - .Machine$double.eps # ensure lowest value falls in first bin

# create bins
data <- data |>
  mutate(
    temp_bins = cut(TEMP, breaks = temp_brks, include.lowest = TRUE)
  )

temp_SI <- data |>
  group_by(temp_bins) |>
  summarize(abundance = mean(NperMsq, na.rm = TRUE), .groups = "drop") |>
  mutate(
    SI_temp = {
      mn <- min(abundance, na.rm = TRUE)
      mx <- max(abundance, na.rm = TRUE)
      if (mx == mn) {
        rep(0, length(abundance))
      } else {
        (abundance - mn) / (mx - mn)
      }
    }
  ) |>
  mutate(
    mid = (temp_brks[-1] + temp_brks[-length(temp_brks)])[as.integer(
      temp_bins
    )] /
      2
  )

# smooth SI across bin midpoints using GAM
temp_gam <- gam(SI_temp ~ s(mid, k = k), data = temp_SI)

# predict SI_temp across bins and re-scale to 0-1
temp_SI <- temp_SI |>
  mutate(
    SI_temp = predict(temp_gam, newdata = temp_SI)
  ) |>
  mutate(
    SI_temp = {
      mn <- min(SI_temp, na.rm = TRUE)
      mx <- max(SI_temp, na.rm = TRUE)
      if (mx == mn) {
        rep(0, length(SI_temp))
      } else {
        (SI_temp - mn) / (mx - mn)
      }
    }
  )

data <- left_join(
  data,
  dplyr::select(temp_SI, temp_bins, SI_temp),
  by = "temp_bins"
)

saveRDS(temp_SI, "../results/HSI/scallop/scallop_temp_SI.rds")

remove(temp, temp_brks, temp_SI, temp_gam)

## ------------------------------------------------------------- ##

# ---------------------------------------------------------------------------------- #

# 2. Plot Individual SI graphs  ----

sed_levels <- levels(data$sediment_bins)
dep_levels <- levels(data$depth_bins)
tmp_levels <- levels(data$temp_bins)

sed_plot <- data %>%
  distinct(sediment_bins, SI_sediment) %>%
  rename(bin = sediment_bins, SI = SI_sediment) %>%
  mutate(var = "Sediment (phi)", bin = factor(bin, levels = sed_levels))

depth_plot <- data %>%
  distinct(depth_bins, SI_depth) %>%
  rename(bin = depth_bins, SI = SI_depth) %>%
  mutate(var = "Depth (m)", bin = factor(bin, levels = dep_levels))

temp_plot <- data %>%
  distinct(temp_bins, SI_temp) %>%
  rename(bin = temp_bins, SI = SI_temp) %>%
  mutate(var = "Temperature (°C)", bin = factor(bin, levels = tmp_levels))

plot_df <- bind_rows(sed_plot, depth_plot, temp_plot)

p <- ggplot(plot_df, aes(x = bin, y = SI, group = 1)) +
  geom_line() +
  geom_point(size = 1.5, alpha = 0.8) +
  facet_wrap(~var, scales = "free_x", ncol = 1) +
  labs(x = NULL, y = "Suitability index (0–1)") +
  # ggtitle(paste0("k=", k)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


remove(
  sed_levels,
  dep_levels,
  tmp_levels,
  sed_plot,
  depth_plot,
  temp_plot,
  plot_df
)

png(
  "../plot/appendix/scallop_SI.png",
  width = 6,
  height = 8,
  units = 'in',
  res = 800
)
print(p)
dev.off()

# quick save for now
write.csv(
  data,
  "../results/HSI/scallop/SS_tows_with_separate_SI.csv",
  row.names = FALSE
)

remove(data, k, p)

# ---------------------------------------------------------------------------------- #

# 3. Composit HSI  ----

library(gbm)
library(dismo)

# fit BRT using original environmental variables
# Gaussian is appropriate here because response is continuous abundance
data <- read.csv("../results/HSI/scallop/SS_tows_with_separate_SI.csv") # load seperate SI

# fit BRT using original environmental variables
# Gaussian is appropriate here because response is continuous abundance

set.seed(123)

brt_fit <- gbm.step(
  data = data,
  gbm.x = c("SED_PHI_MEAN", "DEPTH", "TEMP"),
  gbm.y = "NperMsq",
  family = "gaussian",
  tree.complexity = 3,
  learning.rate = 0.008,
  bag.fraction = 0.75,
)

# inspect variable relative influence
plot(brt_fit)
brt_fit$contributions # 0.01 - 37, 37, 25
brt_fit$n.trees

# extract and standardize weights
brt_weights <- brt_fit$contributions |>
  as.data.frame() |>
  transmute(
    var = var,
    weight = rel.inf / sum(rel.inf)
  )

brt_weights

# pull weights into named objects
w_sed <- brt_weights$weight[brt_weights$var == "SED_PHI_MEAN"]
w_dep <- brt_weights$weight[brt_weights$var == "DEPTH"]
w_tmp <- brt_weights$weight[brt_weights$var == "TEMP"]

# compute weighted HSI
data <- data |>
  mutate(
    HSI_weighted = w_sed * SI_sediment + w_dep * SI_depth + w_tmp * SI_temp
  )

summary(data$HSI_weighted)

write.csv(
  brt_weights,
  "../results/HSI/scallop/scallop_BRT_weights.csv",
  row.names = FALSE
)
write.csv(
  data,
  "../results/HSI/scallop/SS_tows_with_composite_SI.csv",
  row.names = FALSE
)

# ---------------------------------------------------------------------------------- #
