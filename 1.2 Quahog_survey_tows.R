library(tidyverse)

# 1. prepare data for GIS ----

catch_by_tow_df <- read.csv("../data/quahog data/OQtows_original.csv") |>
  filter(STRATUM %in% paste0(1:6, "Q")) |> # filter out the SVAtoSNE region
  mutate(ID.temp = paste(CRUISE6, STATION, sep = "."), .before = 1) |>
  mutate(LON = -as.double(LON)) |>
  filter(!is.na(NperMsq))

write.csv(
  catch_by_tow_df,
  "../data/quahog data/OQ_tows_MA.csv",
  row.names = FALSE
)

# ----------------------------------------------------------------- #

# 2. GIS intersection analyses ----

# done in ArcGIS

# ----------------------------------------------------------------- #

# 3. survey effort loss ----

overlay.df <- read.csv("../data/quahog data/Quahog_tow_overlap.csv") |>
  mutate(ID.temp = paste(CRUISE6, STATION, sep = ".")) # ID.temp was messed up in ArcGIS, so re-create it here

## 3.1 get the WEE dataset ----
catch_by_tow_WEE_df <- catch_by_tow_df |>
  filter(!ID.temp %in% overlay.df$ID.temp)

# quick plot to check the spatial range
ggplot(catch_by_tow_WEE_df, aes(x = LON, y = LAT)) +
  geom_point()

write.csv(
  catch_by_tow_WEE_df,
  "../data/quahog data/OQ_tows_MA_WEE.csv",
  row.names = FALSE
)


## 3.2 survey effort loss by year in retrospective ----
overlay.tow.df <- overlay.df %>%
  group_by(YR) %>%
  summarise(overlay.n = length(STATION))

total.tow.df <- catch_by_tow_df %>%
  group_by(YR) %>%
  summarise(total.n = length(STATION))

overlay.prop.df <- merge(overlay.tow.df, total.tow.df)
remove(overlay.tow.df, total.tow.df)
overlay.prop.df$ratio <- overlay.prop.df$overlay.n / overlay.prop.df$total.n

se.plot <- ggplot(overlay.prop.df) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_line(aes(x = YR, y = ratio)) +
  geom_point(aes(x = YR, y = ratio)) +
  scale_x_continuous(breaks = unique(overlay.prop.df$YR)) +
  labs(x = "YEAR", y = "proportion of survey effort loss") +
  ylim(c(0, 1)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = -0.5))

png(
  "../plot/appendix/quahog_effort_loss.png",
  width = 10,
  height = 3,
  units = 'in',
  res = 800
)
print(se.plot)
dev.off()

# ----------------------------------------------------------------- #

# 4. survey sample loss ----

samp_lost_df <- catch_by_tow_df |>
  transform(overlay = ifelse(ID.temp %in% overlay.df$ID.temp, "y", "n")) |>
  group_by(YR, overlay) |>
  summarize(abun = sum(NperMsq, na.rm = TRUE)) |>
  ungroup() |>
  spread(overlay, abun) |>
  mutate_at(c('y', 'n'), ~ replace_na(., 0)) |>
  mutate(total = y + n, ratio = y / total)

samp_lost_plot <- ggplot(samp_lost_df) +
  geom_line(aes(x = as.numeric(YR), y = ratio)) +
  geom_point(aes(x = as.numeric(YR), y = ratio)) +
  geom_hline(yintercept = 0, color = "grey") +
  scale_x_continuous(breaks = as.numeric(unique(samp_lost_df$YR))) +
  labs(x = "YEAR", y = "proportion of samples loss") +
  ylim(c(0, 1)) +
  theme(axis.text.x = element_text(angle = -45, vjust = -0.5)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = -0.5))

png(
  "../plot/appendix/quahog_sample_loss.png",
  width = 10,
  height = 3,
  units = 'in',
  res = 800
)
print(samp_lost_plot)
dev.off()

# ----------------------------------------------------------------- #
