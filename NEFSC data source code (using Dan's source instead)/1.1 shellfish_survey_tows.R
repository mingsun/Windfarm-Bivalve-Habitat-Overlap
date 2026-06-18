library(tidyverse)
library(stringr)

# 1. surey effort loss ----

## 1.1 acquire full survey station data ----

# trim down the data for the MAB area (LON < -69, LAT < 41.7)
station.df <- read.csv("../data/shellfish survey/22565_UNION_FSCS_SVSTA.csv") |>
  filter(DECDEG_BEGLON <= -69, DECDEG_BEGLAT <= 41.7) |>
  mutate(ID.temp = paste(CRUISE6, STATION, sep = "."), .before = 1) |>
  mutate(YEAR = substr(CRUISE6, 1, 4))

# quick plot to check the spatial range
ggplot(station.df, aes(x = DECDEG_BEGLON, y = DECDEG_BEGLAT)) +
  geom_point()


## 1.2 overlapping stations ----

# this step is done in Arcgis
# dredging tows at 3 knots for 5 minutes, nominal tow distance 154m
# so we create a 154m buffer for all station, just to keep a distance from the OWF leased area

overlay.df <- read.csv(
  "../results/tow data/Shellfish_survey_stations_overlay.csv"
) |>
  select(
    CRUISE6,
    STRATUM,
    TOW,
    STATION,
    ID,
    AREA,
    SVVESSEL,
    SVGEAR,
    EST_YEAR,
    DECDEG_BEGLAT,
    DECDEG_BEGLON
  ) |>
  mutate(ID.temp = paste(CRUISE6, STATION, sep = "."))

# save the tow list outside WEA
station.WEE.df <- station.df |>
  filter(!ID.temp %in% overlay.df$ID.temp)

# quick plot to check the spatial range
ggplot(station.WEE.df, aes(x = DECDEG_BEGLON, y = DECDEG_BEGLAT)) +
  geom_point()

write.csv(
  station.WEE.df,
  "../results/tow data/Shellfish_WEE_stations.csv",
  row.names = FALSE
)


## 1.4 survey effort loss by year in retrospective ----
overlay.no.df <- overlay.df %>%
  group_by(EST_YEAR) %>%
  summarise(overlay.n = length(TOW))

station.no.df <- station.df %>%
  group_by(EST_YEAR) %>%
  summarise(total.n = length(TOW))

overlay.prop.df <- merge(overlay.no.df, station.no.df)
remove(overlay.no.df, station.no.df)
overlay.prop.df$ratio <- overlay.prop.df$overlay.n / overlay.prop.df$total.n

se.plot <- ggplot(overlay.prop.df) +
  geom_line(aes(x = EST_YEAR, y = ratio)) +
  geom_point(aes(x = EST_YEAR, y = ratio)) +
  scale_x_continuous(breaks = unique(overlay.prop.df$EST_YEAR)) +
  labs(x = "YEAR", y = "proportion of survey effort loss") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = -0.5))

png(
  "../plot/appendix/shellfish_survey_effort_loss.png",
  width = 10,
  height = 3,
  units = 'in',
  res = 800
)
print(se.plot)
dev.off()

# ---------------------------------------------------------------------------------------- #

# 2 sample loss (abundance) by year ----

## 2.1 total abundance ranking first ----

# use NEFSC data source here because it
species.abundance.df <- read.csv(
  "../data/shellfish survey/22565_UNION_FSCS_SVCAT.csv"
) |>
  filter(is.na(ID) == FALSE & is.na(SCIENTIFIC_NAME) == FALSE) |>
  group_by(SCIENTIFIC_NAME, SVSPP) |>
  summarize(abun = sum(EXPCATCHNUM, na.rm = TRUE)) |>
  ungroup() |>
  arrange(desc(abun))

# decide to use the following species for this study
# ocean quahog (409), surfclam (403)
# sea scallop save it for the scallop survey

# quick chat to see if these species are spatially balanced
cat.df <- read.csv("../data/shellfish survey/22565_UNION_FSCS_SVCAT.csv") |>
  filter(SVSPP %in% c(409, 403)) |>
  mutate(ID.temp = paste(CRUISE6, STATION, sep = ".")) |>
  select(ID.temp, SVSPP, EXPCATCHNUM, EXPCATCHWT) |>
  mutate(SPECIES = ifelse(SVSPP == 409, "Atlantic Surfclam", "Ocean Quahog")) |>
  filter(ID.temp %in% station.df$ID.temp) |>
  left_join(station.df, by = "ID.temp") |>
  select(
    ID.temp,
    YEAR,
    CRUISE6,
    STRATUM,
    STATION,
    SPECIES,
    EXPCATCHNUM,
    EXPCATCHWT,
    DECDEG_BEGLON,
    DECDEG_BEGLAT,
    BOTTEMP,
    AVGDEPTH
  )

ggplot(cat.df, aes(x = DECDEG_BEGLON, y = DECDEG_BEGLAT, color = EXPCATCHNUM)) +
  geom_point() +
  facet_wrap(. ~ SPECIES)

write.csv(
  cat.df,
  "../results/tow data/Shellfis_catch_by_tow.csv",
  row.names = FALSE
)


## 2.2 calcualte lost samples ----

cat.lost.df <- cat.df |>
  transform(overlay = ifelse(ID.temp %in% overlay.df$ID.temp, "y", "n")) |>
  group_by(SPECIES, YEAR, overlay) |>
  summarize(abun = sum(EXPCATCHNUM, na.rm = TRUE)) |>
  ungroup() |>
  spread(overlay, abun) |>
  mutate_at(c('y', 'n'), ~ replace_na(., 0)) |>
  mutate(total = y + n, ratio = y / total)

# write.csv(cat.lost.df, "results/AS_OQ_sample_loss_all_species.csv", row.names = FALSE)

# quick plot to see the sample loss by year for all species

cat.lost.plot <- ggplot(cat.lost.df) +
  geom_line(aes(x = as.numeric(YEAR), y = ratio)) +
  geom_point(aes(x = as.numeric(YEAR), y = ratio)) +
  geom_hline(yintercept = 0, color = "grey") +
  scale_x_continuous(breaks = as.numeric(unique(cat.df$YEAR))) +
  facet_wrap(. ~ SPECIES, ncol = 1) +
  labs(x = "YEAR", y = "proportion of samples loss") +
  theme(axis.text.x = element_text(angle = -45, vjust = -0.5)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -45, vjust = -0.5))
