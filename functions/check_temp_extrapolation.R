# check extrapolation of temperature values beyond training range
calc_temp_extrapolation <- function(temp_rast, temp_SI, month_id) {
  # training range
  temp_min <- min(temp_SI$mid, na.rm = TRUE)
  temp_max <- max(temp_SI$mid, na.rm = TRUE)

  # extract raster values
  temp_vals <- values(temp_rast)
  temp_vals <- temp_vals[is.finite(temp_vals)]

  # counts
  n_total <- length(temp_vals)
  n_below <- sum(temp_vals < temp_min)
  n_above <- sum(temp_vals > temp_max)
  n_within <- sum(temp_vals >= temp_min & temp_vals <= temp_max)

  # percentages
  pct_below <- n_below / n_total * 100
  pct_above <- n_above / n_total * 100
  pct_within <- n_within / n_total * 100

  # return one-row dataframe
  data.frame(
    month = month_id,
    n_total = n_total,
    n_within = n_within,
    n_below = n_below,
    n_above = n_above,
    pct_within = pct_within,
    pct_below = pct_below,
    pct_above = pct_above
  )
}
