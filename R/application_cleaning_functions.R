# time.lag <- as.Date("2006/01/01") - as.Date("2000/01/01")
# assume X is the cleaned test, this requires the knowledge of outcomes
# default X is the time window

lag_pair_visit_label <- function(X,time.lag, time.window = as.Date("2001/01/01") - as.Date("2000/01/01")){

  X.out <- X %>%
    group_by(id) %>%
    arrange(date, by_group = T) %>%
    mutate(date_diff = date - min(date))

  # only take first day or time lagged day
  X.out <- X.out %>%
    filter((date_diff >= time.lag & date_diff <= time.lag + time.window) | date_diff == 0)

  X.out <- X.out %>%
    group_by(id) %>%
    arrange(date, by_group = T) %>%
    mutate(visit = row_number())

  X.out <- X.out %>%
    filter(visit <= 2)

  X.out <- X.out %>%
    group_by(id) %>%
    mutate(num_visits = max(visit))
  X.out <- X.out %>% filter(num_visits == 2)
  return(X.out)
}



