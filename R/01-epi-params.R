
##
## COVID-19 Corporate Network Model
## Model Parameters
##
## Authors: Samuel M. Jenness
## Date: February 2021
##

# Mortality Rates ---------------------------------------------------------

# Rates per 100,000 for age groups: <1, 1-4, 5-9, 10-14, 15-19, 20-24, 25-29,
#                                   30-34, 35-39, 40-44, 45-49, 50-54, 55-59,
#                                   60-64, 65-69, 70-74, 75-79, 80-84, 85+
# source: https://www.statista.com/statistics/241572/death-rate-by-age-and-sex-in-the-us/
mortality_rate <- c(588.45, 24.8, 11.7, 14.55, 47.85, 88.2, 105.65, 127.2,
                    154.3, 206.5, 309.3, 495.1, 736.85, 1051.15, 1483.45,
                    2294.15, 3642.95, 6139.4, 13938.3)
# rate per person, per day
mr_pp_pd <- mortality_rate / 1e5 / 365

# Build out a mortality rate vector
age_spans <- c(1, 4, rep(5, 16), 1)
mr_vec <- rep(mr_pp_pd, times = age_spans)


# Epi Curve ---------------------------------------------------------------



# CorpMix Data Analysis ---------------------------------------------------

library("tidyverse")
options(tibble.print_max = 1000)

# overall mean degree
md <- readRDS("data/input/df_contact_by_location.RDS")
mix <- readRDS("data/input/df_contact_proportion_by_location.RDS")


# contruct processed mean degree file
md

# categorize all non-work and home as community
md2 <- md
md2$contact_location <- ifelse(md2$contact_location != "Work" &
                                 md2$contact_location != "Home",
                               "Community", md2$contact_location)
md2

# sum up contacts by grouped locations
md3 <- summarise(
  group_by(md2, contact_location, age_participant),
  contacts = sum(count_contacts)
)

# Add counts per age group
md3$n <- rep(md2$n_age[md2$contact_location == "Work"], 3)
md3

# calculate overall counts (for edges terms)
md3.oa <- data.frame(
  contact_location = c("Community", "Home", "Work"),
  age_participant = rep("Overall", 3),
  contacts = c(sum(md3$contacts[md3$contact_location == "Work"]),
               sum(md3$contacts[md3$contact_location == "Home"]),
               sum(md3$contacts[md3$contact_location == "Community"])
  ),
  n = rep(sum(md3$n[md3$contact_location == "Work"]), 3)
)
md3.oa

md3 <- rbind(md3, md3.oa)
print(md3, 100)

# add per-capita rates
md3$rate <- md3$contacts / md3$n
md3

# Add confidence intervals on the rates
md3 <- as.data.frame(md3)
md3$rate.lcl <- rep(NA, nrow(md3))
md3$rate.ucl <- rep(NA, nrow(md3))
for (i in 1:nrow(md3)) {
  md3$rate.lcl[i] <- poisson.test(md3$contacts[i])$conf.int[1] / md3$n[i]
  md3$rate.ucl[i] <- poisson.test(md3$contacts[i])$conf.int[2] / md3$n[i]
}
md3

saveRDS(md3, file = "data/input/cm2.meandegree.rds")


mix



mix

# Population proportions from census data --------------------------------------

library(dplyr)

data <- read.csv("~/git/COVID-Corporate/data/input/ACS Census Data 2019.csv")

data <- data[1:20, ]
data <- data[-2, ]

data <- data %>%
  rename(estimate = United.States..Total..Estimate)

data <- data %>%
  mutate(estimate = gsub(",", "", data$estimate)) %>%
  mutate(estimate = as.numeric(estimate))

total <- data$estimate[1]

data <- data[-1, ]

eightfive <- data %>%
  slice_tail(n = 1) %>%
  mutate(per_year = estimate) %>%
  mutate(age = 85)

data <- data %>%
  mutate(per_year = estimate / 5) %>%
  slice(rep(1:(n()-1), each = 5)) %>%
  mutate(age = 0:84)

data <- rbind(data, eightfive) %>%
  mutate(prop = (per_year / total))

pop_prop <- data %>%
  select(age, prop)

saveRDS(pop_prop, file = "data/input/pop_prop.rds")

