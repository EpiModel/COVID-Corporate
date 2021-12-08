
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
# md <- readRDS("data/input/df_contact_by_location.RDS")
mix <- readRDS("data/input/df_contact_proportion_by_location.RDS")

# ## construct processed mean degree file ##
# md
#
# # categorize all non-work and home as community
# md2 <- md
# md2$contact_location <- ifelse(md2$contact_location != "Work" &
#                                  md2$contact_location != "Home",
#                                "Community", md2$contact_location)
# md2
#
# # recatgorize age
# md2$age_participant <- ifelse(md2$age_participant %in% c("60-69", "70+"),
#                               "60+", md2$age_participant)
# md2
#
# # sum up contacts by grouped locations
# md3 <- summarise(
#   group_by(md2, contact_location, age_participant),
#   contacts = sum(count_contacts)
# )
# md3
#
# # Count up numbers per group
# n <- summarise(
#   group_by(md2[md2$contact_location == "Work", ], age_participant),
#   n = sum(n_age)
# )$n
#
# # Add counts per age group
# md3$n <- rep(n, 3)
# md3
#
# # calculate overall counts (for edges terms)
# md3.oa <- data.frame(
#   contact_location = c("Community", "Home", "Work"),
#   age_participant = rep("Overall", 3),
#   contacts = c(sum(md3$contacts[md3$contact_location == "Work"]),
#                sum(md3$contacts[md3$contact_location == "Home"]),
#                sum(md3$contacts[md3$contact_location == "Community"])
#   ),
#   n = rep(sum(md3$n[md3$contact_location == "Work"]), 3)
# )
# md3.oa
#
# md3 <- rbind(md3, md3.oa)
# md3
#
# # add per-capita rates
# md3$rate <- md3$contacts / md3$n
# md3
#
# # Add confidence intervals on the rates
# # md3 <- as.data.frame(md3)
# md3$rate.lcl <- rep(NA, nrow(md3))
# md3$rate.ucl <- rep(NA, nrow(md3))
# for (i in 1:nrow(md3)) {
#   md3$rate.lcl[i] <- poisson.test(md3$contacts[i])$conf.int[1] / md3$n[i]
#   md3$rate.ucl[i] <- poisson.test(md3$contacts[i])$conf.int[2] / md3$n[i]
# }
# md3
#

## construct processed mixing data file
head(mix)

mix2 <- mix
mix2$contact_location <- ifelse(mix2$contact_location != "Work" &
                                  mix2$contact_location != "Home",
                                "Community", mix2$contact_location)
mix2
str(mix2)

mix2$participant_age <- as.character(mix2$participant_age)
str(mix2)

# recatgorize age
mix2$participant_age <- ifelse(mix2$participant_age %in% c("60-69", "70+"),
                               "60+", mix2$participant_age)
mix2$contact_age <- ifelse(mix2$contact_age %in% c("<1", "1-9", "10-19"),
                           "<20", mix2$contact_age)
mix2$contact_age <- ifelse(mix2$contact_age %in% c("60-69", "70+"),
                           "60+", mix2$contact_age)

mix2

# sum up contacts by grouped locations
mix3 <- summarise(
  group_by(mix2, contact_location, participant_age, contact_age),
  contacts = sum(count_contacts_d1)
)
mix3

mix[mix$contact_location == "Home" & mix$participant_age == "20-29", ]

mix3[mix3$contact_location == "Home" & mix3$participant_age == "20-29", ]
sum(mix3[mix3$contact_location == "Home" & mix3$participant_age == "20-29",]$contacts)

md3
md[md$contact_location == "Home", ]

# unclear why total contacts per age per location are inconsistent
# for example 712 in md versus 576 in mix for Home 20-29


md4 <- md3
md4 <- rename(md4, participant_age = age_participant)

# need to merge with md because need total contacts per age group per location

mix4 <- left_join(mix3, md4, by = c("contact_location", "participant_age"))
mix4

mix4 <- select(mix4, -n, -rate, -rate.lcl, -rate.ucl)
mix4 <- rename(mix4, contacts.grp = contacts.x, contacts.all = contacts.y)
mix4



mix4$grp.prop <- mix4$contacts.grp / mix4$contacts.all
mix4

# Add confidence intervals on the proportions
mix4$grp.prop.lcl <- rep(NA, nrow(mix4))
mix4$grp.prop.ucl <- rep(NA, nrow(mix4))
for (i in 1:nrow(mix4)) {
  mix4$grp.prop.lcl[i] <-
    binom.test(mix4$contacts.grp[i], mix4$contacts.all[i])$conf.int[1]
  mix4$grp.prop.ucl[i] <-
    binom.test(mix4$contacts.grp[i], mix4$contacts.all[i])$conf.int[2]
}
mix4

# saveRDS(mix4, file = "data/input/cm2.mixing.rds")

---

library(tidyverse)

mix <- readRDS("data/input/df_contact_proportion_by_location.RDS")

mix1 <- mix

# Clean up base dataset

mix1$contact_location <- ifelse(mix1$contact_location == "Work", "Work",
                                ifelse(mix1$contact_location == "Home", "Home",
                                       "Community"))

mix1$participant_age <- as.character(mix1$participant_age)

mix1$participant_age <- ifelse(mix1$participant_age %in% c("60-69", "70+"),
                               "60+", mix1$participant_age)
mix1$contact_age <- ifelse(mix1$contact_age %in% c("<1", "1-9", "10-19"),
                           "<20", mix1$contact_age)
mix1$contact_age <- ifelse(mix1$contact_age %in% c("60-69", "70+"),
                           "60+", mix1$contact_age)
mix2 <- mix1 %>%
  select(contact_location, participant_age, contact_age, count_contacts_d1,
         count_contacts_d2, n_age_d1, n_age_d2)

#Create a complete age combination dataset

contact_age <- c("<20", "20-29", "30-39", "40-49", "50-59",
                 "60+")
participant_age <- c(rep("20-29", 6), rep("30-39", 6), rep("40-49", 6),
                     rep("50-59", 6), rep("60+", 6))
contact_location <- c(rep("Community", 30), rep("Work", 30), rep("Home", 30))

mixBase <- data.frame(contact_location = contact_location,
                      participant_age = rep(participant_age, 3),
                      contact_age = rep(contact_age, 5))

mix3 <- full_join(mix2, mixBase) %>%
  arrange(contact_location, participant_age, contact_age)

# Clean up missing data

mix3$count_contacts_d1[is.na(mix3$count_contacts_d1)] <- 0
mix3$count_contacts_d2[is.na(mix3$count_contacts_d2)] <- 0

#Needs to be made reproducible

mix3 <- mix3 %>%
  mutate(n_age_d1 = case_when(
    participant_age == "20-29" ~ 632,
    participant_age == "30-39" ~ 324,
    participant_age == "40-49" ~ 82,
    participant_age == "50-59" ~ 38,
    participant_age == "60+" ~ 20
  )) %>%
  mutate(n_age_d2 = case_when(
    participant_age == "20-29" ~ 497,
    participant_age == "30-39" ~ 292,
    participant_age == "40-49" ~ 72,
    participant_age == "50-59" ~ 34,
    participant_age == "60+" ~ 18
  ))

mix3 <- mix3 %>%
  mutate(n_age_total = n_age_d1 + n_age_d2)

mix4 <- mix3 %>%
  group_by(contact_location, participant_age, contact_age, n_age_total) %>%
  summarise(count_contacts_d1 = sum(count_contacts_d1),
            count_contacts_d2 = sum(count_contacts_d2)) %>%
  mutate(contact_count = count_contacts_d1 + count_contacts_d2) %>%
  mutate(contact_rate = contact_count / n_age_total)

mix5 <- mix4 %>%
  group_by(contact_location, participant_age, n_age_total) %>%
  summarise(contact_count = sum(contact_count),
            contact_rate = sum(contact_rate)) %>%
  mutate(contact_age = "all")


mix6 <- bind_rows(mix5, mix4) %>%
  mutate(contact_rate.lcl = NA,
         contact_rate.ucl = NA) %>%
  arrange(contact_location, participant_age) %>%
  select(contact_location, participant_age, contact_age, n_age_total,
         contact_count, contact_rate, contact_rate.lcl, contact_rate.ucl)

for (i in 1:nrow(mix6)) {
  mix6$contact_rate.lcl[i] <-
    poisson.test(mix6$contact_count[i])$conf.int[1] /
    mix6$n_age_total[i]

  mix6$contact_rate.ucl[i] <-
    poisson.test(mix6$contact_count[i])$conf.int[2] /
    mix6$n_age_total[i]
}

# Final data set with overall contact rate defined for each location-specific
# participant/contact combo and overall contact rate for age group for each
# location

mix7 <- mix6 %>%
  select(contact_location, participant_age, contact_age,
         n_age_total, contact_count, contact_rate,
         contact_rate.lcl, contact_rate.ucl)


#-------------------------------------------------------------------------------

mix8 <- mix3 %>%
  group_by(participant_age, contact_age, n_age_d1, n_age_d2, n_age_total,
           contact_location) %>%
  summarize(count_contacts_d1_age = sum(count_contacts_d1),
            count_contacts_d2_age = sum(count_contacts_d2)) %>%
  mutate(
    count_contacts_total_age = count_contacts_d1_age + count_contacts_d2_age)

mix9 <- mix8 %>%
  group_by(participant_age, n_age_total, contact_location) %>%
  summarize(count_contacts_by_age = sum(count_contacts_total_age))

mix10 <- full_join(mix8, mix9) %>%
  rename(total_contacts_age = count_contacts_by_age,
         pair_contacts_total = count_contacts_total_age) %>%
  mutate(prop_pair_contacts.age = pair_contacts_total / total_contacts_age)

mix11 <- mix10 %>%
  ungroup() %>%
  select(participant_age, contact_age, prop_pair_contacts.age,
         pair_contacts_total, total_contacts_age, contact_location)

mix11$prop_pair_contacts.age_lcl <- rep(NA, nrow(mix11))
mix11$prop_pair_contacts.age_ucl <- rep(NA, nrow(mix11))

for (i in 1:nrow(mix11)) {
  mix11$prop_pair_contacts.age_lcl[i] <-
    binom.test(mix11$pair_contacts_total[i],
               mix11$total_contacts_age[i])$conf.int[1]
  mix11$prop_pair_contacts.age_ucl[i] <-
    binom.test(mix11$pair_contacts_total[i],
               mix11$total_contacts_age[i])$conf.int[2]
}

# Final data set with proportion of total contacts for each age group that are
# within each contact age group

mix12 <- mix11 %>%
  select(contact_location, participant_age, contact_age, prop_pair_contacts.age,
         prop_pair_contacts.age_lcl, prop_pair_contacts.age_ucl)

mix12 %>%
  subset(contact_location == "Home" & participant_age == "40-49") %>%
  summarize(total = sum(prop_pair_contacts.age))
