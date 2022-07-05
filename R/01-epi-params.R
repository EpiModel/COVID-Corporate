
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

# re-categorize location
mix2$contact_location <- ifelse(mix2$contact_location != "Work" &
                                  mix2$contact_location != "Home",
                                "Community", mix2$contact_location)

# re-categorize age
mix2$participant_age <- ifelse(mix2$participant_age %in% c("60-69", "70+"),
                               "60+", as.character(mix2$participant_age))
mix2$contact_age <- ifelse(mix2$contact_age %in% c("<1", "1-9", "10-19"),
                           "<20", mix2$contact_age)
mix2$contact_age <- ifelse(mix2$contact_age %in% c("60-69", "70+"),
                           "60+", mix2$contact_age)

# pull out number of participants per location, (new) age group, and day
n.age1 <- mix2 %>%
  select(participant_age,n_age_d1) %>%
  distinct() %>%
  filter(!is.na(n_age_d1)) %>%
  group_by(participant_age) %>%
  summarise(n_age_d1 = sum(n_age_d1))

n.age2 <- mix2 %>%
  select(participant_age,n_age_d2) %>%
  distinct() %>%
  filter(!is.na(n_age_d2)) %>%
  group_by(participant_age) %>%
  summarise(n_age_d2 = sum(n_age_d2))

# sum up contacts by grouped locations
# sum unique contacts across both days
mix3 <- mix2 %>%
  select(-starts_with("prop_"),-count_contacts_total_d1,-count_contacts_total_d2) %>%
  group_by(contact_location, participant_age, contact_age) %>%
  summarise(across(starts_with("count"),~ sum(.x, na.rm = TRUE))) %>%
  rowwise() %>%
  mutate(count_contacts = sum(count_contacts_d1,count_contacts_d2,na.rm = T)) %>%
  group_by(contact_location, participant_age) %>%
  mutate(count_contacts_total_d1 = sum(count_contacts_d1, na.rm = T),
         count_contacts_total_d2 = sum(count_contacts_d2, na.rm = T),
         count_contacts_total = sum(count_contacts, na.rm = T)) %>%
  left_join(n.age1,by = "participant_age") %>%
  left_join(n.age2,by = "participant_age") %>%
  rowwise() %>%
  # is this meaningful?? or just use d1?
  mutate(n_age = sum(n_age_d1,n_age_d2))

mix[mix$contact_location == "Home" & mix$participant_age == "20-29", ]

mix3[mix3$contact_location == "Home" & mix3$participant_age == "20-29", ]
sum(mix3[mix3$contact_location == "Home" & mix3$participant_age == "20-29",]$count_contacts)

mix3[mix3$contact_location == "Home", ]

# unclear why total contacts per age per location are inconsistent
# for example 712 in md versus 576 in mix for Home 20-29
# seems like number of participants also varies across datasets - is this supposed to be consistent?


# md4 <- mix3
# md4 <- rename(md4, participant_age = age_participant)
#
# # need to merge with md because need total contacts per age group per location
#
# mix4 <- left_join(mix3, md4, by = c("contact_location", "participant_age"))
# mix4
#
# mix4 <- select(mix4, -n, -rate, -rate.lcl, -rate.ucl)
# mix4 <- rename(mix4, contacts.grp = contacts.x, contacts.all = contacts.y)
# mix4
#
#
#
# mix4$grp.prop <- mix4$contacts.grp / mix4$contacts.all
# mix4
#
# # Add confidence intervals on the proportions
# mix4$grp.prop.lcl <- rep(NA, nrow(mix4))
# mix4$grp.prop.ucl <- rep(NA, nrow(mix4))
# for (i in 1:nrow(mix4)) {
#   mix4$grp.prop.lcl[i] <-
#     binom.test(mix4$contacts.grp[i], mix4$contacts.all[i])$conf.int[1]
#   mix4$grp.prop.ucl[i] <-
#     binom.test(mix4$contacts.grp[i], mix4$contacts.all[i])$conf.int[2]
# }
# mix4

# saveRDS(mix4, file = "data/input/cm2.mixing.rds")



mix4 <- mix3 %>%
  mutate(contact_rate = count_contacts / n_age_d1)

mix5 <- mix4 %>%
  group_by(contact_location, participant_age, n_age_d1) %>%
  summarise(count_contacts = sum(count_contacts),
            contact_rate = sum(contact_rate)) %>%
  mutate(contact_age = "Total")


mix6 <- bind_rows(mix5, mix4) %>%
  mutate(contact_rate.lcl = NA,
         contact_rate.ucl = NA) %>%
  arrange(contact_location, participant_age) %>%
  select(contact_location, participant_age, contact_age, n_age_d1,
         count_contacts, contact_rate, contact_rate.lcl, contact_rate.ucl)

for (i in 1:nrow(mix6)) {
  mix6$contact_rate.lcl[i] <-
    poisson.test(mix6$count_contacts[i])$conf.int[1] /
    mix6$n_age_d1[i]

  mix6$contact_rate.ucl[i] <-
    poisson.test(mix6$count_contacts[i])$conf.int[2] /
    mix6$n_age_d1[i]
}

# Summary statistics for ERGM model

mix.sum <- mix6 %>%
  select(-starts_with("contact_rate")) %>%
  filter(contact_age != "Total") %>%
  mutate(degree = (count_contacts/n_age_d1)/2) %>%
  group_by(contact_location,participant_age) %>%
  mutate(total_degree = sum(degree)) %>%
  mutate(prop_degree = degree/total_degree)


