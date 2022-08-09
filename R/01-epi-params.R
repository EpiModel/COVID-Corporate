
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


# Load libraries and data  ------------------------------------------------

library("tidyverse")

## contact data rounds 2-4
mix2 <- readRDS("data/input/df_merged_rd2.RDS")
mix3 <- readRDS("data/input/df_merged_rd3.RDS")
mix4 <- readRDS("data/input/df_merged_rd4.RDS")

# round 1 (match rounds 2-4 format)
contacts1 <- readRDS("data/input/df_contact_clean_rd1.RDS") %>%
  rename(id = part_id)
part1 <- readRDS("data/input/df_participant_clean_rd1.RDS") %>%
  rename(id = part_id)
mix1 <- part1 %>%
  select(id,age_cat) %>%
  rename(participant_age_2 = age_cat) %>%
  right_join(contacts1,by = "id") %>%
  rename(contact_age_3 = contact_age) %>%
  mutate(contact_location2 = case_when(cont_home == 1 ~ "Home",
                                       work == 1 ~ "Work",
                                       street == 1 |
                                         store == 1 |
                                         cont_otherhome == 1 |
                                         school == 1 |
                                         worship == 1 |
                                         leisure == 1 |
                                         transport == 1 |
                                         healthcare == 1 |
                                         gym == 1 |
                                         playground == 1 |
                                         loc_other == 1 ~ "Community")) %>%
  filter(!is.na(contact_location2)) %>%
  mutate(contact_location2 = factor(contact_location2,
                                    levels = c("Home","Work","Community"),
                                    ordered = TRUE)) %>%
  mutate(participant_age_2 = factor(participant_age_2,
                                    levels = c("20-29","30-39","40-49",
                                               "50-59","60+"),
                                    ordered = TRUE))

# create list with all rounds
mix.list <- list(mix1,mix2,mix3,mix4)

# Format data  ------------------------------------------------------------

## function to collapse data and estimate mean degree
summariseContacts <- function(df) {

  df <- df %>%
    select(id,participant_age_2,contact_age_3,contact_location2) %>%
    rename(p_age = participant_age_2,
           location = contact_location2) %>%
    mutate(contact_age = ifelse(contact_age_3 %in% c("0-9","10-19"),
                                "<20",as.character(contact_age_3))) %>%
    mutate(contact_age = factor(contact_age)) %>%
    mutate(p_age = factor(p_age)) %>%
    mutate(location = factor(location)) %>%
    select(-contact_age_3) %>%
    group_by(p_age) %>%
    mutate(n_participants = n_distinct(id)) %>%  # total number of participants by age group
    group_by(location,p_age,contact_age,n_participants,.drop = FALSE) %>%
    summarise(n_contacts = n()) %>%  # total number of contacts reported over 2 days
    ungroup() %>%
    complete(location,p_age,contact_age) %>%
    group_by(location,p_age) %>%
    fill(n_participants,.direction = "downup") %>%
    mutate(degree = (n_contacts/n_participants)/2) %>%   # average mean degree for 1 day
    mutate(total_degree = sum(degree)) %>%  # total degree by location and participant age
    mutate(prop_degree = degree/total_degree)  # age distribution of contacts by location and participant age

}

## apply function to list of contact data
mix.list <- lapply(mix.list, summariseContacts)

## split up data sets for each round
mix1.sum <- mix.list[[1]] %>%
  rename(n_participants = n_participants,
         n_contacts = n_contacts,
         degree = degree,
         total_degree = total_degree,
         prop_degree = prop_degree) %>%
  mutate(round = 1) %>%
  ungroup()

mix2.sum <- mix.list[[2]] %>%
  rename(n_participants = n_participants,
         n_contacts = n_contacts,
         degree = degree,
         total_degree = total_degree,
         prop_degree = prop_degree) %>%
  mutate(round = 2) %>%
  ungroup()

mix3.sum <- mix.list[[3]] %>%
  rename(n_participants = n_participants,
         n_contacts = n_contacts,
         degree = degree,
         total_degree = total_degree,
         prop_degree = prop_degree) %>%
  mutate(round = 3) %>%
  ungroup()

mix4.sum <- mix.list[[4]] %>%
  rename(n_participants = n_participants,
         n_contacts = n_contacts,
         degree = degree,
         total_degree = total_degree,
         prop_degree = prop_degree) %>%
  mutate(round = 4) %>%
  ungroup()

## join all rounds to compare
mix.sum <- mix1.sum %>%
  add_row(mix2.sum) %>%
  add_row(mix3.sum) %>%
  add_row(mix4.sum)

## plot degrees to compare across rounds -----------------------------------
library(ggplot2)
library(NatParksPalettes)

# Total degree by participant age group and location
mix.sum.total <- mix.sum %>%
  select(location,p_age,total_degree,round) %>%
  distinct()

ggplot(data = mix.sum.total, mapping = aes(x = as.factor(round),
                                           y = total_degree,
                                           group = p_age,
                                           color = p_age))+
  geom_line(size = 1)+
  facet_grid(cols = vars(location))+
  scale_color_manual(values=natparks.pals("Yellowstone", 5))+
  theme_bw()+
  labs(x = "Round", y = "Total degree", color = "Participant age")
ggsave(filename = "out/total_degree.png",height = 3,width = 6)

# Degree by participant age, contact age, and location
ggplot(data = mix.sum %>%
         filter(contact_age != "40-59") %>%
         mutate(contact_age = factor(contact_age,
                                     levels = c("<20","20-29","30-39",
                                                "40-49","50-59","60+"))),
       mapping = aes(x = round,
                     y = degree,
                     color = contact_age,
                     group = contact_age))+
  geom_line(size = 0.8)+
  facet_grid(cols = vars(location),
             rows = vars(p_age))+
  scale_color_manual(values=natparks.pals("Yellowstone", 6))+
  theme_bw()+
  labs(x = "Round", y = "Degree", color = "Contact age")
ggsave(filename = "out/degree_by_age.png",height = 7,width = 7)

# Distribution of contact age groups by participant age and location
ggplot(data = mix.sum %>%
         filter(round != 1) %>%
         filter(contact_age != "40-59") %>%
         mutate(contact_age = factor(contact_age,
                                     levels = c("<20","20-29","30-39",
                                                "40-49","50-59","60+"))),
       mapping = aes(x = round,
                     y = prop_degree,
                     fill = contact_age,
                     group = contact_age))+
  geom_area()+
  facet_grid(cols = vars(location),
             rows = vars(p_age))+
  scale_fill_manual(values=natparks.pals("Yellowstone", 6))+
  theme_bw()+
  labs(x = "Round", y = "Proportion of contacts", fill = "Contact age")
ggsave(filename = "out/prop_degree_by_age.png",height = 7,width = 7)
