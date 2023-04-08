
##
## COVID-19 Corporate Network Model
## Network Parameters
##
## Authors: Maria Garcia Quesada
## Date: October 2022
##

# Load libraries and data  ------------------------------------------------

library("tidyverse")
library("socialmixr")

## contact data rounds 2-5
mix2 <- readRDS("data/input/df_merged_rd2.RDS")
mix3 <- readRDS("data/input/df_merged_rd3.RDS")
mix4 <- readRDS("data/input/df_merged_rd4.RDS")
mix5 <- readRDS("data/input/df_merged_rd5.RDS")

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
                                    ordered = TRUE)) %>%
  mutate(contact_duration = case_when(totaltime == "Less than 5 minutes" ~ "<5 mins",
                                      totaltime == "Between 5 to 15 minutes" ~ "5-15 mins",
                                      totaltime == "Between 15 minutes to 1 hour" ~ "15 mins-1 hr",
                                      totaltime == "Between 1 hour to 4 hours" ~ "1-4 hrs",
                                      totaltime == "More than 4 hours" ~ "4+ hrs"))

# participants with 0 contacts in round 1
rd1.no.contacts <- part1 %>%
  filter(!(id %in% contacts1$id)) %>%
  group_by(age_cat) %>%
  tally(name = "n_participants0")

# create list with all rounds
mix5x <- mix5 %>%
  mutate(participant_age_2 = participant_age)
mix.list.raw <- list(mix1,mix2,mix3,mix4,mix5x)

## function to collapse data and estimate mean degree------
summariseContacts <- function(df,filter.length = c("")) {

  df <- df %>%
    select(id,participant_age_2,contact_age_3,contact_location2,contact_duration) %>%
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
    ungroup() %>%
    filter(!(contact_duration %in% filter.length)) %>%
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

############################## ALL CONTACTS ############################################################
# Format data  ------------------------------------------------------------

## apply function to list of contact data
mix.list <- lapply(mix.list.raw, summariseContacts)

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

mix5.sum <- mix.list[[5]] %>%
  rename(n_participants = n_participants,
         n_contacts = n_contacts,
         degree = degree,
         total_degree = total_degree,
         prop_degree = prop_degree) %>%
  mutate(round = 5) %>%
  ungroup()

## join all rounds to compare
mix.sum <- mix1.sum %>%
  add_row(mix2.sum) %>%
  add_row(mix3.sum) %>%
  add_row(mix4.sum) %>%
  add_row(mix5.sum)

## apply ratio of 40-49-50:59 contacts from round 2 to round 1 -------------
r40.50 <- mix2.sum %>%
  filter(contact_age %in% c("40-49","50-59")) %>%
  group_by(location,p_age) %>%
  mutate(ratio = n_contacts/sum(n_contacts)) %>%
  select(location,p_age,contact_age,ratio) %>%
  mutate(contact_age = paste0("age",substr(contact_age,1,2),"s")) %>%
  pivot_wider(names_from = contact_age,values_from = ratio)

mix1.40.59 <- mix1.sum %>%
  filter(contact_age == "40-59") %>%
  left_join(r40.50,by = c("location","p_age")) %>%
  mutate(contacts40s = age40s*n_contacts,
         contacts50s = age50s*n_contacts) %>%
  select(location,p_age,n_participants,contacts40s,contacts50s) %>%
  pivot_longer(cols = c("contacts40s","contacts50s"),
               names_to = "contact_age",values_to = "n_contacts") %>%
  mutate(contact_age = case_when(contact_age == "contacts40s" ~ "40-49",
                                 contact_age == "contacts50s" ~ "50-59")) %>%
  select(location,p_age,contact_age,n_participants,n_contacts)

mix1.sum.edit <- mix1.sum %>%
  select(names(mix1.40.59)) %>%
  add_row(mix1.40.59) %>%
  filter(contact_age != "40-59") %>%
  arrange(location,p_age,contact_age)

## add participants with 0 contacts in round 1 -----------------------------
mix1.sum.edit <- mix1.sum.edit %>%
  left_join(rd1.no.contacts,by = c("p_age" = "age_cat")) %>%
  mutate(n_participants = n_participants + n_participants0) %>%
  select(-n_participants0) %>%
  group_by(location,p_age) %>%
  mutate(degree = (n_contacts/n_participants)/2) %>%   # average mean degree for 1 day
  mutate(total_degree = sum(degree)) %>%  # total degree by location and participant age
  mutate(prop_degree = degree/total_degree) %>%   # age distribution of contacts by location and participant age
  mutate(round = 1) %>%
  ungroup()


## join all rounds with edited round 1 data --------------------------------
mix.sum <- mix1.sum.edit %>%
  add_row(mix2.sum) %>%
  add_row(mix3.sum) %>%
  add_row(mix4.sum) %>%
  add_row(mix5.sum)

############################## EXCLUDING CONTACTS < 15 MIN #############################################
# Format data: Excluding contacts <15 mins  ------------------------------------------------------------



## apply function to list of contact data
mix.list.no15 <- lapply(mix.list.raw, summariseContacts, filter.length = c("<5 mins","5-15 mins"))

## split up data sets for each round
mix1.sum.no15 <- mix.list.no15[[1]] %>%
  rename(n_participants = n_participants,
         n_contacts = n_contacts,
         degree = degree,
         total_degree = total_degree,
         prop_degree = prop_degree) %>%
  mutate(round = 1) %>%
  ungroup()

mix2.sum.no15 <- mix.list.no15[[2]] %>%
  rename(n_participants = n_participants,
         n_contacts = n_contacts,
         degree = degree,
         total_degree = total_degree,
         prop_degree = prop_degree) %>%
  mutate(round = 2) %>%
  ungroup()

mix3.sum.no15 <- mix.list.no15[[3]] %>%
  rename(n_participants = n_participants,
         n_contacts = n_contacts,
         degree = degree,
         total_degree = total_degree,
         prop_degree = prop_degree) %>%
  mutate(round = 3) %>%
  ungroup()

mix4.sum.no15 <- mix.list.no15[[4]] %>%
  rename(n_participants = n_participants,
         n_contacts = n_contacts,
         degree = degree,
         total_degree = total_degree,
         prop_degree = prop_degree) %>%
  mutate(round = 4) %>%
  ungroup()

mix5.sum.no15 <- mix.list.no15[[5]] %>%
  rename(n_participants = n_participants,
         n_contacts = n_contacts,
         degree = degree,
         total_degree = total_degree,
         prop_degree = prop_degree) %>%
  mutate(round = 5) %>%
  ungroup()

## join all rounds to compare
mix.sum.no15 <- mix1.sum.no15 %>%
  add_row(mix2.sum.no15) %>%
  add_row(mix3.sum.no15) %>%
  add_row(mix4.sum.no15) %>%
  add_row(mix5.sum.no15)

## apply ratio of 40-49-50:59 contacts from round 2 to round 1 -------------
r40.50.no15 <- mix2.sum.no15 %>%
  filter(contact_age %in% c("40-49","50-59")) %>%
  group_by(location,p_age) %>%
  mutate(ratio = n_contacts/sum(n_contacts)) %>%
  select(location,p_age,contact_age,ratio) %>%
  mutate(contact_age = paste0("age",substr(contact_age,1,2),"s")) %>%
  pivot_wider(names_from = contact_age,values_from = ratio)

mix1.40.59.no15 <- mix1.sum.no15 %>%
  filter(contact_age == "40-59") %>%
  left_join(r40.50.no15,by = c("location","p_age")) %>%
  mutate(contacts40s = age40s*n_contacts,
         contacts50s = age50s*n_contacts) %>%
  select(location,p_age,n_participants,contacts40s,contacts50s) %>%
  pivot_longer(cols = c("contacts40s","contacts50s"),
               names_to = "contact_age",values_to = "n_contacts") %>%
  mutate(contact_age = case_when(contact_age == "contacts40s" ~ "40-49",
                                 contact_age == "contacts50s" ~ "50-59")) %>%
  select(location,p_age,contact_age,n_participants,n_contacts)

mix1.sum.edit.no15 <- mix1.sum.no15 %>%
  select(names(mix1.40.59)) %>%
  add_row(mix1.40.59) %>%
  filter(contact_age != "40-59") %>%
  arrange(location,p_age,contact_age)

## add participants with 0 contacts in round 1 -----------------------------
mix1.sum.edit.no15 <- mix1.sum.edit.no15 %>%
  left_join(rd1.no.contacts,by = c("p_age" = "age_cat")) %>%
  mutate(n_participants = n_participants + n_participants0) %>%
  select(-n_participants0) %>%
  group_by(location,p_age) %>%
  mutate(degree = (n_contacts/n_participants)/2) %>%   # average mean degree for 1 day
  mutate(total_degree = sum(degree)) %>%  # total degree by location and participant age
  mutate(prop_degree = degree/total_degree) %>%   # age distribution of contacts by location and participant age
  mutate(round = 1) %>%
  ungroup()


## join all rounds with edited round 1 data --------------------------------
mix.sum.no15 <- mix1.sum.edit.no15 %>%
  add_row(mix2.sum.no15) %>%
  add_row(mix3.sum.no15) %>%
  add_row(mix4.sum.no15) %>%
  add_row(mix5.sum.no15)

## save round 4 and 5 contact data ------
saveRDS(mix4.sum.no15, file = "data/input/rd4_mean_degree.rds")
saveRDS(mix5.sum.no15, file = "data/input/rd5_mean_degree.rds")

## Round 4: create symmetrical social mixing matrix by taking average ----

avg.mix <- mix4 %>%
  mutate(p_age = as.character(participant_age),
         contact_age = as.character(contact_age_2),
         location = as.character(contact_location2)) %>%
  select(id,p_age,contact_age,location,contact_duration) %>%
  mutate(contact_age = factor(contact_age)) %>%
  mutate(p_age = factor(p_age)) %>%
  mutate(location = factor(location)) %>%
  group_by(p_age) %>%
  mutate(n_participants = n_distinct(id)) %>%  # total number of participants by age group
  ungroup() %>%
  filter(!(contact_duration %in% c("<5 mins","5-15 mins"))) %>%
  group_by(location,p_age,contact_age,n_participants,.drop = FALSE) %>%
  summarise(n_contacts = n()) %>%  # total number of contacts reported over 2 days
  ungroup() %>%
  complete(location,p_age,contact_age) %>%
  group_by(p_age) %>%
  fill(n_participants,.direction = "downup") %>%
  group_by(location,p_age) %>%
  mutate(degree = (n_contacts/n_participants)/2) %>%   # average mean degree for 1 day
  mutate(mix = case_when(p_age == "20-29" & contact_age == "20-29" ~ "20s-20s",
                         (p_age == "20-29" & contact_age == "30-39") |
                           (p_age == "30-39" & contact_age == "20-29") ~ "20s-30s",
                         (p_age == "20-29" & contact_age == "40-49") |
                           (p_age == "40-49" & contact_age == "20-29") ~ "20s-40s",
                         (p_age == "20-29" & contact_age == "50-59") |
                           (p_age == "50-59" & contact_age == "20-29") ~ "20s-50s",
                         (p_age == "20-29" & contact_age == "60-69") |
                           (p_age == "60-69" & contact_age == "20-29") ~ "20s-60s",
                         (p_age == "20-29" & contact_age == "70+") |
                           (p_age == "70+" & contact_age == "20-29") ~ "20s-70s",
                         p_age == "30-39" & contact_age == "30-39" ~ "30s-30s",
                         (p_age == "30-39" & contact_age == "40-49") |
                           (p_age == "40-49" & contact_age == "30-39") ~ "30s-40s",
                         (p_age == "30-39" & contact_age == "50-59") |
                           (p_age == "50-59" & contact_age == "30-39") ~ "30s-50s",
                         (p_age == "30-39" & contact_age == "60-69") |
                           (p_age == "60-69" & contact_age == "30-39") ~ "30s-60s",
                         (p_age == "30-39" & contact_age == "70+") |
                           (p_age == "70+" & contact_age == "30-39") ~ "30s-70s",
                         p_age == "40-49" & contact_age == "40-49" ~ "40s-40s",
                         (p_age == "40-49" & contact_age == "50-59") |
                           (p_age == "50-59" & contact_age == "40-49") ~ "40s-50s",
                         (p_age == "40-49" & contact_age == "60-69") |
                           (p_age == "60-69" & contact_age == "40-49") ~ "40s-60s",
                         (p_age == "40-49" & contact_age == "70+") |
                           (p_age == "70+" & contact_age == "40-49") ~ "40s-70s",
                         p_age == "50-59" & contact_age == "50-59" ~ "50s-50s",
                         (p_age == "50-59" & contact_age == "60-69") |
                           (p_age == "60-69" & contact_age == "50-59") ~ "50s-60s",
                         (p_age == "50-59" & contact_age == "70+") |
                           (p_age == "70+" & contact_age == "50-59") ~ "50s-70s",
                         p_age == "60-69" & contact_age == "60-69" ~ "60s-60s",
                         (p_age == "60-69" & contact_age == "70+") |
                           (p_age == "70+" & contact_age == "60-69") ~ "60s-70s",
                         p_age == "70+" & contact_age == "70+" ~ "70s-70s",
                         p_age == "20-29" & contact_age == "0-9" ~ "0s-20s",
                         p_age == "30-39" & contact_age == "0-9" ~ "0s-30s",
                         p_age == "40-49" & contact_age == "0-9" ~ "0s-40s",
                         p_age == "50-59" & contact_age == "0-9" ~ "0s-50s",
                         p_age == "60-69" & contact_age == "0-9" ~ "0s-60s",
                         p_age == "70+" & contact_age == "0-9" ~ "0s-70s",
                         p_age == "20-29" & contact_age == "10-19" ~ "10s-20s",
                         p_age == "30-39" & contact_age == "10-19" ~ "10s-30s",
                         p_age == "40-49" & contact_age == "10-19" ~ "10s-40s",
                         p_age == "50-59" & contact_age == "10-19" ~ "10s-50s",
                         p_age == "60-69" & contact_age == "10-19" ~ "10s-60s",
                         p_age == "70+" & contact_age == "10-19" ~ "10s-70s")) %>%
  ungroup() %>%
  group_by(location,mix) %>%
  summarise(degree = mean(degree))


avg.matrix <- avg.mix %>%
  separate(mix,c("age1","age2")) %>%
  pivot_wider(names_from = "age1",values_from = "degree")

avg.total <- avg.mix %>%
  ungroup() %>%
  mutate(age20s = ifelse(grepl("20s",mix),degree,NA),
         age30s = ifelse(grepl("30s",mix),degree,NA),
         age40s = ifelse(grepl("40s",mix),degree,NA),
         age50s = ifelse(grepl("50s",mix),degree,NA),
         age60s = ifelse(grepl("60s",mix),degree,NA),
         age70s = ifelse(grepl("70s",mix),degree,NA)) %>%
  select(-mix,-degree) %>%
  group_by(location) %>%
  summarise(age20s = sum(age20s, na.rm = T),
            age30s = sum(age30s, na.rm = T),
            age40s = sum(age40s, na.rm = T),
            age50s = sum(age50s, na.rm = T),
            age60s = sum(age60s, na.rm = T),
            age70s = sum(age70s, na.rm = T)) %>%
  pivot_longer(cols = -location,names_to = "age2", values_to = "total_degree") %>%
  mutate(age2 = substr(age2,4,6))

avg.matrix <- avg.matrix %>%
  left_join(avg.total, by = c("location","age2")) %>%
  rename(age = age2) %>%
  mutate(prop_match = case_when(age == "20s" ~ `20s`/total_degree,
                                age == "30s" ~ `30s`/total_degree,
                                age == "40s" ~ `40s`/total_degree,
                                age == "50s" ~ `50s`/total_degree,
                                age == "60s" ~ `60s`/total_degree,
                                age == "70s" ~ `70s`/total_degree))

saveRDS(avg.matrix, file = "data/input/rd4_matrix.rds") #matrix with mean degree in each cell
write.csv(avg.matrix,"cm_mixing_matrix.csv",row.names = F,na = "")

# complete matrix with proportions instead of mean degree, incorporating comix data for children
avg.mix2 <- avg.mix %>%
  separate(mix,c("age1","age2")) %>%
  rename(age2 = age1,
         age1 = age2) %>%
  filter(age1 != age2)

avg.matrix2 <- avg.mix %>%
  separate(mix,c("age1","age2")) %>%
  ungroup() %>%
  add_row(avg.mix2) %>%
  mutate(location = factor(location),
         age1 = factor(age1),
         age2 = factor(age2)) %>%
  complete(location,age1,age2) %>%
  pivot_wider(names_from = "age2",values_from = "degree")

avg.total2 <- avg.total %>%
  add_row(data.frame(location = c("Home","Home","Community","Community"),
                     age2 = c("0s","10s","0s","10s"),
                     #total_degree = c(2.25,2.25,4,6) # numbers from comix, includes contacts <15min
                     total_degree = c(1.36,1.14,2.08,3.42) # numbers from polymod, excludes contacts <15min
                     ))

avg.matrix.prop <- avg.matrix2 %>%
  left_join(avg.total2, by = c("location","age1" = "age2")) %>%
  rename(age = age1) %>%
  mutate(`20s` = `20s`/total_degree,
         `30s` = `30s`/total_degree,
         `40s` = `40s`/total_degree,
         `50s` = `50s`/total_degree,
         `60s` = `60s`/total_degree,
         `70s` = `70s`/total_degree)


saveRDS(avg.matrix.prop, file = "data/input/rd4_matrix_prop.rds") # matrix with proportion of total mean degree in each cell
write.csv(avg.matrix.prop,"cm_mixing_matrix_prop.csv",row.names = F,na = "")

## Round 5: create symmetrical social mixing matrix by taking average ----

avg.mix5 <- mix5 %>%
  mutate(p_age = as.character(participant_age),
         contact_age = as.character(contact_age_3),
         location = as.character(contact_location2)) %>%
  select(id,p_age,contact_age,location,contact_duration) %>%
  mutate(contact_age = factor(contact_age)) %>%
  mutate(p_age = factor(p_age)) %>%
  mutate(location = factor(location)) %>%
  group_by(p_age) %>%
  mutate(n_participants = n_distinct(id)) %>%  # total number of participants by age group
  ungroup() %>%
  filter(!(contact_duration %in% c("<5 mins","5-15 mins"))) %>%
  group_by(location,p_age,contact_age,n_participants,.drop = FALSE) %>%
  summarise(n_contacts = n()) %>%  # total number of contacts reported over 2 days
  ungroup() %>%
  complete(location,p_age,contact_age) %>%
  group_by(p_age) %>%
  fill(n_participants,.direction = "downup") %>%
  group_by(location,p_age) %>%
  mutate(degree = (n_contacts/n_participants)/2) %>%   # average mean degree for 1 day
  mutate(mix = case_when(p_age == "20-29" & contact_age == "20-29" ~ "20s-20s",
                         (p_age == "20-29" & contact_age == "30-39") |
                           (p_age == "30-39" & contact_age == "20-29") ~ "20s-30s",
                         (p_age == "20-29" & contact_age == "40-49") |
                           (p_age == "40-49" & contact_age == "20-29") ~ "20s-40s",
                         (p_age == "20-29" & contact_age == "50-59") |
                           (p_age == "50-59" & contact_age == "20-29") ~ "20s-50s",
                         (p_age == "20-29" & contact_age == "60+") |
                           (p_age == "60+" & contact_age == "20-29") ~ "20s-60s",
                         p_age == "30-39" & contact_age == "30-39" ~ "30s-30s",
                         (p_age == "30-39" & contact_age == "40-49") |
                           (p_age == "40-49" & contact_age == "30-39") ~ "30s-40s",
                         (p_age == "30-39" & contact_age == "50-59") |
                           (p_age == "50-59" & contact_age == "30-39") ~ "30s-50s",
                         (p_age == "30-39" & contact_age == "60+") |
                           (p_age == "60+" & contact_age == "30-39") ~ "30s-60s",
                         p_age == "40-49" & contact_age == "40-49" ~ "40s-40s",
                         (p_age == "40-49" & contact_age == "50-59") |
                           (p_age == "50-59" & contact_age == "40-49") ~ "40s-50s",
                         (p_age == "40-49" & contact_age == "60+") |
                           (p_age == "60+" & contact_age == "40-49") ~ "40s-60s",
                         p_age == "50-59" & contact_age == "50-59" ~ "50s-50s",
                         (p_age == "50-59" & contact_age == "60+") |
                           (p_age == "60+" & contact_age == "50-59") ~ "50s-60s",
                         p_age == "60+" & contact_age == "60+" ~ "60s-60s",
                         p_age == "20-29" & contact_age == "0-9" ~ "0s-20s",
                         p_age == "30-39" & contact_age == "0-9" ~ "0s-30s",
                         p_age == "40-49" & contact_age == "0-9" ~ "0s-40s",
                         p_age == "50-59" & contact_age == "0-9" ~ "0s-50s",
                         p_age == "60+" & contact_age == "0-9" ~ "0s-60s",
                         p_age == "20-29" & contact_age == "10-19" ~ "10s-20s",
                         p_age == "30-39" & contact_age == "10-19" ~ "10s-30s",
                         p_age == "40-49" & contact_age == "10-19" ~ "10s-40s",
                         p_age == "50-59" & contact_age == "10-19" ~ "10s-50s",
                         p_age == "60+" & contact_age == "10-19" ~ "10s-60s")) %>%
  ungroup() %>%
  group_by(location,mix) %>%
  summarise(degree = mean(degree))


avg.matrix5 <- avg.mix5 %>%
  separate(mix,c("age1","age2")) %>%
  pivot_wider(names_from = "age1",values_from = "degree")

avg.total5 <- avg.mix5 %>%
  ungroup() %>%
  mutate(age20s = ifelse(grepl("20s",mix),degree,NA),
         age30s = ifelse(grepl("30s",mix),degree,NA),
         age40s = ifelse(grepl("40s",mix),degree,NA),
         age50s = ifelse(grepl("50s",mix),degree,NA),
         age60s = ifelse(grepl("60s",mix),degree,NA),
         age70s = ifelse(grepl("70s",mix),degree,NA)) %>%
  select(-mix,-degree) %>%
  group_by(location) %>%
  summarise(age20s = sum(age20s, na.rm = T),
            age30s = sum(age30s, na.rm = T),
            age40s = sum(age40s, na.rm = T),
            age50s = sum(age50s, na.rm = T),
            age60s = sum(age60s, na.rm = T),
            age70s = sum(age70s, na.rm = T)) %>%
  pivot_longer(cols = -location,names_to = "age2", values_to = "total_degree") %>%
  mutate(age2 = substr(age2,4,6))

avg.matrix5 <- avg.matrix5 %>%
  left_join(avg.total5, by = c("location","age2")) %>%
  rename(age = age2) %>%
  mutate(prop_match = case_when(age == "20s" ~ `20s`/total_degree,
                                age == "30s" ~ `30s`/total_degree,
                                age == "40s" ~ `40s`/total_degree,
                                age == "50s" ~ `50s`/total_degree,
                                age == "60s" ~ `60s`/total_degree))

saveRDS(avg.matrix5, file = "data/input/rd5_matrix.rds") #matrix with mean degree in each cell
write.csv(avg.matrix5,"rd5_cm_mixing_matrix.csv",row.names = F,na = "")

# complete matrix with proportions instead of mean degree, incorporating comix data for children
avg.mix2 <- avg.mix5 %>%
  separate(mix,c("age1","age2")) %>%
  rename(age2 = age1,
         age1 = age2) %>%
  filter(age1 != age2)

avg.matrix2 <- avg.mix5 %>%
  separate(mix,c("age1","age2")) %>%
  ungroup() %>%
  add_row(avg.mix2) %>%
  mutate(location = factor(location),
         age1 = factor(age1),
         age2 = factor(age2)) %>%
  complete(location,age1,age2) %>%
  pivot_wider(names_from = "age2",values_from = "degree")

avg.total2 <- avg.total5 %>%
  add_row(data.frame(location = c("Home","Home","Community","Community"),
                     age2 = c("0s","10s","0s","10s"),
                     #total_degree = c(2.25,2.25,4,6) # numbers from comix, includes contacts <15min
                     total_degree = c(1.36,1.14,2.08,3.42) # numbers from polymod, excludes contacts <15min
  ))

avg.matrix.prop5 <- avg.matrix2 %>%
  left_join(avg.total2, by = c("location","age1" = "age2")) %>%
  rename(age = age1) %>%
  mutate(`20s` = `20s`/total_degree,
         `30s` = `30s`/total_degree,
         `40s` = `40s`/total_degree,
         `50s` = `50s`/total_degree,
         `60s` = `60s`/total_degree)


saveRDS(avg.matrix.prop5, file = "data/input/rd5_matrix_prop.rds") # matrix with proportion of total mean degree in each cell
write.csv(avg.matrix.prop5,"rd5_cm_mixing_matrix_prop.csv",row.names = F,na = "")

############################## DESCRIPTIVE FIGURES #####################################################
############################## ALL CONTACTS ############################################################
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

# summary of other contact level variables --------------------------------

library(gtsummary)
library(flextable)

mix1 %>%
  select(cont_attr,totaltime,contact_location2) %>%
  tbl_summary(by = contact_location2)

mix2 %>%
  select(contact_type,contact_relation,
         contact_freq,contact_duration,repeat_contact,
         mask_wearing,contact_location2) %>%
  tbl_summary(by = contact_location2) -> tb2

mix3 %>%
  select(contact_type,contact_relation,
         contact_freq,contact_duration,repeat_contact,
         mask_wearing,contact_location2) %>%
  tbl_summary(by = contact_location2) -> tb3

mix4 %>%
  select(contact_type,contact_relation,
         contact_freq,contact_duration,repeat_contact,
         mask_wearing,contact_location2) %>%
  tbl_summary(by = contact_location2) -> tb4

tbl_merge(list(tb2,tb3,tb4),
          tab_spanner = c("**Round2**","**Round 3**","**Round 4**")) %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = "cm_summary_tbls.docx")



############################## DESCRIPTIVE FIGURES #####################################################
############################## EXCLUDING CONTACTS < 15 MIN #############################################
## plot degrees to compare across rounds -----------------------------------
library(ggplot2)
library(NatParksPalettes)

# Total degree by participant age group and location
mix.sum.total.no15 <- mix.sum.no15 %>%
  select(location,p_age,total_degree,round) %>%
  distinct()

ggplot(data = mix.sum.total.no15, mapping = aes(x = as.factor(round),
                                                y = total_degree,
                                                group = p_age,
                                                color = p_age))+
  geom_line(size = 1)+
  facet_grid(cols = vars(location))+
  scale_color_manual(values=natparks.pals("Yellowstone", 5))+
  theme_bw()+
  labs(x = "Round", y = "Total degree", color = "Participant age")
ggsave(filename = "out/total_degree_no15.png",height = 3,width = 6)

# Degree by participant age, contact age, and location
ggplot(data = mix.sum.no15 %>%
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
ggsave(filename = "out/degree_by_age_no15.png",height = 7,width = 7)

# Distribution of contact age groups by participant age and location
ggplot(data = mix.sum.no15 %>%
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
ggsave(filename = "out/prop_degree_by_age_no15.png",height = 7,width = 7)

# summary of other contact level variables --------------------------------

library(gtsummary)
library(flextable)

mix1 %>%
  filter(!(contact_duration %in% c("<5 mins","5-15 mins"))) %>%
  select(cont_attr,totaltime,contact_location2) %>%
  tbl_summary(by = contact_location2)

mix2 %>%
  filter(!(contact_duration %in% c("<5 mins","5-15 mins"))) %>%
  select(contact_type,contact_relation,
         contact_freq,contact_duration,repeat_contact,
         mask_wearing,contact_location2) %>%
  tbl_summary(by = contact_location2) -> tb2

mix3 %>%
  filter(!(contact_duration %in% c("<5 mins","5-15 mins"))) %>%
  select(contact_type,contact_relation,
         contact_freq,contact_duration,repeat_contact,
         mask_wearing,contact_location2) %>%
  tbl_summary(by = contact_location2) -> tb3

mix4 %>%
  filter(!(contact_duration %in% c("<5 mins","5-15 mins"))) %>%
  select(contact_type,contact_relation,
         contact_freq,contact_duration,repeat_contact,
         mask_wearing,contact_location2) %>%
  tbl_summary(by = contact_location2) -> tb4

tbl_merge(list(tb2,tb3,tb4),
          tab_spanner = c("**Round2**","**Round 3**","**Round 4**")) %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = "cm_summary_tbls_no15.docx")



# summary of one-off, persistent, and other contacts (round 4) ----
persist <- mix4 %>%
  select(id,participant_age_2,contact_age_3,contact_location2,contact_duration,contact_freq) %>%
  mutate(contact_freq = ifelse(contact_freq %in% c("Once/Twice a week",
                                                   "Once/Twice a month",
                                                   "Less than once per month"),
                               "Sometimes",as.character(contact_freq))) %>%
  mutate(contact_freq = factor(contact_freq,levels = c("Daily/Almost daily","Sometimes","Never met before"))) %>%
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
  ungroup() %>%
  filter(!(contact_duration %in% c("<5 mins","5-15 mins"))) %>%
  group_by(location,p_age,contact_age,contact_freq,n_participants,.drop = FALSE) %>%
  summarise(n_contacts = n()) %>%  # total number of contacts reported over 2 days
  ungroup() %>%
  complete(location,p_age,contact_age,contact_freq) %>%
  group_by(location,p_age) %>%
  fill(n_participants,.direction = "downup") %>%
  mutate(degree = (n_contacts/n_participants)/2) %>%   # average mean degree for 1 day
  group_by(location,p_age,contact_age) %>%
  mutate(total_degree = sum(degree)) %>%  # total degree by location and participant age
  mutate(prop_degree = degree/total_degree)  # age distribution of contacts by location and participant age

# distribution of contact rate across participants (round 4) ----
cont.dist <- mix4 %>%
  select(id,participant_age_2,contact_duration) %>%
  mutate(short_duration = ifelse(contact_duration %in% c("<5 mins","5-15 mins"),"Yes","No")) %>%
  mutate(short_duration = factor(short_duration)) %>%
  mutate(p_age = as.character(participant_age_2)) %>%
  group_by(id,p_age,short_duration,.drop = FALSE) %>%
  summarise(n_contacts = n()/2) %>%
  filter(short_duration == "No")

ggplot(data = cont.dist, aes(x = n_contacts))+
  geom_histogram(stat = "count")+
  facet_grid(cols = vars(p_age))+
  labs(title="Distribution of contacts across participants by age group",
       x = "Number of contacts", y = "Number of participants")+
  theme_bw()
ggsave(filename = "out/dist_contacts_no15.png",height = 3,width = 6)

cont.dist %>%
  ungroup() %>%
  select(p_age,n_contacts) %>%
  tbl_summary(by = p_age,
              type = all_continuous() ~ "continuous2",
              statistic = all_continuous() ~ c("{median} ({p25}, {p75})",
                                               "{mean}",
                                               "{min}, {max}"))
