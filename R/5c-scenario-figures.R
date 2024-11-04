# libraries ----
pacman::p_load(tidyr,
               dplyr,
               ggplot2,
               readxl,
               stringr,
               scales)

# data from tables ----
d.iso <- read_excel("data/corpmix_table_data.xlsx",
                    sheet = "Isolation")
d <- read_excel("data/corpmix_table_data.xlsx",
                sheet = "All outcomes")

# format data ----
d2 <- d %>%
  mutate(
    p.avert.lo = as.numeric(str_extract(p.avert, "(?<=\\().*?(?=,)"))/100, # Lower interval
    p.avert.hi = as.numeric(str_extract(p.avert, "(?<=, ).*?(?=\\))"))/100, # Upper interval
    p.avert = as.numeric(str_extract(p.avert, "^[^ ]+"))/100 # Point estimate
    ) %>%
  mutate(
    scenario2 = paste0("Office: ",p.office*100,"%\n",
                       "Others: ",p.other*100,"%"),
    scenario2 = factor(scenario2, levels = unique(scenario2)),
    scenario3 = case_when(change == "office workers" ~ paste0(p.office*100,"%"),
                          change == "others" ~ paste0(p.other*100,"%"),
                          change == "total population" ~ paste0(p.office*100,"%")),
    scenario3 = factor(scenario3,
                       levels = c("0%","25%","50%","75%","100%")),
    change2 = paste0("Varying adherence among\n",change),
    change2 = factor(change2),
    type = factor(type,
                  levels = c("total","office"),
                  labels = c("Total population",
                             "Office workers")),
    outcome = factor(outcome,
                     levels = c("Infections","Cases","Deaths"))
    )

# figure of % averted ----

ggplot(d2 %>%
         filter(change != "Reference"),
       aes(x = scenario3, y = p.avert, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = p.avert.lo, ymax = p.avert.hi),
                position = position_dodge(width = 0.9), width = 0.25) +
  scale_y_continuous(labels = percent,
                     limits = c(-0.5,0.5),
                     breaks = seq(-0.5,0.5,0.1))+
  facet_grid(outcome ~ change2, scales = "free") +
  labs(x = "Adherence", y = "Percent difference",
       fill = "Population") +
  theme_bw() +
    scale_fill_grey()

ggsave(here::here("results","figures","Figure_outcomes.tiff"),
       height = 7,width = 8)
