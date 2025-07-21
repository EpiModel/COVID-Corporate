# Load libraries ---------------------------------------------------------------
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(patchwork)))
suppressWarnings(suppressMessages(library(viridis)))
suppressWarnings(suppressMessages(library(cowplot)))
################################################################################
# Calibration Figures ----------------------------------------------------------
baseline_sims <- readRDS(here::here("data","output","results",
                                   "baseline_sim.rds"))

num.sim <- 256

times <- data.frame(DayNum = c(608, 577, 546, 516, 485,
                               455, 424, 396, 365, 334,
                               304, 273, 243, 212, 181,
                               151, 120, 90, 59, 31),
                    Month = c("Aug-22", "Jul-22", "Jun-22", "May-22", "Apr-22",
                              "Mar-22", "Feb-22", "Jan-22", "Dec-21", "Nov-21",
                              "Oct-21", "Sep-21", "Aug-21", "Jul-21", "Jun-21",
                              "May-21", "Apr-21", "Mar-21", "Feb-21", "Jan-21"))
times <- times %>% map_df(rev)

colors1 <- c("Result" = "dodgerblue4", "Target (Empirical)" = "deepskyblue")
colors <- c("Dose 1: Result" = "dodgerblue4", "Dose 2: Result" = "darkgreen",
            "Dose 3: Result" = "darkred", "Dose 4: Result" = "darkorchid4",
            "Dose 1: Target" = "deepskyblue", "Dose 2: Target" = "chartreuse",
            "Dose 3: Target" = "deeppink", "Dose 4: Target" = "purple")

# Cases -------------------------------------------------------------------

case_targets <- c(8373, 3165, 1590, 1301, 784, 382, 1439, 7256, 6079,
                  1981, 958, 6259, 19544, 3194, 1278, 701, 1749, 2988,
                  4000, 3504)

for (i in 1:num.sim){
  cases <- baseline_sims[["epi"]][["se.flow"]][, i]
  cases <- data.frame(ts = 1:length(cases), cases = cases)
  cases$cumCases <- cumsum(cases$cases)

  cases <- merge(times, cases, by.x = "DayNum", by.y = "ts")
  cases$cases <- cases$cumCases - lag(cases$cumCases)
  cases$cases[1] <- cases$cumCases[1]
  cases <- cases %>% rowid_to_column("ROW")
  assign(paste0("cases_", i), cases)
}

names <- paste0("cases_", 1:num.sim)
cases_all_runs <- bind_rows(mget(names))
cases_all_runs <- cases_all_runs %>% group_by(ROW) %>%
  summarise(Month = min(Month), minResult = min(cases), meanResult = mean(cases), maxResult = max(cases))

cases_all_runs$target <- case_targets

f1 <- ggplot(data = cases_all_runs, aes(x = ROW, y = minResult)) +
  geom_ribbon(data = cases_all_runs, aes(ymin = minResult, ymax = maxResult), fill="dodgerblue4", alpha=0.5) +
  geom_line(aes(y = meanResult, color = "Result")) +
  geom_line(aes(y = target, color = "Target (Empirical)")) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(1, 20),
                     breaks = seq(1,19,2),
                     labels = c("1" = "Jan 21", "3" = "Mar 21", "5" = "May 21",
                                "7" = "Jul 21", "9" = "Sep 21", "11" = "Nov 21",
                                "13" = "Jan 22", "15" = "Mar 22", "17" = "May 22",
                                "19" = "Jul 22")) +
  scale_y_continuous(expand = c(0, 0), lim = c(0, max(cases_all_runs$maxResult)+5)) +
  xlab("Month") +
  ylab("Incident Infections") +
  scale_color_manual(values = colors1, name = "") +
  ggtitle("Infections") + theme_classic()

# Hospitalizations --------------------------------------------------------

hosp_targets <- c(113, 47, 27, 24, 16, 9, 27, 98, 66, 24, 12, 44, 88, 31, 9, 5,
                  9, 19, 29, 26)

for (i in 1:num.sim){
  hosp <- baseline_sims[["epi"]][["ich.flow"]][, i]
  hosp <- data.frame(ts = 1:length(hosp), hosp = hosp)
  hosp$cumHosp <- cumsum(hosp$hosp)

  hosp <- merge(times, hosp, by.x = "DayNum", by.y = "ts")
  hosp$hosp <- hosp$cumHosp - lag(hosp$cumHosp)
  hosp$hosp[1] <- hosp$cumHosp[1]
  hosp <- hosp %>% rowid_to_column("ROW")
  assign(paste0("hosp_", i), hosp)
}

names <- paste0("hosp_", 1:num.sim)
hosp_all_runs <- bind_rows(mget(names))
hosp_all_runs <- hosp_all_runs %>% group_by(ROW) %>%
  summarise(Month = min(Month), minResult = min(hosp), meanResult = mean(hosp), maxResult = max(hosp))

hosp_all_runs$target <- hosp_targets

f2 <- ggplot(data = hosp_all_runs, aes(x = ROW, y = minResult)) +
  geom_ribbon(data = hosp_all_runs, aes(ymin = minResult, ymax = maxResult), fill="dodgerblue4", alpha=0.5) +
  geom_line(aes(y = meanResult, color = "Result")) +
  geom_line(aes(y = target, color = "Target (Empirical)")) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(1, 20),
                     breaks = seq(1,19,2),
                     labels = c("1" = "Jan 21", "3" = "Mar 21", "5" = "May 21",
                                "7" = "Jul 21", "9" = "Sep 21", "11" = "Nov 21",
                                "13" = "Jan 22", "15" = "Mar 22", "17" = "May 22",
                                "19" = "Jul 22")) +
  scale_y_continuous(expand = c(0, 0), lim = c(0, max(hosp_all_runs$maxResult)+5)) +
  xlab("Month") +
  ylab("Incident Hospitalizations") +
  scale_color_manual(values = colors1, name = "") +
  ggtitle("Hospitalizations") + theme_classic()

# Deaths ------------------------------------------------------------------

death_targets <- c(25, 23, 14, 9, 5, 4, 2, 9, 26, 22, 8, 7, 11, 20, 13, 5, 2,
                   2, 3, 5 )

for (i in 1:num.sim){
  death <- baseline_sims[["epi"]][["d.h.flow"]][, i]
  death <- data.frame(ts = 1:length(death), death = death)
  death$cumDeath <- cumsum(death$death)

  death <- merge(times, death, by.x = "DayNum", by.y = "ts")
  death$death <- death$cumDeath - lag(death$cumDeath)
  death$death[1] <- death$cumDeath[1]
  death <- death %>% rowid_to_column("ROW")
  assign(paste0("death_", i), death)
}

names <- paste0("death_", 1:num.sim)
death_all_runs <- bind_rows(mget(names))
death_all_runs <- death_all_runs %>% group_by(ROW) %>%
  summarise(Month = min(Month), minResult = min(death), meanResult = mean(death), maxResult = max(death))

death_all_runs$target <- death_targets

f3 <- ggplot(data = death_all_runs, aes(x = ROW, y = minResult)) +
  geom_ribbon(data = death_all_runs, aes(ymin = minResult, ymax = maxResult), fill="dodgerblue4", alpha=0.5) +
  geom_line(aes(y = meanResult, color = "Result")) +
  geom_line(aes(y = target, color = "Target (Empirical)")) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(1, 20),
                     breaks = seq(1,19,2),
                     labels = c("1" = "Jan 21", "3" = "Mar 21", "5" = "May 21",
                                "7" = "Jul 21", "9" = "Sep 21", "11" = "Nov 21",
                                "13" = "Jan 22", "15" = "Mar 22", "17" = "May 22",
                                "19" = "Jul 22")) +
  scale_y_continuous(expand = c(0, 0), lim = c(0, max(death_all_runs$maxResult)+5)) +
  xlab("Month") +
  ylab("Incident Deaths") +
  scale_color_manual(values = colors1, name = "") +
  ggtitle("Deaths") + theme_classic()

combined <- f1 + f2 + f3 & theme(legend.position = "bottom")
combined + plot_layout(ncol = 1, guides = "collect")
ggsave(here::here("results","figures","Figure1.tiff"),
       dpi=300, compression = 'lzw', width = 3000, height = 2400, units = "px")


# Vaccine Coverage --------------------------------------------------------
source("R/00-calib-targets.R")

for (i in 1:num.sim) {
  vax_cov <- data.frame(cov.vax1.0to4 = baseline_sims[["epi"]][["cov_vax1_0to4"]][, i],
                        cov.vax1.5to17 = baseline_sims[["epi"]][["cov_vax1_5to17"]][, i],
                        cov.vax1.18to64 = baseline_sims[["epi"]][["cov_vax1_18to64"]][, i],
                        cov.vax1.65p = baseline_sims[["epi"]][["cov_vax1_65p"]][, i],
                        cov.vax2.0to4 = baseline_sims[["epi"]][["cov_vax2_0to4"]][, i],
                        cov.vax2.5to17 = baseline_sims[["epi"]][["cov_vax2_5to17"]][, i],
                        cov.vax2.18to64 = baseline_sims[["epi"]][["cov_vax2_18to64"]][, i],
                        cov.vax2.65p = baseline_sims[["epi"]][["cov_vax2_65p"]][, i],
                        cov.vax3.5to17 = baseline_sims[["epi"]][["cov_vax3_5to17"]][, i],
                        cov.vax3.18to49 = baseline_sims[["epi"]][["cov_vax3_18to49"]][, i],
                        cov.vax3.50to64 = baseline_sims[["epi"]][["cov_vax3_50to64"]][, i],
                        cov.vax3.65p = baseline_sims[["epi"]][["cov_vax3_65p"]][, i],
                        cov.vax4.50to64 = baseline_sims[["epi"]][["cov_vax4_50to64"]][, i],
                        cov.vax4.65p = baseline_sims[["epi"]][["cov_vax4_65p"]][, i])
  vax_cov[is.na(vax_cov)] <- 0
  colnames(vax_cov) <- c("cov.vax1.0to4", "cov.vax1.5to17", "cov.vax1.18to64",
                         "cov.vax1.65p", "cov.vax2.0to4", "cov.vax2.5to17",
                         "cov.vax2.18to64", "cov.vax2.65p", "cov.vax3.5to17",
                         "cov.vax3.18to49", "cov.vax3.50to64", "cov.vax3.65p",
                         "cov.vax4.50to64", "cov.vax4.65p")
  vax_cov <- as.data.frame(sapply(vax_cov, function(x) round(x, 3)*100))
  vax_cov$ts <- 1:nrow(vax_cov)
  vax_cov <- merge(times, vax_cov, by.x = "DayNum", by.y = "ts")
  vax_cov <- vax_cov %>% rowid_to_column("ROW")
  assign(paste0("vax_cov_summary_", i), vax_cov)
}

names <- paste0("vax_cov_summary_", 1:num.sim)
vax_all_runs <- bind_rows(mget(names))

vax_all_runs <- vax_all_runs %>% group_by(ROW) %>%
  summarise(mean.vax1.0to4 = mean(cov.vax1.0to4),
            mean.vax1.5to17 = mean(cov.vax1.5to17),
            mean.vax1.18to64 = mean(cov.vax1.18to64),
            mean.vax1.65p = mean(cov.vax1.65p),
            mean.vax2.0to4 = mean(cov.vax2.0to4),
            mean.vax2.5to17 = mean(cov.vax2.5to17),
            mean.vax2.18to64 = mean(cov.vax2.18to64),
            mean.vax2.65p = mean(cov.vax2.65p),
            mean.vax3.5to17 = mean(cov.vax3.5to17),
            mean.vax3.18to49 = mean(cov.vax3.18to49),
            mean.vax3.50to64 = mean(cov.vax3.50to64),
            mean.vax3.65p = mean(cov.vax3.65p),
            mean.vax4.50to64 = mean(cov.vax4.50to64),
            mean.vax4.65p = mean(cov.vax4.65p),
            min.vax1.0to4 = min(cov.vax1.0to4),
            min.vax1.5to17 = min(cov.vax1.5to17),
            min.vax1.18to64 = min(cov.vax1.18to64),
            min.vax1.65p = min(cov.vax1.65p),
            min.vax2.0to4 = min(cov.vax2.0to4),
            min.vax2.5to17 = min(cov.vax2.5to17),
            min.vax2.18to64 = min(cov.vax2.18to64),
            min.vax2.65p = min(cov.vax2.65p),
            min.vax3.5to17 = min(cov.vax3.5to17),
            min.vax3.18to49 = min(cov.vax3.18to49),
            min.vax3.50to64 = min(cov.vax3.50to64),
            min.vax3.65p = min(cov.vax3.65p),
            min.vax4.50to64 = min(cov.vax4.50to64),
            min.vax4.65p = min(cov.vax4.65p),
            max.vax1.0to4 = max(cov.vax1.0to4),
            max.vax1.5to17 = max(cov.vax1.5to17),
            max.vax1.18to64 = max(cov.vax1.18to64),
            max.vax1.65p = max(cov.vax1.65p),
            max.vax2.0to4 = max(cov.vax2.0to4),
            max.vax2.5to17 = max(cov.vax2.5to17),
            max.vax2.18to64 = max(cov.vax2.18to64),
            max.vax2.65p = max(cov.vax2.65p),
            max.vax3.5to17 = max(cov.vax3.5to17),
            max.vax3.18to49 = max(cov.vax3.18to49),
            max.vax3.50to64 = max(cov.vax3.50to64),
            max.vax3.65p = max(cov.vax3.65p),
            max.vax4.50to64 = max(cov.vax4.50to64),
            max.vax4.65p = max(cov.vax4.65p))

vax_targets <- rbind(vax_targets, c("February 2021", rep(NA, 14)))
vax_targets <- rbind(vax_targets, c("January 2021", rep(NA, 14)))
names(vax_targets) <- c("Month", "target.vax1.0to4", "target.vax1.5to17",
                        "target.vax1.18to64", "target.vax1.65p", "target.vax2.0to4",
                        "target.vax2.5to17", "target.vax2.18to64", "target.vax2.65p",
                        "target.vax3.5to17", "target.vax3.18to49", "target.vax3.50to64",
                        "target.vax3.65p", "target.vax4.50to64", "target.vax4.65p")
vax_targets <- vax_targets %>% mutate_at(c("target.vax1.0to4", "target.vax1.5to17",
                                           "target.vax1.18to64", "target.vax1.65p",
                                           "target.vax2.0to4", "target.vax2.5to17",
                                           "target.vax2.18to64", "target.vax2.65p",
                                           "target.vax3.5to17", "target.vax3.18to49",
                                           "target.vax3.50to64", "target.vax3.65p",
                                           "target.vax4.50to64", "target.vax4.65p"), as.numeric)
vax_targets <- vax_targets %>% map_df(rev)
vax_all_runs <- cbind(vax_all_runs, vax_targets)

# Age Group 1 -------------------------------------------------------------

vax_all_runs_1 <- vax_all_runs[, c("ROW", "mean.vax1.0to4", "mean.vax2.0to4",
                                   "min.vax1.0to4",  "min.vax2.0to4",
                                   "max.vax1.0to4", "max.vax2.0to4",
                                   "target.vax1.0to4", "target.vax2.0to4")]

v1 <- ggplot(data = vax_all_runs_1, aes(x = ROW, y = mean.vax1.0to4)) +
  geom_ribbon(data=vax_all_runs_1, aes(ymin = min.vax1.0to4, ymax = max.vax1.0to4), fill="dodgerblue4", alpha=0.5) +
  geom_ribbon(data=vax_all_runs_1, aes(ymin = min.vax2.0to4, ymax = max.vax2.0to4), fill="darkgreen", alpha=0.5) +
  geom_line(aes(y = mean.vax1.0to4, color = "Dose 1: Result")) +
  geom_line(aes(y = mean.vax2.0to4, color = "Dose 2: Result")) +
  geom_point(aes(y = target.vax1.0to4, color = "Dose 1: Target"),size = 1) +
  geom_point(aes(y = target.vax2.0to4, color = "Dose 2: Target"),size = 1) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(1, 20),
                     breaks = c(1, 7, 13, 19),
                     labels = c("1" = "Jan 21", "7" = "Jul 21",
                                "13" = "Jan 22", "19" = "Jul 22")) +
  scale_y_continuous(expand = c(0, 0.05)) +
  xlab("Month") +
  ylab("Vaccine Coverage (%)") +
  ggtitle("Ages 0 to 4") +
  scale_color_manual(values = colors, guide = "none") +
  theme_classic()
v1

# Age Group 2 -------------------------------------------------------------

vax_all_runs_2 <- vax_all_runs[, c("ROW", "mean.vax1.5to17", "mean.vax2.5to17",
                                   "mean.vax3.5to17", "min.vax1.5to17", "min.vax2.5to17",
                                   "min.vax3.5to17", "max.vax1.5to17", "max.vax2.5to17",
                                   "max.vax3.5to17", "target.vax1.5to17",
                                   "target.vax2.5to17", "target.vax3.5to17")]

v2 <- ggplot(data = vax_all_runs_2, aes(x = ROW, y = mean.vax1.5to17)) +
  geom_ribbon(data=vax_all_runs_2, aes(ymin = min.vax1.5to17, ymax = max.vax1.5to17), fill="dodgerblue4", alpha=0.5) +
  geom_ribbon(data=vax_all_runs_2, aes(ymin = min.vax2.5to17, ymax = max.vax2.5to17), fill="darkgreen", alpha=0.5) +
  geom_ribbon(data=vax_all_runs_2, aes(ymin = min.vax3.5to17, ymax = max.vax3.5to17), fill="darkred", alpha=0.5) +
  geom_line(aes(y = mean.vax1.5to17, color = "Dose 1: Result")) +
  geom_line(aes(y = mean.vax2.5to17, color = "Dose 2: Result")) +
  geom_line(aes(y = mean.vax3.5to17, color = "Dose 3: Result")) +
  geom_point(aes(y = target.vax1.5to17, color = "Dose 1: Target"),size = 1) +
  geom_point(aes(y = target.vax2.5to17, color = "Dose 2: Target"),size = 1) +
  geom_point(aes(y = target.vax3.5to17, color = "Dose 3: Target"),size = 1) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(1, 20),
                     breaks = c(1, 7, 13, 19),
                     labels = c("1" = "Jan 21", "7" = "Jul 21",
                                "13" = "Jan 22", "19" = "Jul 22")) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Month") +
  ylab("Vaccine Coverage (%)") +
  ggtitle("Ages 5 to 17") +
  scale_color_manual(values = colors, guide = "none") +
  theme_classic()



# Age Group 3 -------------------------------------------------------------

vax_all_runs_3 <- vax_all_runs[, c("ROW", "mean.vax1.18to64", "mean.vax2.18to64",
                                   "mean.vax3.18to49", "min.vax1.18to64",
                                   "min.vax2.18to64", "min.vax3.18to49",
                                   "max.vax1.18to64", "max.vax2.18to64",
                                   "max.vax3.18to49", "target.vax1.18to64",
                                   "target.vax2.18to64", "target.vax3.18to49")]

v3 <- ggplot(data = vax_all_runs_3, aes(x = ROW, y = mean.vax1.18to64)) +
  geom_ribbon(data=vax_all_runs_3, aes(ymin = min.vax1.18to64, ymax = max.vax1.18to64), fill="dodgerblue4", alpha=0.5) +
  geom_ribbon(data=vax_all_runs_3, aes(ymin = min.vax2.18to64, ymax = max.vax2.18to64), fill="darkgreen", alpha=0.5) +
  geom_ribbon(data=vax_all_runs_3, aes(ymin = min.vax3.18to49, ymax = max.vax3.18to49), fill="darkred", alpha=0.5) +
  geom_line(aes(y = mean.vax1.18to64, color = "Dose 1: Result")) +
  geom_line(aes(y = mean.vax2.18to64, color = "Dose 2: Result")) +
  geom_line(aes(y = mean.vax3.18to49, color = "Dose 3: Result")) +
  geom_point(aes(y = target.vax1.18to64, color = "Dose 1: Target"),size = 1) +
  geom_point(aes(y = target.vax2.18to64, color = "Dose 2: Target"),size = 1) +
  geom_point(aes(y = target.vax3.18to49, color = "Dose 3: Target"),size = 1) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(1, 20),
                     breaks = c(1, 7, 13, 19),
                     labels = c("1" = "Jan 21", "7" = "Jul 21",
                                "13" = "Jan 22", "19" = "Jul 22")) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Month") +
  ylab("Vaccine Coverage (%)") +
  ggtitle("Ages 18 to 49") +
  scale_color_manual(values = colors, guide = "none") +
  theme_classic()


# Age Group 4 -------------------------------------------------------------

vax_all_runs_4 <- vax_all_runs[, c("ROW", "mean.vax1.18to64", "mean.vax2.18to64",
                                   "mean.vax3.50to64", "mean.vax4.50to64",
                                   "min.vax1.18to64", "min.vax2.18to64",
                                   "min.vax3.50to64", "min.vax4.50to64",
                                   "max.vax1.18to64", "max.vax2.18to64",
                                   "max.vax3.50to64", "max.vax4.50to64",
                                   "target.vax1.18to64", "target.vax2.18to64",
                                   "target.vax3.50to64", "target.vax4.50to64")]

v4 <- ggplot(data = vax_all_runs_4, aes(x = ROW, y = mean.vax1.18to64)) +
  geom_ribbon(data=vax_all_runs_4, aes(ymin = min.vax1.18to64, ymax = max.vax1.18to64), fill="dodgerblue4", alpha=0.5) +
  geom_ribbon(data=vax_all_runs_4, aes(ymin = min.vax2.18to64, ymax = max.vax2.18to64), fill="darkgreen", alpha=0.5) +
  geom_ribbon(data=vax_all_runs_4, aes(ymin = min.vax3.50to64, ymax = max.vax3.50to64), fill="darkred", alpha=0.5) +
  geom_ribbon(data=vax_all_runs_4, aes(ymin = min.vax4.50to64, ymax = max.vax4.50to64), fill="darkorchid4", alpha=0.5) +
  geom_line(aes(y = mean.vax1.18to64, color = "Dose 1: Result")) +
  geom_line(aes(y = mean.vax2.18to64, color = "Dose 2: Result")) +
  geom_line(aes(y = mean.vax3.50to64, color = "Dose 3: Result")) +
  geom_line(aes(y = mean.vax4.50to64, color = "Dose 4: Result")) +
  geom_point(aes(y = target.vax1.18to64, color = "Dose 1: Target"),size = 1) +
  geom_point(aes(y = target.vax2.18to64, color = "Dose 2: Target"),size = 1) +
  geom_point(aes(y = target.vax3.50to64, color = "Dose 3: Target"),size = 1) +
  geom_point(aes(y = target.vax4.50to64, color = "Dose 4: Target"),size = 1) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(1, 20),
                     breaks = c(1, 7, 13, 19),
                     labels = c("1" = "Jan 21", "7" = "Jul 21",
                                "13" = "Jan 22", "19" = "Jul 22")) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Month") +
  ylab("Vaccine Coverage (%)") +
  ggtitle("Ages 50 to 64") +
  scale_color_manual(values = colors, guide = "none") +
  theme_classic()


# Age Group 5 -------------------------------------------------------------

vax_all_runs_5 <- vax_all_runs[, c("ROW", "mean.vax1.65p", "mean.vax2.65p",
                                   "mean.vax3.65p", "mean.vax4.65p", "min.vax1.65p",
                                   "min.vax2.65p", "min.vax3.65p", "min.vax4.65p",
                                   "max.vax1.65p", "max.vax2.65p", "max.vax3.65p",
                                   "max.vax4.65p", "target.vax1.65p", "target.vax2.65p",
                                   "target.vax3.65p", "target.vax4.65p")]

v5 <- ggplot(data = vax_all_runs_5, aes(x = ROW, y = mean.vax1.65p)) +
  geom_ribbon(data=vax_all_runs_5, aes(ymin = min.vax1.65p, ymax = max.vax1.65p), fill="navy", alpha=0.5) +
  geom_ribbon(data=vax_all_runs_5, aes(ymin = min.vax2.65p, ymax = max.vax2.65p), fill="darkgreen", alpha=0.5) +
  geom_ribbon(data=vax_all_runs_5, aes(ymin = min.vax3.65p, ymax = max.vax3.65p), fill="darkred", alpha=0.5) +
  geom_ribbon(data=vax_all_runs_5, aes(ymin = min.vax4.65p, ymax = max.vax4.65p), fill="darkorchid4", alpha=0.5) +
  geom_line(aes(y = mean.vax1.65p, color = "Dose 1: Result")) +
  geom_line(aes(y = mean.vax2.65p, color = "Dose 2: Result")) +
  geom_line(aes(y = mean.vax3.65p, color = "Dose 3: Result")) +
  geom_line(aes(y = mean.vax4.65p, color = "Dose 4: Result")) +
  geom_point(aes(y = target.vax1.65p, color = "Dose 1: Target"),size = 1) +
  geom_point(aes(y = target.vax2.65p, color = "Dose 2: Target"),size = 1) +
  geom_point(aes(y = target.vax3.65p, color = "Dose 3: Target"),size = 1) +
  geom_point(aes(y = target.vax4.65p, color = "Dose 4: Target"),size = 1) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(1, 20),
                     breaks = c(1, 7, 13, 19),
                     labels = c("1" = "Jan 21", "7" = "Jul 21",
                                "13" = "Jan 22", "19" = "Jul 22")) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Month") +
  ylab("Vaccine Coverage (%)") +
  ggtitle("Ages 65+") +
  scale_color_manual(values = colors, name = "") + theme_classic() +
  theme(legend.position="bottom")

combined <- v1 + v2 + v3 + v4 + v5
combined + plot_layout(ncol = 3)
ggsave(here::here("results","figures","Figure2.tiff"),
       dpi=300, compression = 'lzw', width = 3000, height = 2400, units = "px")

# Clean up objects in environment ----------------------------------------------
rm(list=setdiff(ls(), c("baseline_sims")))
##################### VARYING ISOLATION ACROSS ENTIRE POPULATION ###############
################################################################################

# Setup data (RUN THIS FIRST) -------------
d.ratio.tot <- readRDS(here::here("results","ratio_df_8000.rds"))

d.ratio.tot <- d.ratio.tot %>%
  filter(scenario2 %in% c(0,2,4,5)) %>%
  mutate(scenario2 = factor(scenario2,
                            levels = c(2,0,4,5),
                            labels = c("Worse adherence (25%)",
                                       "Reference (58%)",
                                       "Improved adherence (75%)",
                                       "Perfect adherence (100%)")))

# Start date
start_date <- as.Date("2021-01-01")

# Convert numeric days to date format by adding to start date
date_variable <- start_date + d.ratio.tot$time - 1

d.ratio.tot <- cbind(d.ratio.tot, date_variable)

# Rename the new column
colnames(d.ratio.tot)[ncol(d.ratio.tot)] <- "date_variable"

# custom palette
custom_palette <- gradient_n_pal(viridis::viridis_pal(option = "magma")(100)[0:90])

# Ratio Plot: Cases to Isolated -------------

p.sum2.tot <- ggplot(d.ratio.tot, aes(x = date_variable, y = med.sum2))+
  geom_line(aes(y = med.sum2, color = scenario2, group = scenario2),
            linewidth = 0.8) +
  geom_ribbon(aes(ymin = p25.sum2, ymax = p75.sum2, group = scenario2,
                  fill = scenario2), alpha = 0.1)+
  scale_x_date(date_breaks = "3 months")+
  scale_color_manual(values = custom_palette(seq(0, 1, length.out = n_distinct(d.ratio.tot$scenario2))))+
  scale_fill_manual(values = custom_palette(seq(0, 1, length.out = n_distinct(d.ratio.tot$scenario2))))+
  ylim(0,1500)+
  labs(x = "Date",
       y = "Sum of cases and persons in isolation",
       color = "Scenario") +
  guides(fill = "none", color = "none")+
  ggtitle("Total number of office workers ill or in isolation") + theme_classic()

p.cases.tot <- ggplot(d.ratio.tot, aes(x = date_variable, y = med.case))+
  geom_line(aes(y = med.case, color = scenario2, group = scenario2),
            linewidth = 0.8) +
  geom_ribbon(aes(ymin = p25.case, ymax = p75.case, group = scenario2,
                  fill = scenario2), alpha = 0.1)+
  scale_x_date(date_breaks = "3 months")+
  scale_color_manual(values = custom_palette(seq(0, 1, length.out = n_distinct(d.ratio.tot$scenario2))))+
  scale_fill_manual(values = custom_palette(seq(0, 1, length.out = n_distinct(d.ratio.tot$scenario2))))+
  ylim(0,1500)+
  labs(x = "Date",
       y = "Number of cases",
       color = "Scenario") +
  guides(fill = "none", color = "none")+
  ggtitle("Total number of office workers ill") + theme_classic()

p.iso.tot <- ggplot(d.ratio.tot, aes(x = date_variable, y = med.iso))+
  geom_line(aes(y = med.iso, color = scenario2, group = scenario2),
            linewidth = 0.8) +
  geom_ribbon(aes(ymin = p25.iso, ymax = p75.iso, group = scenario2,
                  fill = scenario2), alpha = 0.1)+
  scale_x_date(date_breaks = "3 months")+
  scale_color_manual(values = custom_palette(seq(0, 1, length.out = n_distinct(d.ratio.tot$scenario2))))+
  scale_fill_manual(values = custom_palette(seq(0, 1, length.out = n_distinct(d.ratio.tot$scenario2))))+
  ylim(0,1500)+
  labs(x = "Date",
       y = "Number of persons in isolation",
       color = "Scenario") +
  guides(fill = "none", color = "none")+
  ggtitle("Total number of office workers in isolation") + theme_classic()

p.mask.tot <- ggplot(d.ratio.tot, aes(x = date_variable, y = med.mask))+
  geom_line(aes(y = med.mask, color = scenario2, group = scenario2),
            linewidth = 0.8) +
  geom_ribbon(aes(ymin = p25.mask, ymax = p75.mask, group = scenario2,
                  fill = scenario2), alpha = 0.1)+
  scale_x_date(date_breaks = "3 months")+
  scale_color_manual(values = custom_palette(seq(0, 1, length.out = n_distinct(d.ratio.tot$scenario2))))+
  scale_fill_manual(values = custom_palette(seq(0, 1, length.out = n_distinct(d.ratio.tot$scenario2))))+
  ylim(0,1500)+
  labs(x = "Date",
       y = "Number of persons masked",
       color = "Scenario") +
  guides(fill = "none")+
  ggtitle("Total number of office workers masked") + theme_classic()+
  theme(legend.position = c(0.83,0.8))

patch.tot <- (p.sum2.tot | p.mask.tot) / (p.cases.tot | p.iso.tot)
patch.tot <- patch.tot + plot_annotation(
  title = 'Impact of varying adherence to CDC guidelines among entire population',
  caption = 'Figures represent office workers only'
)
patch.tot

ggsave(here::here("results","figures","Figure4A.tiff"),
       width = 15, height = 6, units = "in", dpi=300, compression = 'lzw')

# Stacked area plot: Cases vs isolating ----
d.long <- d.ratio.tot %>%
  select(scenario2,date_variable,med.casenotiso,med.both,med.isonotcase) %>%
  pivot_longer(cols = -c("scenario2","date_variable"), names_to = "cat", values_to = "number") %>%
  mutate(cat = case_when(cat == "med.casenotiso" ~ "Cases not isolating",
                         cat == "med.both" ~ "Cases isolating",
                         cat == "med.isonotcase" ~ "Non-cases isolating")) %>%
  mutate(cat = factor(cat,levels = c("Non-cases isolating",
                                     "Cases isolating",
                                     "Cases not isolating")))

ggplot(d.long, aes(x = date_variable, y = number, fill = cat)) +
  geom_area(position = "stack")+
  facet_wrap(vars(scenario2),nrow = 2)+
  theme_bw()+
  scale_fill_viridis(discrete = TRUE)+
  labs(x = "Date",
       y = "Number of office workers ill and/or isolating",
       fill = "Category",
       title = "Impact of varying adherence among entire population")

ggsave(here::here("results","figures","Figure5A.tiff"),
       width = 8, height = 4, units = "in", dpi=300, compression = 'lzw')

# ggplot(d.long, aes(x = date_variable, y = number, fill = cat)) +
#   geom_area(position = "fill")+
#   facet_wrap(vars(scenario2),nrow = 2)+
#   theme_bw()+
#   scale_fill_viridis(discrete = TRUE)+
#   labs(x = "Date",
#        y = "Proportion of office workers ill and/or isolating",
#        fill = "Category")
#
# ggsave(here::here("results","figures","Figure5A_prop.tiff"),
#        width = 12, height = 6, units = "in", dpi=300, compression = 'lzw')

# Cumulative cases vs isolations in individual simulations ----
d.sim <- readRDS(here::here("results","sim_df.rds")) %>%
  mutate(scenario = case_when(iso.prob.param == 0.25 ~ "Worse adherence (25%)",
                              iso.prob.param == 0.58 ~ "Reference (58%)",
                              iso.prob.param == 0.75 ~ "Improved adherence (75%)",
                              iso.prob.param == 1 ~ "Perfect adherence (100%)")) %>%
  mutate(scenario = factor(scenario,
                           levels = c("Worse adherence (25%)",
                                      "Reference (58%)",
                                      "Improved adherence (75%)",
                                      "Perfect adherence (100%)")))

p.sim1 <- ggplot(d.sim, aes(x = cum.case, y = cum.iso, color = factor(scenario))) +
  geom_point(alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  scale_color_viridis(discrete = TRUE)+
  guides(color = "none")+
  labs(x = "Cumulative cases", y = "Cumulative isolations",
       title = "Total population") +
  theme_bw()

p.sim2 <- ggplot(d.sim, aes(x = cum.case.w, y = cum.iso.w, color = factor(scenario))) +
  geom_point(alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  scale_color_viridis(discrete = TRUE)+
  labs(x = "Cumulative cases", y = "Cumulative isolations", color = "Scenario",
       title = "Office workers") +
  theme_bw()

combined.sim <- p.sim1 + p.sim2 +
  plot_annotation(tag_levels = list(c("A", "B")),
                  title = "Cumulative cases and persons who isolated at the end of the time period",
                  subtitle = "Varying adherence among entire population")
combined.sim

ggsave(here::here("results","figures","Figure6A.tiff"),
       width = 10, height = 5, units = "in", dpi=300)
################################################################################
# Contour plot libraries -------------------------------------------------------
library("viridis")
library("ggplot2")
library("metR")
library("patchwork")

# Load contour plot data (iso prob vs work notif) ------------------------------

cInc <- readRDS(here::here("data","output","results","cInc.rds"))
cDeaths <- readRDS(here::here("data","output","results","cDeaths.rds"))
cIso <- readRDS(here::here("data","output","results","cIso.rds"))
cCases <- readRDS(here::here("data","output","results","cCases.rds"))
cIncWork <- readRDS(here::here("data","output","results","cIncWork.rds"))
cDeathsWork <- readRDS(here::here("data","output","results","cDeathsWork.rds"))
cIsoWork <- readRDS(here::here("data","output","results","cIsoWork.rds"))
cCasesWork <- readRDS(here::here("data","output","results","cCasesWork.rds"))

# Contour Plots: Iso prob vs work notif ----------------------------------------

# Plot cases
minval.i <- floor(min(cInc$cInc))
maxval.i <- ceiling(max(cInc$cInc))
g_inc_cont <- ggplot(cInc, aes(iso.prob, work.notif.prob)) +
  geom_raster(aes(fill = cInc), interpolate = TRUE) +
  geom_contour(aes(z = cInc), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cInc), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability", y = "Notification probability at work") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.i, maxval.i), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot cases - Office Workers
minval.iw <- floor(min(cIncWork$cInc))
maxval.iw <- ceiling(max(cIncWork$cInc))
g_inc_cont_w <- ggplot(cIncWork, aes(iso.prob, work.notif.prob)) +
  geom_raster(aes(fill = cInc), interpolate = TRUE) +
  geom_contour(aes(z = cInc), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cInc), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability", y = "Notification probability at work") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.iw, maxval.iw), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot deaths
minval.d <- min(cDeaths$cDeaths)
maxval.d <- max(cDeaths$cDeaths)
g_death_cont <- ggplot(cDeaths, aes(iso.prob, work.notif.prob)) +
  geom_raster(aes(fill = cDeaths), interpolate = TRUE) +
  geom_contour(aes(z = cDeaths), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cDeaths), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability", y = "Notification probability at work") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.d, maxval.d), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot deaths - Office Workers
minval.dw <- min(cDeathsWork$cDeaths)
maxval.dw <- max(cDeathsWork$cDeaths)
g_death_cont_w <- ggplot(cDeathsWork, aes(iso.prob, work.notif.prob)) +
  geom_raster(aes(fill = cDeaths), interpolate = TRUE) +
  geom_contour(aes(z = cDeaths), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cDeaths), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability", y = "Notification probability at work") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.dw, maxval.dw), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot isolation
minval.v <- min(cIso$cIso)
maxval.v <- max(cIso$cIso)
g_iso_cont <- ggplot(cIso, aes(iso.prob, work.notif.prob)) +
  geom_raster(aes(fill = cIso), interpolate = TRUE) +
  geom_contour(aes(z = cIso), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cIso), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability", y = "Notification probability at work") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.v, maxval.v), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot isolation - Office Workers
minval.vw <- min(cIsoWork$cIso)
maxval.vw <- max(cIsoWork$cIso)
g_iso_cont_w <- ggplot(cIsoWork, aes(iso.prob, work.notif.prob)) +
  geom_raster(aes(fill = cIso), interpolate = TRUE) +
  geom_contour(aes(z = cIso), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cIso), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability", y = "Notification probability at work") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.vw, maxval.vw), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

combined1 <- g_iso_cont + g_iso_cont_w +
  plot_annotation(tag_levels = list(c("A1", "A2")), title = "Total Persons Entering Isolation")
combined2 <- g_inc_cont + g_inc_cont_w +
  plot_annotation(tag_levels = list(c("B1", "B2")), title = "Infections per 100'000 Person-Days")
combined3 <- g_death_cont +
  g_death_cont_w +
  plot_annotation(tag_levels = list(c("C1", "C2")), title = "Deaths per 100'000 Person-Days")

wrap_elements(combined1) / wrap_elements(combined2) / wrap_elements(combined3)

ggsave(here::here("results","figures","Figure3A.tiff"),
       width = 12, height = 12, units = "in", dpi=300, compression = 'lzw')


# Load contour plot data (iso prob vs iso act) ---------------------------------

cInc <- readRDS(here::here("data","output","results","cInc2.rds"))
cDeaths <- readRDS(here::here("data","output","results","cDeaths2.rds"))
cIso <- readRDS(here::here("data","output","results","cIso2.rds"))
cCases <- readRDS(here::here("data","output","results","cCases2.rds"))
cIncWork <- readRDS(here::here("data","output","results","cIncWork2.rds"))
cDeathsWork <- readRDS(here::here("data","output","results","cDeathsWork2.rds"))
cIsoWork <- readRDS(here::here("data","output","results","cIsoWork2.rds"))
cCasesWork <- readRDS(here::here("data","output","results","cCasesWork2.rds"))

# Contour Plots: Iso prob vs iso act -------------------------------------------

# Plot cases
minval.i <- floor(min(cInc$cInc))
maxval.i <- ceiling(max(cInc$cInc))
g_inc_cont <- ggplot(cInc, aes(iso.prob, act.rate.iso)) +
  geom_raster(aes(fill = cInc), interpolate = TRUE) +
  geom_contour(aes(z = cInc), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cInc), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability", y = "Relative contact intensity") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.i, maxval.i), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot cases - Office Workers
minval.iw <- floor(min(cIncWork$cInc))
maxval.iw <- ceiling(max(cIncWork$cInc))
g_inc_cont_w <- ggplot(cIncWork, aes(iso.prob, act.rate.iso)) +
  geom_raster(aes(fill = cInc), interpolate = TRUE) +
  geom_contour(aes(z = cInc), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cInc), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability", y = "Relative contact intensity") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.iw, maxval.iw), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot deaths
minval.d <- min(cDeaths$cDeaths)
maxval.d <- max(cDeaths$cDeaths)
g_death_cont <- ggplot(cDeaths, aes(iso.prob, act.rate.iso)) +
  geom_raster(aes(fill = cDeaths), interpolate = TRUE) +
  geom_contour(aes(z = cDeaths), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cDeaths), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability", y = "Relative contact intensity") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.d, maxval.d), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot deaths - Office Workers
minval.dw <- min(cDeathsWork$cDeaths)
maxval.dw <- max(cDeathsWork$cDeaths)
g_death_cont_w <- ggplot(cDeathsWork, aes(iso.prob, act.rate.iso)) +
  geom_raster(aes(fill = cDeaths), interpolate = TRUE) +
  geom_contour(aes(z = cDeaths), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cDeaths), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability", y = "Relative contact intensity") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.dw, maxval.dw), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot isolation
minval.v <- min(cIso$cIso)
maxval.v <- max(cIso$cIso)
g_iso_cont <- ggplot(cIso, aes(iso.prob, act.rate.iso)) +
  geom_raster(aes(fill = cIso), interpolate = TRUE) +
  geom_contour(aes(z = cIso), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cIso), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability", y = "Relative contact intensity") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.v, maxval.v), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot isolation - Office Workers
minval.vw <- min(cIsoWork$cIso)
maxval.vw <- max(cIsoWork$cIso)
g_iso_cont_w <- ggplot(cIsoWork, aes(iso.prob, act.rate.iso)) +
  geom_raster(aes(fill = cIso), interpolate = TRUE) +
  geom_contour(aes(z = cIso), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cIso), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability", y = "Relative contact intensity") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.vw, maxval.vw), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

combined1 <- g_iso_cont + g_iso_cont_w +
  plot_annotation(tag_levels = list(c("A1", "A2")), title = "Total Persons Entering Isolation")
combined2 <- g_inc_cont + g_inc_cont_w +
  plot_annotation(tag_levels = list(c("B1", "B2")), title = "Infections per 100'000 Person-Days")
combined3 <- g_death_cont +
  g_death_cont_w +
  plot_annotation(tag_levels = list(c("C1", "C2")), title = "Deaths per 100'000 Person-Days")

wrap_elements(combined1) / wrap_elements(combined2) / wrap_elements(combined3)

ggsave(here::here("results","figures","Figure3B.tiff"),
       width = 12, height = 12, units = "in", dpi=300, compression = 'lzw')

##################### VARYING ISOLATION AMONG OFFICE WORKERS ONLY ##############
################################################################################

# Setup data (RUN THIS FIRST) --------

d.ratio.ow <- readRDS(here::here("results","ratio_df_18000.rds"))

d.ratio.ow <- d.ratio.ow %>%
  filter(scenario2 %in% c(0,2,4,5)) %>%
  mutate(scenario2 = factor(scenario2,
                            levels = c(2,0,4,5),
                            labels = c("Worse adherence (25%)",
                                       "Reference (58%)",
                                       "Improved adherence (75%)",
                                       "Perfect adherence (100%)")))

# Start date
start_date <- as.Date("2021-01-01")

# Convert numeric days to date format by adding to start date
date_variable <- start_date + d.ratio.ow$time - 1

d.ratio.ow <- cbind(d.ratio.ow, date_variable)

# Rename the new column
colnames(d.ratio.ow)[ncol(d.ratio.ow)] <- "date_variable"

# custom palette
custom_palette <- gradient_n_pal(viridis::viridis_pal(option = "magma")(100)[0:90])

# Ratio Plot: Cases to Isolated -------------

p.sum2.ow <- ggplot(d.ratio.ow, aes(x = date_variable, y = med.sum2))+
  geom_line(aes(y = med.sum2, color = scenario2, group = scenario2),
            linewidth = 0.8) +
  geom_ribbon(aes(ymin = p25.sum2, ymax = p75.sum2, group = scenario2,
                  fill = scenario2), alpha = 0.1)+
  scale_color_manual(values = custom_palette(seq(0, 1, length.out = n_distinct(d.ratio.tot$scenario2))))+
  scale_fill_manual(values = custom_palette(seq(0, 1, length.out = n_distinct(d.ratio.tot$scenario2))))+
  scale_x_date(date_breaks = "3 months")+
  ylim(0,1500)+
  labs(x = "Date",
       y = "Sum of cases and persons in isolation",
       color = "Scenario") +
  guides(fill = "none", color = "none")+
  ggtitle("Total number of office workers ill or in isolation") + theme_classic()

p.cases.ow <- ggplot(d.ratio.ow, aes(x = date_variable, y = med.case))+
  geom_line(aes(y = med.case, color = scenario2, group = scenario2),
            linewidth = 0.8) +
  geom_ribbon(aes(ymin = p25.case, ymax = p75.case, group = scenario2,
                  fill = scenario2), alpha = 0.1)+
  scale_color_manual(values = custom_palette(seq(0, 1, length.out = n_distinct(d.ratio.tot$scenario2))))+
  scale_fill_manual(values = custom_palette(seq(0, 1, length.out = n_distinct(d.ratio.tot$scenario2))))+
  scale_x_date(date_breaks = "3 months")+
  ylim(0,1500)+
  labs(x = "Date",
       y = "Number of cases",
       color = "Scenario") +
  guides(fill = "none", color = "none")+
  ggtitle("Total number of office workers ill") + theme_classic()

p.iso.ow <- ggplot(d.ratio.ow, aes(x = date_variable, y = med.iso))+
  geom_line(aes(y = med.iso, color = scenario2, group = scenario2),
            linewidth = 0.8) +
  geom_ribbon(aes(ymin = p25.iso, ymax = p75.iso, group = scenario2,
                  fill = scenario2), alpha = 0.1)+
  scale_color_manual(values = custom_palette(seq(0, 1, length.out = n_distinct(d.ratio.tot$scenario2))))+
  scale_fill_manual(values = custom_palette(seq(0, 1, length.out = n_distinct(d.ratio.tot$scenario2))))+
  scale_x_date(date_breaks = "3 months")+
  ylim(0,1500)+
  labs(x = "Date",
       y = "Number of persons in isolation",
       color = "Scenario") +
  guides(fill = "none", color = "none")+
  ggtitle("Total number of office workers in isolation") + theme_classic()

p.mask.ow <- ggplot(d.ratio.ow, aes(x = date_variable, y = med.mask))+
  geom_line(aes(y = med.mask, color = scenario2, group = scenario2),
            linewidth = 0.8) +
  geom_ribbon(aes(ymin = p25.mask, ymax = p75.mask, group = scenario2,
                  fill = scenario2), alpha = 0.1)+
  scale_color_manual(values = custom_palette(seq(0, 1, length.out = n_distinct(d.ratio.tot$scenario2))))+
  scale_fill_manual(values = custom_palette(seq(0, 1, length.out = n_distinct(d.ratio.tot$scenario2))))+
  scale_x_date(date_breaks = "3 months")+
  ylim(0,1500)+
  labs(x = "Date",
       y = "Number of persons masked",
       color = "Scenario") +
  guides(fill = "none")+
  ggtitle("Total number of office workers masked") + theme_classic()+
  theme(legend.position = c(0.83,0.8))

patch.ow <- (p.sum2.ow | p.mask.ow) / (p.cases.ow | p.iso.ow)
patch.ow <- patch.ow + plot_annotation(
  title = 'Impact of varying adherence to CDC guidelines among office workers',
  caption = 'Adherence among non-office workers is fixed at the reference value (58%) across scenarios'
)
patch.ow

ggsave(here::here("results","figures","Figure4B.tiff"),
       width = 15, height = 6, units = "in", dpi=300, compression = 'lzw')

# Stacked area plot: Cases vs isolating ----
d.long <- d.ratio.ow %>%
  select(scenario2,date_variable,med.casenotiso,med.both,med.isonotcase) %>%
  pivot_longer(cols = -c("scenario2","date_variable"), names_to = "cat", values_to = "number") %>%
  mutate(cat = case_when(cat == "med.casenotiso" ~ "Cases not isolating",
                         cat == "med.both" ~ "Cases isolating",
                         cat == "med.isonotcase" ~ "Non-cases isolating")) %>%
  mutate(cat = factor(cat,levels = c("Non-cases isolating",
                                     "Cases isolating",
                                     "Cases not isolating")))

ggplot(d.long, aes(x = date_variable, y = number, fill = cat)) +
  geom_area(position = "stack")+
  facet_wrap(vars(scenario2),nrow = 2)+
  theme_bw()+
  scale_fill_viridis(discrete = TRUE)+
  labs(x = "Date",
       y = "Number of office workers ill and/or isolating",
       fill = "Category",
       title = "Impact of varying adherence among office workers")

ggsave(here::here("results","figures","Figure5B.tiff"),
       width = 8, height = 4, units = "in", dpi=300, compression = 'lzw')

# ggplot(d.long, aes(x = date_variable, y = number, fill = cat)) +
#   geom_area(position = "fill")+
#   facet_wrap(vars(scenario2),nrow = 2)+
#   theme_bw()+
#   scale_fill_viridis(discrete = TRUE)+
#   labs(x = "Date",
#        y = "Proportion of office workers ill and/or isolating",
#        fill = "Category")
#
# ggsave(here::here("results","figures","Figure5B_prop.tiff"),
#        width = 12, height = 6, units = "in", dpi=300, compression = 'lzw')

# Cumulative cases vs isolations in individual simulations ----
d.sim <- readRDS(here::here("results","sim_df2.rds")) %>%
  filter(iso.prob.other.param == 0.58 &
           iso.prob.office.param %in% c(0.25,0.58,0.75,1)) %>%
  mutate(scenario = case_when(iso.prob.office.param == 0.25 ~ "Worse adherence (25%)",
                              iso.prob.office.param == 0.58 ~ "Reference (58%)",
                              iso.prob.office.param == 0.75 ~ "Improved adherence (75%)",
                              iso.prob.office.param == 1 ~ "Perfect adherence (100%)")) %>%
  mutate(scenario = factor(scenario,
                           levels = c("Worse adherence (25%)",
                                      "Reference (58%)",
                                      "Improved adherence (75%)",
                                      "Perfect adherence (100%)")))

p.sim1 <- ggplot(d.sim, aes(x = cum.case, y = cum.iso, color = factor(scenario))) +
  geom_point(alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  scale_color_viridis(discrete = TRUE)+
  guides(color = "none")+
  labs(x = "Cumulative cases", y = "Cumulative isolations",
       title = "Total population") +
  theme_bw()

p.sim2 <- ggplot(d.sim, aes(x = cum.case.w, y = cum.iso.w, color = factor(scenario))) +
  geom_point(alpha = 0.3) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  scale_color_viridis(discrete = TRUE)+
  labs(x = "Cumulative cases", y = "Cumulative isolations", color = "Scenario",
       title = "Office workers") +
  theme_bw()

combined.sim <- p.sim1 + p.sim2 +
  plot_annotation(tag_levels = list(c("A", "B")),
                  title = "Cumulative cases and persons who isolated at the end of the time period",
                  subtitle = "Varying adherence among office workers only (others at 58%)")
combined.sim

ggsave(here::here("results","figures","Figure6B.tiff"),
       width = 10, height = 5, units = "in", dpi=300)
################################################################################
# Contour plot libraries -------------------------------------------------------
library("viridis")
library("ggplot2")
library("metR")
library("patchwork")

# Load contour plot data (iso prob office vs work notif) ------------------------------

cInc <- readRDS(here::here("data","output","results","WcInc.rds"))
cDeaths <- readRDS(here::here("data","output","results","WcDeaths.rds"))
cIso <- readRDS(here::here("data","output","results","WcIso.rds"))
cCases <- readRDS(here::here("data","output","results","WcCases.rds"))
cIncWork <- readRDS(here::here("data","output","results","WcIncWork.rds"))
cDeathsWork <- readRDS(here::here("data","output","results","WcDeathsWork.rds"))
cIsoWork <- readRDS(here::here("data","output","results","WcIsoWork.rds"))
cCasesWork <- readRDS(here::here("data","output","results","WcCasesWork.rds"))

# Contour Plots: Iso prob vs work notif ----------------------------------------

# Plot cases
minval.i <- floor(min(cInc$cInc))
maxval.i <- ceiling(max(cInc$cInc))
g_inc_cont <- ggplot(cInc, aes(iso.prob.off, work.notif.prob)) +
  geom_raster(aes(fill = cInc), interpolate = TRUE) +
  geom_contour(aes(z = cInc), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cInc), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability\namong office workers", y = "Notification probability at work") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.i, maxval.i), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot cases - Office Workers
minval.iw <- floor(min(cIncWork$cInc))
maxval.iw <- ceiling(max(cIncWork$cInc))
g_inc_cont_w <- ggplot(cIncWork, aes(iso.prob.off, work.notif.prob)) +
  geom_raster(aes(fill = cInc), interpolate = TRUE) +
  geom_contour(aes(z = cInc), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cInc), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability\namong office workers", y = "Notification probability at work") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.iw, maxval.iw), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot deaths
minval.d <- min(cDeaths$cDeaths)
maxval.d <- max(cDeaths$cDeaths)
g_death_cont <- ggplot(cDeaths, aes(iso.prob.off, work.notif.prob)) +
  geom_raster(aes(fill = cDeaths), interpolate = TRUE) +
  geom_contour(aes(z = cDeaths), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cDeaths), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability\namong office workers", y = "Notification probability at work") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.d, maxval.d), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot deaths - Office Workers
minval.dw <- min(cDeathsWork$cDeaths)
maxval.dw <- max(cDeathsWork$cDeaths)
g_death_cont_w <- ggplot(cDeathsWork, aes(iso.prob.off, work.notif.prob)) +
  geom_raster(aes(fill = cDeaths), interpolate = TRUE) +
  geom_contour(aes(z = cDeaths), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cDeaths), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability\namong office workers", y = "Notification probability at work") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.dw, maxval.dw), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot isolation
minval.v <- min(cIso$cIso)
maxval.v <- max(cIso$cIso)
g_iso_cont <- ggplot(cIso, aes(iso.prob.off, work.notif.prob)) +
  geom_raster(aes(fill = cIso), interpolate = TRUE) +
  geom_contour(aes(z = cIso), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cIso), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability\namong office workers", y = "Notification probability at work") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.v, maxval.v), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot isolation - Office Workers
minval.vw <- min(cIsoWork$cIso)
maxval.vw <- max(cIsoWork$cIso)
g_iso_cont_w <- ggplot(cIsoWork, aes(iso.prob.off, work.notif.prob)) +
  geom_raster(aes(fill = cIso), interpolate = TRUE) +
  geom_contour(aes(z = cIso), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cIso), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability\namong office workers", y = "Notification probability at work") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.vw, maxval.vw), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

combined1 <- g_iso_cont + g_iso_cont_w +
  plot_annotation(tag_levels = list(c("A1", "A2")), title = "Total Persons Entering Isolation")
combined2 <- g_inc_cont + g_inc_cont_w +
  plot_annotation(tag_levels = list(c("B1", "B2")), title = "Infections per 100'000 Person-Days")
combined3 <- g_death_cont +
  g_death_cont_w +
  plot_annotation(tag_levels = list(c("C1", "C2")), title = "Deaths per 100'000 Person-Days")

wrap_elements(combined1) / wrap_elements(combined2) / wrap_elements(combined3)

ggsave(here::here("results","figures","Figure3A_office.tiff"),
       width = 12, height = 12, units = "in", dpi=300, compression = 'lzw')


# Load contour plot data (iso prob office vs iso act) ---------------------------------

cInc <- readRDS(here::here("data","output","results","WcInc2.rds"))
cDeaths <- readRDS(here::here("data","output","results","WcDeaths2.rds"))
cIso <- readRDS(here::here("data","output","results","WcIso2.rds"))
cCases <- readRDS(here::here("data","output","results","WcCases2.rds"))
cIncWork <- readRDS(here::here("data","output","results","WcIncWork2.rds"))
cDeathsWork <- readRDS(here::here("data","output","results","WcDeathsWork2.rds"))
cIsoWork <- readRDS(here::here("data","output","results","WcIsoWork2.rds"))
cCasesWork <- readRDS(here::here("data","output","results","WcCasesWork2.rds"))

# Contour Plots: Iso prob vs iso act -------------------------------------------

# Plot cases
minval.i <- floor(min(cInc$cInc))
maxval.i <- ceiling(max(cInc$cInc))
g_inc_cont <- ggplot(cInc, aes(iso.prob.off, act.rate.iso)) +
  geom_raster(aes(fill = cInc), interpolate = TRUE) +
  geom_contour(aes(z = cInc), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cInc), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability\namong office workers", y = "Relative contact intensity") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.i, maxval.i), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot cases - Office Workers
minval.iw <- floor(min(cIncWork$cInc))
maxval.iw <- ceiling(max(cIncWork$cInc))
g_inc_cont_w <- ggplot(cIncWork, aes(iso.prob.off, act.rate.iso)) +
  geom_raster(aes(fill = cInc), interpolate = TRUE) +
  geom_contour(aes(z = cInc), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cInc), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability\namong office workers", y = "Relative contact intensity") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.iw, maxval.iw), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot deaths
minval.d <- min(cDeaths$cDeaths)
maxval.d <- max(cDeaths$cDeaths)
g_death_cont <- ggplot(cDeaths, aes(iso.prob.off, act.rate.iso)) +
  geom_raster(aes(fill = cDeaths), interpolate = TRUE) +
  geom_contour(aes(z = cDeaths), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cDeaths), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability\namong office workers", y = "Relative contact intensity") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.d, maxval.d), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot deaths - Office Workers
minval.dw <- min(cDeathsWork$cDeaths)
maxval.dw <- max(cDeathsWork$cDeaths)
g_death_cont_w <- ggplot(cDeathsWork, aes(iso.prob.off, act.rate.iso)) +
  geom_raster(aes(fill = cDeaths), interpolate = TRUE) +
  geom_contour(aes(z = cDeaths), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cDeaths), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability\namong office workers", y = "Relative contact intensity") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.dw, maxval.dw), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot isolation
minval.v <- min(cIso$cIso)
maxval.v <- max(cIso$cIso)
g_iso_cont <- ggplot(cIso, aes(iso.prob.off, act.rate.iso)) +
  geom_raster(aes(fill = cIso), interpolate = TRUE) +
  geom_contour(aes(z = cIso), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cIso), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability\namong office workers", y = "Relative contact intensity") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.v, maxval.v), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot isolation - Office Workers
minval.vw <- min(cIsoWork$cIso)
maxval.vw <- max(cIsoWork$cIso)
g_iso_cont_w <- ggplot(cIsoWork, aes(iso.prob.off, act.rate.iso)) +
  geom_raster(aes(fill = cIso), interpolate = TRUE) +
  geom_contour(aes(z = cIso), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cIso), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability\namong office workers", y = "Relative contact intensity") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.vw, maxval.vw), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

combined1 <- g_iso_cont + g_iso_cont_w +
  plot_annotation(tag_levels = list(c("A1", "A2")), title = "Total Persons Entering Isolation")
combined2 <- g_inc_cont + g_inc_cont_w +
  plot_annotation(tag_levels = list(c("B1", "B2")), title = "Infections per 100'000 Person-Days")
combined3 <- g_death_cont +
  g_death_cont_w +
  plot_annotation(tag_levels = list(c("C1", "C2")), title = "Deaths per 100'000 Person-Days")

wrap_elements(combined1) / wrap_elements(combined2) / wrap_elements(combined3)

ggsave(here::here("results","figures","Figure3B_office.tiff"),
       width = 12, height = 12, units = "in", dpi=300, compression = 'lzw')

# Load contour plot data (iso prob office vs iso prob other) ---------------------------------

cInc <- readRDS(here::here("data","output","results","WcInc3.rds"))
cDeaths <- readRDS(here::here("data","output","results","WcDeaths3.rds"))
cIso <- readRDS(here::here("data","output","results","WcIso3.rds"))
cCases <- readRDS(here::here("data","output","results","WcCases3.rds"))
cIncWork <- readRDS(here::here("data","output","results","WcIncWork3.rds"))
cDeathsWork <- readRDS(here::here("data","output","results","WcDeathsWork3.rds"))
cIsoWork <- readRDS(here::here("data","output","results","WcIsoWork3.rds"))
cCasesWork <- readRDS(here::here("data","output","results","WcCasesWork3.rds"))

# Contour Plots: Iso prob office vs iso prob other -------------------------------------------

# Plot cases
minval.i <- floor(min(cInc$cInc))
maxval.i <- ceiling(max(cInc$cInc))
g_inc_cont <- ggplot(cInc, aes(iso.prob.off, iso.prob.non)) +
  geom_raster(aes(fill = cInc), interpolate = TRUE) +
  geom_contour(aes(z = cInc), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cInc), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability\namong office workers", y = "Adherence probability\namong non-office workers") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.i, maxval.i), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot cases - Office Workers
minval.iw <- floor(min(cIncWork$cInc))
maxval.iw <- ceiling(max(cIncWork$cInc))
g_inc_cont_w <- ggplot(cIncWork, aes(iso.prob.off, iso.prob.non)) +
  geom_raster(aes(fill = cInc), interpolate = TRUE) +
  geom_contour(aes(z = cInc), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cInc), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability\namong office workers", y = "Adherence probability\namong non-office workers") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.iw, maxval.iw), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot deaths
minval.d <- min(cDeaths$cDeaths)
maxval.d <- max(cDeaths$cDeaths)
g_death_cont <- ggplot(cDeaths, aes(iso.prob.off, iso.prob.non)) +
  geom_raster(aes(fill = cDeaths), interpolate = TRUE) +
  geom_contour(aes(z = cDeaths), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cDeaths), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability\namong office workers", y = "Adherence probability\namong non-office workers") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.d, maxval.d), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot deaths - Office Workers
minval.dw <- min(cDeathsWork$cDeaths)
maxval.dw <- max(cDeathsWork$cDeaths)
g_death_cont_w <- ggplot(cDeathsWork, aes(iso.prob.off, iso.prob.non)) +
  geom_raster(aes(fill = cDeaths), interpolate = TRUE) +
  geom_contour(aes(z = cDeaths), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cDeaths), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability\namong office workers", y = "Adherence probability\namong non-office workers") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.dw, maxval.dw), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot isolation
minval.v <- min(cIso$cIso)
maxval.v <- max(cIso$cIso)
g_iso_cont <- ggplot(cIso, aes(iso.prob.off, iso.prob.non)) +
  geom_raster(aes(fill = cIso), interpolate = TRUE) +
  geom_contour(aes(z = cIso), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cIso), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability\namong office workers", y = "Adherence probability\namong non-office workers") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.v, maxval.v), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot isolation - Office Workers
minval.vw <- min(cIsoWork$cIso)
maxval.vw <- max(cIsoWork$cIso)
g_iso_cont_w <- ggplot(cIsoWork, aes(iso.prob.off, iso.prob.non)) +
  geom_raster(aes(fill = cIso), interpolate = TRUE) +
  geom_contour(aes(z = cIso), col = "white", alpha = 0.5, lwd = 0.5) +
  geom_text_contour(aes(z = cIso), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Adherence probability\namong office workers", y = "Adherence probability\namong non-office workers") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.vw, maxval.vw), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

combined1 <- g_iso_cont + g_iso_cont_w +
  plot_annotation(tag_levels = list(c("A1", "A2")), title = "Total Persons Entering Isolation")
combined2 <- g_inc_cont + g_inc_cont_w +
  plot_annotation(tag_levels = list(c("B1", "B2")), title = "Infections per 100'000 Person-Days")
combined3 <- g_death_cont +
  g_death_cont_w +
  plot_annotation(tag_levels = list(c("C1", "C2")), title = "Deaths per 100'000 Person-Days")

wrap_elements(combined1) / wrap_elements(combined2) / wrap_elements(combined3)

ggsave(here::here("results","figures","Figure3C_office.tiff"),
       width = 12, height = 12, units = "in", dpi=300, compression = 'lzw')

##################### COMBINING SOME PLOTS ##############
################################################################################
# setup ----
my.colors <- c("#641A80", "#DD4969", "#FECD90")

d.ratio.tot.ref <- d.ratio.tot %>%
  filter(scenario2 == "Reference (58%)") %>%
  mutate(scenario2 = as.character(scenario2)) %>%
  distinct()

d.ratio.ow.ref <- d.ratio.ow %>%
  filter(scenario2 == "Reference (58%)") %>%
  mutate(scenario2 = as.character(scenario2)) %>%
  distinct()

# plot for cases/isolation over time -----

# some edits to plots
p.cases.tot <- ggplot(d.ratio.tot %>%
                        filter(scenario2 != "Reference (58%)"),
                      aes(x = date_variable, y = med.case))+
  geom_line(aes(y = med.case, color = scenario2, group = scenario2),
            linewidth = 0.8) +
  geom_ribbon(aes(ymin = p25.case, ymax = p75.case, group = scenario2,
                  fill = scenario2), alpha = 0.1)+
  geom_line(data = d.ratio.tot.ref,
            aes(x = date_variable,y = med.case,linetype = "Reference (58%)"),
            color = "black",linewidth = 0.8)+
  geom_ribbon(data = d.ratio.tot.ref,
              aes(x = date_variable,
                  ymin = p25.case, ymax = p75.case),
              alpha = 0.1, fill = "black")+
  scale_x_date(date_breaks = "3 months")+
  scale_color_manual(values = my.colors)+
  scale_fill_manual(values = my.colors)+
  scale_linetype_manual(name = "", values = c("Reference (58%)" = "dashed"))+
  ylim(0,1500)+
  labs(x = "Date",
       y = "Office worker cases",
       color = "Scenario") +
  guides(fill = "none", color = "none", linetype = "none")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p.cases.tot

p.iso.tot <- ggplot(d.ratio.tot %>%
                      filter(scenario2 != "Reference (58%)"),
                    aes(x = date_variable, y = med.iso))+
  geom_line(aes(y = med.iso, color = scenario2, group = scenario2),
            linewidth = 0.8) +
  geom_ribbon(aes(ymin = p25.iso, ymax = p75.iso, group = scenario2,
                  fill = scenario2), alpha = 0.1)+
  geom_line(data = d.ratio.tot.ref,
            aes(x = date_variable,y = med.iso,linetype = "Reference (58%)"),
            color = "black",linewidth = 0.8)+
  geom_ribbon(data = d.ratio.tot.ref,
              aes(x = date_variable,
                  ymin = p25.iso, ymax = p75.iso),
              alpha = 0.1, fill = "black")+
  scale_x_date(date_breaks = "3 months")+
  scale_color_manual(values = my.colors)+
  scale_fill_manual(values = my.colors)+
  scale_linetype_manual(name = "", values = c("Reference (58%)" = "dashed"))+
  ylim(0,1500)+
  labs(x = "Date",
       y = "Office workers in isolation",
       color = "Scenario") +
  guides(fill = "none", color = "none", linetype = "none")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p.cases.ow <- ggplot(d.ratio.ow %>%
                       filter(scenario2 != "Reference (58%)"),
                     aes(x = date_variable, y = med.case))+
  geom_line(aes(y = med.case, color = scenario2, group = scenario2),
            linewidth = 0.8) +
  geom_ribbon(aes(ymin = p25.case, ymax = p75.case, group = scenario2,
                  fill = scenario2), alpha = 0.1)+
  geom_line(data = d.ratio.ow.ref,
            aes(x = date_variable,y = med.case,linetype = "Reference (58%)"),
            color = "black",linewidth = 0.8)+
  geom_ribbon(data = d.ratio.ow.ref,
              aes(x = date_variable,
                  ymin = p25.case, ymax = p75.case),
              alpha = 0.1, fill = "black")+
  scale_color_manual(values = my.colors)+
  scale_fill_manual(values = my.colors)+
  scale_linetype_manual(name = "", values = c("Reference (58%)" = "dashed"))+
  scale_x_date(date_breaks = "3 months")+
  ylim(0,1500)+
  labs(x = "Date",
       y = "Office worker cases",
       color = "Scenario") +
  guides(fill = "none", color = "none", linetype = "none")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p.iso.ow <- ggplot(d.ratio.ow %>%
                     filter(scenario2 != "Reference (58%)"),
                   aes(x = date_variable, y = med.iso))+
  geom_line(aes(y = med.iso, color = scenario2, group = scenario2),
            linewidth = 0.8) +
  geom_ribbon(aes(ymin = p25.iso, ymax = p75.iso, group = scenario2,
                  fill = scenario2), alpha = 0.1)+
  geom_line(data = d.ratio.ow.ref,
            aes(x = date_variable,y = med.iso,linetype = "Reference (58%)"),
            color = "black",linewidth = 0.8)+
  geom_ribbon(data = d.ratio.ow.ref,
              aes(x = date_variable,
                  ymin = p25.iso, ymax = p75.iso),
              alpha = 0.1, fill = "black")+
  scale_color_manual(values = my.colors)+
  scale_fill_manual(values = my.colors)+
  scale_linetype_manual(name = "", values = c("Reference (58%)" = "dashed"))+
  scale_x_date(date_breaks = "3 months")+
  ylim(0,1500)+
  labs(x = "Date",
       y = "Office workers in isolation",
       color = "Scenario") +
  guides(fill = "none")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.box = "horizontal")
p.iso.ow

# Extract the legend from one plot (e.g., `p.iso.ow`)
legend_plot <- get_legend(p.iso.ow)

# Remove the legend from `p.iso.ow`
p.iso.ow <- p.iso.ow + guides(fill = "none", color = "none", linetype = "none")

# combine plots
patch1 <- (p.iso.tot) / (p.cases.tot ) +
  plot_annotation(tag_levels = list(c("A1", "A2")), title = "Varying adherence among total population")

patch2 <- (p.iso.ow) / (p.cases.ow ) +
  plot_annotation(tag_levels = list(c("B1", "B2")), title = "Varying adherence among office workers")

patch <- (wrap_elements(patch1) | wrap_elements(patch2)) /
  wrap_elements(legend_plot) &
  theme(legend.position = "bottom")

patch <- ((wrap_elements(patch1) | wrap_elements(patch2)) / wrap_elements(legend_plot)) +
  plot_layout(heights = c(1, 0.2))  # Adjust legend height proportion as needed

patch

ggsave(here::here("results","figures","Figure_timevar.tiff"),
       width = 13, height = 9, units = "in", dpi=300, compression = 'lzw')
