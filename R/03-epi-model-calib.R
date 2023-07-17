
##
## COVID-19 Corporate Network Model
## Epidemic Model
##
## Authors: Samuel M. Jenness
## Date: February 2021
##

library("EpiModelCOVID")


# # Read in fitted network models
# est <- readRDS("data/input/est.simple.rds")

source("R/00-calib-targets.R")
source("R/01-epi-params.R")
source("R/02b-household-assign.R")
rm(list=setdiff(ls(), c("est", "hhPairs","mr_vec","vax_targets")))

# Model parameters
param <- param.net(inf.prob = c(0.5 ,0.5, 0.5),
                   inf.add = c(0.15, 0.25), # calibrated
                   inf.sub = c(0.2, 0.2, 0.2),
                   inf.boost.start = 180 + c(-60, 335), # calibrated
                   inf.boost.stop = 180 + c(15, 396), # calibrated
                   inf.supp.start =  180 + c(16, 258, 397), # calibrated
                   inf.supp.stop = 180 + c(180, 334, 640), # calibrated

                   act.rate = c(1, 1, 5),
                   inf.prob.inter.rr = c(0.6, 0.6, 0.6),
                   inf.prob.inter.time = c(Inf, Inf, Inf),
                   act.rate.inter.rr = c(1, 1, 1),
                   act.rate.inter.time = c(Inf, Inf, Inf),
                   inf.prob.a.rr = 0.5,
                   inf.prob.mask.rr = 0.3,

                   # iso.prob = 0.5,
                   # notif.prob = c(0.5,0.3,0.9),
                   iso.prob = 0,
                   notif.prob = c(0,0,0),

                   # act.rate.dx.inter.rr = 0.1,
                   # act.rate.dx.inter.time = Inf,
                   # act.rate.sympt.inter.rr = 0.5,
                   # act.rate.sympt.inter.time = Inf,
                   act.rate.iso.inter.time = 1,
                   act.rate.iso.inter.rr = 0.3,

                   prop.clinical = c(0.573, 0.642, 0.760, 0.800, 0.813, 0.814, 0.769, 0.723, 0.666), # by decade of age
                   prop.hospit = c(0.060, 0.063, 0.081, 0.154, 0.207, 0.268, 0.357, 0.465, 0.539) / 10, # by decade of age

                   hosp.boost.mult = 2.25,
                   hosp.boost.start = 378,
                   hosp.boost.stop = 423,
                   hosp.supp.mult = 0.75,
                   hosp.supp.start = 546,
                   hosp.supp.stop = 604,

                   ea.rate = 1 / 5.5, # from https://academic.oup.com/cid/article/74/9/1678/6359063 (duration of latent period)
                   ar.rate = 1/5.0,
                   eip.rate = 1 / 5.5, # from https://academic.oup.com/cid/article/74/9/1678/6359063 (duration of latent period)
                   ipic.rate = 1/1.5,
                   icr.rate = 1/3.5,
                   ich.rate = 1 / 3,
                   hr.rate = 1 / 15.5, # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7899024/
                   rs.rate = 1 / 300, # from https://www.cdc.gov/coronavirus/2019-ncov/science/science-briefs/vaccine-induced-immunity.html

                   pcr.sens = 0.8,
                   dx.rate.sympt = 0.1,
                   dx.rate.other = 0.01,
                   allow.rescreen = TRUE,

                   vax1.start = 180 + c(534, 132, 84, 74, 11), # age 0-4, 5-17, 18-49, 50 - 64, 65+
                   vax3.start = 180 + c(Inf, 368, 323, 323, 265), # age 0-4, 5-17, 18-49, 50 - 64, 65+
                   vax4.start = 180 + c(Inf, Inf, Inf, 453, 453), # age 0-4, 5-17, 18-49, 50 - 64, 65+
                   vax1.rate = c(0.0005, 0.0012, 0.030, 0.030, 0.022), # manually calibrated; age 0-4, 5-17, 18-49, 50 - 64, 65+
                   vax2.rate = c(0.01, 0.015, 0.40, 0.40, 0.55), # manually calibrated; age 0-4, 5-17, 18-49, 50 - 64, 65+
                   vax3.rate = c(0, 0.003, 0.10, 0.10, 0.05), # manually calibrated; age 0-4, 5-17, 18-49, 50 - 64, 65+
                   vax4.rate = c(0, 0, 0, 0.015, 0.015), # manually calibrated; age 0-4, 5-17, 18-49, 50 - 64, 65+
                   vax2.interval = 21, # 3 weeks
                   vax3.interval = 180, # 6 months
                   vax4.interval = 120, # 4 months
                   vax1.rr.infect = 0.324,
                   vax2.rr.infect = 0.112,
                   vax3.rr.infect = 0.12,
                   vax4.rr.infect = 0.12,
                   vax1.rr.clinical = 0.4,
                   vax2.rr.clinical = 0.09,
                   vax3.rr.clinical = 0.09,
                   vax4.rr.clinical = 0.09,
                   vax1.rr.hosp = 0.3,
                   vax2.rr.hosp = 0.02,
                   vax3.rr.hosp = 0.07,
                   vax4.rr.hosp = 0.07,
                   half.life = 80,

                   a.rate = mean(mr_vec),
                   arrival.age = 0,
                   mort.rates = mr_vec,
                   mort.dis.mult = 180 * 10,

                   hh.pairs = hhPairs)

# Initial conditions
init <- init.net(e.num = 1600)

# Control settings
# pkgload::load_all("~/git/EpiModelCOVID")
pkgload::load_all("C:/Users/mgarc48/git/EpiModelCOVID")
control <- control.net(nsteps = 180 + 608,
                       nsims = 1,
                       ncores = 1,
                       initialize.FUN = init_covid_corporate,
                       aging.FUN = aging_covid,
                       departures.FUN = deaths_covid_corporate,
                       arrivals.FUN = arrival_covid_corporate,
                       resim_nets.FUN = resim_nets_covid_corporate,
                       infection.FUN = infect_covid_corporate,
                       recovery.FUN = progress_covid,
                       dx.FUN = dx_covid,
                       vax.FUN = vax_covid,
                       prevalence.FUN = prevalence_covid_corporate,
                       module.order = c("aging.FUN",
                                        "departures.FUN",
                                        "arrivals.FUN",
                                        "resim_nets.FUN",
                                        "infection.FUN",
                                        "recovery.FUN",
                                        "dx.FUN",
                                        "vax.FUN",
                                        "prevalence.FUN"),
                       resimulate.network = TRUE,
                       save.nwstats = TRUE,
                       nwstats.formula.1 = "formation",
                       nwstats.formula.2 = "formation",
                       nwstats.formula.3 = "formation",
                       tergmLite = TRUE,
                       cumulative.edgelist = TRUE,
                       truncate.el.cuml = 5)

sim <- netsim(est, param, init, control)
# print(sim)

plot(sim, y = c("s.num", "e.num", "a.num"))
plot(sim, y = c("e.num", "a.num","h.num"))
plot(sim, y = c("iso.mild.num","iso.sev.num","mask.post.num","mask.exp.num"))
plot(sim, y = c("v1.num","v2.num","v3.num","v4.num"))

sim <- mutate_epi(sim, se.cuml = cumsum(se.flow),
                       dx.cuml = cumsum(nDx),
                       dx.pos.cuml = cumsum(nDx.pos),
                       totI = e.num + a.num + ip.num + ic.num)

plot(sim, y = c("se.cuml"))
df <- as.data.frame(sim, out = "mean")
names(df)
head(df)

get_nwstats(sim)

df2 <- df %>%
  filter(time >=180)

### CALIBRATION ################################################################
library(patchwork)

times <- data.frame(DayNum = c(608, 577, 546, 516, 485, 455, 424, 396, 365, 334,
                               304, 273, 243, 212, 181, 151, 120, 90, 59, 31),
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

sim180 <- truncate_sim(sim, 181)

# Cases -------------------------------------------------------------------

case_targets <- c(8373, 3165, 1590, 1301, 784, 382, 1439, 7256, 6079,
                  1981, 958, 6259, 19544, 3194, 1278, 701, 1749, 2988,
                  4000, 3504)

cases <- sim180[["epi"]][["se.flow"]]
cases <- data.frame(ts = 1:length(cases), cases = cases)
cases$cumCases <- cumsum(cases$cases)

cases <- merge(times, cases, by.x = "DayNum", by.y = "ts")
cases$cases <- cases$cumCases - lag(cases$cumCases)
cases$cases[1] <- cases$cumCases[1]
cases <- cases %>% rowid_to_column("ROW")

cases$target <- case_targets

f1 <- ggplot(data = cases, aes(x = ROW, y = cases)) +
  geom_line(aes(color = "Result")) +
  geom_line(aes(y = target, color = "Target (Empirical)")) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(1, 20),
                     breaks = c(1, 4, 7, 10, 13, 16, 19),
                     labels = c("1" = "Jan 21", "4" = "Apr 21", "7" = "Jul 21",
                                "10" = "Oct 21", "13" = "Jan 22", "16" = "Apr 22",
                                "19" = "Jul 22")) +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 26500)) +
  xlab("Month") +
  ylab("Incident Infections") +
  scale_color_manual(values = colors1, name = "") +
  ggtitle("Infections") + theme_classic()
f1

# Hospitalizations --------------------------------------------------------

hosp_targets <- c(113, 47, 27, 24, 16, 9, 27, 98, 66, 24, 12, 44, 88, 31, 9, 5,
                  9, 19, 29, 26)

hosp <- sim180[["epi"]][["ich.flow"]]
hosp <- data.frame(ts = 1:length(hosp), hosp = hosp)
hosp$cumHosp <- cumsum(hosp$hosp)

hosp <- merge(times, hosp, by.x = "DayNum", by.y = "ts")
hosp$hosp <- hosp$cumHosp - lag(hosp$cumHosp)
hosp$hosp[1] <- hosp$cumHosp[1]
hosp <- hosp %>% rowid_to_column("ROW")

hosp$target <- hosp_targets

f2 <- ggplot(data = hosp, aes(x = ROW, y = hosp)) +
  geom_line(aes(color = "Result")) +
  geom_line(aes(y = target, color = "Target (Empirical)")) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(1, 20),
                     breaks = c(1, 4, 7, 10, 13, 16, 19),
                     labels = c("1" = "Jan 21", "4" = "Apr 21", "7" = "Jul 21",
                                "10" = "Oct 21", "13" = "Jan 22", "16" = "Apr 22",
                                "19" = "Jul 22")) +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 150)) +
  xlab("Month") +
  ylab("Incident Hospitalizations") +
  scale_color_manual(values = colors1, name = "") +
  ggtitle("Hospitalizations") + theme_classic()
f2

# Deaths ------------------------------------------------------------------

death_targets <- c(25, 23, 14, 9, 5, 4, 2, 9, 26, 22, 8, 7, 11, 20, 13, 5, 2,
                   2, 3, 5 )

death <- sim180[["epi"]][["d.h.flow"]]
death <- data.frame(ts = 1:length(death), death = death)
death$cumDeath <- cumsum(death$death)

death <- merge(times, death, by.x = "DayNum", by.y = "ts")
death$death <- death$cumDeath - lag(death$cumDeath)
death$death[1] <- death$cumDeath[1]
death <- death %>% rowid_to_column("ROW")

death$target <- death_targets

f3 <- ggplot(data = death, aes(x = ROW, y = death)) +
  geom_line(aes(color = "Result")) +
  geom_line(aes(y = target, color = "Target (Empirical)")) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(1, 20),
                     breaks = c(1, 4, 7, 10, 13, 16, 19),
                     labels = c("1" = "Jan 21", "4" = "Apr 21", "7" = "Jul 21",
                                "10" = "Oct 21", "13" = "Jan 22", "16" = "Apr 22",
                                "19" = "Jul 22")) +
  scale_y_continuous(expand = c(0, 0), lim = c(0, 55)) +
  xlab("Month") +
  ylab("Incident Deaths") +
  scale_color_manual(values = colors1, name = "") +
  ggtitle("Deaths") + theme_classic()
f3

combined <- f1 + f2 + f3 & theme(legend.position = "bottom")
combined + plot_layout(ncol = 1, guides = "collect")
ggsave("Figure1.tiff", dpi=300, compression = 'lzw', width = 3000, height = 2400, units = "px")

# Vaccine Coverage --------------------------------------------------------

vax_cov <- data.frame(cov.vax1.0to4 = sim180[["epi"]][["cov_vax1_0to4"]],
                      cov.vax1.5to17 = sim180[["epi"]][["cov_vax1_5to17"]],
                      cov.vax1.18to64 = sim180[["epi"]][["cov_vax1_18to64"]],
                      cov.vax1.65p = sim180[["epi"]][["cov_vax1_65p"]],
                      cov.vax2.0to4 = sim180[["epi"]][["cov_vax2_0to4"]],
                      cov.vax2.5to17 = sim180[["epi"]][["cov_vax2_5to17"]],
                      cov.vax2.18to64 = sim180[["epi"]][["cov_vax2_18to64"]],
                      cov.vax2.65p = sim180[["epi"]][["cov_vax2_65p"]],
                      cov.vax3.5to17 = sim180[["epi"]][["cov_vax3_5to17"]],
                      cov.vax3.18to49 = sim180[["epi"]][["cov_vax3_18to49"]],
                      cov.vax3.50to64 = sim180[["epi"]][["cov_vax3_50to64"]],
                      cov.vax3.65p = sim180[["epi"]][["cov_vax3_65p"]],
                      cov.vax4.50to64 = sim180[["epi"]][["cov_vax4_50to64"]],
                      cov.vax4.65p = sim180[["epi"]][["cov_vax4_65p"]])
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

vax_all_runs <- vax_cov

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

vax_all_runs_1 <- vax_all_runs[, c("ROW", "cov.vax1.0to4", "cov.vax2.0to4",
                                   "target.vax1.0to4", "target.vax2.0to4")]

v1 <- ggplot(data = vax_all_runs_1, aes(x = ROW, y = cov.vax1.0to4)) +
  geom_line(aes(y = cov.vax1.0to4, color = "Dose 1: Result")) +
  geom_line(aes(y = cov.vax2.0to4, color = "Dose 2: Result")) +
  geom_line(aes(y = target.vax1.0to4, color = "Dose 1: Target")) +
  geom_line(aes(y = target.vax2.0to4, color = "Dose 2: Target")) +
  scale_x_continuous(expand = c(0,0),
                     limits = c(1, 20),
                     breaks = c(1, 7, 13, 19),
                     labels = c("1" = "Jan 21", "7" = "Jul 21",
                                "13" = "Jan 22", "19" = "Jul 22")) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Month") +
  ylab("Vaccine Coverage (%)") +
  ggtitle("Ages 0 to 4") +
  scale_color_manual(values = colors, guide = "none") +
  theme_classic()
v1


# Age Group 2 -------------------------------------------------------------

vax_all_runs_2 <- vax_all_runs[, c("ROW", "cov.vax1.5to17", "cov.vax2.5to17",
                                   "cov.vax3.5to17", "target.vax1.5to17",
                                   "target.vax2.5to17", "target.vax3.5to17")]

v2 <- ggplot(data = vax_all_runs_2, aes(x = ROW, y = cov.vax1.5to17)) +
  geom_line(aes(y = cov.vax1.5to17, color = "Dose 1: Result")) +
  geom_line(aes(y = cov.vax2.5to17, color = "Dose 2: Result")) +
  geom_line(aes(y = cov.vax3.5to17, color = "Dose 3: Result")) +
  geom_line(aes(y = target.vax1.5to17, color = "Dose 1: Target")) +
  geom_line(aes(y = target.vax2.5to17, color = "Dose 2: Target")) +
  geom_line(aes(y = target.vax3.5to17, color = "Dose 3: Target")) +
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
v2


# Age Group 3 -------------------------------------------------------------

vax_all_runs_3 <- vax_all_runs[, c("ROW", "cov.vax1.18to64", "cov.vax2.18to64",
                                   "cov.vax3.18to49", "target.vax1.18to64",
                                   "target.vax2.18to64", "target.vax3.18to49")]

v3 <- ggplot(data = vax_all_runs_3, aes(x = ROW, y = cov.vax1.18to64)) +
  geom_line(aes(y = cov.vax1.18to64, color = "Dose 1: Result")) +
  geom_line(aes(y = cov.vax2.18to64, color = "Dose 2: Result")) +
  geom_line(aes(y = cov.vax3.18to49, color = "Dose 3: Result")) +
  geom_line(aes(y = target.vax1.18to64, color = "Dose 1: Target")) +
  geom_line(aes(y = target.vax2.18to64, color = "Dose 2: Target")) +
  geom_line(aes(y = target.vax3.18to49, color = "Dose 3: Target")) +
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
v3

# Age Group 4 -------------------------------------------------------------

vax_all_runs_4 <- vax_all_runs[, c("ROW", "cov.vax1.18to64", "cov.vax2.18to64",
                                   "cov.vax3.50to64", "cov.vax4.50to64",
                                   "target.vax1.18to64", "target.vax2.18to64",
                                   "target.vax3.50to64", "target.vax4.50to64")]

v4 <- ggplot(data = vax_all_runs_4, aes(x = ROW, y = cov.vax1.18to64)) +
  geom_line(aes(y = cov.vax1.18to64, color = "Dose 1: Result")) +
  geom_line(aes(y = cov.vax2.18to64, color = "Dose 2: Result")) +
  geom_line(aes(y = cov.vax3.50to64, color = "Dose 3: Result")) +
  geom_line(aes(y = cov.vax4.50to64, color = "Dose 4: Result")) +
  geom_line(aes(y = target.vax1.18to64, color = "Dose 1: Target")) +
  geom_line(aes(y = target.vax2.18to64, color = "Dose 2: Target")) +
  geom_line(aes(y = target.vax3.50to64, color = "Dose 3: Target")) +
  geom_line(aes(y = target.vax4.50to64, color = "Dose 4: Target")) +
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
v4

# Age Group 5 -------------------------------------------------------------

vax_all_runs_5 <- vax_all_runs[, c("ROW", "cov.vax1.65p", "cov.vax2.65p",
                                   "cov.vax3.65p", "cov.vax4.65p",
                                   "target.vax1.65p", "target.vax2.65p",
                                   "target.vax3.65p", "target.vax4.65p")]

v5 <- ggplot(data = vax_all_runs_5, aes(x = ROW, y = cov.vax1.65p)) +
  geom_line(aes(y = cov.vax1.65p, color = "Dose 1: Result")) +
  geom_line(aes(y = cov.vax2.65p, color = "Dose 2: Result")) +
  geom_line(aes(y = cov.vax3.65p, color = "Dose 3: Result")) +
  geom_line(aes(y = cov.vax4.65p, color = "Dose 4: Result")) +
  geom_line(aes(y = target.vax1.65p, color = "Dose 1: Target")) +
  geom_line(aes(y = target.vax2.65p, color = "Dose 2: Target")) +
  geom_line(aes(y = target.vax3.65p, color = "Dose 3: Target")) +
  geom_line(aes(y = target.vax4.65p, color = "Dose 4: Target")) +
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
v5

combined <- v1 + v2 + v3 + v4 + v5
combined + plot_layout(ncol = 3)
ggsave("Figure2.tiff", dpi=300, compression = 'lzw', width = 3000, height = 2400, units = "px")
rm(list=setdiff(ls(), c("est", "vax_targets", "mr_vec", "baseline_sims", "merge_simfiles")))

