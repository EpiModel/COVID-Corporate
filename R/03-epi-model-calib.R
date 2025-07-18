
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
param <- param.net(inf.prob = c(0.38,0.38,0.38),
                   inf.add = c(0.45, 0.10), # calibrated
                   inf.sub = c(0.22, 0.3, 0.23, 0.12),
                   inf.boost.start = 180 + c(-60, 325), # calibrated
                   inf.boost.stop = 180 + c(15, 396), # calibrated
                   inf.supp.start =  180 + c(30, 250, 397, 486), # calibrated
                   inf.supp.stop = 180 + c(165, 324, 485, 608), # calibrated

                   act.rate = c(1, 1, 5),
                   inf.prob.inter.rr = c(0.6, 0.6, 0.6),
                   inf.prob.inter.time = c(Inf, Inf, Inf),
                   act.rate.inter.rr = c(1, 1, 1),
                   act.rate.inter.time = c(Inf, Inf, Inf),
                   inf.prob.a.rr = 0.5, # asymptomatic infection
                   inf.prob.mask.rr = 0.3,

                   iso.prob.office = 0.58, # Fischer et al: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0249891
                                    # Adherence to masking in GA September 2020
                   iso.prob.other = 0.58,
                   notif.prob = c(0.5,0.1,0.9), # assumed

                   # act.rate.dx.inter.rr = 0.1,
                   # act.rate.dx.inter.time = Inf,
                   # act.rate.sympt.inter.rr = 0.5,
                   # act.rate.sympt.inter.time = Inf,
                   act.rate.iso.inter.time = 1,
                   act.rate.iso.inter.rr = 0.3,

                   prop.clinical = c(0.573, 0.642, 0.760, 0.800, 0.813, 0.814, 0.769, 0.723, 0.666), # by decade of age
                   prop.hospit = c(0.060, 0.063, 0.081, 0.154, 0.207, 0.268, 0.357, 0.465, 0.539) / 10, # by decade of age

                   hosp.boost.mult = 5,
                   hosp.boost.start = 180+31,
                   hosp.boost.stop = 180+334,
                   hosp.supp.mult = 1,
                   hosp.supp.start = 180+334,
                   hosp.supp.stop = 180+516,

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
                   vax1.rate = c(0.0005, 0.0009, 0.0019, 0.0019, 0.0039), # manually calibrated; age 0-4, 5-17, 18-49, 50 - 64, 65+
                   vax2.rate = c(0.009, 0.011, 0.010, 0.010, 0.010), # manually calibrated; age 0-4, 5-17, 18-49, 50 - 64, 65+
                   vax3.rate = c(0, 0.0015, 0.0009, 0.0027, 0.0035), # manually calibrated; age 0-4, 5-17, 18-49, 50 - 64, 65+
                   vax4.rate = c(0, 0, 0, 0.002, 0.002), # manually calibrated; age 0-4, 5-17, 18-49, 50 - 64, 65+
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

                   vax1.boost = c(0.02,0.10,0.37,0.37,0.53),
                   vax2.boost = c(0.00,0.60,0.65,0.65,0.65),
                   vax3.boost = c(0.00,0.10,0.18,0.39,0.22),
                   vax4.boost = c(0.00,0.00,0.00,0.06,0.13),
                   vax1.boost.start = 180 + c(30*19,30*6,30*4,30*4,30*2),
                   vax2.boost.start = 180 + c(30*19,30*6,30*4,30*4,30*2) + 1,
                   vax3.boost.start = 180 + c(Inf,30*14,30*12,30*12,30*10),
                   vax4.boost.start = 180 + c(Inf,Inf,Inf,30*16,30*16),

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
control <- control.net(nsteps = 50,
                       #steps = 180 + 608,
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
                       save.network = FALSE,
                       # nwstats.formula.1 = "formation",
                       # nwstats.formula.2 = "formation",
                       # nwstats.formula.3 = "formation",
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

