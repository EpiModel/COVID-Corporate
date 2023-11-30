library("EpiModelCOVID")
library("EpiModelHPC")
library("dplyr")
library("tidyr")
library("gbp")

# Read in fitted network models
est <- readRDS("est.simple.rds")

# Retrieve attributes
age.grp <- est[[2]]$newnetwork$attr$age.grp
age <- est[[2]]$newnetwork$attr$age
n <- est[[2]]$newnetwork$gal$n

# Set number of households
hh.size <- 2.7 #from https://www.census.gov/quickfacts/GA
n.hh <- round(n / hh.size)
ids <- 1:n
hh.ids <- 1:n.hh

# Create empty data frames to track hh assignment
hh.by.age <- data.frame(hh.ids = hh.ids, mem.u18 = NA, mem.18t64 = NA, mem.65p = NA)
persons.by.hh <- data.frame(ids = ids, age = age, hh = NA)

# Set proportions
prop.hh.with.child <- 0.292 # from Census Table H2 2020
prop.hh.with.adult <- 0.791 # from Census Table H2 2020
prop.hh.with.elderly <- 0.314 # from Census Table H2 2020
prop.children.with.adult <- 0.979 # from Census Table C4 2021

# Household Assignment ----------------------------------------------------

# Determine which households will have a member under 18
hh.u18 <- sample(x = hh.ids, size = round(prop.hh.with.child * n.hh))
hh.by.age$mem.u18 = ifelse(hh.by.age$hh.ids %in% hh.u18, TRUE, FALSE)

# Assign children under 18 to the selected households
persons.by.hh[persons.by.hh$ids %in% which(age<18)[1:length(hh.u18)], ]$hh <- hh.u18
children.hh.assign <- sample(x = hh.u18, size = length(which(age<18)) - length(hh.u18), replace = TRUE)
persons.by.hh[persons.by.hh$ids %in% which(age<18)[(length(hh.u18) + 1):length(which(age<18))], ]$hh <- children.hh.assign

# Determine which hh with children will have at least one 18 - 64 member
num.children.wo.adult <- round((1-prop.children.with.adult) * (sum(age<18)))
num.children <- persons.by.hh[persons.by.hh$age<18, ] %>% group_by(hh) %>% summarize(num = n())
hh.select <- gbp1d_solver_dpp(p = num.children$num, w = num.children$num, c = num.children.wo.adult)
hh.wo.adult <- num.children$hh[as.logical(hh.select$k)]
hh.with.adult <- setdiff(num.children$hh, hh.wo.adult)
hh.by.age$mem.18t64 = ifelse(hh.by.age$hh.ids %in% hh.with.adult, TRUE, NA)

# Determine which hh without children will have at least one 18 - 64 member
num.hh.add.adult <- round(n.hh * prop.hh.with.adult - sum(hh.by.age$mem.18t64 == TRUE, na.rm = TRUE))
hh.add.adult <- sample(x = hh.by.age[is.na(hh.by.age$mem.18t64) & hh.by.age$mem.u18 == FALSE, ]$hh.ids,
                       size = num.hh.add.adult)
hh.by.age[hh.by.age$hh.ids %in% hh.add.adult, ]$mem.18t64 <- TRUE
hh.by.age[is.na(hh.by.age$mem.18t64), ]$mem.18t64 <- FALSE
hh.18t64 <- hh.by.age[hh.by.age$mem.18t64 == TRUE, ]$hh.ids

# Assign adults to the selected households
persons.by.hh[persons.by.hh$ids %in% which(age >= 18 & age < 65)[1:length(hh.18t64)], ]$hh <- hh.18t64
adults.hh.assign <- sample(x = hh.18t64, size = length(which(age >= 18 & age < 65)) - length(hh.18t64), replace = TRUE)
persons.by.hh[persons.by.hh$ids %in% which(age >= 18 & age < 65)[(length(hh.18t64) + 1):length(which(age >= 18 & age < 65))], ]$hh <- adults.hh.assign

# Determine which hh will have a 65+ member
hh.must.elderly <- hh.by.age[hh.by.age$mem.18t64 == FALSE, ]$hh.ids
hh.by.age[hh.by.age$mem.18t64 == FALSE, ]$mem.65p <- TRUE
num.hh.add.elderly <- round(n.hh * prop.hh.with.elderly - length(hh.must.elderly))
hh.add.elderly <- sample(x = hh.by.age[is.na(hh.by.age$mem.65p), ]$hh.ids, size = num.hh.add.elderly)
hh.by.age[hh.by.age$hh.ids %in% hh.add.elderly, ]$mem.65p <- TRUE
hh.by.age[is.na(hh.by.age$mem.65p), ]$mem.65p <- FALSE
hh.65p <- hh.by.age[hh.by.age$mem.65p == TRUE, ]$hh.ids

# Assign elderly to selected households
persons.by.hh[persons.by.hh$ids %in% which(age >= 65)[1:length(hh.65p)], ]$hh <- hh.65p
elderly.hh.assign <- sample(x = hh.65p, size = length(which(age >= 65)) - length(hh.65p), replace = TRUE)
persons.by.hh[persons.by.hh$ids %in% which(age >= 65)[(length(hh.65p) + 1):length(which(age >= 65))], ]$hh <- elderly.hh.assign

# Prepare household network ----------------------------------------------

# Set household attribute
est[[1]]$newnetwork <- set.vertex.attribute(est[[1]]$newnetwork, "household", persons.by.hh$hh)
est[[2]]$newnetwork <- set.vertex.attribute(est[[2]]$newnetwork, "household", persons.by.hh$hh)

# Save household edge list
hhPairs <- merge(persons.by.hh, persons.by.hh, by = "hh")
hhPairs <- subset(hhPairs, (ids.x < ids.y))
hhPairs <- hhPairs[, c(2, 4)]
names(hhPairs) <- c(".head", ".tail")

mortality_rate <- c(
  607.6, 29.6, 12.8, 21.6, 62.8, 116.1, 142.8, 186.5, 228.3, 300.4, 416.1,
  600.4, 945.0, 1453.1, 1952.3, 2817.0, 4368.9, 7158.8, 15626.4
)

mr_pp_pd <- mortality_rate / 1e5 / 365

age_spans <- c(1, 4, rep(5, 16), 1)
mr_vec <- rep(mr_pp_pd, times = age_spans)

# est <- readRDS("est_hh.rds")
# hhPairs <- readRDS("hhPairs.rds")

# Model parameters
pull_env_vars(num.vars = c("inf.p",
                           "boost.start1","boost.start2",
                           "boost.stop1","boost.stop2",
                           "inf.add1","inf.add2",
                           "supp.start1","supp.start2","supp.start3","supp.start4",
                           "supp.stop1","supp.stop2","supp.stop3","supp.stop4",
                           "inf.sub1","inf.sub2","inf.sub3","inf.sub4",
                           "hosp.b.mult","hosp.b.start","hosp.b.stop",
                           "hosp.s.mult","hosp.s.start","hosp.s.stop",
                           "hosp.div","mort.mult"))
param <- param.net(inf.prob = c(inf.p,inf.p,inf.p),
                   inf.add = c(inf.add1, inf.add2), # calibrated
                   inf.sub = c(inf.sub1, inf.sub2, inf.sub3, inf.sub4),
                   inf.boost.start = 180 + c(boost.start1, boost.start2), # calibrated
                   inf.boost.stop = 180 + c(boost.stop1, boost.stop2), # calibrated
                   inf.supp.start =  180 + c(supp.start1, supp.start2, supp.start3, supp.start4), # calibrated
                   inf.supp.stop = 180 + c(supp.stop1, supp.stop2, supp.stop3, supp.stop4), # calibrated

                   act.rate = c(1, 1, 5),
                   inf.prob.inter.rr = c(0.6, 0.6, 0.6),
                   inf.prob.inter.time = c(Inf, Inf, Inf),
                   act.rate.inter.rr = c(1, 1, 1),
                   act.rate.inter.time = c(Inf, Inf, Inf),
                   inf.prob.a.rr = 0.5, # asymptomatic infection
                   inf.prob.mask.rr = 0.3,

                   iso.prob = 0.58, # Fischer et al: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0249891
                   # Adherence to masking in GA September 2020
                   notif.prob = c(0.5,0.1,0.9), # assumed

                   # act.rate.dx.inter.rr = 0.1,
                   # act.rate.dx.inter.time = Inf,
                   # act.rate.sympt.inter.rr = 0.5,
                   # act.rate.sympt.inter.time = Inf,
                   act.rate.iso.inter.time = 1,
                   act.rate.iso.inter.rr = 0.3,

                   prop.clinical = c(0.573, 0.642, 0.760, 0.800, 0.813, 0.814, 0.769, 0.723, 0.666), # by decade of age
                   prop.hospit = c(0.060, 0.063, 0.081, 0.154, 0.207, 0.268, 0.357, 0.465, 0.539) / c(hosp.div), # by decade of age

                   hosp.boost.mult = c(hosp.b.mult),
                   hosp.boost.start = c(hosp.b.start),
                   hosp.boost.stop = c(hosp.b.stop),
                   hosp.supp.mult = c(hosp.s.mult),
                   hosp.supp.start = c(hosp.s.start),
                   hosp.supp.stop = c(hosp.s.stop),

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
                   mort.dis.mult = c(mort.mult),

                   hh.pairs = hhPairs)

# Initial conditions
init <- init.net(e.num = 1600)

# Control settings
control <- control.net(nsteps = 180 + 608,
                       nsims = 32,
                       ncores = 32,
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
fn <- paste0("sim.0824.", fsimno, ".rds")
saveRDS(sim, file = fn)
