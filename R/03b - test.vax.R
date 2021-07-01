##Testing vax.FUN ##

library("EpiModelCOVID")

# Read in fitted network models
est <- readRDS("data/input/est.simple.rds")

# Model parameters
source("R/01-epi-params.R")
param <- param.net(inf.prob = c(0.11, 0.11, 0.11),
                   act.rate = c(5, 1, 1),
                   inf.prob.inter.rr = c(0.6, 0.6, 0.6),
                   inf.prob.inter.time = c(Inf, 15, 15),
                   act.rate.inter.rr = c(1, 1, 1),
                   act.rate.inter.time = c(Inf, Inf, Inf),
                   inf.prob.a.rr = 0.5,

                   act.rate.dx.inter.rr = 0.1,
                   act.rate.dx.inter.time = 1,
                   act.rate.sympt.inter.rr = 0.5,
                   act.rate.sympt.inter.time = 1,

                   prop.clinical = c(0.40, 0.25, 0.37, 0.42, 0.51, 0.59, 0.72,
                                     0.76),
                   prop.hospit = c(0, 0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),

                   ea.rate = 1/4.0,
                   ar.rate = 1/5.0,
                   eip.rate = 1/4.0,
                   ipic.rate = 1/1.5,
                   icr.rate = 1/3.5,
                   ich.rate = 1/3.5,
                   hr.rate = 1/7,

                   pcr.sens = 0.8,
                   dx.rate.sympt = 0.1,
                   dx.rate.other = 0.01,
                   allow.rescreen = TRUE,

                   vax.start = Inf,
                   vax1.rate = c(0, 0, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01),
                   vax2.interval = 21,
                   vax1.rr.infect = 0.75,
                   vax2.rr.infect = 0.25,
                   vax.rr.clinical = 0.05,
                   vax1.immune = 7,
                   vax2.immune = 14,

                   a.rate = mean(mr_vec),
                   arrival.age = 0,
                   mort.rates = mr_vec,
                   mort.dis.mult = 180)

# Initial conditions
init <- init.net(e.num = 100)

# Set control settings with nsims = 25
pkgload::load_all("~/git/EpiModelCOVID")
control <- control.net(nsteps = 100,
                       nsims = 5,
                       ncores = 1,
                       initialize.FUN = init_covid_corporate,
                       aging.FUN = aging_covid,
                       departures.FUN = deaths_covid_corporate,
                       arrivals.FUN = arrival_covid_corporate,
                       edges_correct.FUN = NULL,
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
                       resimulate.networks = TRUE,
                       skip.check = TRUE,
                       tergmLite = TRUE)

sim_vaxInf <- netsim(est, param, init, control)
sim_vaxInf <- mutate_epi(sim_vaxInf, se.cuml = cumsum(se.flow),
                         dx.cuml = cumsum(nDx),
                         dx.pos.cuml = cumsum(nDx.pos),
                         totI = e.num + a.num + ip.num + ic.num)
df_vaxInf <- as.data.frame(sim_vaxInf, out = "mean")

## Look at how vax.FUN works when vax.start ≠ Inf -----------------------------

# Model parameters when vax.start = 1
param_vax1 <- param.net(inf.prob = c(0.11, 0.11, 0.11),
                        act.rate = c(5, 1, 1),
                        inf.prob.inter.rr = c(0.6, 0.6, 0.6),
                        inf.prob.inter.time = c(Inf, 15, 15),
                        act.rate.inter.rr = c(1, 1, 1),
                        act.rate.inter.time = c(Inf, Inf, Inf),
                        inf.prob.a.rr = 0.5,

                        act.rate.dx.inter.rr = 0.1,
                        act.rate.dx.inter.time = 1,
                        act.rate.sympt.inter.rr = 0.5,
                        act.rate.sympt.inter.time = 1,

                        prop.clinical = c(0.40, 0.25, 0.37, 0.42, 0.51, 0.59,
                                          0.72, 0.76),
                        prop.hospit = c(0, 0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),

                        ea.rate = 1/4.0,
                        ar.rate = 1/5.0,
                        eip.rate = 1/4.0,
                        ipic.rate = 1/1.5,
                        icr.rate = 1/3.5,
                        ich.rate = 1/3.5,
                        hr.rate = 1/7,

                        pcr.sens = 0.8,
                        dx.rate.sympt = 0.1,
                        dx.rate.other = 0.01,
                        allow.rescreen = TRUE,

                        vax.start = 1,
                        vax1.rate = c(0, 0, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01),
                        vax2.interval = 21,
                        vax1.rr.infect = 0.75,
                        vax2.rr.infect = 0.25,
                        vax.rr.clinical = 0.05,
                        vax1.immune = 7,
                        vax2.immune = 14,

                        a.rate = mean(mr_vec),
                        arrival.age = 0,
                        mort.rates = mr_vec,
                        mort.dis.mult = 180)

sim_vax1 <- netsim(est, param_vax1, init, control)
sim_vax1 <- mutate_epi(sim_vax1,
                       se.cuml = cumsum(se.flow),
                       dx.cuml = cumsum(nDx),
                       dx.pos.cuml = cumsum(nDx.pos),
                       totI = e.num + a.num + ip.num + ic.num)
df_vax1 <- as.data.frame(sim_vax1, out = "mean")

# Model parameters when vax.start = 50
param_vax50 <- param.net(inf.prob = c(0.11, 0.11, 0.11),
                         act.rate = c(5, 1, 1),
                         inf.prob.inter.rr = c(0.6, 0.6, 0.6),
                         inf.prob.inter.time = c(Inf, 15, 15),
                         act.rate.inter.rr = c(1, 1, 1),
                         act.rate.inter.time = c(Inf, Inf, Inf),
                         inf.prob.a.rr = 0.5,

                         act.rate.dx.inter.rr = 0.1,
                         act.rate.dx.inter.time = 1,
                         act.rate.sympt.inter.rr = 0.5,
                         act.rate.sympt.inter.time = 1,

                         prop.clinical = c(0.40, 0.25, 0.37, 0.42, 0.51, 0.59,
                                           0.72, 0.76),
                         prop.hospit = c(0, 0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),

                         ea.rate = 1/4.0,
                         ar.rate = 1/5.0,
                         eip.rate = 1/4.0,
                         ipic.rate = 1/1.5,
                         icr.rate = 1/3.5,
                         ich.rate = 1/3.5,
                         hr.rate = 1/7,

                         pcr.sens = 0.8,
                         dx.rate.sympt = 0.1,
                         dx.rate.other = 0.01,
                         allow.rescreen = TRUE,

                         vax.start = 50,
                         vax1.rate = c(0, 0, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01),
                         vax2.interval = 21,
                         vax1.rr.infect = 0.75,
                         vax2.rr.infect = 0.25,
                         vax.rr.clinical = 0.05,
                         vax1.immune = 7,
                         vax2.immune = 14,

                         a.rate = mean(mr_vec),
                         arrival.age = 0,
                         mort.rates = mr_vec,
                         mort.dis.mult = 180)
sim_vax50 <- netsim(est, param_vax50, init, control)
sim_vax50 <- mutate_epi(sim_vax50,
                        se.cuml = cumsum(se.flow),
                        dx.cuml = cumsum(nDx),
                        dx.pos.cuml = cumsum(nDx.pos),
                        totI = e.num + a.num + ip.num + ic.num)
df_vax50 <- as.data.frame(sim_vax50, out = "mean")

#Checking effect of changing vax.start works as expected
par(mfrow = c(1,1))

plot(sim_vax1, y = "v1.num", mean.col = "blue", qnts.col = "blue")
plot(sim_vaxInf, y = "v1.num", add = TRUE, mean.col = "black",
     qnts.col = "black" )
plot(sim_vax50, y = "v1.num", add = TRUE, mean.col = "red", qnts.col = "red")
legend("topright", legend = c("Inf", "50", "1"),
       col = c("black", "red", "blue"), lty =1, lwd = 2,
       title = "Time Step of Introduction")
abline(v = 1, col = "blue", lty = 2)
abline(v = 50, col = "red", lty = 2)

plot(sim_vax1, y = "v2.num", mean.col = "blue", qnts.col = "blue")
plot(sim_vaxInf, y = "v2.num", add = TRUE, mean.col = "black",
     qnts.col = "black" )
plot(sim_vax50, y = "v2.num", add = TRUE, mean.col = "red", qnts.col = "red")
abline(v = 1 + 21, col = "blue", lty = 2)
abline(v = 50 + 21, col = "red", lty = 2)
legend("topleft", legend = c("Inf", "50", "1"),
       col = c("black", "red", "blue"), lty =1, lwd = 2,
       title = "Time Step of Introduction")

# Hospitalization

plot(sim_vaxInf, y = "h.num", mean.col = "black", qnts.col = "black",
     main = "Hospitalized")
plot(sim_vax1, y = "h.num", add = TRUE, mean.col = "blue", qnts.col = "blue")
plot(sim_vax50, y = "h.num", add = TRUE, mean.col = "red", qnts.col = "red")
abline(v = 1 + 14 + 21, col = "blue", lty = 3)
abline(v = 50 + 21 + 14, col = "red", lty = 3)
abline(v = 1 + 14, col = "blue", lty = 2)
abline(v = 50 + 14, col = "red", lty = 2)
legend("topright", legend = c("Inf", "50", "1"),
       col = c("black", "red", "blue"), lty = 1, lwd = 2,
       title = "Time Step of Introduction")

# Cumulative Incidence

start.date.cum.inc <- data.frame(
  vax.start.inf = sum(df_vaxInf$se.flow, na.rm = TRUE),
  vax.start.1 = sum(df_vax1$se.flow, na.rm = TRUE),
  vax.start.50 = sum(df_vax50$se.flow, na.rm = TRUE))

start.date.cum.inc

# Total infections

plot(sim_vaxInf, y = "totI", mean.col = "black", qnts.col = "black")
plot(sim_vax50, y = "totI", mean.col = "red", qnts.col = "red", add = TRUE)
plot(sim_vax1, y = "totI", mean.col = "blue", qnts.col = "blue", add = TRUE)

# Asymptomatic infections

plot(sim_vaxInf, y = "a.num", mean.col = "black", qnts.col = "black",
     main = "Asymptomatic Infections")
plot(sim_vax50, y = "a.num", mean.col = "red", qnts.col = "red", add = TRUE)
plot(sim_vax1, y = "a.num", mean.col = "blue", qnts.col = "blue", add = TRUE)
legend("topright", legend = c("Inf", "50", "1"),
       col = c("black", "red", "blue"), lty = 1, lwd = 2,
       title = "Time Step of Introduction")

# Infectious Pre-Clinical infections

plot(sim_vaxInf, y = "ip.num", mean.col = "black", qnts.col = "black",
     main = "Infectious Pre-Clinical Infections")
plot(sim_vax50, y = "ip.num", mean.col = "red", qnts.col = "red", add = TRUE)
plot(sim_vax1, y = "ip.num", mean.col = "blue", qnts.col = "blue", add = TRUE)
legend("topright", legend = c("Inf", "50", "1"),
       col = c("black", "red", "blue"), lty = 1, lwd = 2,
       title = "Time Step of Introduction")

# Infectious Clinical infections

plot(sim_vaxInf, y = "ic.num", mean.col = "black", qnts.col = "black",
     main = "Infectious Clinical Infections")
plot(sim_vax50, y = "ic.num", mean.col = "red", qnts.col = "red", add = TRUE)
plot(sim_vax1, y = "ic.num", mean.col = "blue", qnts.col = "blue", add = TRUE)
abline(v = 1 + 14, col = "blue", lty = 2)
abline(v = 50 + 14, col = "red", lty = 2)
abline(v = 1 + 14 + 21, col = "blue", lty = 4)
abline(v = 50 + 14 + 21, col = "red", lty = 4)
legend("topright", legend = c("Inf", "50", "1"),
       col = c("black", "red", "blue"), lty = 1, lwd = 2,
       title = "Time Step of Introduction")

## Changing vaccination uptake rate --------------------------------------------
param_vax1_uptake10 <- param.net(inf.prob = c(0.11, 0.11, 0.11),
                                 act.rate = c(5, 1, 1),
                                 inf.prob.inter.rr = c(0.6, 0.6, 0.6),
                                 inf.prob.inter.time = c(Inf, 15, 15),
                                 act.rate.inter.rr = c(1, 1, 1),
                                 act.rate.inter.time = c(Inf, Inf, Inf),
                                 inf.prob.a.rr = 0.5,

                                 act.rate.dx.inter.rr = 0.1,
                                 act.rate.dx.inter.time = 1,
                                 act.rate.sympt.inter.rr = 0.5,
                                 act.rate.sympt.inter.time = 1,

                                 prop.clinical = c(0.40, 0.25, 0.37, 0.42, 0.51,
                                                   0.59, 0.72, 0.76),
                                 prop.hospit = c(0, 0, 0.1, 0.1, 0.1, 0.1, 0.1,
                                                 0.1),

                                 ea.rate = 1/4.0,
                                 ar.rate = 1/5.0,
                                 eip.rate = 1/4.0,
                                 ipic.rate = 1/1.5,
                                 icr.rate = 1/3.5,
                                 ich.rate = 1/3.5,
                                 hr.rate = 1/7,

                                 pcr.sens = 0.8,
                                 dx.rate.sympt = 0.1,
                                 dx.rate.other = 0.01,
                                 allow.rescreen = TRUE,

                                 vax.start = 1,
                                 vax1.rate = c(0, 0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
                                 vax2.interval = 21,
                                 vax1.rr.infect = 0.75,
                                 vax2.rr.infect = 0.25,
                                 vax.rr.clinical = 0.05,
                                 vax1.immune = 7,
                                 vax2.immune = 14,

                                 a.rate = mean(mr_vec),
                                 arrival.age = 0,
                                 mort.rates = mr_vec,
                                 mort.dis.mult = 180)
sim_vax1_uptake10 <- netsim(est, param_vax1_uptake10, init, control)
sim_vax1_uptake10 <- mutate_epi(sim_vax1_uptake10,
                                se.cuml = cumsum(se.flow),
                                dx.cuml = cumsum(nDx),
                                dx.pos.cuml = cumsum(nDx.pos),
                                totI = e.num + a.num + ip.num + ic.num)
df_vax1_uptake10 <- as.data.frame(sim_vax1_uptake10, out = "mean")

param_vax1_uptake100 <- param.net(inf.prob = c(0.11, 0.11, 0.11),
                                 act.rate = c(5, 1, 1),
                                 inf.prob.inter.rr = c(0.6, 0.6, 0.6),
                                 inf.prob.inter.time = c(Inf, 15, 15),
                                 act.rate.inter.rr = c(1, 1, 1),
                                 act.rate.inter.time = c(Inf, Inf, Inf),
                                 inf.prob.a.rr = 0.5,

                                 act.rate.dx.inter.rr = 0.1,
                                 act.rate.dx.inter.time = 1,
                                 act.rate.sympt.inter.rr = 0.5,
                                 act.rate.sympt.inter.time = 1,

                                 prop.clinical = c(0.40, 0.25, 0.37, 0.42, 0.51,
                                                   0.59, 0.72, 0.76),
                                 prop.hospit = c(0, 0, 0.1, 0.1, 0.1, 0.1, 0.1,
                                                 0.1),

                                 ea.rate = 1/4.0,
                                 ar.rate = 1/5.0,
                                 eip.rate = 1/4.0,
                                 ipic.rate = 1/1.5,
                                 icr.rate = 1/3.5,
                                 ich.rate = 1/3.5,
                                 hr.rate = 1/7,

                                 pcr.sens = 0.8,
                                 dx.rate.sympt = 0.1,
                                 dx.rate.other = 0.01,
                                 allow.rescreen = TRUE,

                                 vax.start = 1,
                                 vax1.rate = c(0, 0, 1, 1, 1, 1, 1, 1),
                                 vax2.interval = 21,
                                 vax1.rr.infect = 0.75,
                                 vax2.rr.infect = 0.25,
                                 vax.rr.clinical = 0.05,
                                 vax1.immune = 7,
                                 vax2.immune = 14,

                                 a.rate = mean(mr_vec),
                                 arrival.age = 0,
                                 mort.rates = mr_vec,
                                 mort.dis.mult = 180)
sim_vax1_uptake100 <- netsim(est, param_vax1_uptake100, init, control)
sim_vax1_uptake100 <- mutate_epi(sim_vax1_uptake100,
                                se.cuml = cumsum(se.flow),
                                dx.cuml = cumsum(nDx),
                                dx.pos.cuml = cumsum(nDx.pos),
                                totI = e.num + a.num + ip.num + ic.num)
df_vax1_uptake100 <- as.data.frame(sim_vax1_uptake100, out = "mean")

param_vax1_uptake30 <- param.net(inf.prob = c(0.11, 0.11, 0.11),
                                  act.rate = c(5, 1, 1),
                                  inf.prob.inter.rr = c(0.6, 0.6, 0.6),
                                  inf.prob.inter.time = c(Inf, 15, 15),
                                  act.rate.inter.rr = c(1, 1, 1),
                                  act.rate.inter.time = c(Inf, Inf, Inf),
                                  inf.prob.a.rr = 0.5,

                                  act.rate.dx.inter.rr = 0.1,
                                  act.rate.dx.inter.time = 1,
                                  act.rate.sympt.inter.rr = 0.5,
                                  act.rate.sympt.inter.time = 1,

                                  prop.clinical = c(0.40, 0.25, 0.37, 0.42,
                                                    0.51, 0.59, 0.72, 0.76),
                                  prop.hospit = c(0, 0, 0.1, 0.1, 0.1, 0.1, 0.1,
                                                  0.1),

                                  ea.rate = 1/4.0,
                                  ar.rate = 1/5.0,
                                  eip.rate = 1/4.0,
                                  ipic.rate = 1/1.5,
                                  icr.rate = 1/3.5,
                                  ich.rate = 1/3.5,
                                  hr.rate = 1/7,

                                  pcr.sens = 0.8,
                                  dx.rate.sympt = 0.1,
                                  dx.rate.other = 0.01,
                                  allow.rescreen = TRUE,

                                  vax.start = 1,
                                  vax1.rate = c(0, 0, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3),
                                  vax2.interval = 21,
                                  vax1.rr.infect = 0.75,
                                  vax2.rr.infect = 0.25,
                                  vax.rr.clinical = 0.05,
                                  vax1.immune = 7,
                                  vax2.immune = 14,

                                  a.rate = mean(mr_vec),
                                  arrival.age = 0,
                                  mort.rates = mr_vec,
                                  mort.dis.mult = 180)

sim_vax1_uptake30 <- netsim(est, param_vax1_uptake30, init, control)
sim_vax1_uptake30 <- mutate_epi(sim_vax1_uptake30,
                                 se.cuml = cumsum(se.flow),
                                 dx.cuml = cumsum(nDx),
                                 dx.pos.cuml = cumsum(nDx.pos),
                                 totI = e.num + a.num + ip.num + ic.num)
df_vax1_uptake30 <- as.data.frame(sim_vax1_uptake30, out = "mean")

# Visualizations

# Exposed

plot(sim_vax1, y = "e.num", main = "Exposed", mean.col = "blue",
     qnts.col = "blue")
plot(sim_vax1_uptake100, y = "e.num", mean.col = "red", qnts.col = "red",
     add = TRUE)
plot(sim_vax1_uptake30, y = "e.num", mean.col = "yellow", qnts.col = "yellow",
     add = TRUE)
plot(sim_vax1_uptake10, y = "e.num", mean.col = "darkgreen",
     qnts.col = "darkgreen", add = TRUE)
#abline(v = 1 + 14, col = "gray", lty = 2)
legend("topright", legend = c("0.01", "0.3","0.1", "1"), col = c("blue",
                                                                 "yellow",
                                                                 "darkgreen",
                                                                 "red"),
       lwd = 2, title = "vax1.rate")

# Total Infected

plot(sim_vax1, y = "totI", main = "Total Infected", mean.col = "blue",
     qnts.col = "blue")
plot(sim_vax1_uptake100, y = "totI", mean.col = "red", qnts.col = "red",
     add = TRUE)
plot(sim_vax1_uptake30, y = "totI", mean.col = "yellow", qnts.col = "yellow",
     add = TRUE)
plot(sim_vax1_uptake10, y = "totI", mean.col = "darkgreen",
     qnts.col = "darkgreen", add = TRUE)
#abline(v = 1 + 14, col = "gray", lty = 2)
legend("topright", legend = c("0.01", "0.3","0.1", "1"), col = c("blue",
                                                                 "yellow",
                                                                 "darkgreen",
                                                                 "red"),
       lwd = 2, title = "vax1.rate")

# Asymptomatic

plot(sim_vax1, y = "a.num", main = "Asymptomatic", mean.col = "blue",
     qnts.col = "blue")
plot(sim_vax1_uptake100, y = "a.num", mean.col = "red", qnts.col = "red",
     add = TRUE)
plot(sim_vax1_uptake30, y = "a.num", mean.col = "yellow", qnts.col = "yellow",
     add = TRUE)
plot(sim_vax1_uptake10, y = "a.num", mean.col = "darkgreen",
     qnts.col = "darkgreen", add = TRUE)
#abline(v = 1 + 14, col = "gray", lty = 2)
legend("topright", legend = c("0.01", "0.3","0.1", "1"), col = c("blue",
                                                                 "yellow",
                                                                 "darkgreen",
                                                                 "red"),
       lwd = 2, title = "vax1.rate")

# Infectious Clinical

plot(sim_vax1, y = "ic.num", main = "Infectious Clinical", mean.col = "blue",
     qnts.col = "blue")
plot(sim_vax1_uptake100, y = "ic.num", mean.col = "red", qnts.col = "red",
     add = TRUE)
plot(sim_vax1_uptake30, y = "ic.num", mean.col = "yellow", qnts.col = "yellow",
     add = TRUE)
plot(sim_vax1_uptake10, y = "ic.num", mean.col = "darkgreen",
     qnts.col = "darkgreen", add = TRUE)
#abline(v = 1 + 14, col = "gray", lty = 2)
legend("topright", legend = c("0.01", "0.3","0.1", "1"), col = c("blue",
                                                                 "yellow",
                                                                 "darkgreen",
                                                                 "red"),
       lwd = 2, title = "vax1.rate")

# Infectious Pre-Clinical

plot(sim_vax1, y = "ip.num", main = "Infectious Pre-Clinical",
     mean.col = "blue", qnts.col = "blue")
plot(sim_vax1_uptake100, y = "ip.num", mean.col = "red", qnts.col = "red",
     add = TRUE)
plot(sim_vax1_uptake30, y = "ip.num", mean.col = "yellow", qnts.col = "yellow",
     add = TRUE)
plot(sim_vax1_uptake10, y = "ip.num", mean.col = "darkgreen",
     qnts.col = "darkgreen", add = TRUE)
#abline(v = 1 + 14, col = "gray", lty = 2)
legend("topright", legend = c("0.01", "0.3","0.1", "1"), col = c("blue",
                                                                 "yellow",
                                                                 "darkgreen",
                                                                 "red"),
       lwd = 2, title = "vax1.rate")

# Hospitalized

plot(sim_vax1_uptake100, y = "h.num", mean.col = "red", qnts.col = "red")
plot(sim_vax1, y = "h.num", main = "Hospitalized", mean.col = "blue",
     qnts.col = "blue", add = TRUE)
plot(sim_vax1_uptake30, y = "h.num", mean.col = "yellow", qnts.col = "yellow",
     add = TRUE)
plot(sim_vax1_uptake10, y = "h.num", mean.col = "darkgreen",
     qnts.col = "darkgreen", add = TRUE)
#abline(v = 1 + 14, col = "gray", lty = 2)
legend("topright", legend = c("0.01", "0.3","0.1", "1"), col = c("blue",
                                                                 "yellow",
                                                                 "darkgreen",
                                                                 "red"),
       lwd = 2, title = "vax1.rate")

# Recovered

plot(sim_vax1, y = "r.num", main = "Recovered", mean.col = "blue",
     qnts.col = "blue")
plot(sim_vax1_uptake100, y = "r.num", mean.col = "red", qnts.col = "red",
     add = TRUE)
plot(sim_vax1_uptake30, y = "r.num", mean.col = "yellow", qnts.col = "yellow",
     add = TRUE)
plot(sim_vax1_uptake10, y = "r.num", mean.col = "darkgreen",
     qnts.col = "darkgreen", add = TRUE)
#abline(v = 1 + 14, col = "gray", lty = 2)
legend("topright", legend = c("0.01", "0.3","0.1", "1"), col = c("blue",
                                                                 "yellow",
                                                                 "darkgreen",
                                                                 "red"),
       lwd = 2, title = "vax1.rate")

# diff_vax1_uptake30 <- data.frame(diff_exposed = max(df_vax1$e.num) -
#                                  max(df_vax1_uptake30$e.num),
#                                diff_infected = max(df_vax1$totI) -
#                                  max(df_vax1_uptake30$totI),
#                                diff_asymp = max(df_vax1$a.num) -
#                                  max(df_vax1_uptake30$a.num),
#                                diff_ic = max(df_vax1$ic.num) -
#                                  max(df_vax1_uptake30$ic.num),
#                                diff_ip = max(df_vax1$ip.num) -
#                                  max(df_vax1_uptake30$ip.num),
#                                diff_hospitalized = max(df_vax1$h.num) -
#                                  max(df_vax1_uptake30$h.num),
#                                diff_recovered = max(df_vax1$r.num) -
#                                  max(df_vax1_uptake30$r.num))
# diff_vax1_uptake30

#Changing vaccine effectiveness ------------------------------------------------

param_vax1_uptake100_rr0 <- param.net(inf.prob = c(0.11, 0.11, 0.11),
                        act.rate = c(5, 1, 1),
                        inf.prob.inter.rr = c(0.6, 0.6, 0.6),
                        inf.prob.inter.time = c(Inf, 15, 15),
                        act.rate.inter.rr = c(1, 1, 1),
                        act.rate.inter.time = c(Inf, Inf, Inf),
                        inf.prob.a.rr = 0.5,

                        act.rate.dx.inter.rr = 0.1,
                        act.rate.dx.inter.time = 1,
                        act.rate.sympt.inter.rr = 0.5,
                        act.rate.sympt.inter.time = 1,

                        prop.clinical = c(0.40, 0.25, 0.37, 0.42, 0.51, 0.59,
                                          0.72, 0.76),
                        prop.hospit = c(0, 0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),

                        ea.rate = 1/4.0,
                        ar.rate = 1/5.0,
                        eip.rate = 1/4.0,
                        ipic.rate = 1/1.5,
                        icr.rate = 1/3.5,
                        ich.rate = 1/3.5,
                        hr.rate = 1/7,

                        pcr.sens = 0.8,
                        dx.rate.sympt = 0.1,
                        dx.rate.other = 0.01,
                        allow.rescreen = TRUE,

                        vax.start = 1,
                        vax1.rate = 1,
                        vax2.interval = 21,
                        vax1.rr.infect = 0,
                        vax2.rr.infect = 0,
                        vax.rr.clinical = 0,
                        vax1.immune = 7,
                        vax2.immune = 14,

                        a.rate = mean(mr_vec),
                        arrival.age = 0,
                        mort.rates = mr_vec,
                        mort.dis.mult = 180)
sim_vax1_uptake100_rr0 <- netsim(est, param_vax1_uptake100_rr0, init, control)
sim_vax1_uptake100_rr0 <- mutate_epi(sim_vax1_uptake100_rr0,
                                 se.cuml = cumsum(se.flow),
                                 dx.cuml = cumsum(nDx),
                                 dx.pos.cuml = cumsum(nDx.pos),
                                 totI = e.num + a.num + ip.num + ic.num)
df_vax1_uptake100_rr0 <- as.data.frame(sim_vax1_uptake100_rr0, out = "mean")

param_vax1_uptake100_rr1 <- param.net(inf.prob = c(0.11, 0.11, 0.11),
                                      act.rate = c(5, 1, 1),
                                      inf.prob.inter.rr = c(0.6, 0.6, 0.6),
                                      inf.prob.inter.time = c(Inf, 15, 15),
                                      act.rate.inter.rr = c(1, 1, 1),
                                      act.rate.inter.time = c(Inf, Inf, Inf),
                                      inf.prob.a.rr = 0.5,

                                      act.rate.dx.inter.rr = 0.1,
                                      act.rate.dx.inter.time = 1,
                                      act.rate.sympt.inter.rr = 0.5,
                                      act.rate.sympt.inter.time = 1,

                                      prop.clinical = c(0.40, 0.25, 0.37, 0.42,
                                                        0.51, 0.59, 0.72, 0.76),
                                      prop.hospit = c(0, 0, 0.1, 0.1, 0.1, 0.1,
                                                      0.1, 0.1),

                                      ea.rate = 1/4.0,
                                      ar.rate = 1/5.0,
                                      eip.rate = 1/4.0,
                                      ipic.rate = 1/1.5,
                                      icr.rate = 1/3.5,
                                      ich.rate = 1/3.5,
                                      hr.rate = 1/7,

                                      pcr.sens = 0.8,
                                      dx.rate.sympt = 0.1,
                                      dx.rate.other = 0.01,
                                      allow.rescreen = TRUE,

                                      vax.start = 1,
                                      vax1.rate = c(0, 0, 1, 1, 1, 1, 1, 1),
                                      vax2.interval = 21,
                                      vax1.rr.infect = 1,
                                      vax2.rr.infect = 1,
                                      vax.rr.clinical = 1,
                                      vax1.immune = 7,
                                      vax2.immune = 14,

                                      a.rate = mean(mr_vec),
                                      arrival.age = 0,
                                      mort.rates = mr_vec,
                                      mort.dis.mult = 180)
sim_vax1_uptake100_rr1 <- netsim(est, param_vax1_uptake100_rr1, init, control)
sim_vax1_uptake100_rr1 <- mutate_epi(sim_vax1_uptake100_rr1,
                                     se.cuml = cumsum(se.flow),
                                     dx.cuml = cumsum(nDx),
                                     dx.pos.cuml = cumsum(nDx.pos),
                                     totI = e.num + a.num + ip.num + ic.num)
df_vax1_uptake100_rr1 <- as.data.frame(sim_vax1_uptake100_rr1, out = "mean")

# Exposed

plot(sim_vax1_uptake100, y = "e.num", main = "Exposed", mean.col = "blue",
     qnts.col = "blue")
plot(sim_vax1_uptake100_rr0, y = "e.num", mean.col = "red", qnts.col = "red",
     add = TRUE)
plot(sim_vax1_uptake100_rr1, y = "e.num", mean.col = "darkgreen",
     qnts.col = "darkgreen", add = TRUE)
abline(v = 1 + 14, col = "gray", lty = 2)
legend("topright", legend = c("rr.infect and rr.clinical = 0",
                              "rr.infect and rr.clinical = 1"),
       col = c("red", "darkgreen"),
       lwd = 2)

# Total Infected

plot(sim_vax1_uptake100, y = "totI", main = "Total Infected", mean.col = "blue",
     qnts.col = "blue")
plot(sim_vax1_uptake100_rr0, y = "totI", mean.col = "red", qnts.col = "red",
     add = TRUE)
plot(sim_vax1_uptake100_rr1, y = "totI", mean.col = "darkgreen",
     qnts.col = "darkgreen", add = TRUE)
abline(v = 1 + 14, col = "gray", lty = 2)
legend("topright", legend = c("rr.infect and rr.clinical = 0",
                              "rr.infect and rr.clinical = 1"),
       col = c("red", "darkgreen"),
       lwd = 2)

# Asymptomatic

plot(sim_vax1_uptake100, y = "a.num", main = "Asymptomatic", mean.col = "blue",
     qnts.col = "blue")
plot(sim_vax1_uptake100_rr0, y = "a.num", mean.col = "red", qnts.col = "red",
     add = TRUE)
plot(sim_vax1_uptake100_rr1, y = "a.num", mean.col = "darkgreen",
     qnts.col = "darkgreen", add = TRUE)
abline(v = 1 + 14, col = "gray", lty = 2)
legend("topright", legend = c("rr.infect and rr.clinical = 0",
                              "rr.infect and rr.clinical = 1"),
       col = c("red", "darkgreen"),
       lwd = 2)

# Infectious Clinical

plot(sim_vax1_uptake100, y = "ic.num", main = "Infectious Clinical",
     mean.col = "blue", qnts.col = "blue")
plot(sim_vax1_uptake100_rr0, y = "ic.num", mean.col = "red", qnts.col = "red",
     add = TRUE)
plot(sim_vax1_uptake100_rr1, y = "ic.num", mean.col = "darkgreen",
     qnts.col = "darkgreen", add = TRUE)
abline(v = 1 + 14, col = "gray", lty = 2)
legend("topright", legend = c("rr.infect and rr.clinical = 0",
                              "rr.infect and rr.clinical = 1"),
       col = c("red", "darkgreen"),
       lwd = 2)

# Infectious Pre-Clinical

plot(sim_vax1_uptake100, y = "ip.num", main = "Infectious Pre-Clinical",
     mean.col = "blue", qnts.col = "blue")
plot(sim_vax1_uptake100_rr0, y = "ip.num", mean.col = "red", qnts.col = "red",
     add = TRUE)
plot(sim_vax1_uptake100_rr1, y = "ip.num", mean.col = "darkgreen",
     qnts.col = "darkgreen", add = TRUE)
abline(v = 1 + 14, col = "gray", lty = 2)
legend("topright", legend = c("rr.infect and rr.clinical = 0",
                              "rr.infect and rr.clinical = 1"),
       col = c("red", "darkgreen"),
       lwd = 2)

# Hospitalized

plot(sim_vax1_uptake100, y = "h.num", main = "Hospitalized", mean.col = "blue",
     qnts.col = "blue")
plot(sim_vax1_uptake100_rr0, y = "h.num", mean.col = "red", qnts.col = "red",
     add = TRUE)
plot(sim_vax1_uptake100_rr1, y = "h.num", mean.col = "darkgreen",
     qnts.col = "darkgreen", add = TRUE)
abline(v = 1 + 14, col = "gray", lty = 2)
legend("topright", legend = c("rr.infect and rr.clinical = 0",
                              "rr.infect and rr.clinical = 1"),
       col = c("red", "darkgreen"),
       lwd = 2)

# Recovered

plot(sim_vax1_uptake100, y = "r.num", main = "Recovered", mean.col = "blue",
     qnts.col = "blue")
plot(sim_vax1_uptake100_rr0, y = "r.num", mean.col = "red", qnts.col = "red",
     add = TRUE)
plot(sim_vax1_uptake100_rr1, y = "r.num", mean.col = "darkgreen",
     qnts.col = "darkgreen", add = TRUE)
abline(v = 1 + 14, col = "gray", lty = 2)
legend("bottomright", legend = c("rr.infect and rr.clinical = 0",
                              "rr.infect and rr.clinical = 1"),
       col = c("red", "darkgreen"),
       lwd = 2)

sum(df_vax1_uptake100_rr0$se.flow, na.rm = TRUE)
sum(df_vax1_uptake100$se.flow, na.rm = TRUE)
sum(df_vax1_uptake100_rr1$se.flow, na.rm = TRUE)

#calculate percent of infections averted with perfect vaccine
((sum(df_vaxInf$se.flow, na.rm = TRUE) -
    sum(df_vax1_uptake100_rr0$se.flow, na.rm = TRUE)) /
    sum(df_vaxInf$se.flow, na.rm = TRUE))

#calculate percent of infections averted with realistic vaccine
((sum(df_vaxInf$se.flow, na.rm = TRUE) -
    sum(df_vax1_uptake100$se.flow, na.rm = TRUE)) /
    sum(df_vaxInf$se.flow, na.rm = TRUE))

