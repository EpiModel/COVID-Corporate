
##
## COVID-19 Corporate Network Model
## Network Parameters and Model
##
## Authors: Samuel M. Jenness
## Date: February 2021
##

library("EpiModelCOVID")
library("tidyverse")
library(socialmixr)

# load CorpMix round 5 contact data----------------------------------------

# 'raw' round 4 mixing data
cm <- readRDS("data/input/rd5_mean_degree.rds")

# 'balanced' round 5 mixing matrix
cm2 <- readRDS("data/input/rd5_matrix.rds") %>%
  ungroup()

# round 5 mixing matrix with proportions in each cell
cm.prop <- readRDS("data/input/rd5_matrix_prop.rds") %>%
  ungroup()

# total mean degree by location using raw data
md <- cm %>%
  group_by(location,p_age,n_participants) %>%
  summarise(n_contacts = sum(n_contacts)) %>%
  ungroup() %>%
  group_by(location) %>%
  summarise(n_contacts = sum(n_contacts),
            n_participants = sum(n_participants)) %>%
  mutate(md = (n_contacts/n_participants)/2)

# md.age <- cm %>%
#   group_by(location,p_age,n_participants) %>%
#   mutate(n_contacts = sum(n_contacts)) %>%
#   select(location,p_age,n_participants,n_contacts,total_degree) %>%
#   distinct() %>%
#   ungroup()

# proportion within age group mixing by location
md.match.raw <- cm %>%
  filter(as.character(p_age) == as.character(contact_age)) %>%
  group_by(location) %>%
  summarise(n_contacts = sum(n_contacts),
            n_participants = sum(n_participants)) %>%
  mutate(md.match = (n_contacts/n_participants)/2) %>%
  left_join(md,by = "location") %>%
  mutate(match = md.match/md)

md.match <- cm2 %>%
  mutate(within_age_edges = total_degree*prop_match) %>%
  select(location,age,total_degree,within_age_edges) %>%
  group_by(location) %>%
  summarise(total_degree = sum(total_degree),
            within_age_edges = sum(within_age_edges)) %>%
  mutate(match = within_age_edges/total_degree)

# load POLYMOD data -------------------------------------------------------

data("polymod")
contact_matrix(polymod, age.limits = c(0,10,20,30,40,50,60,70,80,90,100),
               filter = list(cnt_home=1, duration_multi=c(3,4,5)))

# home contactas > 15 mins
poly.home <- contact_matrix(polymod, age.limits = c(0,10,20,30,40,50,60,70,80,90,100),
               filter = list(cnt_home=1, duration_multi=c(3,4,5)))$matrix
rowSums(poly.home)

# home contacts all durations
poly.home.all <- contact_matrix(polymod, age.limits = c(0,10,20,30,40,50,60,70,80,90,100),
                            filter = list(cnt_home=1))$matrix
rowSums(poly.home.all)

# contacts outside the home > 15 mins
poly.nothome <- contact_matrix(polymod, age.limits = c(0,10,20,30,40,50,60,70,80,90,100),
                            filter = list(cnt_home=0, duration_multi=c(3,4,5)))$matrix
rowSums(poly.nothome)

# contacts outside the home all durations
poly.nothome.all <- contact_matrix(polymod, age.limits = c(0,10,20,30,40,50,60,70,80,90,100),
                               filter = list(cnt_home=0))$matrix
rowSums(poly.nothome.all)

# Simple 3-layer Network Model --------------------------------------------

## For use in testing only ##

# Set network
n <- 10000
nw <- network_initialize(n)

age <- round(rnorm(n, 40, 21), 1)
age <- pmax(age, 0)
age <- pmin(age, 99)
hist(age)

age.breaks <- seq(0, 100, 10)
age.grp <- cut(age, age.breaks, labels = FALSE, right = FALSE)

nw <- set_vertex_attribute(nw, "age", age)
nw <- set_vertex_attribute(nw, "age.grp", age.grp)

age.grp.num <- as.numeric(table(age.grp))

## Within household network
nw.loc <- "Home"

# Mean degree from corp mix data
md.hh <- md %>%
  filter(location == nw.loc) %>%
  select(md) %>%
  as.numeric()

# Weighted mean degree
md.hh.age <- cm2 %>%
  filter(location == nw.loc & age != "70s") %>%
  pull(total_degree)
# add polymod data for missing age groups
md.hh.age <- c(1.36,1.14,md.hh.age,0.63,0.68,0.65)
# multiply by weights
w <- age.grp.num/sum(age.grp.num)
md.hh2 <- sum(md.hh.age * w)

edges.hh <- md.hh * n/2

# Number of edges per age group for nodefactor()
# note this only includes age groups 20-60+
md.hh.ag <- cm2 %>%
  filter(location == nw.loc) %>%
  select(total_degree) %>%
  as_vector() %>%
  as.numeric()
# adding missing md for other age groups for now
md.hh.ag <- c(NA,NA,md.hh.ag,NA,NA)

nf.hh.ag <- md.hh.ag * age.grp.num

# Proportion within-group mixing for nodematch()
match.hh.ag <- md.match$match[md.match$location == nw.loc]

nm.hh.ag <- edges.hh * match.hh.ag

formation.hh <- ~edges + nodefactor("age.grp", levels = 3:7) + nodematch("age.grp", diff = F)
# formation.hh <- ~edges + nodefactor("age.grp", levels = 3:7)
# test keeping/dropping nodefactor levels

target.stats.hh <- c(edges.hh, nf.hh.ag[3:7], nm.hh.ag)
# target.stats.hh <- c(edges.hh, nf.hh.ag[3:7])

coef.diss.hh <- dissolution_coefs(dissolution = ~offset(edges), duration = 1e6)

est.hh <- netest(nw, formation.hh, target.stats.hh, coef.diss.hh,
                 set.control.ergm = control.ergm(MCMLE.maxit = 500))
summary(est.hh)

# est <- netest(nw, formation, target.stats, coef.diss, edapprox = TRUE,
#               set.control.ergm = control.ergm(MCMC.burnin=1e7, MCMC.interval=1e7,
#                                               MCMC.samplesize = 20000,
#                                               init.MPLE.samplesize = 1e8,
#                                               SAN.control = control.san(SAN.nsteps = 1e8)))

dx.hh <- netdx(est.hh, nsims = 1, ncores = 1, nsteps = 1000, dynamic = TRUE,
               set.control.ergm = control.simulate.formula(MCMC.burnin = 1e6),
               nwstats.formula = ~edges + nodefactor("age.grp",levels = NULL) + nodematch("age.grp",diff = T) +
                 nodemix("age.grp",levels2 = NULL), keep.tedgelist = TRUE)
print(dx.hh)

#list of edges over time steps
#tail and head are positional ids, could look at ages for these
#can calculate my own statistics from this
as.data.frame(dx.hh)

# pull positional ids and corresponding age groups
vertex.names <- get_vertex_attribute(nw, "vertex.names")
age.grp <- get_vertex_attribute(nw, "age.grp")
id.age <- as.data.frame(cbind(vertex.names,age.grp))
num.age <- id.age %>%
  group_by(age.grp) %>%
  tally(name = "num.participants")

# generate network summary statistics
dx.hh.sum <- as.data.frame(dx.hh) %>%
  select(tail,head) %>%
  left_join(id.age, by = c("head" = "vertex.names")) %>%
  rename(head.age = age.grp) %>%
  left_join(id.age, by = c("tail" = "vertex.names")) %>%
  rename(tail.age = age.grp) %>%
  left_join(num.age, by = c("head.age" = "age.grp")) %>%
  group_by(head.age,tail.age,num.participants) %>%
  tally(name = "num.contacts") %>%
  ungroup() %>%
  group_by(head.age,num.participants) %>%
  mutate(total.contacts = sum(num.contacts)) %>%
  mutate(age.degree = num.contacts/num.participants,
         total.degree = total.contacts/num.participants)

# second version with one row for each head AND tail contact
dx.hh.sum <- as.data.frame(dx.hh) %>%
  rename(head = "tail", tail = "head") %>%
  add_row(as.data.frame(dx.hh)) %>%
  select(tail,head) %>%
  left_join(id.age, by = c("head" = "vertex.names")) %>%
  rename(head.age = age.grp) %>%
  left_join(id.age, by = c("tail" = "vertex.names")) %>%
  rename(tail.age = age.grp) %>%
  left_join(num.age, by = c("head.age" = "age.grp")) %>%
  group_by(head.age,tail.age,num.participants) %>%
  tally(name = "num.contacts") %>%
  ungroup() %>%
  group_by(head.age,num.participants) %>%
  mutate(total.contacts = sum(num.contacts)) %>%
  mutate(age.degree = num.contacts/num.participants,
         total.degree = total.contacts/num.participants)

dx.hh.match <- dx.hh.sum %>%
  filter(head.age == tail.age) %>%
  mutate(prop.match =  age.degree/total.degree)

dx.hh.mix <- dx.hh.sum %>%
  ungroup() %>%
  select(head.age,tail.age,age.degree) %>%
  pivot_wider(values_from = "age.degree",names_from = "tail.age")

# dx <- netdx(est, nsims = 1, nsteps = 1000, keep.tedgelist = TRUE,
#             set.control.ergm = control.simulate.stergm(MCMC.init.maxchanges = 1e8,
#                                                        MCMC.burnin.min = 1.5e5,
#                                                        MCMC.burnin.max = 1.5e5))


## Within office network
md.oo <- md %>%
  filter(location == "Work") %>%
  select(md) %>%
  as.numeric()
target.stats.oo <- md.oo * n/2

formation.oo <- ~edges

coef.diss.oo <- dissolution_coefs(dissolution = ~offset(edges), duration = 1e5)

est.oo <- netest(nw, formation.oo, target.stats.oo, coef.diss.oo,
                 set.control.ergm = control.ergm(MCMLE.maxit = 500))
summary(est.oo)

dx.oo <- netdx(est.oo, nsims = 25, ncores = 5, nsteps = 1000, dynamic = TRUE,
               set.control.ergm = control.simulate.formula(MCMC.burnin = 1e6))

print(dx.oo)


## Community network
md.cc <- md %>%
  filter(location == "Community") %>%
  select(md) %>%
  as.numeric()
target.stats.cc <- md.cc * n/2

formation.cc <- ~edges

coef.diss.cc <- dissolution_coefs(dissolution = ~offset(edges), duration = 1)

est.cc <- netest(nw, formation.cc, target.stats.cc, coef.diss.cc,
                 set.control.ergm = control.ergm(MCMLE.maxit = 500))
summary(est.cc)

dx.cc <- netdx(est.cc, nsims = 1000, dynamic = FALSE,
               set.control.ergm = control.simulate.formula(MCMC.burnin = 1e6))
print(dx.cc)

est.hh <- trim_netest(est.hh)
est.oo <- trim_netest(est.oo)
est.cc <- trim_netest(est.cc)

est.simple <- list(est.hh, est.oo, est.cc)
saveRDS(est.simple, file = "data/input/est.simple.rds")



# # Varying Parameters (Table 3) --------------------------------------------
#
# n.sectors <- 10
# prop.within.sector.post <- 0.98
# n.sectors.pass <- NULL  # If NULL, only within-cabin pass/pass mixing
#
#
# # Network Setup -----------------------------------------------------------
#
# n.crew <- 1045
# n.pass <- 2666
# n <- n.crew + n.pass
#
# n.rooms <- 1337
# n.pass.per.room <- n.pass/n.rooms
# n.pass.per.room
#
# pass.ids <- 1:n.pass
#
# room.ids <- 1:n.rooms
# room.ids.pass <- apportion_lr(n.pass, room.ids, rep(1/n.rooms, n.rooms))
# table(room.ids.pass)
#
# type.attr <- rep(c("p", "c"), times = c(n.pass, n.crew))
#
# crew.ids <- which(type.attr == "c")
# crew.ids
# room.ids
#
# # n.sectors <- 10
# room.sectors <- apportion_lr(length(room.ids), 1:n.sectors,
#                              rep(1/n.sectors, n.sectors))
# # head(data.frame(room.ids, room.sectors), 200)
# df.match <- data.frame(room.ids, room.sectors)
#
# sector <- rep(NA, n)
# for (id in pass.ids) {
#   sector[id] <- df.match[df.match$room.ids == room.ids.pass[id], "room.sectors"]
# }
# sector[pass.ids]
#
# room.sectors.c <- apportion_lr(length(crew.ids), 1:n.sectors,
#                                rep(1/n.sectors, n.sectors), shuffled = TRUE)
# sector[crew.ids] <- room.sectors.c
# sector[crew.ids]
# sector
#
# table(sector, type.attr)
#
# # median age for crew was 36 (IQR:29-43) and
# # the median age of the passengers was 69 (IQR: 62-73)
# ages.pass <- round(rnorm(n.pass, 69, 6), 1)
# summary(ages.pass)
# table(ages.pass)
# ages.crew <- pmax(round(rnorm(n.crew, 36, 10), 1), 18)
# summary(ages.crew)
# table(ages.crew)
#
# age <- rep(0, n)
# age[type.attr == "p"] <- ages.pass
# age[type.attr == "c"] <- ages.crew
# summary(age)
#
# # Initialize the network
# nw <- network.initialize(n, directed = FALSE)
# nw <- set.vertex.attribute(nw, "type", type.attr)
# nw <- set.vertex.attribute(nw, "pass.room", room.ids.pass, pass.ids)
# nw <- set.vertex.attribute(nw, "age", age)
# nw <- set.vertex.attribute(nw, "sector", sector)
#
#
# # Model 1. Pass/Pass Contacts ---------------------------------------------
#
# ## 1a. Pre-Lockdown Model
# md.pre <- 5
# edges.pre <- n.pass * md.pre/2
# target.stats1.pre <- c(edges.pre, edges.pre/md.pre)
#
# formation1.pre <- ~edges +
#                    nodematch("pass.room") +
#                    offset(nodefactor("type", levels = -2))
#
# coef.diss1 <- dissolution_coefs(dissolution = ~offset(edges), duration = 1)
# est1.pre <- netest(nw, formation1.pre, target.stats1.pre, coef.diss1, coef.form = -Inf,
#                    set.control.ergm = control.ergm(MCMLE.maxit = 500))
# summary(est1.pre)
#
# mcmc.diagnostics(est1.pre$fit)
# dx1.pre <- netdx(est1.pre, nsims = 1e4, dynamic = FALSE,
#                  nwstats.formula = ~edges + nodefactor("type", levels = NULL),
#                  set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e6))
# print(dx1.pre)
#
#
# ## 1b. Post-Lockdown Model: 20% of contact rate, all within cabin
# md.post <- 1
# edges.post <- n.pass * md.post/2
# target.stats1.post <- c(edges.post, edges.post, edges.post/md.post)
#
# formation1.post <- ~edges +
#                     nodematch("sector") +
#                     nodematch("pass.room") +
#                     offset(nodefactor("type", levels = -2))
#
# est1.post <- netest(nw, formation1.post, target.stats1.post, coef.diss1, coef.form = -Inf,
#                     set.control.ergm = control.ergm(MCMLE.maxit = 500))
# summary(est1.post)
#
# mcmc.diagnostics(est1.post$fit)
# dx1.post <- netdx(est1.post, nsims = 1e4, dynamic = FALSE,
#                   nwstats.formula = ~edges + nodefactor("type", levels = NULL) +
#                                      nodematch("sector", diff = FALSE) +
#                                      nodematch("pass.room", diff = FALSE),
#                   set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e6))
# print(dx1.post)
#
#
#
#
# # Model 2: Crew/Crew Contacts ---------------------------------------------
#
# ## 2a. Pre-Lockdown Model
# md.pre <- 10
# edges.pre <- n.crew * md.pre/2
# prop.within.sector.pre <- 0.50
# target.stats2.pre <- c(edges.pre, edges.pre*prop.within.sector.pre)
#
# formation2.pre <- ~edges +
#                    nodematch("sector") +
#                    offset(nodefactor("type", levels = -1))
# coef.diss2 <- dissolution_coefs(dissolution = ~offset(edges), duration = 1)
#
# est2.pre <- netest(nw, formation2.pre, target.stats2.pre, coef.diss2, coef.form = -Inf,
#                    set.control.ergm = control.ergm(MCMLE.maxit = 500))
# summary(est2.pre)
#
# mcmc.diagnostics(est2.pre$fit)
# dx2.pre <- netdx(est2.pre, nsims = 1e4, dynamic = FALSE,
#                   nwstats.formula = ~edges + nodefactor("type", levels = NULL),
#                   set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e6))
# print(dx2.pre)
#
#
# ## 2b. Post-Lockdown Model
# md.post <- 2
# edges.post <- n.crew * md.post/2
# target.stats2.post <- c(edges.post, edges.post*prop.within.sector.post)
#
# formation2.post <- ~edges +
#                     nodematch("sector") +
#                     offset(nodefactor("type", levels = -1))
#
# est2.post <- netest(nw, formation2.post, target.stats2.post, coef.diss2, coef.form = -Inf,
#                     set.control.ergm = control.ergm(MCMLE.maxit = 500))
# summary(est2.post)
#
# mcmc.diagnostics(est2.post$fit)
# dx2.post <- netdx(est2.post, nsims = 1e4, dynamic = FALSE,
#                   nwstats.formula = ~edges + nodefactor("type", levels = NULL),
#                   set.control.ergm = control.simulate.ergm(MCMC.burnin = 1e6))
# print(dx2.post)
#
#
# # Model 3: Crew/Pass Contacts ---------------------------------------------
#
# ## 3a. Pre-Lockdown Model
# md.pre <- 8
# edges.pre <- md.pre * n.rooms
# target.stats3.pre <- c(edges.pre, 0)
#
# formation3.pre <- ~edges + nodematch("type")
# coef.diss3 <- dissolution_coefs(dissolution = ~offset(edges), duration = 1)
#
# est3.pre <- netest(nw, formation3.pre, target.stats3.pre, coef.diss3,
#                    set.control.ergm = control.ergm(MCMLE.maxit = 500))
# summary(est3.pre)
#
# mcmc.diagnostics(est3.pre$fit)
# dx3.pre <- netdx(est3.pre, nsims = 1e4, dynamic = FALSE,
#                  nwstats.formula = ~edges + nodemix("type", levels = NULL) +
#                                     nodematch("sector"))
# print(dx3.pre)
#
#
# ## 3b. Post-Lockdown Model
# edges.post <- 3*n.rooms
# target.stats3 <- c(edges.post, edges.post*prop.within.sector.post, 0)
#
# formation3.post <- ~edges + nodematch("sector") + nodematch("type")
#
# est3.post <- netest(nw, formation3.post, target.stats3, coef.diss3,
#                     set.control.ergm = control.ergm(MCMLE.maxit = 500))
# summary(est3.post)
#
# mcmc.diagnostics(est3.post$fit)
# dx3.post <- netdx(est3.post, nsims = 1e4, dynamic = FALSE,
#                   nwstats.formula = ~edges + nodemix("type", levels = NULL))
# print(dx3.post)
#
#
# # Save Data ---------------------------------------------------------------
#
# # Pre-Lockdown
# est.pre <- list(est1.pre, est2.pre, est3.pre)
# saveRDS(est.pre, file = "est/est.pre.rds")
#
# # Post-Lockdown
# est.post <- list(est1.post, est2.post, est3.post)
# saveRDS(est.post, file = "est/est.post.base.rds")
