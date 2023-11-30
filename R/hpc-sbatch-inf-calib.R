
# sbatch_master test script

# Set up a list of vectors associated with your free parameters.
# Below is an example of setting up four parameters with lower and upper bounds with each parameter.
# There are different ways to set up this list, but make sure each parameter length is the same
# with this approach.
vax.vars <- list(
  inf.p = rep(0.4,10),
  boost.start1 = rep(-60,10),
  boost.start2 = rep(325,10), # changed from 335
  boost.stop1 = rep(15,10),
  boost.stop2 = rep(396,10),
  inf.add1 = rep(0.45,10),
  inf.add2 = rep(0.1,10),
  supp.start1 = rep(30,10),
  supp.start2 = rep(250,10), #changed this from 258
  supp.start3 = rep(397,10),
  supp.start4 = rep(486,10),
  supp.stop1 = rep(165,10), #changed this from 180
  supp.stop2 = rep(324,10), # changed from 334
  supp.stop3 = rep(485,10),
  supp.stop4 = rep(608,10),
  inf.sub1 = rep(0.22,10),
  inf.sub2 = rep(0.3,10),
  inf.sub3 = rep(0.23,10),
  inf.sub4 = rep(0.12,10),
  hosp.b.mult = rep(1,10),
  hosp.b.start = rep(320,10), #changed from 378
  hosp.b.stop = rep(396,10), #changed from 423
  hosp.s.mult = rep(0.75,10),
  hosp.s.start = rep(397,10),
  hosp.s.stop = rep(604,10),
  hosp.div = seq(5,10,length.out=10),
  mort.mult = rep(180*10,10)
)

### Probably need to try further decreasing inf.p (maybe 0.3 or 0.35)
## this will make the first case spike match better
## then will need to adjust the inf.add/inf.sub

# then that list gets passed into this sbatch_function
EpiModelHPC::sbatch_master(
  vars = vax.vars, # the list above
  expand.vars = FALSE, # do not create a run for each combination of parameters
  master.file = "master_inf_calib.sh", # function creates a new file with this name in your working directory
  runsim.file = "runsim_inf_calib.sh", # this is the name of your runsim file
  simno.start = 2190, # to name each parameter set, a starting number for your runs
  nsims = 32, # number of total simulations per parameter set (32 is fine for testing now)
  ncores = 32, # number of cores per parameter set simulation
  append = TRUE, # If TRUE, you can write more parameter set lines at the bottom of an existing master.sh file
  mem = "0", # Good default for memory use
  walltime = "1:00:00", # Good default for run time
  partition.main = "epimodel", # Name of partition for runs
  account.main = "epimodel", # Name of account for runs
)
