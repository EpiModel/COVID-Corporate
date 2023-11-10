
# sbatch_master test script

# Set up a list of vectors associated with your free parameters.
# Below is an example of setting up four parameters with lower and upper bounds with each parameter.
# There are different ways to set up this list, but make sure each parameter length is the same
# with this approach.
vax.vars <- list(
  boost.start1 = rep(-60,10),
  boost.start2 = rep(335,10),
  boost.stop1 = rep(15,10),
  boost.stop2 = rep(396,10),
  inf.add1 = seq(0.15,0.35,length.out=10),
  inf.add2 = seq(0.30,0.35,length.out=10),
  supp.start1 = seq(16,60,length.out=10),
  supp.start2 = rep(258,10),
  supp.start3 = rep(397,10),
  supp.stop1 = rep(180,10),
  supp.stop2 = rep(334,10),
  supp.stop3 = rep(Inf,10),
  inf.sub1 = rep(0.2,10),
  inf.sub2 = rep(0.2,10),
  inf.sub3 = rep(0.2,10),
  hosp.b.mult = rep(2.25,10),
  hosp.b.start = seq(378,450,length.out=10),
  hosp.b.stop = rep(423,10),
  hosp.s.mult = seq(0.5,0.75,length.out=10),
  hosp.s.start = rep(546,10),
  hosp.s.stop = rep(604,10)
)

# then that list gets passed into this sbatch_function
EpiModelHPC::sbatch_master(
  vars = vax.vars, # the list above
  expand.vars = FALSE, # do not create a run for each combination of parameters
  master.file = "master_inf_calib.sh", # function creates a new file with this name in your working directory
  runsim.file = "runsim_inf_calib.sh", # this is the name of your runsim file
  simno.start = 2000, # to name each parameter set, a starting number for your runs
  nsims = 32, # number of total simulations per parameter set (32 is fine for testing now)
  ncores = 32, # number of cores per parameter set simulation
  append = TRUE, # If TRUE, you can write more parameter set lines at the bottom of an existing master.sh file
  mem = "0", # Good default for memory use
  walltime = "1:00:00", # Good default for run time
  partition.main = "epimodel", # Name of partition for runs
  account.main = "epimodel", # Name of account for runs
)
