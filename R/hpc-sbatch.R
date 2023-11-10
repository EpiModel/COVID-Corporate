
# sbatch_master test script

# Set up a list of vectors associated with your free parameters.
# Below is an example of setting up four parameters with lower and upper bounds with each parameter.
# There are different ways to set up this list, but make sure each parameter length is the same
# with this approach.
iso.vars <- list(
  iso.prob.x =    seq(0, 1, length.out = 11),
  notif.work =    seq(0, 1, length.out = 11),
  notif.com =     0.1,
  notif.home =    0.9,
  act.rate.iso =  0.3,
  inf.prob.mask = 0.3
)

# then that list gets passed into this sbatch_function
EpiModelHPC::sbatch_master(
  vars = iso.vars, # the list above
  expand.vars = TRUE, # do not create a run for each combination of parameters
  master.file = "master.sh", # function creates a new file with this name in your working directory
  runsim.file = "runsim_corp.sh", # this is the name of your runsim file
  simno.start = 6000, # to name each parameter set, a starting number for your runs
  nsims = 32, # number of total simulations per parameter set (32 is fine for testing now)
  ncores = 32, # number of cores per parameter set simulation
  append = TRUE, # If TRUE, you can write more parameter set lines at the bottom of an existing master.sh file
  mem = "0", # Good default for memory use
  walltime = "1:00:00", # Good default for run time
  partition.main = "epimodel", # Name of partition for runs
  account.main = "epimodel", # Name of account for runs
)

# set scenarios
iso.vars <- list(
  iso.prob.x =    c(0.58, 1.0, rep(0.58,3), rep(1.0,3), rep(0.58,3) ),
  notif.work =    c(rep(0.5,2), 1.0, rep(0.5,2), 1.0, rep(0.5,5) ),
  notif.com =     c(rep(0.1,3), 1.0, rep(0.1,2), 1.0, rep(0.1,4)),
  notif.home =    c(rep(0.9,4), 1.0, rep(0.9,2), 1.0, rep(0.9,3)),
  act.rate.iso =  c(rep(0.3,8), 0, 0.3, 0),
  inf.prob.mask = c(rep(0.3,9), 0, 0)
)

# then that list gets passed into this sbatch_function
EpiModelHPC::sbatch_master(
  vars = iso.vars, # the list above
  expand.vars = FALSE, # do not create a run for each combination of parameters
  master.file = "master.sh", # function creates a new file with this name in your working directory
  runsim.file = "runsim_corp.sh", # this is the name of your runsim file
  simno.start = 5000, # to name each parameter set, a starting number for your runs
  nsims = 32, # number of total simulations per parameter set (32 is fine for testing now)
  ncores = 32, # number of cores per parameter set simulation
  append = TRUE, # If TRUE, you can write more parameter set lines at the bottom of an existing master.sh file
  mem = "0", # Good default for memory use
  walltime = "1:00:00", # Good default for run time
  partition.main = "epimodel", # Name of partition for runs
  account.main = "epimodel", # Name of account for runs
)
