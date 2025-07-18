
# sbatch_master test script

# Set up a list of vectors associated with your free parameters.
# Below is an example of setting up four parameters with lower and upper bounds with each parameter.
# There are different ways to set up this list, but make sure each parameter length is the same
# with this approach.
vax.vars <- list(
  vax1.rate1 = seq(0.0003, 0.0007, length.out = 10),
  vax1.rate2 = seq(0.0009, 0.0009, length.out = 10), #final
  vax1.rate3 = seq(0.0019, 0.0019, length.out = 10), #final
  vax1.rate4 = seq(0.0019, 0.0019, length.out = 10), #final
  vax1.rate5 = seq(0.003, 0.005, length.out = 10),

  vax2.rate1 = seq(0.005, 0.010, length.out = 10),
  vax2.rate2 = seq(0.011, 0.011, length.out = 10), #final
  vax2.rate3 = seq(0.005, 0.02, length.out = 10),
  vax2.rate4 = seq(0.005, 0.02, length.out = 10),
  vax2.rate5 = seq(0.005, 0.02, length.out = 10),

  vax3.rate1 = seq(0, 0, length.out = 10),
  vax3.rate2 = seq(0.0025, 0.0035, length.out = 10),
  vax3.rate3 = seq(0.0006, 0.0008, length.out = 10),
  vax3.rate4 = seq(0.002, 0.003, length.out = 10),
  vax3.rate5 = seq(0.002, 0.003, length.out = 10),

  vax4.rate1 = seq(0, 0, length.out = 10),
  vax4.rate2 = seq(0, 0, length.out = 10),
  vax4.rate3 = seq(0, 0, length.out = 10),
  vax4.rate4 = seq(0.0016, 0.0020, length.out = 10),
  vax4.rate5 = seq(0.0013, 0.0018, length.out = 10)
)

# then that list gets passed into this sbatch_function
EpiModelHPC::sbatch_master(
  vars = vax.vars, # the list above
  expand.vars = FALSE, # do not create a run for each combination of parameters
  master.file = "master_vax_calib.sh", # function creates a new file with this name in your working directory
  runsim.file = "runsim_vax_calib.sh", # this is the name of your runsim file
  simno.start = 1040, # to name each parameter set, a starting number for your runs
  nsims = 32, # number of total simulations per parameter set (32 is fine for testing now)
  ncores = 32, # number of cores per parameter set simulation
  append = TRUE, # If TRUE, you can write more parameter set lines at the bottom of an existing master.sh file
  mem = "0", # Good default for memory use
  walltime = "1:00:00", # Good default for run time
  partition.main = "epimodel", # Name of partition for runs
  account.main = "epimodel", # Name of account for runs
)
