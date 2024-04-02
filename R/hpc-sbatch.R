
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
  expand.vars = TRUE, # create a run for each combination of parameters
  master.file = "master.sh", # function creates a new file with this name in your working directory
  runsim.file = "runsim_corp.sh", # this is the name of your runsim file
  simno.start = 6000, # to name each parameter set, a starting number for your runs
  nsims = 250, # number of total simulations per parameter set (32 is fine for testing now)
  ncores = 32, # number of cores per parameter set simulation
  #narray = 1,
  append = TRUE, # If TRUE, you can write more parameter set lines at the bottom of an existing master.sh file
  mem = "0", # Good default for memory use
  walltime = "1:00:00", # Good default for run time
  partition.main = "epimodel", # Name of partition for runs
  account.main = "epimodel", # Name of account for runs
)

iso.vars <- list(
  iso.prob.x =    seq(0, 1, length.out = 11),
  notif.work =    0.5,
  notif.com =     0.1,
  notif.home =    0.9,
  act.rate.iso =  seq(0, 1, length.out = 11),
  inf.prob.mask = 0.3
)

# then that list gets passed into this sbatch_function
EpiModelHPC::sbatch_master(
  vars = iso.vars, # the list above
  expand.vars = TRUE, # create a run for each combination of parameters
  master.file = "master.sh", # function creates a new file with this name in your working directory
  runsim.file = "runsim_corp.sh", # this is the name of your runsim file
  simno.start = 7000, # to name each parameter set, a starting number for your runs
  nsims = 250, # number of total simulations per parameter set (32 is fine for testing now)
  ncores = 32, # number of cores per parameter set simulation
  #narray = 1,
  append = TRUE, # If TRUE, you can write more parameter set lines at the bottom of an existing master.sh file
  mem = "0", # Good default for memory use
  walltime = "1:00:00", # Good default for run time
  partition.main = "epimodel", # Name of partition for runs
  account.main = "epimodel", # Name of account for runs
)

# # set scenarios
# iso.vars <- list(
#   iso.prob.x =    c(0.58, 1.0, rep(0.58,3), rep(1.0,3), rep(0.58,3) ),
#   notif.work =    c(rep(0.5,2), 1.0, rep(0.5,2), 1.0, rep(0.5,5) ),
#   notif.com =     c(rep(0.1,3), 1.0, rep(0.1,2), 1.0, rep(0.1,4)),
#   notif.home =    c(rep(0.9,4), 1.0, rep(0.9,2), 1.0, rep(0.9,3)),
#   act.rate.iso =  c(rep(0.3,8), 0, 0.3, 0),
#   inf.prob.mask = c(rep(0.3,9), 0, 0)
# )

iso.vars <- list(
  iso.prob.x =    c(0.58,0,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.9,1,
                    rep(0.58,10),
                    rep(0.58,10),
                    rep(0.58,10),
                    rep(0.58,10),
                    rep(0.58,10),
                    0,1,0.75,0.75,0.75,0.75),
  notif.work =    c(rep(0.5,14),
                    0,0.1,0.2,0.3,0.4,0.6,0.7,0.8,0.9,1,
                    rep(0.5,10),
                    rep(0.5,10),
                    rep(0.5,10),
                    rep(0.5,10),
                    0,1,0.75,0.5,0.5,0.5),
  notif.com =     c(rep(0.1,14),
                    rep(0.1,10),
                    0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,
                    rep(0.1,10),
                    rep(0.1,10),
                    rep(0.1,10),
                    0,1,0.1,0.1,0.1,0.1),
  notif.home =    c(rep(0.9,14),
                    rep(0.9,10),
                    rep(0.9,10),
                    0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1,
                    rep(0.9,10),
                    rep(0.9,10),
                    0,1,0.9,0.9,0.9,0.9),
  act.rate.iso =  c(rep(0.3,14),
                    rep(0.3,10),
                    rep(0.3,10),
                    rep(0.3,10),
                    0,0.1,0.2,0.4,0.5,0.6,0.7,0.8,0.9,1,
                    rep(0.3,10),
                    1,0,0.3,0.1,0.3,0.1),
  inf.prob.mask = c(rep(0.3,14),
                    rep(0.3,10),
                    rep(0.3,10),
                    rep(0.3,10),
                    rep(0.3,10),
                    0,0.1,0.2,0.4,0.5,0.6,0.7,0.8,0.9,1,
                    1,0,0.3,0.3,0.1,0.1)
)

# then that list gets passed into this sbatch_function
EpiModelHPC::sbatch_master(
  vars = iso.vars, # the list above
  expand.vars = FALSE, # do not create a run for each combination of parameters
  master.file = "master.sh", # function creates a new file with this name in your working directory
  runsim.file = "runsim_corp.sh", # this is the name of your runsim file
  simno.start = 5000, # to name each parameter set, a starting number for your runs
  nsims = 250, # number of total simulations per parameter set (32 is fine for testing now)
  ncores = 32, # number of cores per parameter set simulation
  #narray = 1,
  append = TRUE, # If TRUE, you can write more parameter set lines at the bottom of an existing master.sh file
  mem = "0", # Good default for memory use
  walltime = "1:00:00", # Good default for run time
  partition.main = "epimodel", # Name of partition for runs
  account.main = "epimodel", # Name of account for runs
)


# Final set of scenarios for results
iso.vars <- list(
  iso.prob.x =    c(0.58,0,0.25,0.5,0.75,1),
  notif.work =    c(rep(0.5,6)),
  notif.com =     c(rep(0.1,6)),
  notif.home =    c(rep(0.9,6)),
  act.rate.iso =  c(rep(0.3,6)),
  inf.prob.mask = c(rep(0.3,6))
)

# then that list gets passed into this sbatch_function
EpiModelHPC::sbatch_master(
  vars = iso.vars, # the list above
  expand.vars = FALSE, # do not create a run for each combination of parameters
  master.file = "master.sh", # function creates a new file with this name in your working directory
  runsim.file = "runsim_corp.sh", # this is the name of your runsim file
  simno.start =8000, # to name each parameter set, a starting number for your runs
  nsims = 250, # number of total simulations per parameter set (32 is fine for testing now)
  ncores = 32, # number of cores per parameter set simulation
  #narray = 1,
  append = TRUE, # If TRUE, you can write more parameter set lines at the bottom of an existing master.sh file
  mem = "0", # Good default for memory use
  walltime = "1:00:00", # Good default for run time
  partition.main = "epimodel", # Name of partition for runs
  account.main = "epimodel", # Name of account for runs
)
