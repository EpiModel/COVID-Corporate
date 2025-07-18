
### NEW SBATCH FOR SPLIT ISO.PROB BTWN OFFICE WORKERS AND OTHERS

# comparing isolation prob in office workers vs notification prob among office workers ----
iso.vars <- list(
  iso.prob.off =  seq(0, 1, length.out = 11),
  iso.prob.non =  0.58,
  notif.work =    seq(0, 1, length.out = 11),
  notif.com =     0.1,
  notif.home =    0.9,
  act.rate.iso =  0.3,
  inf.prob.mask = 0.3
)

EpiModelHPC::sbatch_master(
  vars = iso.vars, # the list above
  expand.vars = TRUE, # create a run for each combination of parameters
  master.file = "master2.sh", # function creates a new file with this name in your working directory
  runsim.file = "runsim_corp.sh", # this is the name of your runsim file
  simno.start = 16000, # to name each parameter set, a starting number for your runs
  nsims = 250, # number of total simulations per parameter set (32 is fine for testing now)
  ncores = 32, # number of cores per parameter set simulation
  #narray = 1,
  append = TRUE, # If TRUE, you can write more parameter set lines at the bottom of an existing master.sh file
  mem = "0", # Good default for memory use
  walltime = "1:00:00", # Good default for run time
  partition.main = "epimodel", # Name of partition for runs
  account.main = "epimodel", # Name of account for runs
)

# comparing isolation prob in office workers vs contact intensity ----
iso.vars <- list(
  iso.prob.off =  seq(0, 1, length.out = 11),
  iso.prob.non =  0.58,
  notif.work =    0.5,
  notif.com =     0.1,
  notif.home =    0.9,
  act.rate.iso =  seq(0, 1, length.out = 11),
  inf.prob.mask = 0.3
)

EpiModelHPC::sbatch_master(
  vars = iso.vars, # the list above
  expand.vars = TRUE, # create a run for each combination of parameters
  master.file = "master2.sh", # function creates a new file with this name in your working directory
  runsim.file = "runsim_corp.sh", # this is the name of your runsim file
  simno.start = 17000, # to name each parameter set, a starting number for your runs
  nsims = 250, # number of total simulations per parameter set (32 is fine for testing now)
  ncores = 32, # number of cores per parameter set simulation
  #narray = 1,
  append = TRUE, # If TRUE, you can write more parameter set lines at the bottom of an existing master.sh file
  mem = "0", # Good default for memory use
  walltime = "1:00:00", # Good default for run time
  partition.main = "epimodel", # Name of partition for runs
  account.main = "epimodel", # Name of account for runs
)

# comparing isolation prob in office workers vs non-office-workers ----
iso.vars <- list(
  iso.prob.off =  seq(0, 1, length.out = 11),
  iso.prob.non =  seq(0, 1, length.out = 11),
  notif.work =    0.5,
  notif.com =     0.1,
  notif.home =    0.9,
  act.rate.iso =  0.3,
  inf.prob.mask = 0.3
)

EpiModelHPC::sbatch_master(
  vars = iso.vars, # the list above
  expand.vars = TRUE, # create a run for each combination of parameters
  master.file = "master2.sh", # function creates a new file with this name in your working directory
  runsim.file = "runsim_corp.sh", # this is the name of your runsim file
  simno.start = 11000, # to name each parameter set, a starting number for your runs
  nsims = 250, # number of total simulations per parameter set (32 is fine for testing now)
  ncores = 32, # number of cores per parameter set simulation
  #narray = 1,
  append = TRUE, # If TRUE, you can write more parameter set lines at the bottom of an existing master.sh file
  mem = "0", # Good default for memory use
  walltime = "1:00:00", # Good default for run time
  partition.main = "epimodel", # Name of partition for runs
  account.main = "epimodel", # Name of account for runs
)

# Final set of scenarios for results ----
iso.vars <- list(
  iso.prob.off =    c(0.58,0,0.25,0.5,0.75,1,rep(0.58,5)),
  iso.prob.non =    c(rep(0.58,6),0,0.25,0.5,0.75,1),
  notif.work =    c(rep(0.5,11)),
  notif.com =     c(rep(0.1,11)),
  notif.home =    c(rep(0.9,11)),
  act.rate.iso =  c(rep(0.3,11)),
  inf.prob.mask = c(rep(0.3,11))
)

# then that list gets passed into this sbatch_function
EpiModelHPC::sbatch_master(
  vars = iso.vars, # the list above
  expand.vars = FALSE, # do not create a run for each combination of parameters
  master.file = "master2.sh", # function creates a new file with this name in your working directory
  runsim.file = "runsim_corp.sh", # this is the name of your runsim file
  simno.start = 18000, # to name each parameter set, a starting number for your runs
  nsims = 250, # number of total simulations per parameter set (32 is fine for testing now)
  ncores = 32, # number of cores per parameter set simulation
  #narray = 1,
  append = TRUE, # If TRUE, you can write more parameter set lines at the bottom of an existing master.sh file
  mem = "0", # Good default for memory use
  walltime = "1:00:00", # Good default for run time
  partition.main = "epimodel", # Name of partition for runs
  account.main = "epimodel", # Name of account for runs
)
