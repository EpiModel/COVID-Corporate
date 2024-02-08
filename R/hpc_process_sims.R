# Load libraries ---------------------------------------------------------------
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(EpiModelCOVID)))
library(future.apply)

# Set baseline simulation number and directory ---------------------------------
base.sim <- 5000
my.dir <- "sims/"
s.num <- 69 # number of scenarios, excluding reference scenario

# Function to merge sim files (same as in 04-calib-figures) --------------------
merge_simfiles <- function(simno, indir = "sims/", vars = NULL,  truncate.at = NULL, verbose = TRUE) {

  fn <- list.files(indir, pattern = paste0("sim.*", simno, ".[0-9].rds"), full.names = TRUE)

  if (length(fn) == 0) {
    stop("No files of that simno in the specified indir", call. = FALSE)
  }

  for (i in seq_along(fn)) {
    sim <- readRDS(fn[i])

    if (!is.null(truncate.at)) {
      sim <- truncate_sim(sim, truncate.at)
    }

    sim$network <- NULL
    sim$attr <- NULL
    sim$temp <- NULL
    sim$el <- NULL
    sim$p <- NULL

    if (inherits(sim, "list") && all(c("epi", "param", "control") %in% names(sim))) {
      class(sim) <- "netsim"
    }

    if (!is.null(vars)) {
      sim$epi <- sim$epi[vars]
      sim$stats <- NULL
      if (!is.null(sim$riskh)) {
        sim$riskh <- NULL
      }
    }

    if (i == 1) {
      out <- sim
    } else {
      out <- merge(out, sim, param.error = FALSE, keep.other = FALSE)
    }

    if (verbose == TRUE) {
      cat("File ", i, "/", length(fn), " Loaded ... \n", sep = "")
    }
  }

  return(out)
}

# Tables: FOR NEW OUTPUT Merge & save files ------------------------------------

baseline_sims <- merge_simfiles(base.sim, indir = "ref_sim/",
                                truncate.at = 181, verbose = TRUE)
# scenario0 <- baseline_sims

saveRDS(baseline_sims, file = here::here("processed_sims",
                                         "baseline_sim.rds"))

# for (i in 0:s.num) {
#   scenario_list[[i]] <- merge_simfiles(base.sim + i, indir = my.dir,
#                                        vars = NULL,  truncate.at = 181, verbose = TRUE)
#
# }

future::plan("multisession", workers = 4)

scenario_list <- future.apply::future_lapply(
  0:s.num, function(i) {
    merge_simfiles(base.sim + i, indir = my.dir, vars = NULL,
                   truncate.at = 181, verbose = TRUE)
  }
)

saveRDS(scenario_list, file = here::here("processed_sims",
                                         "scenario_list.rds"))


# Contour Plots: Load and merge files ------------------------------------------

##Isolation prob vs work notification prob
cInc <- data.frame(iso.prob = rep(NA, 121), work.notif.prob = rep(NA, 121), cInc = rep(NA, 121))
cDeaths <- data.frame(iso.prob = rep(NA, 121), work.notif.prob = rep(NA, 121), cDeaths = rep(NA, 121))
cIso <- data.frame(iso.prob = rep(NA, 121), work.notif.prob = rep(NA, 121), cIso = rep(NA, 121))
cCases <- data.frame(iso.prob = rep(NA, 121), work.notif.prob = rep(NA, 121), cCases = rep(NA, 121))

cIncWork <- data.frame(iso.prob = rep(NA, 121), work.notif.prob = rep(NA, 121), cInc = rep(NA, 121))
cDeathsWork <- data.frame(iso.prob = rep(NA, 121), work.notif.prob = rep(NA, 121), cDeaths = rep(NA, 121))
cIsoWork <- data.frame(iso.prob = rep(NA, 121), work.notif.prob = rep(NA, 121), cIso = rep(NA, 121))
cCasesWork <- data.frame(iso.prob = rep(NA, 121), work.notif.prob = rep(NA, 121), cCases = rep(NA, 121))

for (i in 6000:6120) {

  print(i)

  temp <- merge_simfiles(i, indir = my.dir, vars = NULL,  truncate.at = 181, verbose = TRUE)

  scale.1 <- temp$param$iso.prob
  scale.2 <- temp$param$notif.prob[1]

  scenario.cInc <- median(colSums(temp$epi$se.flow) / colSums(temp$epi$num) * 100000)
  scenario.cDeaths <- median(colSums(temp$epi$d.h.flow) / colSums(temp$epi$num) * 100000)
  scenario.cIso <- median(colSums(temp$epi$iso1.flow))
  scenario.cCases <- median(colSums(temp$epi$ipic.flow) / colSums(temp$epi$num) * 100000)

  scenario.cIncWork <- median(colSums(temp$epi$se.flow.l1) / colSums(temp$epi$num.w) * 100000)
  # scenario.cDeathsWork <- median(colSums(temp$epi$d.h.flow.w) / colSums(temp$epi$num.w) * 100000)
  scenario.cIsoWork <- median(colSums(temp$epi$iso1.flow.w))
  scenario.cCasesWork <- median(colSums(temp$epi$ipic.flow.w) / colSums(temp$epi$num.w) * 100000)

  cInc[(i - 5999), ] <- c(scale.1, scale.2, scenario.cInc)
  cDeaths[(i - 5999), ] <- c(scale.1, scale.2, scenario.cDeaths)
  cIso[(i - 5999), ] <- c(scale.1, scale.2, scenario.cIso)
  cCases[(i - 5999), ] <- c(scale.1, scale.2, scenario.cCases)

  cIncWork[(i - 5999), ] <- c(scale.1, scale.2, scenario.cIncWork)
  # cDeathsWork[(i - 5999), ] <- c(scale.1, scale.2, scenario.cDeathsWork)
  cIsoWork[(i - 5999), ] <- c(scale.1, scale.2, scenario.cIsoWork)
  cCasesWork[(i - 5999), ] <- c(scale.1, scale.2, scenario.cCasesWork)

  remove(temp)
}

saveRDS(cInc, file = here::here("processed_sims","cInc.rds"))
saveRDS(cDeaths, file = here::here("processed_sims","cDeaths.rds"))
saveRDS(cIso, file = here::here("processed_sims","cIso.rds"))
saveRDS(cCases, file = here::here("processed_sims","cCases.rds"))
saveRDS(cIncWork, file = here::here("processed_sims","cIncWork.rds"))
# saveRDS(cDeathsWork, file = here::here("processed_sims","cDeathsWork.rds"))
saveRDS(cIsoWork, file = here::here("processed_sims","cIsoWork.rds"))
saveRDS(cCasesWork, file = here::here("processed_sims","cCasesWork.rds"))
