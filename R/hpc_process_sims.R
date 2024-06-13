# Load libraries ---------------------------------------------------------------
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(EpiModelCOVID)))
library(future.apply)

# Set baseline simulation number and directory ---------------------------------
base.sim <- 18000
my.dir <- "sims/"
s.num <- 10 # number of scenarios, excluding reference scenario

# Function to merge sim files (same as in 04-calib-figures) --------------------
merge_simfiles <- function(simno, indir = "sims/", vars = NULL,  truncate.at = NULL, verbose = TRUE) {

  fn <- list.files(indir, pattern = paste0("sim.*", simno, ".[0-9].rds"), full.names = TRUE)

  if (length(fn) == 0) {
    stop(paste0("No files with simno ",simno," in the specified indir"),
         call. = FALSE)
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

baseline_sims <- merge_simfiles(base.sim, indir = "sims/",
                                truncate.at = 181, verbose = TRUE)
# scenario0 <- baseline_sims

saveRDS(baseline_sims, file = here::here("processed_sims",
                                         "baseline_sim2.rds"))

for (i in 0:s.num) {
  scenario_list[[i]] <- merge_simfiles(base.sim + i, indir = my.dir,
                                       vars = NULL,  truncate.at = 181, verbose = TRUE)

}

future::plan("multisession", workers = 4)

scenario_list <- future.apply::future_lapply(
  0:s.num, function(i) {
    merge_simfiles(base.sim + i, indir = my.dir, vars = NULL,
                   truncate.at = 181, verbose = TRUE)
  }
)

saveRDS(scenario_list, file = here::here("processed_sims",
                                         "scenario_list2.rds"))


# Contour Plots: Isolation prob vs work notification prob ----------------------

WcInc <- data.frame(iso.prob.off = rep(NA, 121), work.notif.prob = rep(NA, 121), cInc = rep(NA, 121))
WcDeaths <- data.frame(iso.prob.off = rep(NA, 121), work.notif.prob = rep(NA, 121), cDeaths = rep(NA, 121))
WcIso <- data.frame(iso.prob.off = rep(NA, 121), work.notif.prob = rep(NA, 121), cIso = rep(NA, 121))
WcCases <- data.frame(iso.prob.off = rep(NA, 121), work.notif.prob = rep(NA, 121), cCases = rep(NA, 121))

WcIncWork <- data.frame(iso.prob.off = rep(NA, 121), work.notif.prob = rep(NA, 121), cInc = rep(NA, 121))
WcDeathsWork <- data.frame(iso.prob.off = rep(NA, 121), work.notif.prob = rep(NA, 121), cDeaths = rep(NA, 121))
WcIsoWork <- data.frame(iso.prob.off = rep(NA, 121), work.notif.prob = rep(NA, 121), cIso = rep(NA, 121))
WcCasesWork <- data.frame(iso.prob.off = rep(NA, 121), work.notif.prob = rep(NA, 121), cCases = rep(NA, 121))

for (i in 16000:16120) {

  print(i)

  temp <- merge_simfiles(i, indir = my.dir, vars = NULL,  truncate.at = 181, verbose = TRUE)

  scale.1 <- temp$param$iso.prob.office
  scale.2 <- temp$param$notif.prob[1]

  scenario.cInc <- median(colSums(temp$epi$se.flow) / colSums(temp$epi$num) * 100000)
  scenario.cDeaths <- median(colSums(temp$epi$d.h.flow) / colSums(temp$epi$num) * 100000)
  scenario.cIso <- median(colSums(temp$epi$iso1.flow))
  scenario.cCases <- median(colSums(temp$epi$ipic.flow) / colSums(temp$epi$num) * 100000)

  scenario.cIncWork <- median(colSums(temp$epi$se.flow.l1) / colSums(temp$epi$num.w) * 100000)
  scenario.cDeathsWork <- median(colSums(temp$epi$d.h.flow.w) / colSums(temp$epi$num.w) * 100000)
  scenario.cIsoWork <- median(colSums(temp$epi$iso1.flow.w))
  scenario.cCasesWork <- median(colSums(temp$epi$ipic.flow.w) / colSums(temp$epi$num.w) * 100000)

  WcInc[(i - 15999), ] <- c(scale.1, scale.2, scenario.cInc)
  WcDeaths[(i - 15999), ] <- c(scale.1, scale.2, scenario.cDeaths)
  WcIso[(i - 15999), ] <- c(scale.1, scale.2, scenario.cIso)
  WcCases[(i - 15999), ] <- c(scale.1, scale.2, scenario.cCases)

  WcIncWork[(i - 15999), ] <- c(scale.1, scale.2, scenario.cIncWork)
  WcDeathsWork[(i - 15999), ] <- c(scale.1, scale.2, scenario.cDeathsWork)
  WcIsoWork[(i - 15999), ] <- c(scale.1, scale.2, scenario.cIsoWork)
  WcCasesWork[(i - 15999), ] <- c(scale.1, scale.2, scenario.cCasesWork)

  remove(temp)
}

saveRDS(WcInc, file = here::here("processed_sims","WcInc.rds"))
saveRDS(WcDeaths, file = here::here("processed_sims","WcDeaths.rds"))
saveRDS(WcIso, file = here::here("processed_sims","WcIso.rds"))
saveRDS(WcCases, file = here::here("processed_sims","WcCases.rds"))
saveRDS(WcIncWork, file = here::here("processed_sims","WcIncWork.rds"))
saveRDS(WcDeathsWork, file = here::here("processed_sims","WcDeathsWork.rds"))
saveRDS(WcIsoWork, file = here::here("processed_sims","WcIsoWork.rds"))
saveRDS(WcCasesWork, file = here::here("processed_sims","WcCasesWork.rds"))

# Contour Plots: Isolation prob vs iso act rate --------------------------------

WcInc2 <- data.frame(iso.prob.off = rep(NA, 121), act.rate.iso = rep(NA, 121), cInc = rep(NA, 121))
WcDeaths2 <- data.frame(iso.prob.off = rep(NA, 121), act.rate.iso = rep(NA, 121), cDeaths = rep(NA, 121))
WcIso2 <- data.frame(iso.prob.off = rep(NA, 121), act.rate.iso = rep(NA, 121), cIso = rep(NA, 121))
WcCases2 <- data.frame(iso.prob.off = rep(NA, 121), act.rate.iso = rep(NA, 121), cCases = rep(NA, 121))

WcIncWork2 <- data.frame(iso.prob.off = rep(NA, 121), act.rate.iso = rep(NA, 121), cInc = rep(NA, 121))
WcDeathsWork2 <- data.frame(iso.prob.off = rep(NA, 121), act.rate.iso = rep(NA, 121), cDeaths = rep(NA, 121))
WcIsoWork2 <- data.frame(iso.prob.off = rep(NA, 121), act.rate.iso = rep(NA, 121), cIso = rep(NA, 121))
WcCasesWork2 <- data.frame(iso.prob.off = rep(NA, 121), act.rate.iso = rep(NA, 121), cCases = rep(NA, 121))

for (i in 17000:17120) {

  print(i)

  temp <- merge_simfiles(i, indir = my.dir, vars = NULL,  truncate.at = 181, verbose = TRUE)

  scale.1 <- temp$param$iso.prob.office
  scale.2 <- temp$param$act.rate.iso.inter.rr

  scenario.cInc <- median(colSums(temp$epi$se.flow) / colSums(temp$epi$num) * 100000)
  scenario.cDeaths <- median(colSums(temp$epi$d.h.flow) / colSums(temp$epi$num) * 100000)
  scenario.cIso <- median(colSums(temp$epi$iso1.flow))
  scenario.cCases <- median(colSums(temp$epi$ipic.flow) / colSums(temp$epi$num) * 100000)

  scenario.cIncWork <- median(colSums(temp$epi$se.flow.l1) / colSums(temp$epi$num.w) * 100000)
  scenario.cDeathsWork <- median(colSums(temp$epi$d.h.flow.w) / colSums(temp$epi$num.w) * 100000)
  scenario.cIsoWork <- median(colSums(temp$epi$iso1.flow.w))
  scenario.cCasesWork <- median(colSums(temp$epi$ipic.flow.w) / colSums(temp$epi$num.w) * 100000)

  WcInc2[(i - 16999), ] <- c(scale.1, scale.2, scenario.cInc)
  WcDeaths2[(i - 16999), ] <- c(scale.1, scale.2, scenario.cDeaths)
  WcIso2[(i - 16999), ] <- c(scale.1, scale.2, scenario.cIso)
  WcCases2[(i - 16999), ] <- c(scale.1, scale.2, scenario.cCases)

  WcIncWork2[(i - 16999), ] <- c(scale.1, scale.2, scenario.cIncWork)
  WcDeathsWork2[(i - 16999), ] <- c(scale.1, scale.2, scenario.cDeathsWork)
  WcIsoWork2[(i - 16999), ] <- c(scale.1, scale.2, scenario.cIsoWork)
  WcCasesWork2[(i - 16999), ] <- c(scale.1, scale.2, scenario.cCasesWork)

  remove(temp)
}

saveRDS(WcInc2, file = here::here("processed_sims","WcInc2.rds"))
saveRDS(WcDeaths2, file = here::here("processed_sims","WcDeaths2.rds"))
saveRDS(WcIso2, file = here::here("processed_sims","WcIso2.rds"))
saveRDS(WcCases2, file = here::here("processed_sims","WcCases2.rds"))
saveRDS(WcIncWork2, file = here::here("processed_sims","WcIncWork2.rds"))
saveRDS(WcDeathsWork2, file = here::here("processed_sims","WcDeathsWork2.rds"))
saveRDS(WcIsoWork2, file = here::here("processed_sims","WcIsoWork2.rds"))
saveRDS(WcCasesWork2, file = here::here("processed_sims","WcCasesWork2.rds"))

# Contour Plots: Office isolation prob vs Other isolation prob --------------------------------

WcInc3 <- data.frame(iso.prob.off = rep(NA, 121), iso.prob.non = rep(NA, 121), cInc = rep(NA, 121))
WcDeaths3 <- data.frame(iso.prob.off = rep(NA, 121), iso.prob.non = rep(NA, 121), cDeaths = rep(NA, 121))
WcIso3 <- data.frame(iso.prob.off = rep(NA, 121), iso.prob.non = rep(NA, 121), cIso = rep(NA, 121))
WcCases3 <- data.frame(iso.prob.off = rep(NA, 121), iso.prob.non = rep(NA, 121), cCases = rep(NA, 121))

WcIncWork3 <- data.frame(iso.prob.off = rep(NA, 121), iso.prob.non = rep(NA, 121), cInc = rep(NA, 121))
WcDeathsWork3 <- data.frame(iso.prob.off = rep(NA, 121), iso.prob.non = rep(NA, 121), cDeaths = rep(NA, 121))
WcIsoWork3 <- data.frame(iso.prob.off = rep(NA, 121), iso.prob.non = rep(NA, 121), cIso = rep(NA, 121))
WcCasesWork3 <- data.frame(iso.prob.off = rep(NA, 121), iso.prob.non = rep(NA, 121), cCases = rep(NA, 121))

for (i in 11000:11120) {

  print(i)

  temp <- merge_simfiles(i, indir = my.dir, vars = NULL,  truncate.at = 181, verbose = TRUE)

  scale.1 <- temp$param$iso.prob.office
  scale.2 <- temp$param$iso.prob.other

  scenario.cInc <- median(colSums(temp$epi$se.flow) / colSums(temp$epi$num) * 100000)
  scenario.cDeaths <- median(colSums(temp$epi$d.h.flow) / colSums(temp$epi$num) * 100000)
  scenario.cIso <- median(colSums(temp$epi$iso1.flow))
  scenario.cCases <- median(colSums(temp$epi$ipic.flow) / colSums(temp$epi$num) * 100000)

  scenario.cIncWork <- median(colSums(temp$epi$se.flow.l1) / colSums(temp$epi$num.w) * 100000)
  scenario.cDeathsWork <- median(colSums(temp$epi$d.h.flow.w) / colSums(temp$epi$num.w) * 100000)
  scenario.cIsoWork <- median(colSums(temp$epi$iso1.flow.w))
  scenario.cCasesWork <- median(colSums(temp$epi$ipic.flow.w) / colSums(temp$epi$num.w) * 100000)

  WcInc3[(i - 10999), ] <- c(scale.1, scale.2, scenario.cInc)
  WcDeaths3[(i - 10999), ] <- c(scale.1, scale.2, scenario.cDeaths)
  WcIso3[(i - 10999), ] <- c(scale.1, scale.2, scenario.cIso)
  WcCases3[(i - 10999), ] <- c(scale.1, scale.2, scenario.cCases)

  WcIncWork3[(i - 10999), ] <- c(scale.1, scale.2, scenario.cIncWork)
  WcDeathsWork3[(i - 10999), ] <- c(scale.1, scale.2, scenario.cDeathsWork)
  WcIsoWork3[(i - 10999), ] <- c(scale.1, scale.2, scenario.cIsoWork)
  WcCasesWork3[(i - 10999), ] <- c(scale.1, scale.2, scenario.cCasesWork)

  remove(temp)
}

saveRDS(WcInc3, file = here::here("processed_sims","WcInc3.rds"))
saveRDS(WcDeaths3, file = here::here("processed_sims","WcDeaths3.rds"))
saveRDS(WcIso3, file = here::here("processed_sims","WcIso3.rds"))
saveRDS(WcCases3, file = here::here("processed_sims","WcCases3.rds"))
saveRDS(WcIncWork3, file = here::here("processed_sims","WcIncWork3.rds"))
saveRDS(WcDeathsWork3, file = here::here("processed_sims","WcDeathsWork3.rds"))
saveRDS(WcIsoWork3, file = here::here("processed_sims","WcIsoWork3.rds"))
saveRDS(WcCasesWork3, file = here::here("processed_sims","WcCasesWork3.rds"))
