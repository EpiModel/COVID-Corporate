# Load libraries ---------------------------------------------------------------
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(patchwork)))
suppressWarnings(suppressMessages(library(EpiModelCOVID)))

# Set baseline simulation number and directory ---------------------------------
base.sim <- 5000
my.dir <- "data/output/prelim-results"

# Function to merge sim files (same as in 04-calib-figures) --------------------
merge_simfiles <- function(simno, indir = "data/", vars = NULL,  truncate.at = NULL, verbose = TRUE) {

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

baseline_sims <- merge_simfiles(base.sim, indir = my.dir,
                                truncate.at = 181, verbose = TRUE)

# Table shells -----------------------------------------------------------------

scenarios = c("Baseline",
              "Max isolation",
              "Max isolation + work notif.",
              "Max isolation + community notif.",
              "Max isolation + home notif.",
              "Perfect isolation",
              "Perfect masking",
              "Perfect isolation + masking")

table2 <- data.frame(scenario = scenarios,
                     IsoMild = NA, IsoSevere = NA,
                     PostIsoMask = NA, ExpMask = NA,
                     IsoEnd = NA)

table3 <- data.frame(scenario = scenarios,
                     CInc = NA, IncAverted = NA, PctAverted = NA, AvertedPerDose = NA)

table4 <- data.frame(scenario = scenarios,
                     CDeaths = NA, DeathsAverted = NA, PctAverted = NA, AvertedPerDose = NA)

# Table fill functions ---------------------------------------------------------
table2_fill <- function(sim, rownum, table2) {
  iso.mild <- colSums(sim$epi$iso1.flow)
  iso.severe <- colSums(sim$epi$iso2.flow)
  iso.postmask <- colSums(sim$epi$iso3.flow)
  iso.expmask <- colSums(sim$epi$iso4.flow)
  iso.end <- colSums(sim$epi$isoend.flow)
  pds <- colSums(sim$epi$num)
  iso1_stats <- quantile(iso.mild / pds * 100000, probs = c(.25, 0.5, .75))
  iso2_stats <- quantile(iso.severe / pds * 100000, probs = c(.25, 0.5, .75))
  iso3_stats <- quantile(iso.postmask / pds * 100000, probs = c(.25, 0.5, .75))
  iso4_stats <- quantile(iso.expmask / pds * 100000, probs = c(.25, 0.5, .75))
  iso5_stats <- quantile(iso.end / pds * 100000, probs = c(.25, 0.5, .75))

  table2[rownum, 2] <- paste0(round(iso1_stats[2], 1), " (", round(iso1_stats[1], 1), ", ", round(iso1_stats[3], 1), ")")
  table2[rownum, 3] <- paste0(round(iso2_stats[2], 1), " (", round(iso2_stats[1], 1), ", ", round(iso2_stats[3], 1), ")")
  table2[rownum, 4] <- paste0(round(iso3_stats[2], 1), " (", round(iso3_stats[1], 1), ", ", round(iso3_stats[3], 1), ")")
  table2[rownum, 5] <- paste0(round(iso4_stats[2], 1), " (", round(iso4_stats[1], 1), ", ", round(iso4_stats[3], 1), ")")
  table2[rownum, 6] <- paste0(round(iso5_stats[2], 1), " (", round(iso5_stats[1], 1), ", ", round(iso5_stats[3], 1), ")")

  return(table2)
}


### EDIT AVERTED PER DOSE TO BE AVERTED PER ISOLATION
table3_fill <- function(sim, rownum, table3, ref_inf_rate = NA, ref_inf_count = NA, ref_doses = NA){
  inc <- colSums(sim$epi$se.flow)
  pd <- colSums(sim$epi$num)
  inc_rate <- inc / pd * 100000
  inc_stats <- quantile(inc_rate, probs = c(.25, 0.5, .75))
  table3[rownum, 2] <- paste0(round(inc_stats[2], 1), " (", round(inc_stats[1], 1), ", ", round(inc_stats[3], 1), ")")

  if (!is.na(ref_inf_rate)){
    inf_averted.1 <- round(ref_inf_rate - inc_stats[1], 1)
    inf_averted.2 <- round(ref_inf_rate - inc_stats[2], 1)
    inf_averted.3 <- round(ref_inf_rate - inc_stats[3], 1)

    table3[rownum, 3] <- paste0(inf_averted.2, " (", inf_averted.3, ", ", inf_averted.1, ")")

    pct_averted.1 <- round((ref_inf_rate - inc_stats[1]) / ref_inf_rate, 3) * 100
    pct_averted.2 <- round((ref_inf_rate - inc_stats[2]) / ref_inf_rate, 3) * 100
    pct_averted.3 <- round((ref_inf_rate - inc_stats[3]) / ref_inf_rate, 3) * 100

    table3[rownum, 4] <- paste0(pct_averted.2, " (", pct_averted.3, ", ", pct_averted.1, ")")

    total_doses <- colSums(sim$epi$nVax1) + colSums(sim$epi$nVax2) + colSums(sim$epi$nVax3) + colSums(sim$epi$nVax4)
    addtl_doses <- total_doses - ref_doses
    total_inf <- colSums(sim$epi$se.flow)
    total_inf_averted <- ref_inf_count - total_inf
    averted_per_dose <- total_inf_averted / addtl_doses * 1000

    averted_per_dose <- round(quantile(averted_per_dose, probs = c(.25, 0.5, .75)), 1)

    table3[rownum, 5] <- paste0(averted_per_dose[2], " (", averted_per_dose[1], ", ", averted_per_dose[3], ")")
  }

  return(table3)
}

### EDIT AVERTED PER DOSE TO BE AVERTED PER ISOLATION
table4_fill <- function(sim, rownum, table4, ref_death_rate = NA, ref_death_count = NA, ref_doses = NA){
  deaths <- colSums(sim$epi$d.h.flow)
  pds <- colSums(sim$epi$num)
  death_rates <- deaths / pds * 100000
  death_stats <- quantile(death_rates, probs = c(.25, 0.5, .75))
  table4[rownum, 2] <- paste0(round(death_stats[2], 3), " (", round(death_stats[1], 3), ", ", round(death_stats[3], 3), ")")

  if (!is.na(ref_death_rate)){
    deaths_averted.1 <- round(ref_death_rate - death_stats[1], 3)
    deaths_averted.2 <- round(ref_death_rate - death_stats[2], 3)
    deaths_averted.3 <- round(ref_death_rate - death_stats[3], 3)

    table4[rownum, 3] <- paste0(deaths_averted.2, " (", deaths_averted.3, ", ", deaths_averted.1, ")")

    pct_averted.1 <- round((ref_death_rate - death_stats[1]) / ref_death_rate, 3) * 100
    pct_averted.2 <- round((ref_death_rate - death_stats[2]) / ref_death_rate, 3) * 100
    pct_averted.3 <- round((ref_death_rate - death_stats[3]) / ref_death_rate, 3) * 100

    table4[rownum, 4] <- paste0(pct_averted.2, " (", pct_averted.3, ", ", pct_averted.1, ")")

    total_doses <- colSums(sim$epi$nVax1) + colSums(sim$epi$nVax2) + colSums(sim$epi$nVax3) + colSums(sim$epi$nVax4)
    addtl_doses <- total_doses - ref_doses
    total_deaths <- colSums(sim$epi$d.h.flow)
    total_deaths_averted <- ref_death_count - total_deaths
    averted_per_dose <- total_deaths_averted / addtl_doses * 1000

    averted_per_dose <- round(quantile(averted_per_dose, probs = c(.25, 0.5, .75)), 1)

    table4[rownum, 5] <- paste0(averted_per_dose[2], " (", averted_per_dose[1], ", ", averted_per_dose[3], ")")
  }
  return(table4)
}

# Load and merge files ---------------------------------------------------------

scenario1 <- merge_simfiles(base.sim + 1, indir = my.dir,
                            vars = NULL,  truncate.at = 181, verbose = TRUE)
scenario2 <- merge_simfiles(base.sim + 2, indir = my.dir,
                            vars = NULL,  truncate.at = 181, verbose = TRUE)
scenario3 <- merge_simfiles(base.sim + 3, indir = my.dir,
                            vars = NULL,  truncate.at = 181, verbose = TRUE)
scenario4 <- merge_simfiles(base.sim + 4, indir = my.dir,
                            vars = NULL,  truncate.at = 181, verbose = TRUE)
scenario5 <- merge_simfiles(base.sim + 5, indir = my.dir,
                            vars = NULL,  truncate.at = 181, verbose = TRUE)
scenario6 <- merge_simfiles(base.sim + 6, indir = my.dir,
                            vars = NULL,  truncate.at = 181, verbose = TRUE)
scenario7 <- merge_simfiles(base.sim + 7, indir = my.dir,
                            vars = NULL,  truncate.at = 181, verbose = TRUE)

# Fill in tables ---------------------------------------------------------------

ref_inf_rate <- quantile(colSums(baseline_sims$epi$se.flow) / colSums(baseline_sims$epi$num) * 100000, probs = 0.5)
ref_inf_count <- quantile(colSums(baseline_sims$epi$se.flow), probs = 0.5)
ref_case_rate <- quantile(colSums(baseline_sims$epi$ipic.flow) / colSums(baseline_sims$epi$num) * 100000, probs = 0.5)
ref_doses <- median(colSums(baseline_sims$epi$nVax1) + colSums(baseline_sims$epi$nVax2) + colSums(baseline_sims$epi$nVax3) + colSums(baseline_sims$epi$nVax4))
ref_death_rate <- quantile(colSums(baseline_sims$epi$d.h.flow) / colSums(baseline_sims$epi$num) * 100000, probs = 0.5)
ref_death_count <- quantile(colSums(baseline_sims$epi$d.h.flow), probs = 0.5)


# Table 2
table2 <- table2_fill(baseline_sims,1,table2)
for (i in seq(1,7)) {
  table2 <- table2_fill(eval(parse(text = paste0("scenario", i))),1+i,table2)
}

# Table 3

table3 <- table3_fill(baseline_sims, 1, table3,
                      ref_inf_rate = ref_inf_rate, ref_inf_count = ref_inf_count,
                      ref_doses = ref_doses)
for (i in seq(1,7)) {
  table3 <- table3_fill(eval(parse(text = paste0("scenario", i))),1+i,table3,
                        ref_inf_rate = ref_inf_rate, ref_inf_count = ref_inf_count,
                        ref_doses = ref_doses)
}

# Table 4
table4 <- table4_fill(baseline_sims, 1, table4,
                      ref_death_rate = ref_death_rate, ref_death_count = ref_death_count,
                      ref_doses = ref_doses)
for (i in seq(1,7)) {
  table4 <- table3_fill(eval(parse(text = paste0("scenario", i))),1+i,table4,
                        ref_inf_rate = ref_inf_rate, ref_inf_count = ref_inf_count,
                        ref_doses = ref_doses)
}

