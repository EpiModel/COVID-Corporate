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

# Tables: Load and merge files -------------------------------------------------

baseline_sims <- merge_simfiles(base.sim, indir = my.dir,
                                truncate.at = 181, verbose = TRUE)
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
                     CInc = NA, IncAverted = NA, PctAverted = NA, AvertedPerIso = NA)

table4 <- data.frame(scenario = scenarios,
                     CDeaths = NA, DeathsAverted = NA, PctAverted = NA, AvertedPerIso = NA)

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


table3_fill <- function(sim, rownum, table3, ref_inf_rate = NA, ref_inf_count = NA, ref_iso = NA){
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

    total_iso <- colSums(sim$epi$iso1.flow) + colSums(sim$epi$iso4.flow)
    addtl_iso <- total_iso - ref_iso
    total_inf <- colSums(sim$epi$se.flow)
    total_inf_averted <- ref_inf_count - total_inf
    averted_per_iso <- total_inf_averted / addtl_iso

    averted_per_iso <- round(quantile(averted_per_iso, probs = c(.25, 0.5, .75)), 1)

    table3[rownum, 5] <- paste0(averted_per_iso[2], " (", averted_per_iso[1], ", ", averted_per_iso[3], ")")
  }

  return(table3)
}

table4_fill <- function(sim, rownum, table4, ref_death_rate = NA, ref_death_count = NA, ref_iso = NA){
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

    total_iso <- colSums(sim$epi$iso1.flow) + colSums(sim$epi$iso4.flow)
    addtl_iso <- total_iso - ref_iso
    total_deaths <- colSums(sim$epi$d.h.flow)
    total_deaths_averted <- ref_death_count - total_deaths
    averted_per_iso <- total_deaths_averted / addtl_iso

    averted_per_iso <- round(quantile(averted_per_iso, probs = c(.25, 0.5, .75)), 1)

    table4[rownum, 5] <- paste0(averted_per_iso[2], " (", averted_per_iso[1], ", ", averted_per_iso[3], ")")
  }
  return(table4)
}

# Fill in tables ---------------------------------------------------------------

ref_inf_rate <- quantile(colSums(baseline_sims$epi$se.flow) / colSums(baseline_sims$epi$num) * 100000, probs = 0.5)
ref_inf_count <- quantile(colSums(baseline_sims$epi$se.flow), probs = 0.5)
ref_case_rate <- quantile(colSums(baseline_sims$epi$ipic.flow) / colSums(baseline_sims$epi$num) * 100000, probs = 0.5)
ref_death_rate <- quantile(colSums(baseline_sims$epi$d.h.flow) / colSums(baseline_sims$epi$num) * 100000, probs = 0.5)
ref_death_count <- quantile(colSums(baseline_sims$epi$d.h.flow), probs = 0.5)
ref_iso <- median(colSums(baseline_sims$epi$iso1.flow) + colSums(baseline_sims$epi$iso4.flow))

# Table 2
table2 <- table2_fill(baseline_sims,1,table2)
for (i in seq(1,7)) {
  table2 <- table2_fill(eval(parse(text = paste0("scenario", i))),1+i,table2)
}

# Table 3

table3 <- table3_fill(baseline_sims, 1, table3,
                      ref_inf_rate = ref_inf_rate, ref_inf_count = ref_inf_count,
                      ref_iso = ref_iso)
for (i in seq(1,7)) {
  table3 <- table3_fill(eval(parse(text = paste0("scenario", i))),1+i,table3,
                        ref_inf_rate = ref_inf_rate, ref_inf_count = ref_inf_count,
                        ref_iso = ref_iso)
}

# Table 4
table4 <- table4_fill(baseline_sims, 1, table4,
                      ref_death_rate = ref_death_rate, ref_death_count = ref_death_count,
                      ref_iso = ref_iso)
for (i in seq(1,7)) {
  table4 <- table3_fill(eval(parse(text = paste0("scenario", i))),1+i,table4,
                        ref_inf_rate = ref_inf_rate, ref_inf_count = ref_inf_count,
                        ref_iso = ref_iso)
}

# Contour Plots: Load and merge files ------------------------------------------

##Isolation prob vs work notification prob
cInc <- data.frame(iso.prob = rep(NA, 100), work.notif.prob = rep(NA, 100), cInc = rep(NA, 100))
cDeaths <- data.frame(iso.prob = rep(NA, 100), work.notif.prob = rep(NA, 100), cDeaths = rep(NA, 100))
cIso <- data.frame(iso.prob = rep(NA, 100), work.notif.prob = rep(NA, 100), cIso = rep(NA, 100))
cCases <- data.frame(iso.prob = rep(NA, 100), work.notif.prob = rep(NA, 100), cCases = rep(NA, 100))

for (i in 6000:6099) {

  print(i)

  temp <- merge_simfiles(i, indir = my.dir, vars = NULL,  truncate.at = 181, verbose = TRUE)

  scale.1 <- temp$param$iso.prob
  scale.2 <- temp$param$notif.prob[1]

  scenario.cInc <- median(colSums(temp$epi$se.flow) / colSums(temp$epi$num) * 100000)
  scenario.cDeaths <- median(colSums(temp$epi$d.h.flow) / colSums(temp$epi$num) * 100000)
  scenario.cIso <- median(colSums(baseline_sims$epi$iso1.flow) + colSums(baseline_sims$epi$iso4.flow))
  scenario.cCases <- median(colSums(temp$epi$ipic.flow) / colSums(temp$epi$num) * 100000)

  cInc[(i - 5999), ] <- c(scale.1, scale.2, scenario.cInc)
  cDeaths[(i - 5999), ] <- c(scale.1, scale.2, scenario.cDeaths)
  cIso[(i - 5999), ] <- c(scale.1, scale.2, scenario.cIso)
  cCases[(i - 5999), ] <- c(scale.1, scale.2, scenario.cCases)

  remove(temp)
}

# Contour Plots: Plots ---------------------------------------------------------

library("viridis")
library("ggplot2")
library("metR")
library("patchwork")

# Plot cases
minval.i <- floor(min(cInc$cInc))
maxval.i <- ceiling(max(cInc$cInc))
g_inc_cont <- ggplot(cInc, aes(iso.prob, work.notif.prob)) +
  geom_raster(aes(fill = cInc), interpolate = TRUE) +
  geom_contour(aes(z = cInc), col = "white", alpha = 0.5, lwd = 0.5) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Isolation Probability", y = "Notification Probability at Work") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.i, maxval.i), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot deaths
minval.d <- min(cDeaths$cDeaths)
maxval.d <- max(cDeaths$cDeaths)
g_death_cont <- ggplot(cDeaths, aes(iso.prob, work.notif.prob)) +
  geom_raster(aes(fill = cDeaths), interpolate = TRUE) +
  geom_contour(aes(z = cDeaths), col = "white", alpha = 0.5, lwd = 0.5) +
  #geom_text_contour(aes(z = cDeaths), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Isolation Probability", y = "Notification Probability at Work") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.d, maxval.d), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

# Plot doses
minval.v <- min(cIso$cIso)
maxval.v <- max(cIso$cIso)
g_iso_cont <- ggplot(cIso, aes(iso.prob, work.notif.prob)) +
  geom_raster(aes(fill = cIso), interpolate = TRUE) +
  geom_contour(aes(z = cIso), col = "white", alpha = 0.5, lwd = 0.5) +
  #geom_text_contour(aes(z = cIso), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Isolation Probability", y = "Notification Probability at Work") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.v, maxval.v), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

combined1 <- g_iso_cont + plot_annotation(tag_levels = list(c("A")), title = "Total Persons Entering Isolation")
combined2 <- g_inc_cont + plot_annotation(tag_levels = list(c("B")), title = "Infections per 100'000 Person-Days")
combined3 <- g_death_cont + plot_annotation(tag_levels = list(c("C")), title = "Deaths per 100'000 Person-Days")

wrap_elements(combined1) / wrap_elements(combined2) / wrap_elements(combined3)

ggsave(here::here("results","figures","Figure3.tiff"),
       width = 12, height = 12, units = "in", dpi=300, compression = 'lzw')
