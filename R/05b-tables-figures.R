# Load libraries ---------------------------------------------------------------
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(patchwork)))
suppressWarnings(suppressMessages(library(EpiModelCOVID)))

# Set baseline simulation number and directory ---------------------------------
base.sim <- 5000
my.dir <- "data/output/results"
s.num <- 74 # number of scenarios, excluding reference scenario

# Tables: FOR EXISTING OUTPUT Load files ---------------------------------------

scenario_list <- readRDS(here::here("data","output","results",
                                "scenario_list.rds"))

baseline_sims <- readRDS(here::here("data","output","results",
                                    "baseline_sims.rds"))

# scenario0 <- baseline_sims

# Table shells -----------------------------------------------------------------

## Need to edit
# scenarios = c("Baseline",
#               "Max isolation",
#               "Max work notification",
#               "Max community notification",
#               "Max home notification",
#               "Max isolation + work notif.",
#               "Max isolation + community notif.",
#               "Max isolation + home notif.",
#               "Perfect isolation",
#               "Perfect masking",
#               "Perfect isolation + masking")
scenarios = seq(0,s.num)

table2 <- data.frame(scenario = scenarios,
                     IsoMild = NA, IsoSevere = NA,
                     PostIsoMask = NA, ExpMask = NA,
                     IsoEnd = NA,
                     IsoMildWork = NA, IsoSevereWork = NA,
                     PostIsoMaskWork = NA, ExpMaskWork = NA,
                     IsoEndWork = NA)

table3 <- data.frame(scenario = scenarios,
                     CInc = NA, IncAverted = NA, PctAverted = NA,
                     AvertedPerIso = NA,
                     CIncWork = NA, IncAvertedWork = NA, PctAvertedWork = NA,
                     AvertedPerIsoWork = NA)

table4 <- data.frame(scenario = scenarios,
                     CDeaths = NA, DeathsAverted = NA, PctAverted = NA,
                     AvertedPerIso = NA,
                     CDeathsWork = NA, DeathsAvertedWork = NA, PctAvertedWork = NA,
                     AvertedPerIsoWork = NA)

table5 <- data.frame(scenario = scenarios,
                     CCase = NA, CaseAverted = NA, PctAverted = NA,
                     AvertedPerIso = NA,
                     CCaseWork = NA, CaseAvertedWork = NA, PctAvertedWork = NA,
                     AvertedPerIsoWork = NA)

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

  iso.mild.w <- colSums(sim$epi$iso1.flow.w)
  iso.severe.w <- colSums(sim$epi$iso2.flow.w)
  iso.postmask.w <- colSums(sim$epi$iso3.flow.w)
  iso.expmask.w <- colSums(sim$epi$iso4.flow.w)
  iso.end.w <- colSums(sim$epi$isoend.flow.w)
  pds.w <- colSums(sim$epi$num.w)
  iso1_stats.w <- quantile(iso.mild.w / pds.w * 100000, probs = c(.25, 0.5, .75))
  iso2_stats.w <- quantile(iso.severe.w / pds.w * 100000, probs = c(.25, 0.5, .75))
  iso3_stats.w <- quantile(iso.postmask.w / pds.w * 100000, probs = c(.25, 0.5, .75))
  iso4_stats.w <- quantile(iso.expmask.w / pds.w * 100000, probs = c(.25, 0.5, .75))
  iso5_stats.w <- quantile(iso.end.w / pds.w * 100000, probs = c(.25, 0.5, .75))

  table2[rownum, 2] <- paste0(round(iso1_stats[2], 1), " (", round(iso1_stats[1], 1), ", ", round(iso1_stats[3], 1), ")")
  table2[rownum, 3] <- paste0(round(iso2_stats[2], 1), " (", round(iso2_stats[1], 1), ", ", round(iso2_stats[3], 1), ")")
  table2[rownum, 4] <- paste0(round(iso3_stats[2], 1), " (", round(iso3_stats[1], 1), ", ", round(iso3_stats[3], 1), ")")
  table2[rownum, 5] <- paste0(round(iso4_stats[2], 1), " (", round(iso4_stats[1], 1), ", ", round(iso4_stats[3], 1), ")")
  table2[rownum, 6] <- paste0(round(iso5_stats[2], 1), " (", round(iso5_stats[1], 1), ", ", round(iso5_stats[3], 1), ")")
  table2[rownum, 7] <- paste0(round(iso1_stats.w[2], 1), " (", round(iso1_stats.w[1], 1), ", ", round(iso1_stats.w[3], 1), ")")
  table2[rownum, 8] <- paste0(round(iso2_stats.w[2], 1), " (", round(iso2_stats.w[1], 1), ", ", round(iso2_stats.w[3], 1), ")")
  table2[rownum, 9] <- paste0(round(iso3_stats.w[2], 1), " (", round(iso3_stats.w[1], 1), ", ", round(iso3_stats.w[3], 1), ")")
  table2[rownum, 10] <- paste0(round(iso4_stats.w[2], 1), " (", round(iso4_stats.w[1], 1), ", ", round(iso4_stats.w[3], 1), ")")
  table2[rownum, 11] <- paste0(round(iso5_stats.w[2], 1), " (", round(iso5_stats.w[1], 1), ", ", round(iso5_stats.w[3], 1), ")")

  return(table2)
}


table3_fill <- function(sim, rownum, table3,
                        ref_inf_rate = NA, ref_inf_rate_w = NA,
                        ref_inf_count = NA, ref_inf_count_w = NA,
                        ref_iso = NA, ref_iso_w = NA){
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

    total_iso <- colSums(sim$epi$iso1.flow)
    addtl_iso <- total_iso - ref_iso
    total_inf <- colSums(sim$epi$se.flow)
    total_inf_averted <- ref_inf_count - total_inf
    averted_per_iso <- total_inf_averted / addtl_iso

    averted_per_iso <- round(quantile(averted_per_iso, probs = c(.25, 0.5, .75)), 1)

    table3[rownum, 5] <- paste0(averted_per_iso[2], " (", averted_per_iso[1], ", ", averted_per_iso[3], ")")
  }

  inc.w <- colSums(sim$epi$se.flow.l1)
  pd.w <- colSums(sim$epi$num.w)
  inc_rate.w <- inc.w / pd.w * 100000
  inc_stats.w <- quantile(inc_rate.w, probs = c(.25, 0.5, .75))
  table3[rownum, 6] <- paste0(round(inc_stats.w[2], 1), " (", round(inc_stats.w[1], 1), ", ", round(inc_stats.w[3], 1), ")")

  if (!is.na(ref_inf_rate_w)){
    inf_averted.1.w <- round(ref_inf_rate_w - inc_stats.w[1], 1)
    inf_averted.2.w <- round(ref_inf_rate_w - inc_stats.w[2], 1)
    inf_averted.3.w <- round(ref_inf_rate_w - inc_stats.w[3], 1)

    table3[rownum, 7] <- paste0(inf_averted.2.w, " (", inf_averted.3.w, ", ", inf_averted.1.w, ")")

    pct_averted.1.w <- round((ref_inf_rate_w - inc_stats.w[1]) / ref_inf_rate_w, 3) * 100
    pct_averted.2.w <- round((ref_inf_rate_w - inc_stats.w[2]) / ref_inf_rate_w, 3) * 100
    pct_averted.3.w <- round((ref_inf_rate_w - inc_stats.w[3]) / ref_inf_rate_w, 3) * 100

    table3[rownum, 8] <- paste0(pct_averted.2.w, " (", pct_averted.3.w, ", ", pct_averted.1.w, ")")

    total_iso.w <- colSums(sim$epi$iso1.flow.w)
    addtl_iso.w <- total_iso.w - ref_iso_w
    total_inf.w <- colSums(sim$epi$se.flow.l1)
    total_inf_averted.w <- ref_inf_count_w - total_inf.w
    averted_per_iso.w <- total_inf_averted.w / addtl_iso.w

    averted_per_iso.w <- round(quantile(averted_per_iso.w, probs = c(.25, 0.5, .75)), 1)

    table3[rownum, 9] <- paste0(averted_per_iso.w[2], " (", averted_per_iso.w[1], ", ", averted_per_iso.w[3], ")")
  }

  return(table3)
}

table4_fill <- function(sim, rownum, table4,
                        ref_death_rate = NA,
                        # ref_death_rate_w = NA,
                        ref_death_count = NA,
                        # ref_death_count_w = NA,
                        ref_iso = NA
                        # ref_iso_w = NA
                        ){

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

    total_iso <- colSums(sim$epi$iso1.flow)
    addtl_iso <- total_iso - ref_iso
    total_deaths <- colSums(sim$epi$d.h.flow)
    total_deaths_averted <- ref_death_count - total_deaths
    averted_per_iso <- total_deaths_averted / addtl_iso

    averted_per_iso <- round(quantile(averted_per_iso, probs = c(.25, 0.5, .75)), 1)

    table4[rownum, 5] <- paste0(averted_per_iso[2], " (", averted_per_iso[1], ", ", averted_per_iso[3], ")")
  }

  # deaths.w <- colSums(sim$epi$d.h.flow.w)
  # pds.w <- colSums(sim$epi$num.w)
  # death_rates.w <- deaths.w / pds.w * 100000
  # death_stats.w <- quantile(death_rates.w, probs = c(.25, 0.5, .75))
  # table4[rownum, 2] <- paste0(round(death_stats.w[2], 3), " (", round(death_stats.w[1], 3), ", ", round(death_stats.w[3], 3), ")")
  #
  # if (!is.na(ref_death_rate_w_w)){
  #   deaths_averted.1.w <- round(ref_death_rate_w - death_stats.w[1], 3)
  #   deaths_averted.2.w <- round(ref_death_rate_w - death_stats.w[2], 3)
  #   deaths_averted.3.w <- round(ref_death_rate_w - death_stats.w[3], 3)
  #
  #   table4[rownum, 3] <- paste0(deaths_averted.2.w, " (", deaths_averted.3.w, ", ", deaths_averted.1.w, ")")
  #
  #   pct_averted.1.w <- round((ref_death_rate_w - death_stats.w[1]) / ref_death_rate_w, 3) * 100
  #   pct_averted.2.w <- round((ref_death_rate_w - death_stats.w[2]) / ref_death_rate_w, 3) * 100
  #   pct_averted.3.w <- round((ref_death_rate_w - death_stats.w[3]) / ref_death_rate_w, 3) * 100
  #
  #   table4[rownum, 4] <- paste0(pct_averted.2.w, " (", pct_averted.3.w, ", ", pct_averted.1.w, ")")
  #
  #   total_iso.w <- colSums(sim$epi$iso1.flow.w)
  #   addtl_iso.w <- total_iso.w - ref_iso_w
  #   total_deaths.w <- colSums(sim$epi$d.h.flow.w)
  #   total_deaths_averted.w <- ref_death_count_w - total_deaths.w
  #   averted_per_iso.w <- total_deaths_averted.w / addtl_iso.w
  #
  #   averted_per_iso.w <- round(quantile(averted_per_iso.w, probs = c(.25, 0.5, .75)), 1)
  #
  #   table4[rownum, 5] <- paste0(averted_per_iso.w[2], " (", averted_per_iso.w[1], ", ", averted_per_iso.w[3], ")")
  # }
  #
  return(table4)
}

table5_fill <- function(sim, rownum, table5,
                        ref_case_rate = NA, ref_case_rate_w = NA,
                        ref_case_count = NA, ref_case_count_w = NA,
                        ref_iso = NA, ref_iso_w = NA){
  inc <- colSums(sim$epi$ipic.flow)
  pd <- colSums(sim$epi$num)
  cinc_rate <- inc / pd * 100000
  cinc_stats <- quantile(cinc_rate, probs = c(.25, 0.5, .75))
  table5[rownum, 2] <- paste0(round(cinc_stats[2], 1), " (", round(cinc_stats[1], 1), ", ", round(cinc_stats[3], 1), ")")

  if (!is.na(ref_case_rate)){
    case_averted.1 <- round(ref_case_rate - cinc_stats[1], 1)
    case_averted.2 <- round(ref_case_rate - cinc_stats[2], 1)
    case_averted.3 <- round(ref_case_rate - cinc_stats[3], 1)

    table5[rownum, 3] <- paste0(case_averted.2, " (", case_averted.3, ", ", case_averted.1, ")")

    pct_averted.1 <- round((ref_case_rate - cinc_stats[1]) / ref_case_rate, 3) * 100
    pct_averted.2 <- round((ref_case_rate - cinc_stats[2]) / ref_case_rate, 3) * 100
    pct_averted.3 <- round((ref_case_rate - cinc_stats[3]) / ref_case_rate, 3) * 100

    table5[rownum, 4] <- paste0(pct_averted.2, " (", pct_averted.3, ", ", pct_averted.1, ")")

    total_iso <- colSums(sim$epi$iso1.flow)
    addtl_iso <- total_iso - ref_iso
    total_inf <- colSums(sim$epi$ipic.flow)
    total_case_averted <- ref_case_count - total_inf
    averted_per_iso <- total_case_averted / addtl_iso

    averted_per_iso <- round(quantile(averted_per_iso, probs = c(.25, 0.5, .75)), 1)

    table5[rownum, 5] <- paste0(averted_per_iso[2], " (", averted_per_iso[1], ", ", averted_per_iso[3], ")")
  }

  inc.w <- colSums(sim$epi$ipic.flow.w)
  pd.w <- colSums(sim$epi$num.w)
  cinc_rate.w <- inc.w / pd.w * 100000
  cinc_stats.w <- quantile(cinc_rate.w, probs = c(.25, 0.5, .75))
  table5[rownum, 6] <- paste0(round(cinc_stats.w[2], 1), " (", round(cinc_stats.w[1], 1), ", ", round(cinc_stats.w[3], 1), ")")

  if (!is.na(ref_case_rate_w)){
    case_averted.1.w <- round(ref_case_rate_w - cinc_stats.w[1], 1)
    case_averted.2.w <- round(ref_case_rate_w - cinc_stats.w[2], 1)
    case_averted.3.w <- round(ref_case_rate_w - cinc_stats.w[3], 1)

    table5[rownum, 7] <- paste0(case_averted.2.w, " (", case_averted.3.w, ", ", case_averted.1.w, ")")

    pct_averted.1.w <- round((ref_case_rate_w - cinc_stats.w[1]) / ref_case_rate_w, 3) * 100
    pct_averted.2.w <- round((ref_case_rate_w - cinc_stats.w[2]) / ref_case_rate_w, 3) * 100
    pct_averted.3.w <- round((ref_case_rate_w - cinc_stats.w[3]) / ref_case_rate_w, 3) * 100

    table5[rownum, 8] <- paste0(pct_averted.2.w, " (", pct_averted.3.w, ", ", pct_averted.1.w, ")")

    total_iso.w <- colSums(sim$epi$iso1.flow.w)
    addtl_iso.w <- total_iso.w - ref_iso_w
    total_inf.w <- colSums(sim$epi$ipic.flow.w)
    total_case_averted.w <- ref_case_count_w - total_inf.w
    averted_per_iso.w <- total_case_averted.w / addtl_iso.w

    averted_per_iso.w <- round(quantile(averted_per_iso.w, probs = c(.25, 0.5, .75)), 1)

    table5[rownum, 9] <- paste0(averted_per_iso.w[2], " (", averted_per_iso.w[1], ", ", averted_per_iso.w[3], ")")
  }

  return(table5)
}

# Fill in tables ---------------------------------------------------------------

ref_inf_rate <- quantile(colSums(baseline_sims$epi$se.flow) / colSums(baseline_sims$epi$num) * 100000, probs = 0.5)
ref_inf_count <- quantile(colSums(baseline_sims$epi$se.flow), probs = 0.5)
ref_case_rate <- quantile(colSums(baseline_sims$epi$ipic.flow) / colSums(baseline_sims$epi$num) * 100000, probs = 0.5)
ref_case_count <- quantile(colSums(baseline_sims$epi$ipic.flow), probs = 0.5)
ref_death_rate <- quantile(colSums(baseline_sims$epi$d.h.flow) / colSums(baseline_sims$epi$num) * 100000, probs = 0.5)
ref_death_count <- quantile(colSums(baseline_sims$epi$d.h.flow), probs = 0.5)
ref_iso <- median(colSums(baseline_sims$epi$iso1.flow))

# reference values among office workers
ref_inf_rate_w <- quantile(colSums(baseline_sims$epi$se.flow.l1) / colSums(baseline_sims$epi$num.w) * 100000, probs = 0.5)
ref_inf_count_w <- quantile(colSums(baseline_sims$epi$se.flow.l1), probs = 0.5)
ref_case_rate_w <- quantile(colSums(baseline_sims$epi$ipic.flow.w) / colSums(baseline_sims$epi$num.w) * 100000, probs = 0.5)
ref_case_count_w <- quantile(colSums(baseline_sims$epi$ipic.flow.w), probs = 0.5)
## would need to add these to death module
# ref_death_rate_w <- quantile(colSums(baseline_sims$epi$d.h.flow.w) / colSums(baseline_sims$epi$num.w) * 100000, probs = 0.5)
# ref_death_count_w <- quantile(colSums(baseline_sims$epi$d.h.flow.w), probs = 0.5)
ref_iso_w <- median(colSums(baseline_sims$epi$iso1.flow.w))

# Table 2
table2 <- table2_fill(baseline_sims,1,table2)
for (i in seq(1,s.num)) {
  table2 <- table2_fill(scenario_list[[i]],1+i,table2)
}

# Table 3
table3 <- table3_fill(baseline_sims, 1, table3,
                      ref_inf_rate = ref_inf_rate,
                      ref_inf_rate_w = ref_inf_rate_w,
                      ref_inf_count = ref_inf_count,
                      ref_inf_count_w = ref_inf_count_w,
                      ref_iso = ref_iso,
                      ref_iso_w = ref_iso_w)
for (i in seq(1,s.num)) {
  table3 <- table3_fill(scenario_list[[i]],1+i,table3,
                        ref_inf_rate = ref_inf_rate,
                        ref_inf_rate_w = ref_inf_rate_w,
                        ref_inf_count = ref_inf_count,
                        ref_inf_count_w = ref_inf_count_w,
                        ref_iso = ref_iso,
                        ref_iso_w = ref_iso_w)
}

# Table 4
table4 <- table4_fill(baseline_sims, 1, table4,
                      ref_death_rate = ref_death_rate,
                      # ref_death_rate_w = ref_death_rate_w,
                      ref_death_count = ref_death_count,
                      # ref_death_count_w = ref_death_count_w,
                      ref_iso = ref_iso
                      # ,
                      # ref_iso_w = ref_iso_w
                      )
for (i in seq(1,s.num)) {
  table4 <- table4_fill(scenario_list[[i]],1+i,table4,
                        ref_death_rate = ref_death_rate,
                        # ref_death_rate_w = ref_death_rate_w,
                        ref_death_count = ref_death_count,
                        # ref_death_count_w = ref_death_count_w,
                        ref_iso = ref_iso
                        # ,
                        # ref_iso_w = ref_iso_w
                        )
}

# Table 5
table5 <- table5_fill(baseline_sims, 1, table5,
                      ref_case_rate = ref_case_rate,
                      ref_case_rate_w = ref_case_rate_w,
                      ref_case_count = ref_case_count,
                      ref_case_count_w = ref_case_count_w,
                      ref_iso = ref_iso,
                      ref_iso_w = ref_iso_w)
for (i in seq(1,s.num)) {
  table5 <- table5_fill(scenario_list[[i]],1+i,table5,
                        ref_case_rate = ref_case_rate,
                        ref_case_rate_w = ref_case_rate_w,
                        ref_case_count = ref_case_count,
                        ref_case_count_w = ref_case_count_w,
                        ref_iso = ref_iso,
                        ref_iso_w = ref_iso_w)
}

library(openxlsx)

# Create an Excel workbook
wb <- createWorkbook()

# Add data frames to different sheets
addWorksheet(wb, "Table2")
writeData(wb, sheet = 1, x = table2)

addWorksheet(wb, "Table3")
writeData(wb, sheet = 2, x = table3)

# addWorksheet(wb, "Table4")
# writeData(wb, sheet = 3, x = table4)

addWorksheet(wb, "Table5")
writeData(wb, sheet = 3, x = table5)

# Save the workbook to a file
saveWorkbook(wb, here::here("results","prelim-tables.xlsx"))

# Ratio Plot: Cases to Isolated ------------------------------------------------

ratio.plot <- function(sim,i) {
  sim <- mutate_epi(sim,
                    totCase = ic.num.w + h.num.w,
                    totIso = iso.mild.num.w + iso.sev.num.w,
                    totMask = mask.post.num.w + mask.exp.num.w)

  sim <- mutate_epi(sim,
                    CaseIsoRatio = totCase/totIso,
                    CaseIsoSum = totCase+totIso)

  df <- as.data.frame(sim) %>%
    ungroup() %>%
    group_by(time) %>%
    summarise(med.case = median(totCase),
           p25.case = quantile(totCase, probs = 0.25),
           p75.case = quantile(totCase, probs = 0.75),
           med.iso = median(totIso),
           p25.iso = quantile(totIso, probs = 0.25),
           p75.iso = quantile(totIso, probs = 0.75),
           med.mask = median(totMask),
           p25.mask = quantile(totMask, probs = 0.25),
           p75.mask = quantile(totMask, probs = 0.75),
           med.ratio = median(CaseIsoRatio),
           p25.ratio = quantile(CaseIsoRatio, probs = 0.25),
           p75.ratio = quantile(CaseIsoRatio, probs = 0.75),
           med.sum = median(CaseIsoSum),
           p25.sum = quantile(CaseIsoSum, probs = 0.25),
           p75.sum = quantile(CaseIsoSum, probs = 0.75)) %>%
    mutate(scenario = i)

  return(df)
}

d.ratio <- data.frame(time = NA,
                      med.case = NA, p25.case = NA, p75.case = NA,
                      med.iso = NA, p25.iso = NA, p75.iso = NA,
                      med.mask = NA, p25.mask = NA, p75.mask = NA,
                      med.ratio = NA, p25.ratio = NA, p75.ratio = NA,
                      med.sum = NA, p25.sum = NA, p75.sum = NA,
                      scenario = i)

for (i in seq(0,s.num,1)) {

  temp <- ratio.plot(scenario_list[[i]],i)

  d.ratio <- d.ratio %>%
    add_row(temp) %>%
    filter(!is.na(time))

  rm(temp)

}

d.ratio <- d.ratio %>%
  mutate(scenario2 = factor(scenario,
                            levels = seq(0,s.num),
                            labels = c("Baseline",
                                       "Max isolation",
                                       "Max work notification",
                                       "Max community notification",
                                       "Max home notification",
                                       "Max isolation + work notif.",
                                       "Max isolation + community notif.",
                                       "Max isolation + home notif.",
                                       "Perfect isolation",
                                       "Perfect masking",
                                       "Perfect isolation + masking")))

p.ratio <- ggplot(d.ratio, aes(x = time, y = med.ratio))+
  geom_line(aes(y = med.ratio, color = scenario2, group = scenario2)) +
  # scale_x_continuous(expand = c(0,0),
  #                    limits = c(1, 20),
  #                    breaks = c(1, 4, 7, 10, 13, 16, 19),
  #                    labels = c("1" = "Jan 21", "4" = "Apr 21", "7" = "Jul 21",
  #                               "10" = "Oct 21", "13" = "Jan 22", "16" = "Apr 22",
  #                               "19" = "Jul 22")) +
  # scale_y_continuous(expand = c(0, 0), lim = c(0, 26500)) +
  labs(x = "Time",
       y = "Number cases/Number in isolation",
       color = "Scenario") +
  guides(color = "none")+
  ggtitle("Ratio of clinical cases to persons in isolation") + theme_classic()

p.sum <- ggplot(d.ratio, aes(x = time, y = med.sum))+
  geom_line(aes(y = med.sum, color = scenario2, group = scenario2)) +
  labs(x = "Time",
       y = "Sum of cases and persons in isolation",
       color = "Scenario") +
  guides(color = "none")+
  ggtitle("Total number of people ill or in isolation") + theme_classic()

p.cases <- ggplot(d.ratio, aes(x = time, y = med.case))+
  geom_line(aes(y = med.case, color = scenario2, group = scenario2)) +
  labs(x = "Time",
       y = "Number of cases",
       color = "Scenario") +
  guides(color = "none")+
  ggtitle("Total number of people ill") + theme_classic()

p.iso <- ggplot(d.ratio, aes(x = time, y = med.iso))+
  geom_line(aes(y = med.iso, color = scenario2, group = scenario2)) +
  labs(x = "Time",
       y = "Number of persons in isolation",
       color = "Scenario") +
  guides(color = "none")+
  ggtitle("Total number of people in isolation") + theme_classic()

p.mask <- ggplot(d.ratio, aes(x = time, y = med.mask))+
  geom_line(aes(y = med.mask, color = scenario2, group = scenario2)) +
  labs(x = "Time",
       y = "Number of persons masked",
       color = "Scenario") +
  ggtitle("Total number of people masked") + theme_classic()

(p.ratio | p.sum) / (p.cases | p.iso) / p.mask

ggsave(here::here("results","figures","Figure4.tiff"),
       width = 20, height = 12, units = "in", dpi=300, compression = 'lzw')

################################################################################
# Load contour plot data -------------------------------------------------------

cInc <- readRDS(here::here("data","output","results","cInc.rds"))
cDeaths <- readRDS(here::here("data","output","results","cDeaths.rds"))
cIso <- readRDS(here::here("data","output","results","cIso.rds"))
cCases <- readRDS(here::here("data","output","results","cCases.rds"))
cIncWork <- readRDS(here::here("data","output","results","cIncWork.rds"))
# cDeathsWork <- readRDS(here::here("data","output","results","cDeathsWork.rds"))
cIsoWork <- readRDS(here::here("data","output","results","cIsoWork.rds"))
cCasesWork <- readRDS(here::here("data","output","results","cCasesWork.rds"))

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

# Plot cases - Office Workers
minval.iw <- floor(min(cIncWork$cInc))
maxval.iw <- ceiling(max(cIncWork$cInc))
g_inc_cont_w <- ggplot(cIncWork, aes(iso.prob, work.notif.prob)) +
  geom_raster(aes(fill = cInc), interpolate = TRUE) +
  geom_contour(aes(z = cInc), col = "white", alpha = 0.5, lwd = 0.5) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Isolation Probability", y = "Notification Probability at Work") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.iw, maxval.iw), guide = "none") +
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

# # Plot deaths - Office Workers
# minval.dw <- min(cDeathsWork$cDeathsWork)
# maxval.dw <- max(cDeathsWork$cDeathsWork)
# g_death_cont_w <- ggplot(cDeathsWork, aes(iso.prob, work.notif.prob)) +
#   geom_raster(aes(fill = cDeathsWork), interpolate = TRUE) +
#   geom_contour(aes(z = cDeathsWork), col = "white", alpha = 0.5, lwd = 0.5) +
#   #geom_text_contour(aes(z = cDeathsWork), stroke = 0.1) +
#   theme_minimal() +
#   scale_y_continuous(expand = c(0, 0)) +
#   scale_x_continuous(expand = c(0, 0)) +
#   labs(x = "Isolation Probability", y = "Notification Probability at Work") +
#   scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
#                      name = "", limits = c(minval.dw, maxval.dw), guide = "none") +
#   theme(axis.text = element_text(size = rel(1.2)),
#         axis.title = element_text(size = rel(1.2)),
#         legend.text = element_text(size = rel(1.2))) +
#   coord_fixed()

# Plot isolation
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

# Plot isolation - Office Workers
minval.vw <- min(cIsoWork$cIso)
maxval.vw <- max(cIsoWork$cIso)
g_iso_cont_w <- ggplot(cIsoWork, aes(iso.prob, work.notif.prob)) +
  geom_raster(aes(fill = cIso), interpolate = TRUE) +
  geom_contour(aes(z = cIso), col = "white", alpha = 0.5, lwd = 0.5) +
  #geom_text_contour(aes(z = cIsoWork), stroke = 0.1) +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Isolation Probability", y = "Notification Probability at Work") +
  scale_fill_viridis(discrete = FALSE, alpha = 1, option = "D", direction = 1,
                     name = "", limits = c(minval.vw, maxval.vw), guide = "none") +
  theme(axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2))) +
  coord_fixed()

combined1 <- g_iso_cont + g_iso_cont_w +
  plot_annotation(tag_levels = list(c("A1", "A2")), title = "Total Persons Entering Isolation")
combined2 <- g_inc_cont + g_inc_cont_w +
  plot_annotation(tag_levels = list(c("B1", "B2")), title = "Infections per 100'000 Person-Days")
combined3 <- g_death_cont +
  # g_death_cont_w +
  plot_annotation(tag_levels = list(c("C")), title = "Deaths per 100'000 Person-Days")

wrap_elements(combined1) / wrap_elements(combined2) / wrap_elements(combined3)

ggsave(here::here("results","figures","Figure3.tiff"),
       width = 12, height = 12, units = "in", dpi=300, compression = 'lzw')
