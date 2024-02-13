# Load libraries ---------------------------------------------------------------
# pacman::p_load(ggplot2,
#                tidyverse,
#                tidyr,
#                dplyr,
#                patchwork,
#                EpiModelCOVID,
#                openxlsx)

suppressMessages({
  library("tidyverse")
  library("patchwork")
  library("tidyr")
  library("dplyr")
  library("ggplot2")
  library("EpiModelCOVID")
  library("openxlsx")
})

# Set baseline simulation number and directory ---------------------------------
base.sim <- 5000
my.dir <- "sims/"
s.num <- 69 # number of scenarios

# Tables: FOR EXISTING OUTPUT Load files ---------------------------------------

scenario_list <- readRDS(here::here("processed_sims",
                                    "scenario_list.rds"))

baseline_sims <- readRDS(here::here("processed_sims",
                                    "baseline_sim.rds"))

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
                     AvertedPerIso = NA, NNT = NA,
                     CIncWork = NA, IncAvertedWork = NA, PctAvertedWork = NA,
                     AvertedPerIsoWork = NA, NNTWork= NA)

table4 <- data.frame(scenario = scenarios,
                     CDeaths = NA, DeathsAverted = NA, PctAverted = NA,
                     AvertedPerIso = NA,
                     CDeathsWork = NA, DeathsAvertedWork = NA, PctAvertedWork = NA,
                     AvertedPerIsoWork = NA)

table5 <- data.frame(scenario = scenarios,
                     CCase = NA, CaseAverted = NA, PctAverted = NA,
                     AvertedPerIso = NA, NNT = NA,
                     CCaseWork = NA, CaseAvertedWork = NA, PctAvertedWork = NA,
                     AvertedPerIsoWork = NA, NNTWork = NA)

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

    nnt <- addtl_iso/total_inf_averted
    nnt <- round(quantile(nnt, probs = c(.25, 0.5, .75)), 1)

    table3[rownum, 5] <- paste0(averted_per_iso[2], " (", averted_per_iso[1], ", ", averted_per_iso[3], ")")
    table3[rownum, 6] <- paste0(nnt[2], " (", nnt[1], ", ", nnt[3], ")")
  }

  inc.w <- colSums(sim$epi$se.flow.l1)
  pd.w <- colSums(sim$epi$num.w)
  inc_rate.w <- inc.w / pd.w * 100000
  inc_stats.w <- quantile(inc_rate.w, probs = c(.25, 0.5, .75))
  table3[rownum, 7] <- paste0(round(inc_stats.w[2], 1), " (", round(inc_stats.w[1], 1), ", ", round(inc_stats.w[3], 1), ")")

  if (!is.na(ref_inf_rate_w)){
    inf_averted.1.w <- round(ref_inf_rate_w - inc_stats.w[1], 1)
    inf_averted.2.w <- round(ref_inf_rate_w - inc_stats.w[2], 1)
    inf_averted.3.w <- round(ref_inf_rate_w - inc_stats.w[3], 1)

    table3[rownum, 8] <- paste0(inf_averted.2.w, " (", inf_averted.3.w, ", ", inf_averted.1.w, ")")

    pct_averted.1.w <- round((ref_inf_rate_w - inc_stats.w[1]) / ref_inf_rate_w, 3) * 100
    pct_averted.2.w <- round((ref_inf_rate_w - inc_stats.w[2]) / ref_inf_rate_w, 3) * 100
    pct_averted.3.w <- round((ref_inf_rate_w - inc_stats.w[3]) / ref_inf_rate_w, 3) * 100

    table3[rownum, 9] <- paste0(pct_averted.2.w, " (", pct_averted.3.w, ", ", pct_averted.1.w, ")")

    total_iso.w <- colSums(sim$epi$iso1.flow.w)
    addtl_iso.w <- total_iso.w - ref_iso_w
    total_inf.w <- colSums(sim$epi$se.flow.l1)
    total_inf_averted.w <- ref_inf_count_w - total_inf.w

    averted_per_iso.w <- total_inf_averted.w / addtl_iso.w
    averted_per_iso.w <- round(quantile(averted_per_iso.w, probs = c(.25, 0.5, .75)), 1)

    nnt.w <- addtl_iso.w/total_inf_averted.w
    nnt.w <- round(quantile(nnt.w, probs = c(.25, 0.5, .75)), 1)

    table3[rownum, 10] <- paste0(averted_per_iso.w[2], " (", averted_per_iso.w[1], ", ", averted_per_iso.w[3], ")")
    table3[rownum, 11] <- paste0(nnt.w[2], " (", nnt.w[1], ", ", nnt.w[3], ")")
  }

  return(table3)
}

table4_fill <- function(sim, rownum, table4,
                        ref_death_rate = NA,
                        ref_death_rate_w = NA,
                        ref_death_count = NA,
                        ref_death_count_w = NA,
                        ref_iso = NA,
                        ref_iso_w = NA){

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

  deaths.w <- colSums(sim$epi$d.h.flow.w)
  pds.w <- colSums(sim$epi$num.w)
  death_rates.w <- deaths.w / pds.w * 100000
  death_stats.w <- quantile(death_rates.w, probs = c(.25, 0.5, .75))
  table4[rownum, 6] <- paste0(round(death_stats.w[2], 3), " (", round(death_stats.w[1], 3), ", ", round(death_stats.w[3], 3), ")")

  if (!is.na(ref_death_rate_w)){
    deaths_averted.1.w <- round(ref_death_rate_w - death_stats.w[1], 3)
    deaths_averted.2.w <- round(ref_death_rate_w - death_stats.w[2], 3)
    deaths_averted.3.w <- round(ref_death_rate_w - death_stats.w[3], 3)

    table4[rownum, 7] <- paste0(deaths_averted.2.w, " (", deaths_averted.3.w, ", ", deaths_averted.1.w, ")")

    pct_averted.1.w <- round((ref_death_rate_w - death_stats.w[1]) / ref_death_rate_w, 3) * 100
    pct_averted.2.w <- round((ref_death_rate_w - death_stats.w[2]) / ref_death_rate_w, 3) * 100
    pct_averted.3.w <- round((ref_death_rate_w - death_stats.w[3]) / ref_death_rate_w, 3) * 100

    table4[rownum, 8] <- paste0(pct_averted.2.w, " (", pct_averted.3.w, ", ", pct_averted.1.w, ")")

    total_iso.w <- colSums(sim$epi$iso1.flow.w)
    addtl_iso.w <- total_iso.w - ref_iso_w
    total_deaths.w <- colSums(sim$epi$d.h.flow.w)
    total_deaths_averted.w <- ref_death_count_w - total_deaths.w
    averted_per_iso.w <- total_deaths_averted.w / addtl_iso.w

    averted_per_iso.w <- round(quantile(averted_per_iso.w, probs = c(.25, 0.5, .75)), 1)

    table4[rownum, 9] <- paste0(averted_per_iso.w[2], " (", averted_per_iso.w[1], ", ", averted_per_iso.w[3], ")")
  }

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

    nnt <- addtl_iso/total_case_averted
    nnt <- round(quantile(nnt, probs = c(.25, 0.5, .75)), 1)

    table5[rownum, 5] <- paste0(averted_per_iso[2], " (", averted_per_iso[1], ", ", averted_per_iso[3], ")")
    table5[rownum, 6] <- paste0(nnt[2], " (", nnt[1], ", ", nnt[3], ")")
  }

  inc.w <- colSums(sim$epi$ipic.flow.w)
  pd.w <- colSums(sim$epi$num.w)
  cinc_rate.w <- inc.w / pd.w * 100000
  cinc_stats.w <- quantile(cinc_rate.w, probs = c(.25, 0.5, .75))
  table5[rownum, 7] <- paste0(round(cinc_stats.w[2], 1), " (", round(cinc_stats.w[1], 1), ", ", round(cinc_stats.w[3], 1), ")")

  if (!is.na(ref_case_rate_w)){
    case_averted.1.w <- round(ref_case_rate_w - cinc_stats.w[1], 1)
    case_averted.2.w <- round(ref_case_rate_w - cinc_stats.w[2], 1)
    case_averted.3.w <- round(ref_case_rate_w - cinc_stats.w[3], 1)

    table5[rownum, 8] <- paste0(case_averted.2.w, " (", case_averted.3.w, ", ", case_averted.1.w, ")")

    pct_averted.1.w <- round((ref_case_rate_w - cinc_stats.w[1]) / ref_case_rate_w, 3) * 100
    pct_averted.2.w <- round((ref_case_rate_w - cinc_stats.w[2]) / ref_case_rate_w, 3) * 100
    pct_averted.3.w <- round((ref_case_rate_w - cinc_stats.w[3]) / ref_case_rate_w, 3) * 100

    table5[rownum, 9] <- paste0(pct_averted.2.w, " (", pct_averted.3.w, ", ", pct_averted.1.w, ")")

    total_iso.w <- colSums(sim$epi$iso1.flow.w)
    addtl_iso.w <- total_iso.w - ref_iso_w
    total_inf.w <- colSums(sim$epi$ipic.flow.w)
    total_case_averted.w <- ref_case_count_w - total_inf.w

    averted_per_iso.w <- total_case_averted.w / addtl_iso.w
    averted_per_iso.w <- round(quantile(averted_per_iso.w, probs = c(.25, 0.5, .75)), 1)

    nnt.w <- addtl_iso.w/total_case_averted.w
    nnt.w <- round(quantile(nnt.w, probs = c(.25, 0.5, .75)), 1)

    table5[rownum, 10] <- paste0(averted_per_iso.w[2], " (", averted_per_iso.w[1], ", ", averted_per_iso.w[3], ")")
    table5[rownum, 11] <- paste0(nnt.w[2], " (", nnt.w[1], ", ", nnt.w[3], ")")
  }

  return(table5)
}

# Fill in tables ---------------------------------------------------------------

# ref_inf_rate <- quantile(colSums(baseline_sims$epi$se.flow) / colSums(baseline_sims$epi$num) * 100000, probs = 0.5)
# ref_inf_count <- quantile(colSums(baseline_sims$epi$se.flow), probs = 0.5)
# ref_case_rate <- quantile(colSums(baseline_sims$epi$ipic.flow) / colSums(baseline_sims$epi$num) * 100000, probs = 0.5)
# ref_case_count <- quantile(colSums(baseline_sims$epi$ipic.flow), probs = 0.5)
# ref_death_rate <- quantile(colSums(baseline_sims$epi$d.h.flow) / colSums(baseline_sims$epi$num) * 100000, probs = 0.5)
# ref_death_count <- quantile(colSums(baseline_sims$epi$d.h.flow), probs = 0.5)
# ref_iso <- median(colSums(baseline_sims$epi$iso1.flow))
#
# # reference values among office workers
# ref_inf_rate_w <- quantile(colSums(baseline_sims$epi$se.flow.l1) / colSums(baseline_sims$epi$num.w) * 100000, probs = 0.5)
# ref_inf_count_w <- quantile(colSums(baseline_sims$epi$se.flow.l1), probs = 0.5)
# ref_case_rate_w <- quantile(colSums(baseline_sims$epi$ipic.flow.w) / colSums(baseline_sims$epi$num.w) * 100000, probs = 0.5)
# ref_case_count_w <- quantile(colSums(baseline_sims$epi$ipic.flow.w), probs = 0.5)
# ref_death_rate_w <- quantile(colSums(baseline_sims$epi$d.h.flow.w) / colSums(baseline_sims$epi$num.w) * 100000, probs = 0.5)
# ref_death_count_w <- quantile(colSums(baseline_sims$epi$d.h.flow.w), probs = 0.5)
# ref_iso_w <- median(colSums(baseline_sims$epi$iso1.flow.w))
#
# # Table 2
# for (i in seq(1,s.num+1)) {
#   table2 <- table2_fill(scenario_list[[i]],i,table2)
# }
#
# print("Table 2 done!")
#
# # Table 3
# for (i in seq(1,s.num+1)) {
#   table3 <- table3_fill(scenario_list[[i]],i,table3,
#                         ref_inf_rate = ref_inf_rate,
#                         ref_inf_rate_w = ref_inf_rate_w,
#                         ref_inf_count = ref_inf_count,
#                         ref_inf_count_w = ref_inf_count_w,
#                         ref_iso = ref_iso,
#                         ref_iso_w = ref_iso_w)
# }
#
# print("Table 3 done!")
#
# # Table 4
# for (i in seq(1,s.num+1)) {
#   table4 <- table4_fill(scenario_list[[i]],i,table4,
#                         ref_death_rate = ref_death_rate,
#                         ref_death_rate_w = ref_death_rate_w,
#                         ref_death_count = ref_death_count,
#                         ref_death_count_w = ref_death_count_w,
#                         ref_iso = ref_iso,
#                         ref_iso_w = ref_iso_w
#   )
# }
#
# print("Table 4 done!")
#
# # Table 5
# for (i in seq(1,s.num+1)) {
#   table5 <- table5_fill(scenario_list[[i]],i,table5,
#                         ref_case_rate = ref_case_rate,
#                         ref_case_rate_w = ref_case_rate_w,
#                         ref_case_count = ref_case_count,
#                         ref_case_count_w = ref_case_count_w,
#                         ref_iso = ref_iso,
#                         ref_iso_w = ref_iso_w)
# }
#
# print("Table 5 done!")
#
# # Create an Excel workbook
# wb <- createWorkbook()
#
# # Add data frames to different sheets
# addWorksheet(wb, "Table2")
# writeData(wb, sheet = 1, x = table2)
#
# addWorksheet(wb, "Table3")
# writeData(wb, sheet = 2, x = table3)
#
# addWorksheet(wb, "Table4")
# writeData(wb, sheet = 3, x = table4)
#
# addWorksheet(wb, "Table5")
# writeData(wb, sheet = 4, x = table5)
#
# # Save the workbook to a file
# saveWorkbook(wb, here::here("tbl_fig","result-tables.xlsx"))
#
# print("Excel with tables saved!")

# Ratio Plot: Cases to Isolated ------------------------------------------------
#
# ratio.plot <- function(sim,i) {
#   sim <- mutate_epi(sim,
#                     totCase = ic.num.w + h.num.w,
#                     totIso = iso.mild.num.w + iso.sev.num.w,
#                     totMask = mask.post.num.w + mask.exp.num.w)
#
#   sim <- mutate_epi(sim,
#                     CaseIsoRatio = totCase/totIso,
#                     CaseIsoSum = totCase+totIso)
#
#   df <- as.data.frame(sim) %>%
#     ungroup() %>%
#     group_by(time) %>%
#     summarise(med.case = median(totCase, na.rm = TRUE),
#               p25.case = quantile(totCase, probs = 0.25, na.rm = TRUE),
#               p75.case = quantile(totCase, probs = 0.75, na.rm = TRUE),
#               med.iso = median(totIso, na.rm = TRUE),
#               p25.iso = quantile(totIso, probs = 0.25, na.rm = TRUE),
#               p75.iso = quantile(totIso, probs = 0.75, na.rm = TRUE),
#               med.mask = median(totMask, na.rm = TRUE),
#               p25.mask = quantile(totMask, probs = 0.25, na.rm = TRUE),
#               p75.mask = quantile(totMask, probs = 0.75, na.rm = TRUE),
#               med.ratio = median(CaseIsoRatio, na.rm = TRUE),
#               p25.ratio = quantile(CaseIsoRatio, probs = 0.25, na.rm = TRUE),
#               p75.ratio = quantile(CaseIsoRatio, probs = 0.75, na.rm = TRUE),
#               med.sum = median(CaseIsoSum, na.rm = TRUE),
#               p25.sum = quantile(CaseIsoSum, probs = 0.25, na.rm = TRUE),
#               p75.sum = quantile(CaseIsoSum, probs = 0.75, na.rm = TRUE)) %>%
#     mutate(scenario = i)
#
#   return(df)
# }
#
# d.ratio <- data.frame(time = NA,
#                       med.case = NA, p25.case = NA, p75.case = NA,
#                       med.iso = NA, p25.iso = NA, p75.iso = NA,
#                       med.mask = NA, p25.mask = NA, p75.mask = NA,
#                       med.ratio = NA, p25.ratio = NA, p75.ratio = NA,
#                       med.sum = NA, p25.sum = NA, p75.sum = NA,
#                       scenario = i)
#
# for (i in seq(1,s.num+1,1)) {
#
#   temp <- ratio.plot(scenario_list[[i]],i-1)
#
#   d.ratio <- d.ratio %>%
#     add_row(temp) %>%
#     filter(!is.na(time))
#
#   rm(temp)
#
# }
#
# d.ratio <- d.ratio %>%
#   mutate(scenario2 = factor(scenario,
#                             levels = seq(0,s.num)))
#
# saveRDS(d.ratio, file = here::here("tbl_fig","ratio_df.rds"))
#
# print("Ratio data frame saved!")

# Summary measures by different time periods -----------------------------------

# Time 1 [000-151]: Jan-21 to June-21 (endemic)
# Time 2 [152-304]: June-21 to Nov-21 (Delta wave)
# Time 3 [305-424]: Nov-21 to Mar-22 (Omicron wave)
# Time 4 [425-608]: Mar-22 to Aug-22 (endemic)

table_nnt <- data.frame(scenario = scenarios,
                        time1 = NA,
                        T1_NNT_Inf = NA, T1_NNT_InfW = NA,
                        T1_NNT_Case = NA, T1_NNT_CaseW = NA,
                        time2 = NA,
                        T2_NNT_Inf = NA, T2_NNT_InfW = NA,
                        T2_NNT_Case = NA, T2_NNT_CaseW = NA,
                        time3 = NA,
                        T3_NNT_Inf = NA, T3_NNT_InfW = NA,
                        T3_NNT_Case = NA, T3_NNT_CaseW = NA,
                        time4 = NA,
                        T4_NNT_Inf = NA, T4_NNT_InfW = NA,
                        T4_NNT_Case = NA, T4_NNT_CaseW = NA)

nnt_fill <- function(sim,rownum,table_nnt) {

  colnum <- 2

  # Define row ranges
  row_ranges <- list(c(1, 151),
                     c(152, 304),
                     c(305, 424),
                     c(425, 608))

  # Iterate through row ranges
  for (i in seq_along(row_ranges)) {

    row_range <- row_ranges[[i]]
    table_nnt[rownum,colnum] <- paste0(row_range[[1]],"-",row_range[[2]])
    colnum <- colnum + 1

    ## reference scenario values
    ref_iso <- median(colSums(baseline_sims$epi$iso1.flow[row_range[[1]]:row_range[[2]], ]))
    ref_iso_w <- median(colSums(baseline_sims$epi$iso1.flow.w[row_range[[1]]:row_range[[2]], ]))
    ref_inf_count <- quantile(colSums(baseline_sims$epi$se.flow[row_range[[1]]:row_range[[2]], ]), probs = 0.5)
    ref_inf_count_w <- quantile(colSums(baseline_sims$epi$se.flow.l1[row_range[[1]]:row_range[[2]], ]), probs = 0.5)
    ref_case_count <- quantile(colSums(baseline_sims$epi$ipic.flow[row_range[[1]]:row_range[[2]], ]), probs = 0.5)
    ref_case_count_w <- quantile(colSums(baseline_sims$epi$ipic.flow.w[row_range[[1]]:row_range[[2]], ]), probs = 0.5)

    ## total/additional isolations
    total_iso <- colSums(sim$epi$iso1.flow[row_range[[1]]:row_range[[2]], ])
    addtl_iso <- total_iso - ref_iso

    total_iso.w <- colSums(sim$epi$iso1.flow.w[row_range[[1]]:row_range[[2]], ])
    addtl_iso.w <- total_iso.w - ref_iso_w

    ## infections
    total_inf <- colSums(sim$epi$se.flow[row_range[[1]]:row_range[[2]], ])
    total_inf_averted <- ref_inf_count - total_inf

    nnt.inf <- addtl_iso/total_inf_averted
    nnt.inf <- round(quantile(nnt.inf, probs = c(.25, 0.5, .75)), 1)

    table_nnt[rownum,colnum] <- paste0(nnt.inf[2], " (",
                                       nnt.inf[1], ", ",
                                       nnt.inf[3], ")")
    colnum <- colnum + 1

    total_inf.w <- colSums(sim$epi$se.flow.l1[row_range[[1]]:row_range[[2]], ])
    total_inf_averted.w <- ref_inf_count_w - total_inf.w

    nnt.inf.w <- addtl_iso.w/total_inf_averted.w
    nnt.inf.w <- round(quantile(nnt.inf.w, probs = c(.25, 0.5, .75)), 1)

    table_nnt[rownum,colnum] <- paste0(nnt.inf.w[2], " (",
                                       nnt.inf.w[1], ", ",
                                       nnt.inf.w[3], ")")
    colnum <- colnum + 1

    ## cases
    total_case <- colSums(sim$epi$ipic.flow[row_range[[1]]:row_range[[2]], ])
    total_case_averted <- ref_case_count - total_case

    nnt.case <- addtl_iso/total_case_averted
    nnt.case <- round(quantile(nnt.case, probs = c(.25, 0.5, .75)), 1)

    table_nnt[rownum,colnum] <- paste0(nnt.case[2], " (",
                                       nnt.case[1], ", ",
                                       nnt.case[3], ")")
    colnum <- colnum + 1

    total_case.w <- colSums(sim$epi$ipic.flow.w[row_range[[1]]:row_range[[2]], ])
    total_case_averted.w <- ref_case_count_w - total_case.w

    nnt.case.w <- addtl_iso.w/total_case_averted.w
    nnt.case.w <- round(quantile(nnt.case.w, probs = c(.25, 0.5, .75)), 1)

    table_nnt[rownum,colnum] <- paste0(nnt.case.w[2], " (",
                                       nnt.case.w[1], ", ",
                                       nnt.case.w[3], ")")
    colnum <- colnum + 1

  }

  return(table_nnt)

}

# NNT Table
for (i in seq(1,s.num+1)) {
  table_nnt <- nnt_fill(scenario_list[[i]],i,table_nnt)
}

# Save table

# Create an Excel workbook
wb <- createWorkbook()

# Add data frames to different sheets
addWorksheet(wb, "TableNNT")
writeData(wb, sheet = 1, x = table_nnt)

# Save the workbook to a file
saveWorkbook(wb, here::here("tbl_fig","nnt-table.xlsx"))

print("NNT table done!")
