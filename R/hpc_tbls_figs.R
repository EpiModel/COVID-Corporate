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
base.sim <- 8000
my.dir <- "sims/"
s.num <- 5 # number of scenarios

# Tables: FOR EXISTING OUTPUT Load files ---------------------------------------

scenario_list <- readRDS(here::here("processed_sims",
                                    "scenario_list.rds"))

baseline_sims <- readRDS(here::here("processed_sims",
                                    "baseline_sim.rds"))

# scenario0 <- baseline_sims

# Table shells -----------------------------------------------------------------

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
# saveWorkbook(wb, here::here("tbl_fig","result-tables-8000.xlsx"))
#
# print("Excel with tables saved!")

# Ratio Plot: Cases to Isolated ------------------------------------------------

# ratio.plot <- function(sim,i) {
#   sim <- mutate_epi(sim,
#                     totCase = ic.num.w + h.num.w,
#                     totIso = iso.mild.num.w + iso.sev.num.w,
#                     totMask = mask.post.num.w + mask.exp.num.w,
#                     CaseNotIso = ic.num.w + h.num.w - iso.and.sick.num.w,
#                     IsoNotCase = iso.mild.num.w + iso.sev.num.w - iso.and.sick.num.w)
#
#   sim <- mutate_epi(sim,
#                     #CaseIsoRatio = totCase/totIso,
#                     CaseIsoSum1 = totCase+totIso - iso.and.sick.num.w,
#                     CaseIsoSum2 = iso.or.sick.num.w)
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
#               # med.ratio = median(CaseIsoRatio, na.rm = TRUE),
#               # p25.ratio = quantile(CaseIsoRatio, probs = 0.25, na.rm = TRUE),
#               # p75.ratio = quantile(CaseIsoRatio, probs = 0.75, na.rm = TRUE),
#               med.sum1 = median(CaseIsoSum1, na.rm = TRUE),
#               p25.sum1 = quantile(CaseIsoSum1, probs = 0.25, na.rm = TRUE),
#               p75.sum1 = quantile(CaseIsoSum1, probs = 0.75, na.rm = TRUE),
#               med.sum2 = median(CaseIsoSum2, na.rm = TRUE),
#               p25.sum2 = quantile(CaseIsoSum2, probs = 0.25, na.rm = TRUE),
#               p75.sum2 = quantile(CaseIsoSum2, probs = 0.75, na.rm = TRUE),
#               med.both = median(iso.and.sick.num.w, na.rm = TRUE),
#               p25.both = quantile(iso.and.sick.num.w, probs = 0.25, na.rm = TRUE),
#               p75.both = quantile(iso.and.sick.num.w, probs = 0.75, na.rm = TRUE),
#               med.casenotiso = median(CaseNotIso, na.rm = TRUE),
#               p25.casenotiso = quantile(CaseNotIso, probs = 0.25, na.rm = TRUE),
#               p75.casenotiso = quantile(CaseNotIso, probs = 0.75, na.rm = TRUE),
#               med.isonotcase = median(IsoNotCase, na.rm = TRUE),
#               p25.isonotcase = quantile(IsoNotCase, probs = 0.25, na.rm = TRUE),
#               p75.isonotcase = quantile(IsoNotCase, probs = 0.75, na.rm = TRUE)) %>%
#     mutate(scenario = i)
#
#   return(df)
# }
#
# d.ratio <- data.frame(time = NA,
#                       med.case = NA, p25.case = NA, p75.case = NA,
#                       med.iso = NA, p25.iso = NA, p75.iso = NA,
#                       med.mask = NA, p25.mask = NA, p75.mask = NA,
#                       # med.ratio = NA, p25.ratio = NA, p75.ratio = NA,
#                       med.sum1 = NA, p25.sum1 = NA, p75.sum1 = NA,
#                       med.sum2 = NA, p25.sum2 = NA, p75.sum2 = NA,
#                       med.both = NA, p25.both = NA, p75.both = NA,
#                       med.casenotiso = NA, p25.casenotiso = NA, p75.casenotiso = NA,
#                       med.isonotcase = NA, p25.isonotcase = NA, p75.isonotcase = NA,
#                       scenario = NA)
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
# saveRDS(d.ratio, file = here::here("tbl_fig","ratio_df_8000.rds"))
#
# print("Ratio data frame saved!")

# Summary measures by different time periods -----------------------------------

# Time 1 [000-151]: Jan-21 to June-21 (endemic)
# Time 2 [152-304]: June-21 to Nov-21 (Delta wave)
# Time 3 [305-424]: Nov-21 to Mar-22 (Omicron wave)
# Time 4 [425-608]: Mar-22 to Aug-22 (endemic)

# table_nnt <- data.frame(scenario = scenarios,
#                         time1 = NA,
#                         T1_NNT_Inf = NA, T1_NNT_InfW = NA,
#                         T1_NNT_Case = NA, T1_NNT_CaseW = NA,
#                         time2 = NA,
#                         T2_NNT_Inf = NA, T2_NNT_InfW = NA,
#                         T2_NNT_Case = NA, T2_NNT_CaseW = NA,
#                         time3 = NA,
#                         T3_NNT_Inf = NA, T3_NNT_InfW = NA,
#                         T3_NNT_Case = NA, T3_NNT_CaseW = NA,
#                         time4 = NA,
#                         T4_NNT_Inf = NA, T4_NNT_InfW = NA,
#                         T4_NNT_Case = NA, T4_NNT_CaseW = NA)
#
# table_nnt <- data.frame(scenario = scenarios,
#                         T1_RefIso = NA, T1_TotalIso = NA,
#                         T1_AddtlIso = NA,
#                         T1_RefIsoW = NA, T1_TotalIsoW = NA,
#                         T1_AddtlIsoW = NA,
#
#                         T1_RefInf = NA, T1_TotalInf = NA,
#                         T1_InfAvert = NA,
#                         T1_NNT_Inf = NA,
#
#                         T1_RefInfW = NA, T1_TotalInfW = NA,
#                         T1_InfAvertW = NA,
#                         T1_NNT_InfW = NA,
#
#                         T1_RefCase = NA, T1_TotalCase = NA,
#                         T1_CaseAvert = NA,
#                         T1_NNT_Case = NA,
#
#                         T1_RefCaseW = NA, T1_TotalCaseW = NA,
#                         T1_CaseAvertW = NA,
#                         T1_NNT_CaseW = NA,
#                         ###
#                         T2_RefIso = NA, T2_TotalIso = NA,
#                         T2_AddtlIso = NA,
#                         T2_RefIsoW = NA, T2_TotalIsoW = NA,
#                         T2_AddtlIsoW = NA,
#
#                         T2_RefInf = NA, T2_TotalInf = NA,
#                         T2_InfAvert = NA,
#                         T2_NNT_Inf = NA,
#
#                         T2_RefInfW = NA, T2_TotalInfW = NA,
#                         T2_InfAvertW = NA,
#                         T2_NNT_InfW = NA,
#
#                         T2_RefCase = NA, T2_TotalCase = NA,
#                         T2_CaseAvert = NA,
#                         T2_NNT_Case = NA,
#
#                         T2_RefCaseW = NA, T2_TotalCaseW = NA,
#                         T2_CaseAvertW = NA,
#                         T2_NNT_CaseW = NA,
#                         ###
#                         T3_RefIso = NA, T3_TotalIso = NA,
#                         T3_AddtlIso = NA,
#                         T3_RefIsoW = NA, T3_TotalIsoW = NA,
#                         T3_AddtlIsoW = NA,
#
#                         T3_RefInf = NA, T3_TotalInf = NA,
#                         T3_InfAvert = NA,
#                         T3_NNT_Inf = NA,
#
#                         T3_RefInfW = NA, T3_TotalInfW = NA,
#                         T3_InfAvertW = NA,
#                         T3_NNT_InfW = NA,
#
#                         T3_RefCase = NA, T3_TotalCase = NA,
#                         T3_CaseAvert = NA,
#                         T3_NNT_Case = NA,
#
#                         T3_RefCaseW = NA, T3_TotalCaseW = NA,
#                         T3_CaseAvertW = NA,
#                         T3_NNT_CaseW = NA,
#                         ###
#                         T4_RefIso = NA, T4_TotalIso = NA,
#                         T4_AddtlIso = NA,
#                         T4_RefIsoW = NA, T4_TotalIsoW = NA,
#                         T4_AddtlIsoW = NA,
#
#                         T4_RefInf = NA, T4_TotalInf = NA,
#                         T4_InfAvert = NA,
#                         T4_NNT_Inf = NA,
#
#                         T4_RefInfW = NA, T4_TotalInfW = NA,
#                         T4_InfAvertW = NA,
#                         T4_NNT_InfW = NA,
#
#                         T4_RefCase = NA, T4_TotalCase = NA,
#                         T4_CaseAvert = NA,
#                         T4_NNT_Case = NA,
#
#                         T4_RefCaseW = NA, T4_TotalCaseW = NA,
#                         T4_CaseAvertW = NA,
#                         T4_NNT_CaseW = NA)
#
# nnt_fill <- function(sim,rownum,table_nnt) {
#
#   colnum <- 2
#
#   # Define row ranges
#   row_ranges <- list(c(1, 151),
#                      c(152, 304),
#                      c(305, 424),
#                      c(425, 608))
#
#   # Iterate through row ranges
#   for (i in seq_along(row_ranges)) {
#
#     row_range <- row_ranges[[i]]
#     # table_nnt[rownum,colnum] <- paste0(row_range[[1]],"-",row_range[[2]])
#     # colnum <- colnum + 1
#
#     ## reference scenario values
#     ref_iso <- median(colSums(baseline_sims$epi$iso1.flow[row_range[[1]]:row_range[[2]], ]))
#     ref_iso_w <- median(colSums(baseline_sims$epi$iso1.flow.w[row_range[[1]]:row_range[[2]], ]))
#     ref_inf_count <- quantile(colSums(baseline_sims$epi$se.flow[row_range[[1]]:row_range[[2]], ]), probs = 0.5)
#     ref_inf_count_w <- quantile(colSums(baseline_sims$epi$se.flow.l1[row_range[[1]]:row_range[[2]], ]), probs = 0.5)
#     ref_case_count <- quantile(colSums(baseline_sims$epi$ipic.flow[row_range[[1]]:row_range[[2]], ]), probs = 0.5)
#     ref_case_count_w <- quantile(colSums(baseline_sims$epi$ipic.flow.w[row_range[[1]]:row_range[[2]], ]), probs = 0.5)
#
#     ## total/additional isolations
#     total_iso <- colSums(sim$epi$iso1.flow[row_range[[1]]:row_range[[2]], ])
#     addtl_iso <- total_iso - ref_iso
#
#     total_iso.w <- colSums(sim$epi$iso1.flow.w[row_range[[1]]:row_range[[2]], ])
#     addtl_iso.w <- total_iso.w - ref_iso_w
#
#     variables <- list(ref_iso, total_iso, addtl_iso, ref_iso_w, total_iso.w, addtl_iso.w)
#     for (var in variables) {
#       temp <- round(quantile(var, probs = c(.25, 0.5, .75)), 1)
#       table_nnt[rownum, colnum] <- paste0(temp[2], " (",
#                                           temp[1], ", ",
#                                           temp[3], ")")
#       colnum <- colnum + 1
#     }
#
#     ## infections
#     total_inf <- colSums(sim$epi$se.flow[row_range[[1]]:row_range[[2]], ])
#     total_inf_averted <- ref_inf_count - total_inf
#
#     nnt.inf <- addtl_iso/total_inf_averted
#
#     total_inf.w <- colSums(sim$epi$se.flow.l1[row_range[[1]]:row_range[[2]], ])
#     total_inf_averted.w <- ref_inf_count_w - total_inf.w
#
#     nnt.inf.w <- addtl_iso.w/total_inf_averted.w
#
#     variables <- list(ref_inf_count, total_inf, total_inf_averted, nnt.inf,
#                       ref_inf_count_w, total_inf.w, total_inf_averted.w, nnt.inf.w)
#     for (var in variables) {
#       temp <- round(quantile(var, probs = c(.25, 0.5, .75)), 1)
#       table_nnt[rownum, colnum] <- paste0(temp[2], " (",
#                                           temp[1], ", ",
#                                           temp[3], ")")
#       colnum <- colnum + 1
#     }
#
#     ## cases
#     total_case <- colSums(sim$epi$ipic.flow[row_range[[1]]:row_range[[2]], ])
#     total_case_averted <- ref_case_count - total_case
#
#     nnt.case <- addtl_iso/total_case_averted
#
#     total_case.w <- colSums(sim$epi$ipic.flow.w[row_range[[1]]:row_range[[2]], ])
#     total_case_averted.w <- ref_case_count_w - total_case.w
#
#     nnt.case.w <- addtl_iso.w/total_case_averted.w
#
#     variables <- list(ref_case_count, total_case, total_case_averted, nnt.case,
#                       ref_case_count_w, total_case.w, total_case_averted.w, nnt.case.w)
#     for (var in variables) {
#       temp <- round(quantile(var, probs = c(.25, 0.5, .75)), 1)
#       table_nnt[rownum, colnum] <- paste0(temp[2], " (",
#                                           temp[1], ", ",
#                                           temp[3], ")")
#       colnum <- colnum + 1
#     }
#
#
#   }
#
#   return(table_nnt)
#
# }
#
# # NNT Table
# for (i in seq(1,s.num+1)) {
#   table_nnt <- nnt_fill(scenario_list[[i]],i,table_nnt)
# }
#
# # Save table
#
# # Create an Excel workbook
# wb <- createWorkbook()
#
# # Add data frames to different sheets
# addWorksheet(wb, "TableNNT")
# writeData(wb, sheet = 1, x = table_nnt)
#
# # Save the workbook to a file
# saveWorkbook(wb, here::here("tbl_fig","nnt-table.xlsx"))
#
# print("NNT table done!")


# Cumulative isolation and cases for individual simulations ----
for (i in c(0,2,4,5)) {
  iso.prob <- scenario_list[[i+1]]$param$iso.prob

  sim <- mutate_epi(scenario_list[[i+1]],
                    cum.iso = cumsum(scenario_list[[i+1]]$epi$iso1.flow),
                    cum.case = cumsum(scenario_list[[i+1]]$epi$ipic.flow),
                    cum.iso.w = cumsum(scenario_list[[i+1]]$epi$iso1.flow.w),
                    cum.case.w = cumsum(scenario_list[[i+1]]$epi$ipic.flow.w))

  assign(paste0("d", i),
         as.data.frame(sim) %>%
           filter(time == 608) %>%
           select(sim,time,starts_with("cum.")) %>%
           mutate(iso.prob.param = iso.prob))
}

d <- bind_rows(d0, d2, d4, d5)

saveRDS(d, file = here::here("tbl_fig","sim_df.rds"))
