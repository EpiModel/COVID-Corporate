# Load libraries ---------------------------------------------------------------
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(tidyr)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(patchwork)))

# Ratio Plot: Cases to Isolated ------------------------------------------------

d.ratio <- readRDS(here::here("results","ratio_df.rds"))

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
