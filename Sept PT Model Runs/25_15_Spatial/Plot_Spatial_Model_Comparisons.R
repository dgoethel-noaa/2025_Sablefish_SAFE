# Purpose: To plot and compare among spatial models, as well as with spatially-aggregated models
# Creator: Matthew LH. Cheng
# Date 7/29/25

# Set up ------------------------------------------------------------------

library(here)
library(tidyverse)
library(SPoRC)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(cowplot)

# Plotting dimensions
years <- 1960:2024

# Read in models
spt <- readRDS(here("Sept PT Model Runs", "25_15_Spatial", "Spatial_MltRel", "Spatial_MltRel_model_results.RDS")) # "final" spatial model
sgl <- readRDS(here("Sept PT Model Runs", "25_12_Drop_TS_Upd_M", "25_12_Drop_TS_Upd_M_model_results.RDS")) # "final" spatially aggregated model

# Read in apportionment
apport <- read.csv(here("Data Pulls", "data outputs", "lls_RPW_area_prop.csv"))

## Movement ----------------------------------------------------------------
# Read in maps here
west <- ne_states(c("United States of America", "Russia", "Canada"), returnclass = "sf")
west <- st_shift_longitude(west) # shift ongitude for plotting

# Read in stat areas
nmfs_areas <- read_sf(dsn = here("Sept PT Model Runs", "25_15_Spatial", "NMFS_Stat_Areas", "Sablefish_Longline_Area"), layer = "Sablefish_Longline_Area")
nmfs_areas <- nmfs_areas %>% mutate(GEN_NAME = ifelse(NAME %in% c("East Yakutat / Southeast Alaska", "West Yakutat"), "Eastern Gulf of Alaska", "A")) %>%
  mutate(NAME = case_when(
    NAME == "Aleutian Islands" ~ "AI",
    NAME == "Bering Sea" ~ "BS",
    NAME == "Western Gulf of Alaska" ~ "WGOA",
    NAME == "Central Gulf of Alaska" ~ "CGOA",
    NAME == "West Yakutat" ~ "EGOA",
    NAME == "East Yakutat / Southeast Alaska" ~ "EGOA"
  ), NAME = factor(NAME, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA"))) %>%
  group_by(NAME) %>%
  summarise(geometry = st_union(geometry))

# Coerce longline areas
nmfs_areas <- st_make_valid(nmfs_areas) # make valid so that vertices aren't duplicated
nmfs_areas <- nmfs_areas %>% st_transform(4326) # transform to crs 4326
nmfs_areas <- st_shift_longitude(nmfs_areas) # shift longitude for plotting

# get centroids of the geometry for plotting
centroids <- nmfs_areas %>%
  group_by(NAME) %>%
  summarise(geometry = st_centroid(geometry)) %>%
  ungroup()

# rename variables
move_df <- reshape2::melt(spt$rep$Movement) %>%
  mutate(    Var1 = case_when(
    Var1 == 1 ~ "BS",
    Var1 == 2 ~ "AI",
    Var1 == 3 ~ "WGOA",
    Var1 == 4 ~ "CGOA",
    Var1 == 5 ~ "EGOA"
  ), Var1 = factor(Var1, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA")),
  Var2 = case_when(
    Var2 == 1 ~ "BS",
    Var2 == 2 ~ "AI",
    Var2 == 3 ~ "WGOA",
    Var2 == 4 ~ "CGOA",
    Var2 == 5 ~ "EGOA"
  ), Var2 = factor(Var2, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA"))) %>%
  rename(From = Var1, To = Var2, Year = Var3, Age = Var4, Sex = Var5) %>%
  filter(Sex == 1, Year == 1, Age %in% c(2,8,30)) # filter to appropriate movement parameters

# Left join centroids to dataframe for plotting where from -> to
move_df_with_coords <- move_df %>%
  left_join(centroids, by = c("From" = "NAME")) %>%
  rename(from_geometry = geometry) %>%
  left_join(centroids, by = c("To" = "NAME")) %>%
  rename(to_geometry = geometry) %>%
  mutate(
    from_x = st_coordinates(st_shift_longitude(from_geometry))[, 1], # get from centroids
    from_y = st_coordinates(from_geometry)[, 2],
    to_x = st_coordinates(st_shift_longitude(to_geometry))[, 1], # get to centroids
    to_y = st_coordinates(to_geometry)[, 2],
    Age = case_when(
      Age == 2 ~ "Young (2 - 7)",
      Age == 8 ~ "Intermediate (8 - 15)",
      Age == 30 ~ "Old (16 - 31)"
    ),
    Age = factor(Age, levels = c("Young (2 - 7)", "Intermediate (8 - 15)", "Old (16 - 31)"))
  )

# plot movement
colors <- unname(ggthemes::ggthemes_data[["colorblind"]][["value"]]) # get colors

map_movement_plot <- ggplot() +
  geom_sf(data = nmfs_areas, alpha = 0.55) +
  geom_sf(data = west, lwd = 0.05, color = 'black', alpha = 1) + # World Map
  geom_curve(data = move_df_with_coords %>% filter(From != To),
             aes(x = from_x, y = from_y, xend = to_x, yend = to_y, color = To, size = value),
             curvature = 0.4, alpha = 0.45) + # movement arrows (from != to)
  geom_text(data = move_df_with_coords %>% filter(From != To),
            aes(x = to_x, y = to_y, label = round(value, 2), color = To),
            alpha = 1, size = 13, nudge_y = 0.07, show.legend = F) + # from != to
  scale_color_manual(values = colors[-c(1,5)]) +
  geom_text(data = move_df_with_coords %>% filter(From == To), size = 13, # from == to
            aes(x = from_x, y = from_y, label = round(value, 2)), alpha = 1,
            color = 'black', nudge_y = 0.07, show.legend = F) +
  facet_grid(Age~From) +
  scale_size(guide = "none") +
  guides(color = guide_legend(override.aes = list(
    linetype = 1, linewidth = 5, alpha = 1, size = 3))) +
  coord_sf(ylim = c(45.2, 70.5), xlim = c(165, 230)) + # Restrict Map Area
  theme_bw(base_size = 30) +
  labs(x = "Longitude", y = "Latitude", size = 'P(movement)',
       color = "To", label = "P(residence)") +
  theme(legend.position = 'top',
        legend.box = 'vertical',
        legend.background = element_blank(),
        legend.spacing.y = unit(-1, 'cm'),  # Adjust spacing between items
        plot.background = element_rect(fill = "transparent", colour = NA),
        axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave(map_movement_plot,
       filename = here("Sept PT Model Runs", "25_15_Spatial", "Figures", "Map Movement.png"),
       width = 28, height = 20
)

# Reference Points --------------------------------------------------------

# compute single region reference points
sgl_ref_pts <- SPoRC::Get_Reference_Points(data = sgl$data,
                                           rep = sgl$rep,
                                           SPR_x = 0.4,
                                           t_spwn = 0,
                                           sex_ratio_f = 0.5,
                                           calc_rec_st_yr = 20,
                                           rec_age = 2,
                                           type = "single_region",
                                           what = "SPR")

# compute spatial reference points
spt_ref_pts <- SPoRC::Get_Reference_Points(data = spt$data,
                                           rep = spt$rep,
                                           SPR_x = 0.4,
                                           t_spwn = 0,
                                           sex_ratio_f = 0.5,
                                           calc_rec_st_yr = 20,
                                           rec_age = 2,
                                           type = "multi_region",
                                           what = "global_SPR")

## SSB Comparison ----------------------------------------------------------

### Spatial SSB -------------------------------------------------------------
# munging spatial ssb data
spatial_ssb <- reshape2::melt(spt$rep$SSB) %>%
  rename(Region = Var1, Year = Var2) %>%
  mutate(Type = 'Spatial',
         Region =
           case_when(
             Region == 1 ~ "BS",
             Region == 2 ~ "AI",
             Region == 3 ~ "WGOA",
             Region == 4 ~ "CGOA",
             Region == 5 ~ "EGOA"
           ),
         Region = factor(Region, levels = c("BS", "AI", "WGOA", "CGOA", 'EGOA')),
         Year = Year + 1959)


# Spatial SSB
spt_ssb <- ggplot() +
  geom_line(spatial_ssb, mapping = aes(x = Year, y = value, color = Region), lwd = 1.3) +
  facet_wrap(~Region, nrow = 1) +
  theme_bw(base_size = 15) +
  theme(legend.position = 'none') +
  ggthemes::scale_color_colorblind() +
  coord_cartesian(ylim = c(0, NA)) +
  labs(x = 'Year', y = 'Spawning Stock Biomass (kt)')

ggsave(spt_ssb,
       filename = here("Sept PT Model Runs", "25_15_Spatial", "Figures", "Spatial SSB.png"),
       width = 15, height = 7
)

### Aggregated SSB ----------------------------------------------------------

# Compare aggregated SSB
agg_ssb <- reshape2::melt(colSums(spt$rep$SSB)) %>%
  mutate(Year = years,
         Type = 'Spatial')

# Combine with single region SSB
agg_ssb <- agg_ssb %>%
  bind_rows(
    reshape2::melt(sgl$rep$SSB) %>%
      mutate(Year = years,
             Type = 'Single Region') %>%
      dplyr::select(-c(Var1, Var2))
  ) %>%
  # joint B40 reference points
  left_join(data.frame(Bref = c(sgl_ref_pts$b_ref_pt, sum(spt_ref_pts$b_ref_pt)),
                       B0 = c(sgl_ref_pts$virgin_b_ref_pt, sum(spt_ref_pts$virgin_b_ref_pt)),
                       Type = c("Single Region", "Spatial")), by = 'Type')

ssb_mod_comp <- ggplot() +
  geom_line(agg_ssb, mapping = aes(x = Year, y = value, color = Type)) +
  geom_point(agg_ssb, mapping = aes(x = Year, y = value, color = Type), size = 3) +
  geom_hline(agg_ssb, mapping = aes(yintercept = Bref, color = Type), lty = 2, lwd = 2) +
  theme_bw(base_size = 15) +
  theme(legend.position = 'top') +
  coord_cartesian(ylim = c(0, NA)) +
  labs(x = 'Year', y = 'Spawning Stock Biomass (kt)')

ssb_mod_dep <- ggplot() +
  geom_line(agg_ssb, mapping = aes(x = Year, y = value / B0, color = Type)) +
  geom_point(agg_ssb, mapping = aes(x = Year, y = value / B0, color = Type), size = 3) +
  theme_bw(base_size = 15) +
  theme(legend.position = 'none') +
  coord_cartesian(ylim = c(0, NA)) +
  labs(x = 'Year', y = 'Spawning Stock Biomass (kt) / B0')

ggsave(plot_grid(ssb_mod_comp, ssb_mod_dep, labels = c("A", "B"), ncol = 1, label_size = 25),
       filename = here("Sept PT Model Runs", "25_15_Spatial", "Figures", "SSB Comparison (Spt vs. Agg).png"),
       width = 8, height = 10
)

## Recruitment Comparison ----------------------------------------------------------

### Spatial Recruitment -------------------------------------------------------------
# munging spatial Recruitment data
spatial_rec <- reshape2::melt(spt$rep$Rec) %>%
  rename(Region = Var1, Year = Var2) %>%
  mutate(Type = 'Spatial',
         Region =
           case_when(
             Region == 1 ~ "BS",
             Region == 2 ~ "AI",
             Region == 3 ~ "WGOA",
             Region == 4 ~ "CGOA",
             Region == 5 ~ "EGOA"
           ),
         Region = factor(Region, levels = c("BS", "AI", "WGOA", "CGOA", 'EGOA')),
         Year = Year + 1959)


# Spatial Recruitment
spt_rec <- ggplot() +
  geom_line(spatial_rec, mapping = aes(x = Year, y = value, color = Region), lwd = 1.3) +
  facet_wrap(~Region, nrow = 1) +
  theme_bw(base_size = 15) +
  theme(legend.position = 'none') +
  ggthemes::scale_color_colorblind() +
  coord_cartesian(ylim = c(0, NA)) +
  labs(x = 'Year', y = 'Recruitment (age 2; millions)')

ggsave(spt_rec,
       filename = here("Sept PT Model Runs", "25_15_Spatial", "Figures", "Spatial Recruitment.png"),
       width = 15, height = 7
)

### Aggregated Recruitment ----------------------------------------------------------

# Compare aggregated recruitment
agg_rec <- reshape2::melt(colSums(spt$rep$Rec)) %>%
  mutate(Year = years,
         Type = 'Spatial')

# Combine with single region SSB
agg_rec <- agg_rec %>%
  bind_rows(
    reshape2::melt(sgl$rep$Rec) %>%
      mutate(Year = years,
             Type = 'Single Region') %>%
      dplyr::select(-c(Var1, Var2))
  )

rec_mod_comp <- ggplot() +
  geom_line(agg_rec, mapping = aes(x = Year, y = value, color = Type)) +
  geom_point(agg_rec, mapping = aes(x = Year, y = value, color = Type), size = 3) +
  theme_bw(base_size = 15) +
  theme(legend.position = 'top') +
  coord_cartesian(ylim = c(0, NA)) +
  labs(x = 'Year', y = 'Recruitment (age 2; millions)')

ggsave(rec_mod_comp,
       filename = here("Sept PT Model Runs", "25_15_Spatial", "Figures", "Recruitment Comparison (Spt vs. Agg).png"),
       width = 15, height = 7
)

# Apportionment -----------------------------------------------------------

# Comparison of total biomass from spatial model vs. total biomass apportioned from a single region model
spt_tot_biom <- reshape2::melt(spt$rep$Total_Biom) %>%
  rename(Region = Var1, Year = Var2) %>%
  mutate(Year = Year + 1959,
         Type = 'Spatial',
         Region =
           case_when(
             Region == 1 ~ "BS",
             Region == 2 ~ "AI",
             Region == 3 ~ "WGOA",
             Region == 4 ~ "CGOA",
             Region == 5 ~ "EGOA"
           ),
         Region = factor(Region, levels = c("BS", "AI", "WGOA", "CGOA", 'EGOA')))

# Compute apportioned total biomass
apport_tot_biom <- apport %>%
  select(year, area, five_yr_avg) %>%
  group_by(year, area) %>%
  mutate(area = case_when(
    area == "Aleutians" ~ "AI",
    area == "Bering Sea" ~ "BS",
    area == "Western Gulf of Alaska" ~ "WGOA",
    area == "Central Gulf of Alaska" ~ "CGOA",
    area == "West Yakutat" ~ "EGOA",
    area == "East Yakutat/Southeast" ~ "EGOA"
  ),
  area = factor(area, levels = c("BS", "AI", "WGOA", "CGOA", "EGOA"))) %>%
  summarize(five_yr_avg = sum(as.numeric(five_yr_avg))) %>%
  left_join(data.frame(biom = as.vector(sgl$rep$Total_Biom), year = years), by = "year") %>%
  mutate(apport_biom = biom * five_yr_avg, Type = "Apportionment") %>%
  rename(Year = year, Region = area) %>%
  select(Year, Region, apport_biom)

compare_tot_biom <- spt_tot_biom %>%
  left_join(apport_tot_biom, by = c("Region", "Year")) %>%
  drop_na()

# Time series comparison
ts_apport <- ggplot() +
  geom_line(compare_tot_biom, mapping = aes(x = Year, y = value, color = 'Spatial'), lwd = 1.3) +
  geom_line(compare_tot_biom, mapping = aes(x = Year, y = apport_biom, color = 'Apportionment'), lwd = 1.3) +
  facet_wrap(~Region, nrow = 1) +
  coord_cartesian(ylim = c(0,NA)) +
  scale_color_manual(values = c("#21918c", "#3b528b")) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.05, 0.85)) +
  labs(x = 'Year', y = 'Total Biomass', color = 'Type')

# Regression comparison
reg_apport <- ggplot(compare_tot_biom, aes(x = value, y = apport_biom)) +
  geom_point() +
  ggpmisc::stat_poly_line() +
  ggpmisc::stat_poly_eq(size = 6, label.y = 0.065, label.x = 0.97) +
  facet_wrap(~Region, nrow = 1, scales = 'free') +
  scale_color_manual(values = c("#21918c", "#3b528b")) +
  theme_bw(base_size = 15) +
  coord_cartesian(ylim = c(0,NA)) +
  labs(y = 'Apportioned Total Biomass', x = 'Spatial Model Total Biomass')

ggsave(
  plot_grid(ts_apport, reg_apport,
            ncol = 1, labels = c("A", "B"), label_size = 25),
  filename = here("Sept PT Model Runs", "25_15_Spatial", "Figures", "Apportionment Comparison (Spt vs. Agg).png"),
  width = 20, height = 10
)

ggsave(
  ts_apport + facet_wrap(~Region, ncol = 1) + theme(legend.position = 'top'),
  filename = here("Sept PT Model Runs", "25_15_Spatial", "Figures", "Apportionment TS Comparison (Spt vs. Agg).png"),
  width = 5, height = 10
)

# Comparisons among Spatial Models ----------------------------------------

# List out model files
spt_models <- list.files(here("Sept PT Model Runs", "25_15_Spatial"))
spt_models <- spt_models[!str_detect(spt_models, "Plot|Figures|NMFS|Francis|scratch")]

# read in models
spt_models_rep <- list()
spt_models_data <- list()
spt_models_sd <- list()

for(i in 1:length(spt_models)) {
  # read in temporary rds
  tmp_rds <- readRDS(here("Sept PT Model Runs", "25_15_Spatial", spt_models[i], paste(spt_models[i], "model_results.RDS", sep = "_")))

  # extract report, data, and sdrep
  spt_models_rep[[i]] <- tmp_rds$rep
  spt_models_data[[i]] <- tmp_rds$data
  spt_models_sd[[i]] <- tmp_rds$sd_rep
}


### Index Comparison --------------------------------------------------------

# get index fits
idx_fits_all <- data.frame()
for (i in 1:length(spt_models_rep)) {
  idx_fits <- get_idx_fits(data = spt_models_data[[i]], rep = spt_models_rep[[i]], year_labs = years) %>%
    plyr::mutate(Model = spt_models[i])
  idx_fits_all <- rbind(idx_fits_all, idx_fits)
}

# remove survey 2
idx_fits_all <- idx_fits_all %>%
  filter(Category != "Survey2, Q1")

ggplot2::ggplot() +
  ggplot2::geom_line(idx_fits_all %>%  dplyr::filter(obs != 0),
                     mapping = ggplot2::aes(x = Year, y = value, color = factor(Model)), lwd = 1.3) +
  ggplot2::geom_pointrange(idx_fits_all %>%  dplyr::filter(obs != 0),
                           mapping = ggplot2::aes(x = Year, y = obs, ymin = lci, ymax = uci), color = "black") +
  ggplot2::labs(x = "Year", y = "Index", color = "Model") +
  theme_sablefish() + ggplot2::coord_cartesian(ylim = c(0, NA)) +
  ggplot2::facet_grid(Category ~ Region, scales = "free_y") +
  ggthemes::scale_color_colorblind()


### Time Series Comparison ----------------------------------------------------------
ts_plots <- get_ts_plot(spt_models_rep, spt_models_sd, spt_models, do_ci = F)

#### Recruitment -------------------------------------------------------------
ts_plots[[3]] +
  scale_x_continuous(labels = seq(min(years), max(years), 15), breaks = seq(1, length(years), 15)) +
  ggthemes::scale_color_colorblind() +
  coord_cartesian(ylim = c(0, 130)) +
  facet_grid(Model~Region,
             labeller = labeller(Region = c("Region 1" = "BS", "Region 2" = "AI", "Region 3" = "WGOA",
                                            "Region 4" = "CGOA", "Region 5" = "EGOA"))) +
  theme_bw(base_size = 13) +
  theme(legend.position = 'none')

#### SSB ---------------------------------------------------------------------
ts_plots[[4]] +
  scale_x_continuous(labels = seq(min(years), max(years), 15), breaks = seq(1, length(years), 15)) +
  ggthemes::scale_color_colorblind() +
  # coord_cartesian(ylim = c(0, 250)) +
  facet_grid(Model~Region,
             labeller = labeller(Region = c("Region 1" = "BS", "Region 2" = "AI", "Region 3" = "WGOA",
                                            "Region 4" = "CGOA", "Region 5" = "EGOA"))) +
  theme_bw(base_size = 13) +
  theme(legend.position = 'none')


#### Total Biomass ---------------------------------------------------------------------
ts_plots[[5]] +
  scale_x_continuous(labels = seq(min(years), max(years), 15), breaks = seq(1, length(years), 15)) +
  ggthemes::scale_color_colorblind() +
  coord_cartesian(ylim = c(0, 250)) +
  facet_grid(Model~Region,
             labeller = labeller(Region = c("Region 1" = "BS", "Region 2" = "AI", "Region 3" = "WGOA",
                                            "Region 4" = "CGOA", "Region 5" = "EGOA"))) +
  theme_bw(base_size = 13) +
  theme(legend.position = 'none')

