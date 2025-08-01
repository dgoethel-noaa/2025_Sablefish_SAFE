# Purpose: To plot spatial compositions over time
# Creator: Matthew LH. Cheng
# Date 7/31/25


# Setup -------------------------------------------------------------------

library(here)
library(SPoRC)

# read in model / data
spt <- readRDS(here("Sept PT Model Runs", "25_15_Spatial", "Spatial_MltRel", "Spatial_MltRel_model_results.RDS")) # "final" spatial model

# region mapping
region_map <- data.frame(Region = 1:5, Name = c("BS", "AI", "WGOA", "CGOA", "EGOA"))

# Compositions ------------------------------------------------------------
comp_df <- get_comp_prop(spt$data, spt$rep, age_labels = 2:31, len_labels = seq(41,99,2), year_labels = 1960:2024)


### Fishery -----------------------------------------------------------------

pdf(here("Sept PT Model Runs", "25_15_Spatial", "Figures", "Fishery Ages.pdf"), width = 15, height = 15)
for(i in 1:nrow(region_map)) {

  # get regional comps
  tmp_comps <- comp_df$Fishery_Ages %>% filter(Region == i)

  # adjust comps to get females above, males below
  tmp_comps$obs_adj <- ifelse(tmp_comps$Sex == 2, -tmp_comps$obs, tmp_comps$obs)
  tmp_comps$Sex_Name <- ifelse(tmp_comps$Sex == 2, 'M', 'F')

  # plot
  print(ggplot(tmp_comps, mapping = aes(x = Age, y = obs_adj, fill = Sex_Name)) +
          geom_col(alpha = 0.5, color = 'black') +
          facet_wrap(~Year) +
          theme_sablefish() +
          labs(x = 'Age', y = "Proportions", fill = 'Sex', title = paste(region_map$Name[i], "Fishery Ages"))
        )
}
dev.off()

pdf(here("Sept PT Model Runs", "25_15_Spatial", "Figures", "Fishery Lengths"), width = 15, height = 15)
for(i in 1:nrow(region_map)) {

  # get regional comps
  tmp_comps <- comp_df$Fishery_Lens %>% filter(Region == i, Fleet == 2)

  # adjust comps to get females above, males below
  tmp_comps$obs_adj <- ifelse(tmp_comps$Sex == 2, -tmp_comps$obs, tmp_comps$obs)
  tmp_comps$Sex_Name <- ifelse(tmp_comps$Sex == 2, 'M', 'F')

  # plot
  print(ggplot(tmp_comps, mapping = aes(x = Age, y = obs_adj, fill = Sex_Name)) +
          geom_col(alpha = 0.5, color = 'black') +
          facet_wrap(~Year) +
          theme_sablefish() +
          labs(x = 'Age', y = "Proportions", fill = 'Sex', title = paste(region_map$Name[i], "Trawl Fishery Ages"))
  )
}
dev.off()


### Survey ------------------------------------------------------------------

pdf(here("Sept PT Model Runs", "25_15_Spatial", "Figures", "Survey Ages.pdf"), width = 15, height = 15)
for(i in 1:nrow(region_map)) {

  # get regional comps
  tmp_comps <- comp_df$Survey_Ages %>% filter(Region == i)

  # adjust comps to get females above, males below
  tmp_comps$obs_adj <- ifelse(tmp_comps$Sex == 2, -tmp_comps$obs, tmp_comps$obs)
  tmp_comps$Sex_Name <- ifelse(tmp_comps$Sex == 2, 'M', 'F')

  # plot
  print(ggplot(tmp_comps, mapping = aes(x = Age, y = obs_adj, fill = Sex_Name)) +
          geom_col(alpha = 0.5, color = 'black') +
          facet_wrap(~Year) +
          theme_sablefish() +
          labs(x = 'Age', y = "Proportions", fill = 'Sex', title = paste(region_map$Name[i], "Survey Ages"))
  )
}
dev.off()

