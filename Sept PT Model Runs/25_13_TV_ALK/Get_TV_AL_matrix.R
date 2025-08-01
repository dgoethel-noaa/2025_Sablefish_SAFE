# Purpose: To construct an age-length transition matrix that is time-varying using estiamtes from Cheng et al. 2024
# Creator: Matthew LH. Cheng
# 7/17/25

# Function to construct age length transition
get_al_trans_matrix_ = function(age_bins, len_bins, mean_length, sd) {
  # container array
  age_length = matrix(0.0, nrow = length(age_bins), ncol = length(len_bins))

  for(a in 1:length(age_bins)) {
    # Use actual bin lower limits as per specification
    # Assume len_bins contains the lower limits of each bin
    bin_lower_limits = len_bins

    # Calculate cumulative probabilities at bin lower limits
    AL = pnorm(bin_lower_limits, mean_length[a], sd[a])

    # Calculate bin probabilities according to the specification
    for(l in 1:length(len_bins)) {
      if(l == 1) {
        # First bin: from -infinity to lower limit of bin 2
        age_length[a, l] = pnorm(bin_lower_limits[2], mean_length[a], sd[a])
      } else if(l == length(len_bins)) {
        # Last bin: from lower limit of last bin to +infinity
        age_length[a, l] = 1 - pnorm(bin_lower_limits[l], mean_length[a], sd[a])
      } else {
        # Middle bins: from lower limit of bin l to lower limit of bin l+1
        age_length[a, l] = pnorm(bin_lower_limits[l+1], mean_length[a], sd[a]) -
          pnorm(bin_lower_limits[l], mean_length[a], sd[a])
      }
    }
  }
  return(age_length)
}

# note: estimates start at 1996; using everything from SAFE doc prior to 1996 for consistency

# setup -------------------------------------------------------------------

library(here)
library(tidyverse)

# Read in base data file
base_data <- read_rds(here("Sept PT Model Runs", "Data Files", "SPoRC_one_reg.RDS"))
# setup dimensions
len_bins <- seq(41, 99, 2)
age_bins <- 2:31
n_regions <- 1
n_years <- base_data$nyrs
years <- 1960:2024
srv_years <- length(1960:1996):length(1960:2023)
srv_year_st <- length(1960:1996) - 1

# Create container for WAA
waa_cont <- array(0, dim = dim(base_data$WAA))

# Create container for age length transition
age_len_cont <- array(0, dim = dim(base_data$SizeAgeTrans))


# Populate WAA ------------------------------------------------------------

# Populate WAA with legacy values prior to 1996
waa_cont[,1:length(1960:1995),,] <- base_data$WAA[,1:length(1960:1995),,]

# Extract out quantities for females
growth_f <- read.csv(here("Sept PT Model Runs", "25_13_TV_ALK", "Weight_estimates.csv")) %>%
  filter(re_model == '3dar1', Growth_Model == 'weight', Sex == 'female', peel == 0)

for(i in srv_years) {
  tmp_wt <- growth_f %>% filter(Year == years[i]) # extract mean length
  waa_cont[,i,,1] <- tmp_wt$Mean
}

# Extract out quantities for males
growth_m <- read.csv(here("Sept PT Model Runs", "25_13_TV_ALK", "Weight_estimates.csv")) %>%
  filter(re_model == '3dar1', Growth_Model == 'weight', Sex == 'Male', peel == 0)

for(i in srv_years) {
  tmp_wt <- growth_m %>% filter(Year == years[i]) # extract mean length
  waa_cont[,i,,2] <- tmp_wt$Mean
}

# Fill in terminal year with last year with data
waa_cont[,length(years),,1] <- waa_cont[,length(years)-1,,1]
waa_cont[,length(years),,2] <- waa_cont[,length(years)-1,,2]

write_rds(waa_cont, here("Sept PT Model Runs", "Data Files", "TV_WAA_matrix.rds"))

# Populate AL Matrix ---------------------------------------------------------
# Populate age length transition with legacy values prior to 1996
age_len_cont[,1:length(1960:1995),,,] <- base_data$SizeAgeTrans[,1:length(1960:1995),,,]

# Extract out quantities - Females
load(here("Sept PT Model Runs", "25_13_TV_ALK", paste("length_female_3DAR1_Model.RData", sep = "")))
growth_f <- read.csv(here("Sept PT Model Runs", "25_13_TV_ALK", "Length_estimates.csv")) %>%
  filter(re_model == '3dar1', Growth_Model == 'length', Sex == 'female', peel == 0)

for(i in srv_years) {
  tmp_len <- growth_f %>% filter(Year == years[i]) # extract mean length
  age_len_cont[,i,,,1] <-
    t(get_al_trans_matrix_(age_bins, len_bins,
                           mean_length = tmp_len$Mean,
                           sd = model$rep$obs_sigma_at[,i - srv_year_st]))
}

# Extract out quantities - Males
load(here("Sept PT Model Runs", "25_13_TV_ALK", paste("length_male_3DAR1_Model.RData", sep = "")))
growth_m <- read.csv(here("Sept PT Model Runs", "25_13_TV_ALK", "Length_estimates.csv")) %>%
  filter(re_model == '3dar1', Growth_Model == 'length', Sex == 'Male', peel == 0)

for(i in srv_years) {
  tmp_len <- growth_m %>% filter(Year == years[i]) # extract mean length
  age_len_cont[,i,,,2] <-
    t(get_al_trans_matrix_(age_bins, len_bins,
                           mean_length = tmp_len$Mean,
                           sd = model$rep$obs_sigma_at[,i - srv_year_st]))

}

# Fill in terminal year with last year with data
age_len_cont[,length(years),,,1] <- age_len_cont[,length(years)-1,,,1]
age_len_cont[,length(years),,,2] <- age_len_cont[,length(years)-1,,,2]

write_rds(age_len_cont, here("Sept PT Model Runs", "Data Files", "TV_AL_matrix.rds"))

# Plot to visualize
agelen_df <- reshape2::melt(age_len_cont) %>%
  rename(Region = Var1, Year = Var2,
         Len = Var3, Age = Var4,
         Sex = Var5)

pdf(here("Sept PT Model Runs", "25_13_TV_ALK", "Growth_Summary.pdf"), width = 13, height = 13)
# Female LAA
print(
  ggplot(agelen_df %>% filter(Sex == 1), aes(x = Age, y = Len, fill = value)) +
    geom_tile() +
    scale_fill_distiller(palette = "YlGnBu", direction = 1) +
    facet_wrap(~Year) +
    labs(x = "Age", y = "Length (cm)", fill = "Probability", title = 'Female') +
    theme_bw(base_size = 15)
)

# Male LAA
print(
  ggplot(agelen_df %>% filter(Sex == 2), aes(x = Age, y = Len, fill = value)) +
    geom_tile() +
    scale_fill_distiller(palette = "YlGnBu", direction = 1) +
    facet_wrap(~Year) +
    labs(x = "Age", y = "Length (cm)", fill = "Probability", title = 'Male') +
    theme_bw(base_size = 15)
)

# Plot LAA
laa <- read.csv(here("Sept PT Model Runs", "25_13_TV_ALK", "Length_estimates.csv")) %>%
  filter(re_model == '3dar1', Growth_Model == 'length', peel == 0)

print(
  ggplot(laa %>% filter(Sex == 'female'), aes(x = Age, y = Mean, color = Year, group = Year)) +
    geom_line(lwd = 1.3) +
    scale_color_viridis_c() +
    facet_wrap(~Year) +
    labs(x = 'Year', y = 'Length (cm)', title = 'Female') +
    theme_bw(base_size = 15)
)

print(
  ggplot(laa %>% filter(Sex == 'female'), aes(x = Age, y = Mean, color = Year, group = Year)) +
    geom_line(lwd = 1.3) +
    scale_color_viridis_c() +
    labs(x = 'Year', y = 'Length (cm)', title = 'Female') +
    theme_bw(base_size = 15)
)

print(
  ggplot(laa %>% filter(Sex == 'Male'), aes(x = Age, y = Mean, color = Year, group = Year)) +
    geom_line(lwd = 1.3) +
    scale_color_viridis_c() +
    labs(x = 'Year', y = 'Length (cm)', title = 'Male') +
    theme_bw(base_size = 15)
)

print(
  ggplot(laa %>% filter(Sex == 'Male'), aes(x = Age, y = Mean, color = Year, group = Year)) +
    geom_line(lwd = 1.3) +
    scale_color_viridis_c() +
    facet_wrap(~Year) +
    labs(x = 'Year', y = 'Length (cm)', title = 'Male') +
    theme_bw(base_size = 15)
)

# Plot WAA
waa <- read.csv(here("Sept PT Model Runs", "25_13_TV_ALK", "weight_estimates.csv")) %>%
  filter(re_model == '3dar1', Growth_Model == 'weight', peel == 0)

print(
  ggplot(waa %>% filter(Sex == 'female'), aes(x = Age, y = Mean, color = Year, group = Year)) +
    geom_line(lwd = 1.3) +
    scale_color_viridis_c() +
    facet_wrap(~Year) +
    labs(x = 'Year', y = 'Weight (kg)', title = 'Female') +
    theme_bw(base_size = 15)
)

print(
  ggplot(waa %>% filter(Sex == 'female'), aes(x = Age, y = Mean, color = Year, group = Year)) +
    geom_line(lwd = 1.3) +
    scale_color_viridis_c() +
    labs(x = 'Year', y = 'Weight (kg)', title = 'Female') +
    theme_bw(base_size = 15)
)

print(
  ggplot(waa %>% filter(Sex == 'Male'), aes(x = Age, y = Mean, color = Year, group = Year)) +
    geom_line(lwd = 1.3) +
    scale_color_viridis_c() +
    labs(x = 'Year', y = 'Weight (kg)', title = 'Male') +
    theme_bw(base_size = 15)
)

print(
  ggplot(waa %>% filter(Sex == 'Male'), aes(x = Age, y = Mean, color = Year, group = Year)) +
    geom_line(lwd = 1.3) +
    scale_color_viridis_c() +
    facet_wrap(~Year) +
    labs(x = 'Year', y = 'Weight (kg)', title = 'Male') +
    theme_bw(base_size = 15)
)
dev.off()

