# Purpose: Spatial Sablefish Sensitivity Run - Multinomial Release Conditioned
# Creator: Matthew LH. Cheng
# Date Created: 7/17/25

# Set up ------------------------------------------------------------------
library(here)
library(tidyverse)
library(RTMB)
library(SPoRC)

# Setup folder
root_folder <- here("Sept PT Model Runs", "SR_Spatial")
mod_name <- "Spatial_MltRel"

# read in 2021 assessment data
data("mlt_rg_sable_data")

# read in 2025 assessment data
dat_2025 <- readRDS(here("Data Pulls", "SPoRC_spatial.rds"))

# Setup Model -------------------------------------------------------------

# Initialize model dimensions and data list
input_list <- Setup_Mod_Dim(years = 1:length(dat_2025$yrs),
                            # vector of years (1 - 62)
                            ages = 1:length(dat_2025$ages),
                            # vector of ages (1 - 30)
                            lens = dat_2025$length_bins,
                            # number of lengths (41 - 99)
                            n_regions = dat_2025$nreg,
                            # number of regions (5)
                            n_sexes = dat_2025$nsex,
                            # number of sexes (2)
                            n_fish_fleets = dat_2025$nflt,
                            # number of fishery fleet (2)
                            n_srv_fleets = dat_2025$nsrv_flt,
                            # number of survey fleets (3)
                            verbose = TRUE
)

# Setup recruitment stuff (using defaults for other stuff)
input_list <- Setup_Mod_Rec(input_list = input_list, # input data list from above
                            do_rec_bias_ramp = 0,
                            sigmaR_switch = as.integer(length(dat_2025$styr:1975)),
                            dont_est_recdev_last = 1, # don't estimate last rec dev
                            sexratio = c(0.5, 0.5), # fix sex ratio at 0.5

                            # Model options
                            rec_model = "mean_rec", # recruitment model
                            sigmaR_spec = "fix", # fixing
                            InitDevs_spec = "est_shared_r",
                            # initial deviations are shared across regions,
                            # but recruitment deviations are region specific
                            ln_sigmaR = log(c(0.4, 0.9)),
                            # values to fix sigmaR at, or starting values
                            ln_global_R0 = log(25),
                            Use_Rec_prop_Prior = 1,
                            Rec_prop_prior = 5
)

# Setup biological stuff (using defaults for other stuff)
input_list <- Setup_Mod_Biologicals(input_list = input_list,
                                    WAA = dat_2025$WAA, # weight at age
                                    MatAA = dat_2025$MatAA, # maturity at age
                                    AgeingError = dat_2025$AgeError,
                                    # ageing error matrix
                                    fit_lengths = 1, # fitting lengths
                                    SizeAgeTrans = dat_2025$SizeAgeTrans,
                                    # size age transition matrix
                                    M_spec = "fix", # fix natural mortality
                                    # values to fix natural mortality at
                                    Fixed_natmort = array(0.0988975, dim = c(input_list$data$n_regions,
                                                                            length(input_list$data$years),
                                                                            length(input_list$data$ages),
                                                                            input_list$data$n_sexes))
)

# setting up movement parameterization
input_list <- Setup_Mod_Movement(input_list = input_list,
                                 # Model options
                                 Movement_ageblk_spec = list(c(1:6), c(7:15), c(16:30)),
                                 # estimating movement in 3 age blocks
                                 # (ages 1-6, ages 7-15, ages 16-30)
                                 Movement_yearblk_spec = "constant", # time-invariant movement
                                 Movement_sexblk_spec = "constant", # sex-invariant movement
                                 do_recruits_move = 0, # recruits do not move
                                 use_fixed_movement = 0, # estimating movement
                                 Use_Movement_Prior = 1, # priors used for movement
                                 Movement_prior = 2 # vague prior to penalize movement away from the extremes
)

# setting up tagging parameterization

# setup tagging priors
tag_prior <- data.frame(
  region = 1,
  block = c(1,2,3),
  mu = NA, # no mean, since symmetric beta
  sd = 5, # sd = 5
  type = 0 # symmetric beta
)

input_list <- Setup_Mod_Tagging(input_list = input_list,
                                UseTagging = 1, # using tagging data
                                max_tag_liberty = 15, # maximum number of years to track a cohort

                                # Data Inputs
                                tag_release_indicator = mlt_rg_sable_data$tag_release_indicator,
                                # tag release indicator (first col = tag region,
                                # second col = tag year),
                                # total number of rows = number of tagged cohorts
                                Tagged_Fish = mlt_rg_sable_data$Tagged_Fish, # Released fish
                                # dimensioned by total number of tagged cohorts, (implicitly
                                # tracks the release year and region), age, and sex
                                Obs_Tag_Recap = mlt_rg_sable_data$Obs_Tag_Recap,
                                # dimensioned by max tag liberty, tagged cohorts, regions,
                                # ages, and sexes

                                # Model options
                                Tag_LikeType = "Multinomial_Release", # Multinomial Release Conditioned
                                mixing_period = 2, # Don't fit tagging until release year + 1
                                t_tagging = 0.5, # tagging happens midway through the year,
                                # movement does not occur within that year
                                tag_selex = "SexSp_DomFleet", # tagging recapture selectivity is the dominant fleet (fixed-gear)
                                tag_natmort = "AgeSp_SexSp", # tagging natural mortality is age and sex-specific
                                Use_TagRep_Prior = 1, # tag reporting rate priors are used
                                TagRep_Prior = tag_prior,
                                move_age_tag_pool = as.list(1:30), # whether or not to pool tagging data when fitting (for computational cost)
                                move_sex_tag_pool = list(c(1:2)), # whether or not to pool sex-specific data whezn fitting
                                Init_Tag_Mort_spec = "fix", # fixing initial tag mortality
                                Tag_Shed_spec = "fix", # fixing chronic shedding
                                TagRep_spec = "est_shared_r", # tag reporting rates are not region specific
                                Tag_Reporting_blocks = c(
                                  paste("Block_1_Year_1-36_Region_", c(1:input_list$data$n_regions), sep = ''),
                                  paste("Block_2_Year_36-56_Region_", c(1:input_list$data$n_regions), sep = ''),
                                  paste("Block_3_Year_57-terminal_Region_", c(1:input_list$data$n_regions), sep = '')
                                ),
                                # Specify starting values or fixing values
                                ln_Init_Tag_Mort = log(0.1), # fixing initial tag mortality
                                ln_Tag_Shed = log(0.02),  # fixing tag shedding
                                ln_tag_theta = log(0.5), # starting value for tagging overdispersion
                                Tag_Reporting_Pars = array(log(0.2 / (1-0.2)), dim = c(input_list$data$n_regions, 3))  # starting values for tag reporting pars
)

# setting up catch data
dat_2025$UseCatch[is.na(dat_2025$ObsCatch)] <- 0 # replace w/ 0s

input_list <- Setup_Mod_Catch_and_F(input_list = input_list,
                                    # Data inputs
                                    ObsCatch = dat_2025$ObsCatch,
                                    Catch_Type = array(1, dim = c(length(input_list$data$years), input_list$data$n_fish_fleets)),
                                    UseCatch = dat_2025$UseCatch,
                                    # Model options
                                    Use_F_pen = 1,
                                    # whether to use f penalty, == 0 don't use, == 1 use
                                    sigmaC_spec = 'fix',
                                    ln_sigmaC =
                                      array(log(0.05), dim = c(input_list$data$n_regions, input_list$data$n_fish_fleets)),
                                    # fixing catch sd at small value
                                    ln_F_mean = array(-2, dim = c(input_list$data$n_regions,
                                                                  input_list$data$n_fish_fleets))
                                    # some starting values for fishing mortality
)

# Fishery Indices and Compositions

# Make sure the use indicators for compositions are 0s when NAs
dat_2025$UseFishAgeComps[is.na(apply(dat_2025$ObsFishAgeComps, c(1,2,5), sum))] <- 0
dat_2025$UseFishLenComps[is.na(apply(dat_2025$ObsFishLenComps, c(1,2,5), sum))] <- 0

input_list <- Setup_Mod_FishIdx_and_Comps(input_list = input_list,
                                          # data inputs
                                          ObsFishIdx = dat_2025$ObsFishIdx,
                                          ObsFishIdx_SE = dat_2025$ObsFishIdx_SE,
                                          UseFishIdx =  dat_2025$UseFishIdx,
                                          ObsFishAgeComps = dat_2025$ObsFishAgeComps,
                                          UseFishAgeComps = dat_2025$UseFishAgeComps,
                                          ISS_FishAgeComps = dat_2025$ISS_FishAgeComps,
                                          ObsFishLenComps = dat_2025$ObsFishLenComps,
                                          UseFishLenComps = dat_2025$UseFishLenComps,
                                          ISS_FishLenComps = dat_2025$ISS_FishLenComps,

                                          # Model options
                                          fish_idx_type = c("none", "none"),
                                          # fishery indices not used
                                          FishAgeComps_LikeType =
                                            c("Multinomial", "none"),
                                          # age comp likelihoods for fishery fleet 1 and 2
                                          FishLenComps_LikeType =
                                            c("Multinomial", "Multinomial"),
                                          # length comp likelihoods for fishery fleet 1 and 2
                                          FishAgeComps_Type =
                                            c("spltRjntS_Year_1-terminal_Fleet_1",
                                              "none_Year_1-terminal_Fleet_2"),
                                          # age comp structure for fishery fleet 1 and 2
                                          FishLenComps_Type =
                                            c("spltRjntS_Year_1-terminal_Fleet_1",
                                              "spltRjntS_Year_1-terminal_Fleet_2"),
                                          # length comp structure for fishery fleet 1 and 2
                                          FishAge_comp_agg_type = c(NA,NA),
                                          # ADMB aggregation quirks, ideally get rid of this
                                          FishLen_comp_agg_type = c(NA,NA)
                                          # ADMB aggregation quirks, ideally get rid of this
)

# Survey Indices and Compositions
# Note: These are now sds (i.e., the tmb likelihoods) for the domestic and JP LLS, but not for the trawl survey (TS is still CV * index).
# Not a big issue because the trawl survey is not used here.

# Make sure the use indicators for compositions are 0s when NAs
dat_2025$UseSrvAgeComps[is.na(apply(dat_2025$ObsSrvAgeComps, c(1,2,5), sum))] <- 0
dat_2025$UseSrvLenComps[is.na(apply(dat_2025$ObsSrvLenComps, c(1,2,5), sum))] <- 0

# Don't fit trawl survey
dat_2025$UseSrvLenComps[,,2] <- 0
dat_2025$UseSrvIdx[,,2] <- 0

input_list <- Setup_Mod_SrvIdx_and_Comps(input_list = input_list,
                                         # data inputs
                                         ObsSrvIdx = dat_2025$ObsSrvIdx,
                                         ObsSrvIdx_SE = dat_2025$ObsSrvIdx_SE * 2,
                                         UseSrvIdx =  dat_2025$UseSrvIdx,
                                         ObsSrvAgeComps = dat_2025$ObsSrvAgeComps,
                                         ISS_SrvAgeComps = dat_2025$ISS_SrvAgeComps,
                                         UseSrvAgeComps = dat_2025$UseSrvAgeComps,
                                         ObsSrvLenComps = dat_2025$ObsSrvLenComps,
                                         UseSrvLenComps = dat_2025$UseSrvLenComps,
                                         ISS_SrvLenComps = dat_2025$ISS_SrvLenComps,

                                         # Model options
                                         srv_idx_type = c("abd", 'biom', "abd"),
                                         # abundance and biomass for survey fleet 1 and 2
                                         SrvAgeComps_LikeType =
                                           c("Multinomial", 'none', "Multinomial"),
                                         # survey age composition likelihood for survey fleet
                                         # 1, and 2
                                         SrvLenComps_LikeType =
                                           c("none", "none", "none"),
                                         #  no length compositions used for survey
                                         SrvAgeComps_Type = c("spltRjntS_Year_1-terminal_Fleet_1",
                                                              "none_Year_1-terminal_Fleet_2",
                                                              "spltRjntS_Year_1-terminal_Fleet_3"),
                                         # survey age comp type
                                         SrvLenComps_Type = c("none_Year_1-terminal_Fleet_1",
                                                              "none_Year_1-terminal_Fleet_2",
                                                              "none_Year_1-terminal_Fleet_3"),
                                         # survey length comp type
                                         SrvAge_comp_agg_type = c(NA, NA, NA),
                                         # ADMB aggregation quirks, ideally get rid of this
                                         SrvLen_comp_agg_type = c(NA, NA, NA)
                                         # ADMB aggregation quirks, ideally get rid of this
)

# Fishery Selectivity and Catchability

# defining priors
sex_par <- expand.grid(sex = 1:2, par = 1:2)
fleet_blocks <- data.frame(
  fleet = c(1, 1, 1, 2),
  block = c(1, 2, 3, 1)
)

# merge together (note that unlike the operational assessment, selectivity
# blocks are reduced from 3 to 2)
fish_selex_structure <- merge(fleet_blocks, sex_par)

# Merge to get all valid combinations
fish_selex_structure <- merge(fleet_blocks, sex_par) %>%
  dplyr::filter(!(fleet == 1 & block == 1 & sex == 2 & par == 2)) %>%              # remove priors for any unestimated pars -- par1=a50, par2=delta; NEEDS TO MATCH PARAMETER MAPPING
  dplyr::filter(!(fleet == 2 & block == 1 & sex == 2 & par == 1))                  # remove priors for any unestimated pars -- par1=a50, par2=delta; NEEDS TO MATCH PARAMETER MAPPING

# Add the lognormal prior values - creates a dataframe, each row is a unique parameter combination to apply the prior to
fish_selex_prior <- cbind(
  region = 1,
  fish_selex_structure,
  mu = 1,                                                                      # All selex means = 1 (means should be defined in normal space)
  sd = 5                                                                       # All selex sd = 5
)

fish_selex_prior_tf <- fish_selex_prior %>%                                    # set tighter selex prior for TF
  dplyr::filter((fleet == 2 & par == 1)) %>%
  dplyr::mutate(mu = 2, sd = 1) %>%
  dplyr::full_join(fish_selex_prior %>%  dplyr::filter(!(fleet == 2 & par == 1 )))

fish_selex_prior_tf <- fish_selex_prior_tf %>%                                    # set tighter selex prior for TF
  dplyr::filter((fleet == 2 & par == 2)) %>%
  dplyr::mutate(mu = 1, sd = 2) %>%
  dplyr::full_join(fish_selex_prior_tf %>%  dplyr::filter(!(fleet == 2 & par == 2)))

input_list <- Setup_Mod_Fishsel_and_Q(input_list = input_list,

                                      # Model options
                                      cont_tv_fish_sel = c("none_Fleet_1", "none_Fleet_2"),
                                      # fishery selectivity, whether continuous time-varying

                                      # fishery selectivity blocks
                                      fish_sel_blocks =                        # fishery selectivity time blocks if not TV specified above for a given fleet
                                        c("Block_1_Year_1-35_Fleet_1",         # pre-IFQ time block for fixed gear fishery 1994 and before
                                          "Block_2_Year_36-56_Fleet_1",        # IFQ time block for fixed gear fishery-- 1995 to 2015
                                          "Block_3_Year_57-terminal_Fleet_1",  # Recent time block for fixed gear fishery--2016 to terminal year
                                          "none_Fleet_2"),
                                      # no blocks for trawl fishery

                                      # fishery selectivity form
                                      fish_sel_model =
                                        c("logist1_Fleet_1", "gamma_Fleet_2"),

                                      # fishery catchability blocks
                                      fish_q_blocks =
                                        c("none_Fleet_1", "none_Fleet_2"),
                                      # no blocks since q is not estimated

                                      # sharing fishery selex parameters
                                      fish_fixed_sel_pars =
                                        c("est_shared_r", "est_shared_r"),

                                      # whether to estimate all fixed effects
                                      # for fishery catchability
                                      fish_q_spec =
                                        c("fix", "fix"),
                                      Use_fish_selex_prior = 1,
                                      fish_selex_prior = fish_selex_prior
)


# setup survey selectivity
# Define sex and parameter combinations
sex_par <- expand.grid(sex = 1:2, par = 1:2)

# Define valid fleet-block combinations (only estimating domestic and jp LLS)
fleet_blocks <- data.frame(
  fleet = c(1, 1, 3),
  block = c(1, 2, 1)
)

# Merge to get all valid combinations
srv_selex_structure <- merge(fleet_blocks, sex_par)

# Add the lognormal prior values - creates a dataframe, each row is a unique parameter combination to apply the prior to
srv_selex_prior <- cbind(
  region = 1,
  srv_selex_structure,
  mu = 1,
  sd = 5
) %>%
  filter(!(fleet == 3 & par == 2 & sex == 2))

input_list <- Setup_Mod_Srvsel_and_Q(input_list = input_list,

                                     # Model options
                                     # survey selectivity, whether continuous time-varying
                                     cont_tv_srv_sel =
                                       c("none_Fleet_1",
                                         "none_Fleet_2",
                                         "none_Fleet_3"),

                                     # survey selectivity blocks
                                     srv_sel_blocks =                          # survey selectivity time blocks if not TV specified above for a given fleet
                                       c("Block_1_Year_1-56_Fleet_1",          # Early time block LLS-- 1960 to 2016
                                         "Block_2_Year_57-terminal_Fleet_1",   # Recent time block for LLS--2017 to terminal year
                                         "none_Fleet_2",                       # No blocks for trawl survey
                                         "none_Fleet_3"                        # No blocks for JPN LLS
                                       ),

                                     # survey selectivity form
                                     srv_sel_model =
                                       c("logist1_Fleet_1",
                                         "exponential_Fleet_2",
                                         "logist1_Fleet_3"
                                       ),

                                     # survey catchability blocks
                                     srv_q_blocks =
                                       c("none_Fleet_1",
                                         "none_Fleet_2",
                                         "none_Fleet_3"),

                                     # whether to estiamte all fixed effects
                                     # for survey selectivity and later
                                     # modify to fix/share parameters
                                     srv_fixed_sel_pars_spec =
                                       c("est_shared_r",
                                         "fix",
                                         "est_shared_r"),

                                     # whether to estiamte all
                                     # fixed effects for survey catchability
                                     # spatially-invariant q
                                     srv_q_spec =
                                       c("est_shared_r",
                                         "fix",
                                         "est_shared_r"),
                                     Use_srv_selex_prior = 1,
                                     srv_selex_prior = srv_selex_prior
)

# set up model weighting stuff
input_list <- Setup_Mod_Weighting(input_list = input_list,
                                  sablefish_ADMB = 0,
                                  # don't use sablefish single region ADMB quirks
                                  likelihoods = 1, # using TMB-style likelihoods,
                                  # and weight using sigmas, instead of lambdas
                                  # and sigmas together
                                  Wt_Catch = 1,
                                  Wt_FishIdx = 1,
                                  Wt_SrvIdx = 1,
                                  Wt_Rec = 1,
                                  Wt_F = 1,
                                  Wt_Tagging = 0.5,
                                  # Composition model weighting
                                  Wt_FishAgeComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_fish_fleets)),
                                  Wt_FishLenComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_fish_fleets)),
                                  Wt_SrvAgeComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_srv_fleets)),
                                  Wt_SrvLenComps =
                                    array(1, dim = c(input_list$data$n_regions,
                                                     length(input_list$data$years),
                                                     input_list$data$n_sexes,
                                                     input_list$data$n_srv_fleets))
)

# extract out lists updated with helper functions
data <- input_list$data
parameters <- input_list$par
mapping <- input_list$map


# Additional Model Specifications -----------------------------------------
# Change up ISS
data$ISS_SrvAgeComps[] <- 60
data$ISS_FishAgeComps[] <- 40
data$ISS_FishLenComps[] <- 20
data$ISS_SrvLenComps[] <- 20

# Map off early delta for fishery
map_fish_fixed <- array(mapping$ln_fish_fixed_sel_pars, dim = dim(parameters$ln_fish_fixed_sel_pars))
map_fish_fixed[,2,1,2,1]  <- map_fish_fixed[,2,1,1,1] # share deltas

# Map off bmax for trawl females
map_fish_fixed[,1,1,2,2]  <- map_fish_fixed[,1,1,1,2] # share deltas
mapping$ln_fish_fixed_sel_pars <- factor(map_fish_fixed)

# Map off delta for JP LLS
map_srv_fixed <- array(mapping$ln_srv_fixed_sel_pars, dim = dim(parameters$ln_srv_fixed_sel_pars))
map_srv_fixed[,2,1,2,3]  <- map_srv_fixed[,2,1,1,3] # share deltas
mapping$ln_srv_fixed_sel_pars <- factor(map_srv_fixed)

# Some starting values to help out the model
parameters$ln_srv_fixed_sel_pars[] <- log(3)
parameters$ln_fish_fixed_sel_pars[] <- log(3)

# Fit Model ---------------------------------------------------------------
sabie_rtmb_model <- fit_model(data,
                              parameters,
                              mapping,
                              random = NULL,
                              newton_loops = 3,
                              silent = T
)

reshape2::melt(sabie_rtmb_model$rep$Rec) %>%
  ggplot(aes(x = Var2, y = value)) +
  geom_line() +
  facet_wrap(~Var1)

plot(colSums(sabie_rtmb_model$rep$SSB))

get_idx_fits_plot(list(data), list(sabie_rtmb_model$rep), 1)
get_selex_plot(list(sabie_rtmb_model$rep), 1)


# Get standard error report
sabie_rtmb_model$sd_rep <- RTMB::sdreport(sabie_rtmb_model)

# Save other outputs
sabie_rtmb_model$data <- data
sabie_rtmb_model$mapping <- mapping
sabie_rtmb_model$parameters <- parameters

# Write out RDS file
saveRDS(sabie_rtmb_model, file = here(root_folder, mod_name,paste0(mod_name,"_model_results.RDS")))
