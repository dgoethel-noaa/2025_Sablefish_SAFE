# Purpose: To create an RTMB model for the 2024 sablefish assessment using SPoRC
# Creator: Matthew LH. Cheng (UAF CFOS)
# Date 5/28/25

# Setup -------------------------------------------------------------------
# devtools::install_github("chengmatt/SPoRC", dependencies = TRUE, lib = Sys.getenv("R_LIBS_USER"))
library(here)
library(SPoRC)
library(ggplot2)
library(RTMB)

# Load in 2024 assessment data
data("sgl_rg_sable_data")
years <- 1960:2024 # model years

# Prepare Model Inputs ----------------------------------------------------

# Initialize model dimensions and data list
input_list <- Setup_Mod_Dim(years = 1:length(years), # vector of years
                            ages = 1:30, # vector of ages
                            lens = seq(41,99,2), # number of lengths
                            n_regions = 1, # number of regions
                            n_sexes = 2, # number of sexes == 1, female, == 2 male
                            n_fish_fleets = 2, # number of fishery fleet == 1, fixed gear, == 2 trawl gear
                            n_srv_fleets = 3, # number of survey fleets
                            verbose = TRUE
)

# Setup recruitment stuff (using defaults for other stuff)
input_list <- Setup_Mod_Rec(input_list = input_list, # input data list from above
                            # Model options
                            do_rec_bias_ramp = 1, # do bias ramp (0 == don't do bias ramp, 1 == do bias ramp)
                            # breakpoints for bias ramp (1 == no bias ramp - 1960 - 1980, 2 == ascending limb of bias ramp - 1980 - 1990,
                            # 3 == full bias correction - 1990 - 2022, == 4 no bias correction - terminal year of recruitment estimate)
                            bias_year = c(length(1960:1979), length(1960:1989), (length(1960:2023) - 5), length(1960:2024) - 2) + 1,
                            sigmaR_switch = as.integer(length(1960:1975)), # when to switch from early to late sigmaR
                            dont_est_recdev_last = 1, # don't estimate last recruitment deviate
                            ln_sigmaR = log(c(0.4, 1.2)),
                            rec_model = "mean_rec", # recruitment model
                            sigmaR_spec = "fix_early_est_late", # fix early sigmaR, estiamte late sigmaR
                            InitDevs_spec = NULL, # estimate all initial deviations
                            RecDevs_spec = NULL, # stiamte all recruitment deivations
                            sexratio = as.vector(c(0.5, 0.5)), # recruitment sex ratio
                            init_age_strc = 1,
                            init_F_prop = 0.1
)

input_list <- Setup_Mod_Biologicals(input_list = input_list,

                                    # Data inputs
                                    WAA = sgl_rg_sable_data$WAA,
                                    # weight-at-age
                                    MatAA = sgl_rg_sable_data$MatAA,
                                    # maturity at age
                                    AgeingError = as.matrix(sgl_rg_sable_data$age_error),
                                    # ageing error
                                    SizeAgeTrans = sgl_rg_sable_data$SizeAgeTrans,
                                    # size age transition matrix

                                    # Model options
                                    fit_lengths = 1,
                                    # fitting length compositions
                                    Use_M_prior = 1,
                                    # use natural mortality prior
                                    M_prior = c(0.1, 0.1),
                                    # mean and sd for M prior
                                    M_spec = "est_ln_M_only",
                                    ln_M = log(0.1143034), # starting value for M
                                    M_offset = -0.00819813327864
                                    # starting value /
                                    # fixing Male M offset (accidently jittered
                                    # from original assessment)
)

# setting up movement and tagging (not used)
input_list <- Setup_Mod_Movement(input_list = input_list,
                                 use_fixed_movement = 1, # don't est move
                                 Fixed_Movement = NA, # use identity move
                                 do_recruits_move = 0 # recruits don't move
)

input_list <- Setup_Mod_Tagging(input_list = input_list,
                                UseTagging = 0
)

# setting up catch specifications
input_list <- Setup_Mod_Catch_and_F(input_list = input_list,
                                    # Data inputs
                                    ObsCatch = sgl_rg_sable_data$ObsCatch,
                                    Catch_Type =
                                      array(1, dim = c(length(input_list$data$years),
                                                       input_list$data$n_fish_fleets)),
                                    UseCatch = sgl_rg_sable_data$UseCatch,
                                    # Model options
                                    Use_F_pen = 1,
                                    # whether to use f penalty, == 0 don't use, == 1 use
                                    sigmaC_spec = 'fix',
                                    Catch_Constant = c(0.01, 0.8)
)

# setting up fishery indices and compositions
input_list <- Setup_Mod_FishIdx_and_Comps(input_list = input_list,
                                          # data inputs
                                          ObsFishIdx = sgl_rg_sable_data$ObsFishIdx,
                                          ObsFishIdx_SE = sgl_rg_sable_data$ObsFishIdx_SE,
                                          UseFishIdx =  sgl_rg_sable_data$UseFishIdx,
                                          ObsFishAgeComps = sgl_rg_sable_data$ObsFishAgeComps,
                                          UseFishAgeComps = sgl_rg_sable_data$UseFishAgeComps,
                                          ISS_FishAgeComps = sgl_rg_sable_data$ISS_FishAgeComps,
                                          ObsFishLenComps = sgl_rg_sable_data$ObsFishLenComps,
                                          UseFishLenComps = sgl_rg_sable_data$UseFishLenComps,
                                          ISS_FishLenComps = sgl_rg_sable_data$ISS_FishLenComps,

                                          # Model options
                                          fish_idx_type =
                                            c("biom", "none"),
                                          # biomass indices for fishery fleet 1 and 2
                                          FishAgeComps_LikeType =
                                            c("Multinomial", "none"),
                                          # age comp likelihoods for fishery fleet 1 and 2
                                          FishLenComps_LikeType =
                                            c("Multinomial", "Multinomial"),
                                          # length comp likelihoods for fishery fleet 1 and 2
                                          FishAgeComps_Type =
                                            c("agg_Year_1-terminal_Fleet_1",
                                              "none_Year_1-terminal_Fleet_2"),
                                          # age comp structure for fishery fleet 1 and 2
                                          FishLenComps_Type =
                                            c("spltRspltS_Year_1-terminal_Fleet_1",
                                              "spltRspltS_Year_1-terminal_Fleet_2"),
                                          # length comp structure for fishery fleet 1 and 2
                                          FishAge_comp_agg_type = c(0,NA),
                                          # ADMB aggregation quirks, ideally get rid of this
                                          FishLen_comp_agg_type = c(1,1)
                                          # ADMB aggregation quirks, ideally get rid of this
)

# setting up survey indices and compositions
input_list <- Setup_Mod_SrvIdx_and_Comps(input_list = input_list,
                                         # data inputs
                                         ObsSrvIdx = sgl_rg_sable_data$ObsSrvIdx,
                                         ObsSrvIdx_SE = sgl_rg_sable_data$ObsSrvIdx_SE,
                                         UseSrvIdx =  sgl_rg_sable_data$UseSrvIdx,
                                         ObsSrvAgeComps = sgl_rg_sable_data$ObsSrvAgeComps,
                                         ISS_SrvAgeComps = sgl_rg_sable_data$ISS_SrvAgeComps,
                                         UseSrvAgeComps = sgl_rg_sable_data$UseSrvAgeComps,
                                         ObsSrvLenComps = sgl_rg_sable_data$ObsSrvLenComps,
                                         UseSrvLenComps = sgl_rg_sable_data$UseSrvLenComps,
                                         ISS_SrvLenComps = sgl_rg_sable_data$ISS_SrvLenComps,

                                         # Model options
                                         srv_idx_type = c("abd", "biom", "abd"),
                                         # abundance and biomass for survey fleet 1, 2, and 3
                                         SrvAgeComps_LikeType =
                                           c("Multinomial", "none", "Multinomial"),
                                         # survey age composition likelihood for survey fleet
                                         # 1, 2, and 3
                                         SrvLenComps_LikeType =
                                           c("Multinomial", "Multinomial", "Multinomial"),
                                         #  survey length composition likelihood for survey fleet
                                         # 1, 2, and 3
                                         SrvAgeComps_Type = c("agg_Year_1-terminal_Fleet_1",
                                                              "none_Year_1-terminal_Fleet_2",
                                                              "agg_Year_1-terminal_Fleet_3"),
                                         # survey age comp type
                                         SrvLenComps_Type = c("spltRspltS_Year_1-terminal_Fleet_1",
                                                              "spltRspltS_Year_1-terminal_Fleet_2",
                                                              "spltRspltS_Year_1-terminal_Fleet_3"),
                                         # survey length comp type
                                         SrvAge_comp_agg_type = c(1,NA,1),
                                         # ADMB aggregation quirks, ideally get rid of this
                                         SrvLen_comp_agg_type = c(0,0,0)
                                         # ADMB aggregation quirks, ideally get rid of this
)

# setting up fishery selectivity
input_list <- Setup_Mod_Fishsel_and_Q(input_list = input_list,

                                      # Model options
                                      cont_tv_fish_sel = c("none_Fleet_1", "none_Fleet_2"),
                                      # fishery selectivity, whether continuous time-varying

                                      # fishery selectivity blocks
                                      fish_sel_blocks =
                                        c("Block_1_Year_1-35_Fleet_1",
                                          # block 1, fishery ll selex
                                          "Block_2_Year_36-56_Fleet_1",
                                          # block 2 fishery ll selex
                                          "Block_3_Year_57-terminal_Fleet_1",
                                          # block 3 fishery ll selex
                                          "none_Fleet_2"),
                                      # no blocks for trawl fishery

                                      # fishery selectivity form
                                      fish_sel_model =
                                        c("logist1_Fleet_1",
                                          "gamma_Fleet_2"),

                                      # fishery catchability blocks
                                      fish_q_blocks =
                                        c("Block_1_Year_1-35_Fleet_1",
                                          # block 1, fishery ll selex
                                          "Block_2_Year_36-56_Fleet_1",
                                          # block 2 fishery ll selex
                                          "Block_3_Year_57-terminal_Fleet_1",
                                          # block 3 fishery ll selex
                                          "none_Fleet_2"),
                                      # no blocks for trawl fishery

                                      # whether to estimate all fixed effects
                                      # for fishery selectivity and later modify
                                      # to fix and share parameters
                                      fish_fixed_sel_pars =
                                        c("est_all", "est_all"),

                                      # whether to estimate all fixed effects
                                      # for fishery catchability
                                      fish_q_spec =
                                        c("est_all", "fix")
                                      # estimate fishery q for fleet 1, not for fleet 2
)

# additional mapping for fishery selectivity
# sharing delta across sexes from early domestic fishery (first time block)
input_list$map$ln_fish_fixed_sel_pars <- factor(c(1:7, 2, 8:11, rep(12:13,3), rep(c(14,13),3)))

# setting up survey selectivity
input_list <- Setup_Mod_Srvsel_and_Q(input_list = input_list,

                                     # Model options
                                     # survey selectivity, whether continuous time-varying
                                     cont_tv_srv_sel =
                                       c("none_Fleet_1",
                                         "none_Fleet_2",
                                         "none_Fleet_3"),

                                     # survey selectivity blocks
                                     srv_sel_blocks =
                                       c("Block_1_Year_1-56_Fleet_1",
                                         # block 1 for domestic ll survey
                                         "Block_2_Year_57-terminal_Fleet_1",
                                         # block 2 for domestic ll survey
                                         "none_Fleet_2",
                                         "none_Fleet_3"
                                       ), # no blocks for trawl and jp survey

                                     # survey selectivity form
                                     srv_sel_model =
                                       c("logist1_Fleet_1",
                                         "exponential_Fleet_2",
                                         "logist1_Fleet_3"),

                                     # survey catchability blocks
                                     srv_q_blocks =
                                       c("none_Fleet_1",
                                         "none_Fleet_2",
                                         "none_Fleet_3"),

                                     # whether to estiamte all fixed effects
                                     # for survey selectivity and later
                                     # modify to fix/share parameters
                                     srv_fixed_sel_pars_spec =
                                       c("est_all",
                                         "est_all",
                                         "est_all"),

                                     # whether to estiamte all
                                     # fixed effects for survey catchability
                                     srv_q_spec =
                                       c("est_all",
                                         "est_all",
                                         "est_all")
)

# ll survey, share delta female (index 2) across time blocks and to the coop jp ll survey delta
# ll survey, share delta male (index 5) across time blocks and to the coop jp ll survey delta
# coop jp survey does not estimate parameters and shares deltas with longline survey
# single time block with trawl survey and only one parameter hence,
# only one parameter estimated across blocks (indices 7 and 8)
input_list$map$ln_srv_fixed_sel_pars <-
  factor(c(1:3, 2, 4:6, 5,rep(7,4), rep(8, 4), rep(c(NA,2), 2), rep(c(NA, 5), 2)))

# Coop JP Survey (Logistic) Single time block (some of these estimates are fixed!)
input_list$par$ln_srv_fixed_sel_pars[1,,,1,3] <- c(0.980660760456, 0.9295241)
input_list$par$ln_srv_fixed_sel_pars[1,,,2,3] <- c(1.22224502478, 0.8821623)

# setting up model weighting
# set up data weighting stuff
Wt_FishAgeComps <- array(NA, dim = c(input_list$data$n_regions,
                                     length(input_list$data$years),
                                     input_list$data$n_sexes,
                                     input_list$data$n_fish_fleets))
# weights for fishery age comps
Wt_FishAgeComps[1,,1,1] <- 0.826107286513784
# Weight for fixed gear age comps

# Fishery length comps
Wt_FishLenComps <- array(NA, dim = c(input_list$data$n_regions,
                                     length(input_list$data$years),
                                     input_list$data$n_sexes,
                                     input_list$data$n_fish_fleets))
Wt_FishLenComps[1,,1,1] <- 4.1837057381917
# Weight for fixed gear len comps females
Wt_FishLenComps[1,,2,1] <- 4.26969350917589
# Weight for fixed gear len comps males
Wt_FishLenComps[1,,1,2] <- 0.316485920691651
# Weight for trawl gear len comps females
Wt_FishLenComps[1,,2,2] <- 0.229396580680981
# Weight for trawl gear len comps males

# survey age comps
Wt_SrvAgeComps <- array(NA, dim = c(input_list$data$n_regions,
                                    length(input_list$data$years),
                                    input_list$data$n_sexes,
                                    input_list$data$n_srv_fleets))
# weights for survey age comps
Wt_SrvAgeComps[1,,1,1] <- 3.79224544725927
# Weight for domestic survey ll gear age comps
Wt_SrvAgeComps[1,,1,3] <- 1.31681114024037
# Weight for coop jp survey ll gear age comps

# Survey length comps
Wt_SrvLenComps <- array(0, dim = c(input_list$data$n_regions,
                                   length(input_list$data$years),
                                   input_list$data$n_sexes,
                                   input_list$data$n_srv_fleets))
Wt_SrvLenComps[1,,1,1] <- 1.43792019016567
# Weight for domestic ll survey len comps females
Wt_SrvLenComps[1,,2,1] <- 1.07053763450712
# Weight for domestic ll survey len comps males
Wt_SrvLenComps[1,,1,2] <- 0.670883273592302
# Weight for domestic trawl survey len comps females
Wt_SrvLenComps[1,,2,2] <- 0.465207132450763
# Weight for domestic trawl survey len comps males
Wt_SrvLenComps[1,,1,3] <- 1.27772810174693
# Weight for coop jp ll survey len comps females
Wt_SrvLenComps[1,,2,3] <- 0.857519546948587
# Weight for coop jp ll survey len comps males

input_list <- Setup_Mod_Weighting(input_list = input_list,
                                  sablefish_ADMB = 1,
                                  likelihoods = 0,
                                  Wt_Catch = 50,
                                  Wt_FishIdx = 0.448,
                                  Wt_SrvIdx = 0.448,
                                  Wt_Rec = 1.5,
                                  Wt_F = 0.1,
                                  Wt_Tagging = 0,
                                  Wt_FishAgeComps = Wt_FishAgeComps,
                                  Wt_FishLenComps = Wt_FishLenComps,
                                  Wt_SrvAgeComps = Wt_SrvAgeComps,
                                  Wt_SrvLenComps = Wt_SrvLenComps
)

# extract out lists updated with helper functions
data <- input_list$data
parameters <- input_list$par
mapping <- input_list$map


# Fit Model ---------------------------------------------------------------

sabie_rtmb_model <- fit_model(data,
                              parameters,
                              mapping,
                              random = NULL,
                              newton_loops = 3,
                              silent = TRUE
)

# Get standard error report
sabie_rtmb_model$sd_rep <- RTMB::sdreport(sabie_rtmb_model)


# Reference Points --------------------------------------------------------

ref_pts <- SPoRC::Get_Reference_Points(data = data,
                                       rep = sabie_rtmb_model$rep,
                                       SPR_x = 0.4,
                                       t_spwn = 0,
                                       sex_ratio_f = 0.5,
                                       calc_rec_st_yr = 20,
                                       rec_age = 2, type = "single_region",
                                       what = "SPR")

b40 <- ref_pts$b_ref_pt
f40 <- ref_pts$f_ref_pt
sabie_rtmb_model$b40 <- ref_pts$b_ref_pt
sabie_rtmb_model$f40 <- ref_pts$f_ref_pt

# Catch Advice ------------------------------------------------------------
# Define HCR to use
HCR_function <- function(x, frp, brp, alpha = 0.05) {
  stock_status <- x / brp # define stock status
  # If stock status is > 1
  if(stock_status >= 1) f <- frp
  # If stock status is between brp and alpha
  if(stock_status > alpha && stock_status < 1) f <- frp * (stock_status - alpha) / (1 - alpha)
  # If stock status is less than alpha
  if(stock_status < alpha) f <- 0
  return(f)
}

# Setup necessary inputs
t_spawn <- 0
sexratio <- 0.5
n_proj_yrs <- 15
n_regions <- 1
n_ages <- length(data$ages)
n_sexes <- data$n_sexes
n_fish_fleets <- 2
do_recruits_move <- 0
terminal_NAA <- array(sabie_rtmb_model$rep$NAA[,length(data$years),,], dim = c(n_regions, n_ages, n_sexes))
WAA <- array(rep(data$WAA[,length(data$years),,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes)) # weight at age
WAA_fish <- array(rep(data$WAA[,length(data$years),,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets)) # weight at age
MatAA <- array(rep(data$MatAA[,length(data$years),,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes)) # maturity at age
fish_sel <- array(rep(sabie_rtmb_model$rep$fish_sel[,length(data$years),,,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes, n_fish_fleets)) # selectivity
Movement <- array(rep(sabie_rtmb_model$rep$Movement[,,length(data$years),,], each = n_proj_yrs), dim = c(n_regions, n_regions, n_proj_yrs, n_ages, n_sexes)) # movement - not used
terminal_F <- array(sabie_rtmb_model$rep$Fmort[,length(data$years),], dim = c(n_regions, n_fish_fleets)) # terminal F
natmort <- array(rep(sabie_rtmb_model$rep$natmort[,length(data$years),,], each = n_proj_yrs), dim = c(n_regions, n_proj_yrs, n_ages, n_sexes)) # natural mortaility
recruitment <- array(sabie_rtmb_model$rep$Rec[,20:(length(data$years) - 2)], dim = c(n_regions, length(20:(length(data$years)-2)))) # recruitment from years 20 - terminal (corresponds to 1979)

# Use FABC to project
f_ref_pt <- array(f40, dim = c(n_regions, n_proj_yrs))
b_ref_pt <- array(b40, dim = c(n_regions, n_proj_yrs))

# Do projection to get ABC
get_abc_vals <- Do_Population_Projection(n_proj_yrs = n_proj_yrs, # projection years
                                         n_regions = n_regions, # number of regions
                                         n_ages = n_ages, # number of ages
                                         n_sexes = n_sexes,  # number of sexes
                                         sexratio = sexratio, # sexratio
                                         n_fish_fleets = n_fish_fleets, # number of fishery fleets
                                         do_recruits_move = do_recruits_move, # whether recruits move
                                         recruitment = recruitment, # recruitment evector
                                         terminal_NAA = terminal_NAA, # terminal numbers at age
                                         terminal_F = terminal_F, # terminal F array
                                         natmort = natmort, # natural mortliaty array
                                         WAA = WAA, # weight at age array
                                         WAA_fish = WAA_fish,
                                         MatAA = MatAA, # maturity at age array
                                         fish_sel = fish_sel, # fishery selex array
                                         Movement = Movement, # movement array
                                         f_ref_pt = f_ref_pt, # fishery reference point
                                         b_ref_pt = b_ref_pt, # biological reference point
                                         HCR_function = HCR_function, # HCR function that takes x, frp, and brp
                                         recruitment_opt = "mean_rec", # recruitment assumption
                                         t_spawn = t_spawn # spawn fraction
)

# ABC for subsequent years (dimensions, region, year, fleet)
sum(get_abc_vals$proj_Catch[1,2,]) # year 2 = the actual projected year

# Save Model --------------------------------------------------------------

# Save data, parameters, and mapping in model object
sabie_rtmb_model$data <- data
sabie_rtmb_model$parameters <- parameters
sabie_rtmb_model$mapping <- mapping

# Write out RDS file
saveRDS(sabie_rtmb_model, file = here("2024 Assessment Files", "2024_RTMB_Model.RDS"))
saveRDS(input_list, file = here("2024 Assessment Files", "2024_RTMB_Model_List.RDS"))

