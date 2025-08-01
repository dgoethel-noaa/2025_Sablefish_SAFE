# Purpose: To bridge the 2024 federal sablefish ADMB models to SPoRC (RTMB)
# Creator: Matthew LH. Cheng (UAF CFOS)
# Date 5/28/25

# Workflow:
# 1) Compare the 3 ADMB model versions
# 2) Compare SPoRC model to the 3 ADMB model versions

# Setup -------------------------------------------------------------------
library(here)
library(ggplot2)
library(R2admb)
library(tidyverse)

# Define path to 2024 ADMB assessments
ass_path <- here("2024 Assessment Files")
figs_path <- here("Sept PT Model Runs",   "Bridging")

# Define admb model names
admb_names <- c("23.5", "23.5a", "23.5b")

# Read in RTMB model
rtmb_out <- readRDS(here(ass_path, "2024_RTMB_Model.RDS")) # read in base assessment createad from SPoRC_RTMB_Model.R
data <- rtmb_out$data

# Read in ADMB outputs
admb_out <- list(
  dget(here(ass_path, "2024 Base (23.5)_final model", "tem.rdat")), # base model (23.5)
  dget(here(ass_path, "2024 Base (23.5)_final model_a", "tem.rdat")), # 23.5a
  dget(here(ass_path, "2024 Base (23.5)_final model_b", "tem.rdat")) # 23.5b
)

# Read in ADMB parameters
admb_pars <- list(
  read_pars(here(ass_path, "2024 Base (23.5)_final model", "tem")), # base model (23.5)
  read_pars(here(ass_path, "2024 Base (23.5)_final model_a", "tem")), # 23.5a
  read_pars(here(ass_path, "2024 Base (23.5)_final model_b", "tem")) # 23.5b
)

# Compare ADMB Models -----------------------------------------------------
admb_ts <- data.frame() # empty dataframe for storage

# get admb outputs
for(i in 1:length(admb_out)) {

  out <- admb_out[[i]]
  par <- admb_pars[[i]]

  # Get time-series outputs
  tmp_ssb <- data.frame(Model = admb_names[i], Year = out$t.series$year, Value = out$t.series$spbiom, Type = "SSB")
  tmp_rec <- data.frame(Model = admb_names[i], Year = out$t.series$year, Value = out$t.series$Recr, Type = "Recruitment")
  tmp_f <- data.frame(Model = admb_names[i], Year = out$t.series$year, Value = out$t.series$fmort, Type = "Total F")
  tmp_nf <- data.frame(Model = admb_names[i], Year = out$t.series$year, Value = out$t.series$numbers.f, Type = "Total Females")
  tmp_nm <- data.frame(Model = admb_names[i], Year = out$t.series$year, Value = out$t.series$numbers.m, Type = "Total Males")
  tmp_ssb_se <- data.frame(Model = admb_names[i], Year = out$t.series$year, Value = par$se[str_detect(names(par$se), "ssb")], Type = "SSB (SE)")
  tmp_rec_se <- data.frame(Model = admb_names[i], Year = out$t.series$year, Value = par$se[str_detect(names(par$se), "pred_rec")], Type = "Recruitment (SE)")

  # bind
  admb_ts <- rbind(admb_ts, tmp_ssb, tmp_rec, tmp_f, tmp_nf, tmp_nm, tmp_ssb_se, tmp_rec_se)

} # end i loop

# Plot time series comparisons of ADMB models
png(here(figs_path, "ADMB_TimeSeries_Comparison.png"), width = 1000, height = 800)
print(
  ggplot(admb_ts, aes(x = Year, y = Value, color = Model)) +
    geom_line(lwd = 1.3) +
    coord_cartesian(ylim = c(0, NA)) +
    facet_wrap(~Type, scales = 'free') +
    theme_bw(base_size = 18) +
    theme(legend.position = "top")
)
dev.off()

# Compare SPoRC to ADMB Models --------------------------------------------

### Extract Time Series ---------------------------------------------------------
# Extract RTMB Time Series
rtmb_ssb <- data.frame(Model = "RTMB 24.1", Year = out$t.series$year, Value = as.vector(rtmb_out$rep$SSB), Type = "SSB")
rtmb_rec <- data.frame(Model = "RTMB 24.1", Year = out$t.series$year, Value = as.vector(rtmb_out$rep$Rec), Type = "Recruitment")
rtmb_f <- data.frame(Model = "RTMB 24.1", Year = out$t.series$year, Value = apply(rtmb_out$rep$Fmort,2,sum), Type = "Total F")
rtmb_nf <- data.frame(Model = "RTMB 24.1", Year = out$t.series$year, Value = rowSums(rtmb_out$rep$NAA[1,-66,,1]), Type = "Total Females")
rtmb_nm <- data.frame(Model = "RTMB 24.1", Year = out$t.series$year, Value = rowSums(rtmb_out$rep$NAA[1,-66,,2]), Type = "Total Males")
rtmb_ssb_se <- data.frame(Model = "RTMB 24.1", Year = out$t.series$year, Value = rtmb_out$sd_rep$sd[names(rtmb_out$sd_rep$value) == "SSB"], Type = "SSB (SE)")
rtmb_rec_se <- data.frame(Model = "RTMB 24.1", Year = out$t.series$year, Value = rtmb_out$sd_rep$sd[names(rtmb_out$sd_rep$value) == "Rec"], Type = "Recruitment (SE)")

# bind
rtmb_ts <- rbind(rtmb_ssb, rtmb_rec, rtmb_f, rtmb_nf, rtmb_nm, rtmb_ssb_se, rtmb_rec_se)

# left join to admb time series
all_ts <- rbind(admb_ts, rtmb_ts)

### Extract Selectivities ---------------------------------------------------
# Get selectivities
dom_ll_fish_f1 <- data.frame(Age = 1:30,
                             `RTMB 24.1` = rtmb_out$rep$fish_sel[1,1,,1,1],
                             ADMB = admb_out[[3]]$agesel$fish1sel.f,
                             Type = "Domestic LL Fishery Female Block 1")

dom_ll_fish_m1 <- data.frame(Age = 1:30,
                             `RTMB 24.1` = rtmb_out$rep$fish_sel[1,1,,2,1],
                             ADMB = admb_out[[3]]$agesel$fish1sel.m,
                             Type = "Domestic LL Fishery Male Block 1")

dom_ll_fish_f2 <- data.frame(Age = 1:30,
                             `RTMB 24.1` = rtmb_out$rep$fish_sel[1,40,,1,1],
                             ADMB = admb_out[[3]]$agesel$fish4sel.f,
                             Type = "Domestic LL Fishery Female Block 2")

dom_ll_fish_m2 <- data.frame(Age = 1:30,
                             `RTMB 24.1` = rtmb_out$rep$fish_sel[1,40,,2,1],
                             ADMB = admb_out[[3]]$agesel$fish4sel.m,
                             Type = "Domestic LL Fishery Male Block 2")

dom_ll_fish_f3 <- data.frame(Age = 1:30,
                             `RTMB 24.1` = rtmb_out$rep$fish_sel[1,60,,1,1],
                             ADMB = admb_out[[3]]$agesel$fish5sel.f,
                             Type = "Domestic LL Fishery Female Block 3")

dom_ll_fish_m3 <- data.frame(Age = 1:30,
                             `RTMB 24.1` = rtmb_out$rep$fish_sel[1,60,,2,1],
                             ADMB = admb_out[[3]]$agesel$fish5sel.m,
                             Type = "Domestic LL Fishery Male Block 3")

dom_trwl_fish_f <- data.frame(Age = 1:30,
                              `RTMB 24.1` = rtmb_out$rep$fish_sel[1,1,,1,2],
                              ADMB = admb_out[[3]]$agesel$fish3sel.f,
                              Type = "Domestic Trawl Female")

dom_trwl_fish_m <- data.frame(Age = 1:30,
                              `RTMB 24.1` = rtmb_out$rep$fish_sel[1,1,,2,2],
                              ADMB = admb_out[[3]]$agesel$fish3sel.m,
                              Type = "Domestic Trawl Male")

dom_ll_srv_f1 <- data.frame(Age = 1:30,
                            `RTMB 24.1` = rtmb_out$rep$srv_sel[1,1,,1,1],
                            ADMB = admb_out[[3]]$agesel$srv1sel.f,
                            Type = "Domestic LL Survey Female Block 1")

dom_ll_srv_m1 <- data.frame(Age = 1:30,
                            `RTMB 24.1` = rtmb_out$rep$srv_sel[1,1,,2,1],
                            ADMB = admb_out[[3]]$agesel$srv1sel.m,
                            Type = "Domestic LL Survey Male Block 1")

dom_ll_srv_f2 <- data.frame(Age = 1:30,
                            `RTMB 24.1` = rtmb_out$rep$srv_sel[1,60,,1,1],
                            ADMB = admb_out[[3]]$agesel$srv10sel.f,
                            Type = "Domestic LL Survey Female Block 2")

dom_ll_srv_m2 <- data.frame(Age = 1:30,
                            `RTMB 24.1` = rtmb_out$rep$srv_sel[1,60,,2,1],
                            ADMB = admb_out[[3]]$agesel$srv10sel.m,
                            Type = "Domestic LL Survey Male Block 2")

dom_trwl_srv_f2 <- data.frame(Age = 1:30,
                              `RTMB 24.1` = rtmb_out$rep$srv_sel[1,60,,1,2],
                              ADMB = admb_out[[3]]$agesel$srv7sel.f,
                              Type = "Domestic Trawl Survey Female")

dom_trwl_srv_m2 <- data.frame(Age = 1:30,
                              `RTMB 24.1` = rtmb_out$rep$srv_sel[1,60,,2,2],
                              ADMB = admb_out[[3]]$agesel$srv7sel.m,
                              Type = "Domestic Trawl Survey Male")

coop_ll_srv_f2 <- data.frame(Age = 1:30,
                             `RTMB 24.1` = rtmb_out$rep$srv_sel[1,60,,1,3],
                             ADMB = admb_out[[3]]$agesel$srv2sel.f,
                             Type = "Coop LL Survey Female")

coop_ll_srv_m2 <- data.frame(Age = 1:30,
                             `RTMB 24.1` = rtmb_out$rep$srv_sel[1,60,,2,3],
                             ADMB = admb_out[[3]]$agesel$srv2sel.m,
                             Type = "Coop LL Survey Male")

# bind!
all_sel_df <- rbind(
  dom_ll_fish_m1,
  dom_ll_fish_f2,
  dom_ll_fish_m2,
  dom_ll_fish_f3,
  dom_ll_fish_m3,
  dom_trwl_fish_f,
  dom_trwl_fish_m,
  dom_ll_srv_f1,
  dom_ll_srv_m1,
  dom_ll_srv_f2,
  dom_ll_srv_m2,
  dom_trwl_srv_f2,
  dom_trwl_srv_m2,
  coop_ll_srv_f2,
  coop_ll_srv_m2
)


### Extract Likelihoods -----------------------------------------------------
# Compare likelihoods
like_df <- data.frame(dat_type = c("jnLL", 'Fixed Gear Fishery Age',
                                      "Fixed Gear Fishery Length (F)",
                                      "Fixed Gear Fishery Length (M)",
                                      "Trawl Gear Fishery Length (F)",
                                      "Trawl Gear Fishery Length (M)",
                                      "Domestic Survey LL Age",
                                      "Domestic Survey LL Length (F)",
                                      "Domestic Survey LL Length (M)",
                                      "Domestic Trawl Survey Length (F)",
                                      "Domestic Trawl Survey Length (M)",
                                      "Japanese LL Survey Length (F)",
                                      "Japanese LL Survey Length (M)",
                                      "Catch",
                                      "Domestic LL Survey Index",
                                      "Japanese LL Survey Index",
                                      "Trawl Survey Index",
                                      "Domestic LL Fishery Index",
                                      "Japanese LL Fishery Index",
                                      "FMort Penalty",
                                      "M Prior",
                                      "Rec Penalty"
                                      ),
ADMB = c(
  admb_out[[3]]$likecomp[names(admb_out[[3]]$likecomp) == "obj.fun"],
  admb_out[[3]]$likecomp[names(admb_out[[3]]$likecomp) == "L.fish1age"],
  admb_out[[3]]$likecomp[names(admb_out[[3]]$likecomp) == "L.fish1sizef"],
  admb_out[[3]]$likecomp[names(admb_out[[3]]$likecomp) == "L.fish1sizem"],
  admb_out[[3]]$likecomp[names(admb_out[[3]]$likecomp) == "L.fish3sizef"],
  admb_out[[3]]$likecomp[names(admb_out[[3]]$likecomp) == "L.fish3sizem"],
  admb_out[[3]]$likecomp[names(admb_out[[3]]$likecomp) == "L.surv1age"],
  admb_out[[3]]$likecomp[names(admb_out[[3]]$likecomp) == "L.surv1sizef"],
  admb_out[[3]]$likecomp[names(admb_out[[3]]$likecomp) == "L.surv1sizem"],
  admb_out[[3]]$likecomp[names(admb_out[[3]]$likecomp) == "L.surv7sizef"],
  admb_out[[3]]$likecomp[names(admb_out[[3]]$likecomp) == "L.surv7sizem"],
  admb_out[[3]]$likecomp[names(admb_out[[3]]$likecomp) == "L.surv2sizef"],
  admb_out[[3]]$likecomp[names(admb_out[[3]]$likecomp) == "L.surv2sizem"],
  admb_out[[3]]$likecomp[names(admb_out[[3]]$likecomp) == "Catch"],
  admb_out[[3]]$likecomp[names(admb_out[[3]]$likecomp) == "L.surv3"],
  admb_out[[3]]$likecomp[names(admb_out[[3]]$likecomp) == "L.surv4"],
  admb_out[[3]]$likecomp[names(admb_out[[3]]$likecomp) == "L.surv7"],
  admb_out[[3]]$likecomp[names(admb_out[[3]]$likecomp) == "L.surv5"],
  admb_out[[3]]$likecomp[names(admb_out[[3]]$likecomp) == "L.surv6"],
  admb_out[[3]]$likecomp[names(admb_out[[3]]$likecomp) == "F.reg"],
  admb_out[[3]]$likecomp[names(admb_out[[3]]$likecomp) == "M.prior"],
  admb_out[[3]]$likecomp[names(admb_out[[3]]$likecomp) == "Rec.Pen"]
),
RTMB = c(
  rtmb_out$rep$jnLL,
  sum(rtmb_out$rep$FishAgeComps_nLL),
  sum(rtmb_out$rep$FishLenComps_nLL[,,1,1]),
  sum(rtmb_out$rep$FishLenComps_nLL[,,2,1]),
  sum(rtmb_out$rep$FishLenComps_nLL[,,1,2]),
  sum(rtmb_out$rep$FishLenComps_nLL[,,2,2]),
  sum(rtmb_out$rep$SrvAgeComps_nLL[,,1,1]),
  sum(rtmb_out$rep$SrvLenComps_nLL[,,1,1]),
  sum(rtmb_out$rep$SrvLenComps_nLL[,,2,1]),
  sum(rtmb_out$rep$SrvLenComps_nLL[,,1,2]),
  sum(rtmb_out$rep$SrvLenComps_nLL[,,2,2]),
  sum(rtmb_out$rep$SrvLenComps_nLL[,,1,3]),
  sum(rtmb_out$rep$SrvLenComps_nLL[,,2,3]),
  sum(rtmb_out$rep$Catch_nLL) * data$Wt_Catch,
  sum(rtmb_out$rep$SrvIdx_nLL[1,,1]) * data$Wt_SrvIdx,
  sum(rtmb_out$rep$SrvIdx_nLL[1,,3]) * data$Wt_SrvIdx,
  sum(rtmb_out$rep$SrvIdx_nLL[1,,2]) * data$Wt_SrvIdx,
  sum(rtmb_out$rep$FishIdx_nLL[,36:63,1]) * data$Wt_FishIdx,
  sum(rtmb_out$rep$FishIdx_nLL[1,-c(36:63),1]) * data$Wt_FishIdx,
  sum(rtmb_out$rep$Fmort_nLL) * data$Wt_F,
  rtmb_out$rep$M_nLL,
  sum(rtmb_out$rep$Init_Rec_nLL) * data$Wt_Rec + sum(rtmb_out$rep$Rec_nLL) * data$Wt_Rec
)) %>%
  mutate(ADMB = round(ADMB, 4),
         RTMB = round(RTMB, 4), Difference = ADMB - RTMB,
         Relative_Difference = (ADMB - RTMB) / RTMB * 100)

write.csv(like_df, here("Sept PT Model Runs", "Bridging", "Bridge_nLL.csv"))

### Extract Key Estimates -----------------------------------------
M_df <- data.frame(Par = "M",
                   TMB = exp(rtmb_out$sd_rep$par.fixed[names(rtmb_out$sd_rep$par.fixed) == "ln_M"]),
                   ADMB = exp(admb_pars[[3]]$coefficients[names(admb_pars[[3]]$coefficients) == 'logm']))

R0_df <- data.frame(Par = "Mean Recruitment",
                    TMB = rtmb_out$rep$R0,
                    ADMB = exp(admb_pars[[3]]$coefficients[names(admb_pars[[3]]$coefficients) == 'log_mean_rec']))

srv_q_df <- data.frame(Par = c('srv_domLL_q', 'srv_trwl_q', 'srv_coopLL_q'),
                       TMB = exp(rtmb_out$sd_rep$par.fixed[names(rtmb_out$sd_rep$par.fixed) == "ln_srv_q"]),
                       ADMB = c(exp(admb_pars[[3]]$coefficients[names(admb_pars[[3]]$coefficients) == 'log_q_srv1']),
                                exp(admb_pars[[3]]$coefficients[names(admb_pars[[3]]$coefficients) == 'log_q_srv7']),
                                exp(admb_pars[[3]]$coefficients[names(admb_pars[[3]]$coefficients) == 'log_q_srv2'])))

fish_q_df <- data.frame(Par = c('fish_domLL_q1', 'fish_domLL_q2', 'fish_domLL_q3'),
                        TMB = exp(rtmb_out$sd_rep$par.fixed[names(rtmb_out$sd_rep$par.fixed) == "ln_fish_q"]),
                        ADMB = c(exp(admb_pars[[3]]$coefficients[names(admb_pars[[3]]$coefficients) == 'log_q_srv6']),
                                 exp(admb_pars[[3]]$coefficients[names(admb_pars[[3]]$coefficients) == 'log_q_srv8']),
                                 exp(admb_pars[[3]]$coefficients[names(admb_pars[[3]]$coefficients) == 'log_q_LL_fish_recent'])))

ref_pts_df <- data.frame(Par = c("F40", "B40"),
                         TMB = c(rtmb_out$f40, rtmb_out$b40),
                         ADMB = c(admb_pars[[3]]$coefficients[names(admb_pars[[3]]$coefficients) == 'F40'],
                                  admb_pars[[3]]$coefficients[names(admb_pars[[3]]$coefficients) == 'B40']))

par_df <- rbind(M_df, R0_df, srv_q_df, fish_q_df, ref_pts_df)


### Plots -------------------------------------------------------------------
# Plot time series comparisons of ADMB and RTMB Models
png(here(figs_path, "RTMB_ADMB_TimeSeries_Comparison.png"), width = 1000, height = 800)
ggplot(all_ts %>% filter(Model %in% c("23.5b", "RTMB 24.1")),
       aes(x = Year, y = Value, color = Model, lty = Model)) +
  geom_line(size = 2) +
  facet_wrap(~Type, scales = "free_y") +
  coord_cartesian(ylim = c(0, NA)) +
  theme_bw(base_size = 18) +
  theme(legend.position = "top")
dev.off()

# Plot relative difference between ADMB and RTMB Models
png(here(figs_path, "RTMB_ADMB_TimeSeries_RelDiff.png"), width = 1000, height = 800)
all_ts %>% filter(Model %in% c("23.5b", "RTMB 24.1")) %>%
  pivot_wider(names_from = 'Model', values_from = 'Value') %>%
  ggplot(aes(x = Year, y = (`23.5b` - `RTMB 24.1`) / `23.5b`)) +
  geom_line(size = 1.3) +
  geom_hline(yintercept = 0, lty = 2, size = 1) +
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(ylim = c(-0.005, 0.005)) +
  facet_wrap(~Type, scales = "free_y") +
  labs(y = "Relative difference (%)") +
  theme_bw(base_size = 18)
dev.off()

# Plot Selectivity Comparison of ADMB and RTMB Models
png(here(figs_path, "RTMB_ADMB_Selex_Comparison.png"), width = 1200, height = 800)
all_sel_df %>%
  pivot_longer(cols = c("RTMB.24.1", "ADMB"), names_to = "Model") %>%
  mutate(Model = ifelse(Model == "ADMB", "23.5b", Model)) %>%
  ggplot(aes(x = Age, y = value, color = Model, lty = Model)) +
  geom_line(lwd = 1.3) +
  facet_wrap(~Type) +
  labs(y = "Selectivity") +
  theme_bw(base_size = 18) +
  theme(legend.position = 'top')
dev.off()

# Plot Selectivity Relative Difference of ADMB and RTMB Models
png(here(figs_path, "RTMB_ADMB_Selex_RelDiff.png"), width = 1200, height = 800)
ggplot(all_sel_df, aes(x = Age , y = (ADMB -`RTMB.24.1`) / ADMB)) +
  geom_line(size = 1.3) +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) +
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(ylim = c(-0.005, 0.005)) +
  facet_wrap(~Type) +
  labs(y = "Relative difference (%)") +
  theme_bw(base_size = 18)
dev.off()

# Plot Key Estimates
png(here(figs_path, "RTMB_ADMB_KeyQuantsRelDiff.png"), width = 1200, height = 800)
ggplot(par_df, aes(x = Par, y = (ADMB - TMB) / ADMB)) +
  geom_point(size = 7, position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) +
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(ylim = c(-0.005, 0.005)) +
  labs(y = "Relative difference (%)", x = "Key Estimates") +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# Get key quantites
reference_points_opt <- list(SPR_x = 0.4,
                             t_spwn = 0,
                             sex_ratio_f = 0.5,
                             calc_rec_st_yr = 21,
                             rec_age = 2,
                             type = "single_region",
                             what = "SPR"
)

proj_model_opt <- list(n_proj_yrs = 2,
                       n_avg_yrs = 1,
                       HCR_function = HCR_function <- function(x, frp, brp, alpha = 0.05) {
                         stock_status <- x / brp # define stock status
                         # If stock status is > 1
                         if(stock_status >= 1) f <- frp
                         # If stock status is between brp and alpha
                         if(stock_status > alpha && stock_status < 1) f <- frp * (stock_status - alpha) / (1 - alpha)
                         # If stock status is less than alpha
                         if(stock_status < alpha) f <- 0
                         return(f)
                       },
                       recruitment_opt = 'mean_rec',
                       fmort_opt = 'HCR')

out <- get_key_quants(list(rtmb_out$data), list(rtmb_out$rep), reference_points_opt, proj_model_opt, 1)
