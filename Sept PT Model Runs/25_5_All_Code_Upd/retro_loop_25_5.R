########################################################################################################################

# Purpose: To create an RTMB model for the sablefish assessment using SPoRC
# Creator: Matthew LH. Cheng (UAF CFOS)
# Date 5/28/25

########################################################################################################################
# Setup -------------------------------------------------------------------

library(here)
library(SPoRC)
library(ggplot2)
library(RTMB)
library(tidyr)
library(tidyverse)
library(cowplot)

#############################################################################################################################################################################################################
#---------USER INPUTS--NEED TO UPDATE

mod_name_main <- '25_5_All_Code_Upd'                                                    # Model name to use for saving files and labels
root_folder <- here("Sept PT Model Runs","_True_Retros", mod_name_main)                                           # The root folder where runs will be stored (sub folders named with mod_name)
mod <- readRDS(here(root_folder,paste0(mod_name_main,"_input_list.RDS")))

comp_name <- "Code_Upd"
do_francis <- 1                                                                # whether or not to do francis reweighting for this run, ==0 NO, ==1 YES
peels <- 7
do_mod_comp <- 1                                                               # whether to run comparison graphics with another model retro, if yes (==1) then need the associated rds files in the root folder (data, rep, sdrep, and mod_names)
mod_comp_name <- "Cont_Retro"                                                  # name of model runs to read in

do_high_res_plots <- 1
Rec_Age <- 2
##############################################################################################################################################################################################################
##############################################################################################################################################################################################################
# Make plotting directory
dir.create(here(root_folder,"plots"))
plot_fold <- here(root_folder,"plots")


rtmb_out <- list()
data_list <- list()
rep_list <- list()
sd_rep <- list()
mod_name <- list()
retro_data <- data.frame()



for(i in 0:peels){

  new_input_list <- SPoRC:::truncate_yr(i, mod$data, mod$par, mod$map)           # use truncate function to remove 1 year of data
  data <- new_input_list$retro_data
  parameters <- new_input_list$retro_parameters
  mapping <- new_input_list$retro_mapping

  data$UseSrvAgeComps[,65-i,1] <- 0                                               # no terminal year age comps
  data$UseFishAgeComps[,65-i,1] <- 0                                              # no terminal year age comps
  data$UseFishLenComps[,65-i,1] <- 0                                              # no terminal year length comps fishery
  data$UseFishLenComps[,65-i,2] <- 0                                              # no terminal year length comps fishery

  if(do_francis == 1){

    n_francis_iter <- 10                                                           # number of francis iterations to run

    # Define dimensions
    years <- data$years
    ages <- data$ages
    lens <- data$lens

    # Loop through francis
    for(j in 1:n_francis_iter) {

      if(j == 1) { # reset weights at 1
        data$Wt_FishAgeComps[] <- 1
        data$Wt_FishLenComps[] <- 1
        data$Wt_SrvAgeComps[] <- 1
        data$Wt_SrvLenComps[] <- 1
      } else {
        data$Wt_FishAgeComps[] <- wts$new_fish_age_wts
        data$Wt_FishLenComps[] <- wts$new_fish_len_wts
        data$Wt_SrvAgeComps[] <- wts$new_srv_age_wts
        data$Wt_SrvLenComps[] <- wts$new_srv_len_wts
      }

      # run model
      sabie_rtmb_model <- SPoRC::fit_model(data,
                                           parameters,
                                           mapping,
                                           random = NULL,
                                           newton_loops = 3,
                                           silent = TRUE
      )

      rep <- sabie_rtmb_model$report(sabie_rtmb_model$env$last.par.best) # Get report

      # get francis weights
      wts <- SPoRC::do_francis_reweighting(data = data, rep = rep, age_labels = ages,
                                           len_labels = lens, year_labels = years)

      # get francis weights for first iteration
      if(j == 1) {
        wts1 <- SPoRC::do_francis_reweighting(data = data, rep = rep, age_labels = ages, len_labels = lens, year_labels = years)
      }
    }
  }
  sabie_rtmb_model <- fit_model(data,
                                parameters,
                                mapping,
                                random = NULL,
                                newton_loops = 3,
                                silent = TRUE
  )

  # Get standard error report
  sabie_rtmb_model$sd_rep <- RTMB::sdreport(sabie_rtmb_model)

  # Save Model --------------------------------------------------------------

  # Save data, parameters, and mapping in model object
  sabie_rtmb_model$data <- data
  sabie_rtmb_model$parameters <- parameters
  sabie_rtmb_model$mapping <- mapping
  rep <- sabie_rtmb_model$report(sabie_rtmb_model$env$last.par.best) # Get report

  retro_tmp <- reshape2::melt(rep$SSB) %>% dplyr::rename(Region = Var1, Year = Var2) %>% dplyr::mutate(Type = "SSB") %>%
    bind_rows(reshape2::melt(rep$Rec) %>% dplyr::rename(Region = Var1, Year = Var2) %>% dplyr::mutate(Type = "Recruitment")) %>% dplyr::mutate(peel = i)

  retro_data <- rbind(retro_data, retro_tmp)

  # Check convergence, will let you know if high correlation, large gradient, or high SE for any pars
  # messages <- capture.output({
  #  post_optim_sanity_checks(sabie_rtmb_model$sd_rep,sabie_rtmb_model$rep,gradient_tol = 0.00001,se_tol = 5,corr_tol = 0.95)
  #}, type = "message")

  data_list[[i+1]] <- data
  rep_list[[i+1]] <- rep
  sd_rep[[i+1]] <- sabie_rtmb_model$sd_rep
  mod_name[[i+1]] <- paste0(mod_name_main,"_",2024-i)
}



saveRDS(rep_list, file = here(root_folder,paste0(comp_name,"_rep_list.RDS")))
saveRDS(data_list, file = here(root_folder,paste0(comp_name,"_data_list.RDS")))
saveRDS(sd_rep, file = here(root_folder,paste0(comp_name,"_sd_rep.RDS")))
saveRDS(mod_name, file = here(root_folder,paste0(comp_name,"_mod_name.RDS")))


info <- readRDS(here(root_folder,"SPoRC_one_reg.rds"))                            # read in data file to get years
mod_start_year <- info$yrs[1]
mod_years <- info$yrs

mod_name <- unlist(mod_name)

# get index fits ------------------------------------------------------------------
idx_fit <- get_idx_fits_plot(data_list , rep_list, mod_name)

idx_plot <- idx_fit+
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                              labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))+
  facet_wrap(~Category, scales = "free_y", labeller =
               labeller(Category = c(`Fishery1, Q1` = "JPN Fishery CPUE",
                                     `Fishery1, Q2` = "Fishery CPUE",
                                     `Fishery1, Q3` = "Fishery CPUE",
                                     `Survey1, Q1` = "LL Survey RPNs",
                                     `Survey2, Q1` = "GOA Trawl Survey (kt)",
                                     `Survey3, Q1` = "JPN LL Survey RPNs")))+
  labs(title = "Index Fits")

# get biologicals ------------------------------------------------------------------
bios <- get_biological_plot(data_list, rep_list, mod_name)

M_plot <- bios[[2]]+
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  labs(title = "Natural Mortality")

WAA_plot <- bios[[3]]+
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  labs(title = "Weight-at-Age")

MAA_plot <- bios[[4]]+
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  labs(title = "Maturity")

# get selectivity------------------------------------------------------------------
selex <- get_selex_plot(rep_list, mod_name,Selex_Type = "age",year_indx = c(1:length(mod_years)))

f_selex <- selex[[1]]+
  labs(title = "Fishery Selectivity")

s_selex <- selex[[2]]+
  labs(title = "Survey Selectivity")


# fits to catch --------------------------------------------------------------

catch <- get_catch_fits_plot(data_list, rep_list, mod_name)

catch_plot <- catch + labs(title = "Catch")


# get time series plot ------------------------------------------------------------------
time_series <- get_ts_plot(rep_list, sd_rep, mod_name,do_ci = TRUE)

ts_plot_CI <- time_series[[1]] +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                              labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))+ labs(title = "Time Series Plots")

F_plot_CI <- time_series[[2]]+  facet_wrap(~Type) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                              labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))+ labs(title = "Fishing Mortality")

recr_CI <-      time_series[[3]]+  facet_wrap(~Type) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                              labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))+ labs(title = "Recruitment")

ssb_CI <-      time_series[[4]]+  facet_wrap(~Type) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                              labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))+ labs(title = "SSB")

bio_CI <-      time_series[[5]]+  facet_wrap(~Type) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                              labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))+ labs(title = "Biomass")


# REPEAT PLOTS WITHOUT CIs
time_series <- get_ts_plot(rep_list, sd_rep, mod_name,do_ci = FALSE)

ts_plot <-      time_series[[1]] +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                              labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5)) + labs(title = "Time Series Plots")

F_plot <- time_series[[2]]+  facet_wrap(~Type) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                              labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))+ labs(title = "Fishing Mortality")

recr <-      time_series[[3]]+  facet_wrap(~Type) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                              labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))+ labs(title = "Recruitment")

ssb <-     time_series[[4]]+  facet_wrap(~Type) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                              labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))+ labs(title = "SSB")

bio <-      time_series[[5]]+  facet_wrap(~Type) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind()+
  ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                              labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))+ labs(title = "Biomass")


# get data fitted plot ------------------------------------------------------------------
data_used <- get_data_fitted_plot(data_list, mod_name)

data_plot <- data_used+
  ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                              labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))+ labs(title = "Data")

# get nLL plot ------------------------------------------------------------------
nll <- get_nLL_plot(data_list, rep_list, mod_name)

nll_plot <- nll[[1]]+ labs(title = "Negative Log-Likelihood")


# Add table of terminal year values (SSB, stock status, ABC) and table of likelihood components, convergence status, max grad, etc.

# get key quantities
reference_points_opt <- list(SPR_x = 0.4,
                             t_spwn = 0,
                             sex_ratio_f = 0.5,
                             calc_rec_st_yr = 20,
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
ABC <- get_key_quants(data_list, rep_list, reference_points_opt, proj_model_opt, mod_name)

ABC[[1]]  # key quantities data.frame

tiff(filename=here(plot_fold,paste0(comp_name,"_status+BRPs.png")),units = "in", width=16,height=4, res = 300)
ABC[[2]]  # table plot as cowplot object
dev.off()

saveRDS(ABC, file = here(root_folder,paste0(comp_name,"_ABC.RDS")))

# Recruitment squid plot
squid <- retro_data %>%
  dplyr::mutate(Year = Year, terminal = max(retro_data$Year) - peel, cohort = Year - Rec_Age, years_est = terminal-Year) %>%
  dplyr::filter(Type == 'Recruitment', cohort %in% seq(55, max(retro_data$Year), 1), terminal != Year) %>%                                        # show all cohorts back to 2014
  ggplot2::ggplot(ggplot2::aes(x = years_est - 1, y = value, group = Year, color = factor(cohort))) +
  ggplot2::geom_line(lwd = 1.3) +
  ggplot2::geom_point(size = 4) +
  ggplot2::theme_bw(base_size = 15) +
  ggplot2::labs(x = 'Years since cohort was first estimated', y = 'Recruitment (millions)', color = 'Cohort')

tiff(filename=here(plot_fold,paste0(comp_name,"_Recr_squid.png")),units = "in", width=22,height=16, res = 300)
squid
dev.off()
saveRDS(squid, file = here(root_folder,paste0(comp_name,"_squid.RDS")))

tiff(filename=here(plot_fold,paste0(comp_name,"_ts_F_noCI.png")),units = "in", width=24,height=6, res = 300)
print(F_plot)
dev.off()
tiff(filename=here(plot_fold,paste0(comp_name,"_ts_Recr_noCI.png")),units = "in", width=22,height=10, res = 300)
print(recr)
dev.off()
tiff(filename=here(plot_fold,paste0(comp_name,"_ts_SSB_noCI.png")),units = "in", width=22,height=10, res = 300)
print(ssb)
dev.off()
tiff(filename=here(plot_fold,paste0(comp_name,"_ts_Bio_noCI.png")),units = "in", width=22,height=10, res = 300)
print(bio)
dev.off()


# plot all basic info (basically just a wrapper function for a subset of functions defined above) ------------------------------------------------------------------
all_plots <- plot_all_basic(data_list, rep_list, sd_rep, mod_name, out_path = here(plot_fold))


pdf(file = here(plot_fold,paste0(comp_name,"_Full Comparison.pdf")),  width = 25, height = 13)

print(squid)
print(F_plot)
print(recr)
print(ssb)
print(bio)

print(ABC[[2]])
print(nll[[2]])
print(nll_plot)

print(data_plot)

print(M_plot)
print(WAA_plot)
print(MAA_plot)

print(idx_plot)
print(catch_plot)

print(f_selex)
print(s_selex)
print(ts_plot)
print(ts_plot_CI)
print(F_plot_CI)
print(recr_CI)
print(ssb_CI)
print(bio_CI)

dev.off()


if(do_high_res_plots == 1){
  tiff(filename=here(plot_fold,paste0(comp_name,"_nll.png")),units = "in", width=12,height=10, res = 300)
  print(nll_plot)
  dev.off()
  tiff(filename=here(plot_fold,paste0(comp_name,"_nll_tab.png")),units = "in", width=24,height=2, res = 300)
  print(nll[[2]])
  dev.off()
  tiff(filename=here(plot_fold,paste0(comp_name,"_data.png")),units = "in", width=14,height=12, res = 300)
  print(data_plot)
  dev.off()


  tiff(filename=here(plot_fold,paste0(comp_name,"_M.png")),units = "in", width=12,height=8, res = 300)
  print(M_plot)
  dev.off()
  tiff(filename=here(plot_fold,paste0(comp_name,"_WAA.png")),units = "in", width=12,height=8, res = 300)
  print(WAA_plot)
  tiff(filename=here(plot_fold,paste0(comp_name,"_Mat.png")),units = "in", width=12,height=8, res = 300)
  print(MAA_plot)
  dev.off()

  tiff(filename=here(plot_fold,paste0(comp_name,"_idx_fit.png")),units = "in", width=24,height=14, res = 300)
  print(idx_plot)
  dev.off()
  tiff(filename=here(plot_fold,paste0(comp_name,"_selex_fish.png")),units = "in", width=12,height=10, res = 300)
  print(f_selex)
  dev.off()
  tiff(filename=here(plot_fold,paste0(comp_name,"_selex_srvy.png")),units = "in", width=12,height=10, res = 300)
  print(s_selex)
  dev.off()
  tiff(filename=here(plot_fold,paste0(comp_name,"_llf_catch_fit.png")),units = "in", width=14,height=12, res = 300)
  print(catch_plot)
  dev.off()


  tiff(filename=here(plot_fold,paste0(comp_name,"_ts_noCI.png")),units = "in", width=14,height=12, res = 300)
  print(ts_plot)
  dev.off()


  tiff(filename=here(plot_fold,paste0(comp_name,"_ts.png")),units = "in", width=14,height=12, res = 300)
  print(ts_plot_CI)
  dev.off()
  tiff(filename=here(plot_fold,paste0(comp_name,"_ts_F.png")),units = "in", width=24,height=6, res = 300)
  print(F_plot_CI)
  dev.off()
  tiff(filename=here(plot_fold,paste0(comp_name,"_ts_Recr.png")),units = "in", width=22,height=10, res = 300)
  print(recr_CI)
  dev.off()
  tiff(filename=here(plot_fold,paste0(comp_name,"_ts_SSB.png")),units = "in", width=22,height=10, res = 300)
  print(ssb_CI)
  dev.off()
  tiff(filename=here(plot_fold,paste0(comp_name,"_ts_Bio.png")),units = "in", width=22,height=10, res = 300)
  print(bio_CI)
  dev.off()
}



if(do_mod_comp == 1){                                                        # code for making comparison graphics, if doing comparisons

  rep_list2 <- readRDS(file = here(root_folder,paste0(mod_comp_name,"_rep_list.RDS")))
  data_list2 <-  readRDS(file = here(root_folder,paste0(mod_comp_name,"_data_list.RDS")))
  sd_rep2 <- readRDS(file = here(root_folder,paste0(mod_comp_name,"_sd_rep.RDS")))
  mod_name2 <- unlist(readRDS(file = here(root_folder,paste0(mod_comp_name,"_mod_name.RDS"))))
  ABC2 <-  readRDS(file = here(root_folder,paste0(mod_comp_name,"_ABC.RDS")))
  squid2 <-  readRDS(file = here(root_folder,paste0(mod_comp_name,"_squid.RDS")))

  dir.create(here(root_folder,paste0("plots_",comp_name,"_v_",mod_comp_name)))
  plot_comp_dir <-  here(root_folder,paste0("plots_",comp_name,"_v_",mod_comp_name))


  tiff(filename=here(plot_comp_dir,paste0("status+BRPs.png")),units = "in", width=16,height=10, res = 300)
  print(plot_grid(ABC[[2]], ABC2[[2]], ncol=1, nrow =2))
  dev.off()


  time_series <- get_ts_plot(rep_list, sd_rep, mod_name,do_ci = FALSE)
  ssb = time_series[[4]]+  facet_wrap(~Type) +
    ggthemes::scale_color_colorblind() +
    ggthemes::scale_fill_colorblind()+
    ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                                labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))

  time_series2 <- get_ts_plot(rep_list2, sd_rep2, mod_name2,do_ci = FALSE)
  ssb2 = time_series2[[4]]+  facet_wrap(~Type) +
    ggthemes::scale_color_colorblind() +
    ggthemes::scale_fill_colorblind()+
    ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                                labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))

  tiff(filename=here(plot_comp_dir,paste0("/SSB.png")),units = "in", width=22,height=18, res = 300)
  print(plot_grid(ssb, ssb2, ncol=1, nrow =2))
  dev.off()


  recr <- time_series[[3]]+  facet_wrap(~Type) +
    ggthemes::scale_color_colorblind() +
    ggthemes::scale_fill_colorblind()+
    ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                                labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))

  recr2 <- time_series2[[3]]+  facet_wrap(~Type) +
    ggthemes::scale_color_colorblind() +
    ggthemes::scale_fill_colorblind()+
    ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                                labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))

  tiff(filename=here(plot_comp_dir,paste0("/Recr.png")),units = "in", width=22,height=18, res = 300)
  print(plot_grid(recr, recr2, ncol=1, nrow =2))
  dev.off()


  tiff(filename=here(plot_comp_dir,paste0("/Recr_squid.png")),units = "in", width=22,height=18, res = 300)
  print(plot_grid(squid, squid2, ncol=1, nrow =2))
  dev.off()

  # get index fits ------------------------------------------------------------------

  idx_fit <- get_idx_fits_plot(data_list , rep_list, mod_name)
  idx_plot <- idx_fit+
    ggthemes::scale_color_colorblind() +
    ggthemes::scale_fill_colorblind()+
    ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                                labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))+
    facet_wrap(~Category, scales = "free_y", labeller =
                 labeller(Category = c(`Fishery1, Q1` = "JPN Fishery CPUE",
                                       `Fishery1, Q2` = "Fishery CPUE",
                                       `Fishery1, Q3` = "Fishery CPUE",
                                       `Survey1, Q1` = "LL Survey RPNs",
                                       `Survey2, Q1` = "GOA Trawl Survey (kt)",
                                       `Survey3, Q1` = "JPN LL Survey RPNs")))

  idx_fit2 <- get_idx_fits_plot(data_list2 , rep_list2, mod_name2)
  idx_plot2 <- idx_fit2+
    ggthemes::scale_color_colorblind() +
    ggthemes::scale_fill_colorblind()+
    ggplot2::scale_x_continuous(breaks = seq(1, length(mod_years), 5),
                                labels = (mod_start_year  - 1)+ seq(1, length(mod_years), 5))+
    facet_wrap(~Category, scales = "free_y", labeller =
                 labeller(Category = c(`Fishery1, Q1` = "JPN Fishery CPUE",
                                       `Fishery1, Q2` = "Fishery CPUE",
                                       `Fishery1, Q3` = "Fishery CPUE",
                                       `Survey1, Q1` = "LL Survey RPNs",
                                       `Survey2, Q1` = "GOA Trawl Survey (kt)",
                                       `Survey3, Q1` = "JPN LL Survey RPNs")))


  tiff(filename=here(plot_comp_dir,paste0("idx_fit.png")),units = "in", width=20,height=20, res = 300)
  print(plot_grid(idx_plot, idx_plot2, ncol=1, nrow =2))

  dev.off()
}

