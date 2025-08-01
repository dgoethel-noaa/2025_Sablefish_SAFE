
library(tidyverse)
library(here)
library(fs)


load(here("Sept PT Model Runs","_Survey Analysis","GOA_BSAI_Extrap","Sablefish_Data_2024.RData"))
SPoRC_dat <- readRDS(here("Sept PT Model Runs","_Survey Analysis","GOA_BSAI_Extrap","SPoRC_one_reg_age_joint_no_len_kp_JPN.rds"))


 model_run <- "5yr Extrap"

new_years <- c(2014,2016,2018,2020,2022)
styr <- 1960
filt_reg <- c("Western Gulf of Alaska","Central Gulf of Alaska","East Yakutat/Southeast","West Yakutat")
filt_reg_abb <- c('EG','CG','WG')
filt_reg_abb_alt <- 'AI'

lls_age_filter_sex_SPoRC_temp <- lls_age_effN %>%
  dplyr::rename(area= geographic_area_name) %>%
  dplyr::filter(!area %in% c("Amatuli Gully", "Shelikof Trough", "W-Grounds",
                             "Yakutat Valley", "Southeast Shelf")) %>%          # this is in Dana's code and is essentially same as using exploitable=1 filter (called 'Sigler Way' by Dana)
  dplyr::filter(sex %in% c(1, 2)) %>%
  dplyr::mutate(sex = ifelse(sex == 1, "Male", "Female"))

lls_age_filter_sex_SPoRC_yrs <- lls_age_filter_sex_SPoRC_temp %>%
  dplyr::filter(c(year %in% c(new_years) & !council_sablefish_management_area %in% filt_reg))          # filter out BS in

lls_age_filter_sex_SPoRC <- lls_age_filter_sex_SPoRC_temp %>% filter(!year %in% new_years) %>%
  full_join(lls_age_filter_sex_SPoRC_yrs) %>%
  dplyr::group_by(year,area,Age, sex) %>%
  dplyr::summarize(obs = dplyr::n()) %>%                                        # getting number of observations per age group and sex (numerator)
  dplyr::group_by(year,area) %>%                                            # NOTE: calculating using 'joint' approach  (summed within sex/age)
  dplyr::mutate(tot_obs = sum(obs)) %>%                                         # getting total number of observations per area across sexes/age
  dplyr::group_by(year,area,Age, sex) %>%
  dplyr::mutate(freq = obs/tot_obs)                                             # getting raw (unweighted) age freq per area

lls_age_prop_sex_SPoRC <- left_join(lls_age_filter_sex_SPoRC, lls_area_rpn_filter) %>%   # joining age frequencies with rpns
  dplyr::filter(year > 1995) %>%                                     # filtering out all years before 1995
  dplyr::mutate(wgt_rpn = RPN*freq) %>%                              # area-weighted age freq
  dplyr::group_by(year, Age, sex) %>%
  dplyr::summarize(sum_wgt_rpn = sum(wgt_rpn, na.rm=T)) %>%          # summing the area-weighted age freq across areas within sex (get numerator)
  dplyr::group_by(year) %>%
  dplyr::mutate(Tot_wgt_rpn = sum(sum_wgt_rpn)) %>%                  # summing area-weighted age freq across ages and sexes (get denominator)
  dplyr::group_by(year, Age, sex) %>%
  dplyr::mutate(Prop = sum_wgt_rpn/ Tot_wgt_rpn) %>%                 # dividing area weighted age freq by sex (summed across areas) by summed age freq (summed across sex, ages and areas) to get final age proportion
  dplyr::select(year, sex, Age, Prop) %>%
  vroom::vroom_write(here::here("Sept PT Model Runs","_Survey Analysis","GOA_BSAI_Extrap",model_run, "lls_rpn_age_comp_prop_SPoRC.csv"),delim = ",")

for(i in new_years){
  lls_rpn_dep_area_temp_yrs_repl <- lls_rpn_dep_area_temp %>%
    dplyr::filter(year %in% (i-1) & area %in% filt_reg_abb) %>%
    #full_join(lls_rpn_dep_area_temp %>%
    #            dplyr::filter(year %in% (i-2) & area == filt_reg_abb)) %>%
    dplyr::mutate(year = i) %>%
    full_join(lls_rpn_dep_area_temp %>%
                dplyr::filter(year %in% i & !area %in% filt_reg_abb)) %>%
    full_join(lls_rpn_dep_area_temp %>%
                dplyr::filter(year %in% (i-1) & !area %in% filt_reg_abb)) %>%
    dplyr::summarize(RPN = sum(rpn,na.rm=T))

  SPoRC_dat$ObsSrvIdx[,i-styr+1,1] <- lls_rpn_dep_area_temp_yrs_repl %>% pull(RPN)


  SPoRC_dat$ObsSrvAgeComps[,i-styr+1,,1,1] <- lls_age_prop_sex_SPoRC %>% ungroup() %>% filter(sex == "Female" & year == i) %>%
    full_join(data.frame(Age=rep(2:31),year=rep(i,30),sex=rep('Female',30))) %>%
    mutate(Prop = replace_na(Prop,0)) %>% arrange(Age) %>%
    pull(Prop)


  SPoRC_dat$ObsSrvAgeComps[,i-styr+1,,2,1] <- lls_age_prop_sex_SPoRC %>% ungroup() %>% filter(sex == "Male" & year == i) %>%
    full_join(data.frame(Age=rep(2:31),year=rep(i,30),sex=rep('Male',30))) %>%
    mutate(Prop = replace_na(Prop,0)) %>% arrange(Age) %>%
    pull(Prop)

}

file_copy(here("Sept PT Model Runs","_Survey Analysis","GOA_BSAI_Extrap","srvy_data_cln.r"), here("Sept PT Model Runs","_Survey Analysis","GOA_BSAI_Extrap",model_run,"srvy_data_cln.r"), overwrite = TRUE)

saveRDS(SPoRC_dat, file = here("Sept PT Model Runs","_Survey Analysis","GOA_BSAI_Extrap",model_run,"SPoRC_one_reg_age_joint_no_len_kp_JPN.rds"))

