# NHIS-NDI AAMR
# Rishi Shah

# load necessary libraries
library(data.table)
library(foreign)
library(dplyr)
library(survey)
library(srvyr)
library(data.table)
library(epiDisplay)
library(ggplot2)
library(gtsummary)
library(gt)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# read in NHIS data
df <- fread('nhis_ndi_86-18.csv')
names(df) <- tolower(names(df))

# keeping only those eligible for mortality follow-up
# df <- df[df$astatflg == 1, ] # superfluous
df <- df[df$mortelig == 1, ] # can do this here or within the function
df$mortwtsa <- ifelse(df$year < 1997, df$mortwt, df$mortwtsa) # crosswalk mortality weights from year < 1997

######################################
# 1. Independent Variables           #
######################################

# 1.1 Race/Ethnicity
df$race <- NA
df$race[df$racea == 100 & df$hispeth == 10] <- 1
df$race[df$racea == 200 & df$hispeth == 10] <- 3
df$race[(df$racea > 201 & df$racea < 889) & df$hispeth == 10] <- 7
df$race[(df$racea >= 400 & df$racea < 500) & df$hispeth == 10] <- 5
df$race[df$hispeth > 10 & df$hispeth <= 70] <- 9

levels <- c("NH White", "NH Black", "NH Asian", "NH Other", "Hispanic")
labels <- c(1, 3, 5, 7, 9)
df$race <- factor(df$race, levels = labels, labels = levels)

# 1.2 Sex
# create 'female' variable and set initial value to 0
df$female <- 0

# replace 'female' with 1 if 'sex' is equal to 2
df$female[df$sex == 2] <- 1
levels <- c("Male", "Female")
labels <- c(0, 1)
df$female <- factor(df$female, levels = labels, labels = levels)

######################################
# 2. Dependent Variables             #
######################################

# leading underlying cause of death
# df$mortucodld
# 01	Diseases of heart	
# 02	Malignant neoplasms	
# 03	Chronic lower respiratory diseases	
# 04	Accidents (unintentional injuries)
# 05	Cerebrovascular diseases	
# 06	Alzheimer's disease	
# 07	Diabetes mellitus	
# 08	Influenza and pneumonia	
# 09	Nephritis, nephrotic syndrome, and nephrosis
# 10	All other causes (residual)


######################################
# 3. Analysis                        #
######################################
# calculate survey-adjusted and raw age-adjusted mortality rates

# define population weights for standardization (from the 2000 US census: https://www.cdc.gov/nchs/data/statnt/statnt20.pdf)
age_shares <- c(0.215746, 0.204517, 0.207431, 0.149771, 0.098425, 0.079180, 0.044930)

# function to calculate survey-adjusted and age-adjusted mortality rate for a given year
calculate_age_adjusted_rate_survey <- function(target_year, df, age_shares) {
  # subset data for the mortality year including relevant previous years and exclude individuals who died before the target year
  df_target <- df %>% filter(year <= target_year & mortelig == 1 & !is.na(mortdody) & !(mortdody < target_year))
  
  # adjust age for each respondent to their age at the target year
  df_target <- df_target %>% 
    mutate(age_target = age + (target_year - year)) %>% 
    mutate(age_target_cat = case_when(
      age_target >= 18 & age_target <= 29 ~ 1,
      age_target >= 30 & age_target <= 39 ~ 2,
      age_target >= 40 & age_target <= 49 ~ 3,
      age_target >= 50 & age_target <= 59 ~ 4,
      age_target >= 60 & age_target <= 69 ~ 5,
      age_target >= 70 & age_target <= 79 ~ 6,
      age_target >= 80 & age_target < 100 ~ 7, # NCHS recommends excluding anyone above the age of 100 has probability of survival to that age is low
      age_target >= 100 ~ NA_integer_,
      TRUE ~ NA_integer_
    ))
  # optional: filter for just Black and White races
  df_target <- df_target %>% filter(race %in% c("NH Black", "NH White"))
  
  # define all_cause_death within the function scope
  df_target <- df_target %>% mutate(all_cause_death = ifelse(mortucodld %in% 1:10 & mortdody == target_year & mortstat == 1, 1, 0))
  
  # define survey design for the target year
  # adjust sample weight for pooled analysis
  df_target$pooled_weight <- df_target$mortwtsa / (target_year - 1986 + 1) 
  
  options(survey.lonely.psu="adjust")
  
  design_target <- svydesign(
    ids = ~psu,
    strata = ~strata,
    weights = ~pooled_weight,
    data = df_target,
    nest = TRUE
  )
  
  # calculate crude death rates by age group, race, and gender
  crude_death_rates_by_race_gender <- svyby(
    ~all_cause_death,
    ~age_target_cat + race + female,
    design_target,
    svymean,
    na.rm = TRUE,
    vartype = "ci",
    method = "betaWald" # can also try default Taylor series linearization for CIs for survey means or the Korn-Graubard method; the overall difference between these techniques is negligible
  )
  
  # initialize a data frame to store the results
  results <- data.frame(Year = integer(), Race = character(), Gender = character(), AgeAdjustedRate = numeric(), LowerCI = numeric(), UpperCI = numeric())
  
  # loop through each combination of race and gender
  for (race_category in c("NH Black", "NH White")) {
    for (gender in unique(crude_death_rates_by_race_gender$female)) {
      # filter the crude death rates for the current race and gender combination
      rates_for_race_gender <- crude_death_rates_by_race_gender %>% 
        filter(race == race_category & female == gender)
      
      # apply population weights to get age-adjusted death rate
      weighted_rates <- rates_for_race_gender$all_cause_death * age_shares
      age_adjusted_rate <- sum(weighted_rates, na.rm = TRUE) * 100000
      lower_ci <- sum(rates_for_race_gender$ci_l * age_shares, na.rm = TRUE) * 100000
      upper_ci <- sum(rates_for_race_gender$ci_u * age_shares, na.rm = TRUE) * 100000
      
      # append the results to the data frame
      gender_label <- ifelse(gender == "Male", "Male", "Female")
      results <- rbind(results, data.frame(
        Year = target_year,
        Race = race_category,
        Gender = gender_label,
        AgeAdjustedRate = age_adjusted_rate,
        LowerCI = lower_ci,
        UpperCI = upper_ci))
      }
  }
  
  return(results)
}

# initialize a data frame to store the results for all years
all_years_results <- data.frame(Year = integer(), Race = character(), Gender = character(), AgeAdjustedRate = numeric(), LowerCI = numeric(), UpperCI = numeric())

for (year in 1986:2019) {
  year_results <- calculate_age_adjusted_rate_survey(year, df, age_shares)
  year_results <- year_results %>% mutate(Year = year)
  all_years_results <- rbind(all_years_results, year_results)
}

from_1999 <- subset(all_years_results, Year >= 1999)

# plot survey-adjusted aamr results from 1999-2019
survey_aamr_99_19 <- ggplot(from_1999, aes(x = Year, y = AgeAdjustedRate, color = Race, shape = Gender, group = interaction(Race, Gender))) +
  geom_line(size = 1) +
  geom_point(size = 4) +
  # geom_ribbon(aes(ymin = LowerCI, ymax = UpperCI), alpha = 0.2) + 
  # geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.4) +
  labs(title = "All-Cause Age-Adjusted Mortality Rates (1999-2019)",
       x = "Year",
       y = "Age-Adjusted Mortality Rate per 100,000",
       color = "Race",
       linetype = "Gender") +
  theme_classic(base_size = 15) +   
  scale_x_continuous(breaks = seq(1999, 2019, by = 1)) + 
  scale_y_continuous(breaks = seq(0, 2200, by = 200), limits = c(0, 2200)) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = c(0.75, 0.8),
    legend.direction = "horizontal", 
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) + 
  scale_color_manual(values = c("NH Black" = "#df8f44", "NH White" = "#374e55")) +
  scale_shape_manual(values = c("Male" = 16, "Female" = 18)) +
  guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))

# excess aamr
calculate_excess_aamr <- function(all_years_results) {
  excess_aamr_results <- data.frame(Year = integer(), Gender = character(), ExcessAAMR = numeric())
  
  for (year in unique(all_years_results$Year)) {
    for (gender in c("Male", "Female")) {
      aamr_black <- all_years_results %>% filter(Year == year, Race == "NH Black", Gender == gender) %>% select(AgeAdjustedRate, LowerCI, UpperCI)
      aamr_white <- all_years_results %>% filter(Year == year, Race == "NH White", Gender == gender) %>% select(AgeAdjustedRate, LowerCI, UpperCI)
      
      if (nrow(aamr_black) > 0 & nrow(aamr_white) > 0) {
        excess_aamr <- aamr_black$AgeAdjustedRate - aamr_white$AgeAdjustedRate
        
        excess_aamr_results <- rbind(excess_aamr_results, data.frame(
          Year = year,
          Gender = gender,
          ExcessAAMR = excess_aamr
        ))
      }
    }
  }
  
  return(excess_aamr_results)
}

excess_aamr_results <- calculate_excess_aamr(from_1999)

# plot excess AAMR results from 1999-2019
excess_aamr_99_19 <- ggplot(excess_aamr_results, aes(x = Year, y = ExcessAAMR, color = Gender, group = Gender)) +
  geom_line(size = 1) +
  geom_point(size = 4) +
  labs(title = "Excess Mortality Among the Black Population\nCompared to the White Population (1999-2019)",
       x = "Year",
       y = "Excess Age-Adjusted Mortality Rate per 100,000",
       color = "Gender") +
  theme_classic(base_size = 15) +
  scale_x_continuous(breaks = seq(1999, 2019, by = 1)) +
  scale_y_continuous() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = c(0.75, 0.8),
    legend.direction = "horizontal", 
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  scale_color_manual(values = c("Male" = "#df8f44", "Female" = "#374e55"))
