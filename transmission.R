# This model estimates COVID-19 transmission outcomes within Country A of travel from Country B to Country A, given user-defined national infection and
# vaccine profiles, quarantine and testing policies and traveller stay durations. The output of the model is saved in output_prop_tbl and output_per_million_tbl (tibbles).
# See output.txt for descriptions of the output variables.

###########################################################################

### Required libraries ###
# Please load these libraries BEFORE running the rest of the script!

library(data.table)
library(tidyverse, warn.conflicts = FALSE)
library(mvtnorm)
library(bit)
library(foreach)
library(doParallel)
library(doFuture)
library(doRNG)

###########################################################################

### User-specified parameters ###

## Quarantine and testing parameters
# Option 1 for specifying quarantine and testing policies (a csv file with quarantine length and test days, types and delays)
# (please note that tests included on days after the end of quarantine are ignored)
policy_BA_filename <- "policy_covid.csv"
                                          # Name of csv file containing quarantine and testing policy specifications for travel from Country B to Country A;
                                          # (see policy_covid.csv for default data and expected format; policy.txt for description of format)

# Option 2 for specifying a single quarantine and testing policy
# (a named list whose names are the variables StrategyDescription, PreTestName, PreTestDelay, EntryTest, EntryTestDelay, ExitTest, ExitTestDelay, QDays,
# QTestName1, QTestDelay1, ..., QTestNameN, QTestDelayN, ..., where N denotes a test on the Nth day; the order of the variables does not matter, and variables for tests
# which are not included in a policy should be omitted)
# (please note that tests included on days after the end of quarantine, QTestNameK and QTestDelayK with K > QDaysMax, are ignored)
# E.g.:
# policy_BA <- list(StrategyDescription = "Max 14 Day Quarantine with PCR", PreTest = "PCR", PreTestDelay = 3, QDaysMin = 0, QDaysMax = 14, QTestName2 = "PCR", QTestDelay2 = 2,
# QTestName14 = "PCR", QTestDelay14 = 2) denotes a policy with PCR pre-test with 3-day delay, 14-day quarantine, PCR test on day 2 of quarantine with 2-day delay,
# and PCR test on day 14 of quarantine with 2-day delay; StrategyDescription will be added to the output table for user reference.
# Please see policy.txt for further clarification on variables.
policy_BA <- list(StrategyDescription = "Max 14 Day Quarantine with PCR", PreTest = "PCR", PreTestDelay = 3, QDaysMin = 0, QDaysMax = 14, QTestName2 = "PCR", QTestDelay2 = 2,
                  QTestName14 = "PCR", QTestDelay14 = 2)
                           # Named list containing quarantine and testing variables for travel from Country B to A

# Option 3 for specifying quarantine and testing policies by regular patterns (a csv file with min and max quarantine lengths, test frequencies, types and delays)
policy_reg_BA_filename <- "policy_reg_covid.csv"
                                          # Name of csv file containing quarantine and testing policy specifications by regular patterns for travel from Country B to Country A;
                                          # (see policy_reg_covid.csv for default data and expected format; policy_reg.txt for description of format)

# Option 4 for specifying a single quarantine and testing policy by regular patterns
# (a named list whose names are the variables StrategyDescription, PreTestName, PreTestDelay, EntryTest, EntryTestDelay, ExitTest, ExitTestDelay, QDaysMin, QDaysMax,
# QTestName_1, QTestDelay_1, QTestStart_1, QTestEnd_1, QTestFreq_1, QTestName_2, ..., where N denotes a testing pattern; the order of the variables does not matter, and
# variables for tests which are not included in a policy should be omitted)
# E.g.:
# policy_BA <- list(StrategyDescription = "Max 14 Day Quarantine with PCR", PreTest = "PCR", PreTestDelay = 3, QDaysMin = 0, QDaysMax = 14, QTestName_1 = "PCR", QTestDelay_1 = 2,
# QTestStart_1 = 2, QTestEnd_1 = 11, QTestFreq_1 = 3) denotes a policy with PCR pre-test with 3-day delay, up to 14-day quarantine, PCR tests with 2-day delay every 3 days of
# quarantine, starting on day 2 and stopping by day 11; StrategyDescription will be added to the output table for user reference.
# Please see policy_reg.txt for further clarification on variables.
policy_reg_BA <- list(StrategyDescription = "Max 14 Day Quarantine with PCR", PreTest = "PCR", PreTestDelay = 3, QDaysMin = 0, QDaysMax = 14, QTestName_1 = "PCR", QTestDelay_1 = 2,
                      QTestStart_1 = 2, QTestEnd_1 = 11, QTestFreq_1 = 3)
                           # Named list containing quarantine and testing pattern variables for travel from Country B to Country A

# Select option number 1, 2, 3 or 4 for specifying quarantine and testing policies:
policy_option <- 3

## Viral load parameters
viralload_time_params_sympt_filename <- "viralload_time_params_prealpha_sympt.csv"
                        # Name of csv file containing viral load parameters for symptomatic travellers, by day of infection
                        # (see viralload_time_params_alpha_sympt.csv for default data and expected format;
                        # viralload_time_params.txt for description of format)

viralload_time_covar_sympt_filename <- "viralload_time_covar_prealpha_sympt.csv"
                              # Name of csv file containing covariance matrix for viral load parameters for symptomatic travellers, by day of infection
                              # (see viralload_time_covar_alpha_sympt.csv for default data and expected format;
                              # viralload_time_covar.txt for description of format)

viralload_time_params_asympt_filename <- "viralload_time_params_prealpha_asympt.csv"
                        # Name of csv file containing viral load parameters for asymptomatic travellers, by day of infection
                        # (see viralload_time_params_alpha_sympt.csv for default data and expected format;
                        # viralload_time_params.txt for description of format)

viralload_time_covar_asympt_filename <- "viralload_time_covar_prealpha_asympt.csv"
                              # Name of csv file containing covariance matrix for viral load parameters for asymptomatic travellers, by day of infection
                              # (see viralload_time_covar_alpha_asympt.csv for default data and expected format;
                              # viralload_time_covar.txt for description of format)

sens_viralload_filename <- "sensitivity_viralload_covid.csv"
                        # Name of csv file containing sensitivity of each test, by viral load
                        # (see sensitivity_viralload_covid.csv for default data (for COVID - PCR and Antigen) - and expected format;
                        # sensitivity_viralload.txt for description of format)

viralload_initial_min <- 0
                      # Min initial log viral load generated for each traveller
                      # (0 for COVID, 0? for influenza)

viralload_growth_min <- 0
                     # Min log viral load growth per day generated for each traveller
                     # (0 for COVID, 0? for influenza)

viralload_peak_min <- 0
                   # Min peak log viral load generated for each traveller
                   # (4 for COVID, 0? for influenza)

viralload_decline_min <- 0
                      # Min log viral load decline per day generated for each traveller
                      # (0.3 for COVID, 0? for influenza)

## Transmission parameters
transrisk_vl_filename <- "transrisk_vl_covid.csv"
                      # Name of csv file containing transmission risk by viral load
                      # (see transrisk_vl_covid.csv for expected format; 
                      # transrisk_vl.txt for description of format)

## Output parameters
output_prop_filename <- "covid_outputs/output_prealpha.csv"
                                                 # Name of csv file to write output table (for proportions) to (ensure name is not taken to avoid overwriting)
save_output_prop <- TRUE
                                                 # If TRUE, save output table (for proportions) to output_prop_filename; if FALSE, do not save to csv;
                                                 # outputs are accessible in output_prop_tbl

### End of user-specified parameters ###

###########################################################################

### Statistical parameters ###

## Simulation parameters
no_infected <- 100000
                              # No. of infected travellers to simulate per condition, each condition consisting of day of infection (from 1 to max_inf_day) for travellers
max_inf_day <- 28
                              # Max. day of infection (ranging from 1 to max_inf_day) for infected travellers (on day of travel)
simulants_per_tbl <- no_infected * max_inf_day
                              # Number of infected travellers to simulate per table of travellers
init_travellers_seed <- 50140373
                              # Random seed for init_travellers_tbl
test_mode <- FALSE
                              # If test_mode == TRUE, do not run the main loop (to allow testing of functions)
quick_test <- FALSE
                              # If quick_test == TRUE and test_mode == TRUE, run quick_test loop (to load and process input tables for testing)


### End of statistical parameters ###

###########################################################################

### Pre-processing functions (reading csv files) ###

read_viralload_time_params <- function(filename) {
  # Reads csv file with data for viral load parameters by day of infection, and returns table of data in tidied format, indexed by day and test for fast lookup
  # Arguments: filename (string) - name of csv file containing viral load parameters (see viralload_time_params_alpha_sympt.csv and viralload_time_params.txt for format)
  # Returns: single-row data.table with two columns per parameter ("Parameter_Dist", "Parameter_Mean") as well as "Detect_Min" and "Sympt_Prob"
  
  tbl <- fread(filename)
  
  tbl
}

read_viralload_time_covar <- function(filename) {
  # Reads csv file with data for covariances between viral load parameters by day of infection, and returns table of data in tidied format, indexed by day and test for fast lookup
  # Arguments: filename (string) - name of csv file containing covariance matrix for viral load parameters (see viralload_time_covar_alpha_sympt.csv and 
  # viralload_time_covar.txt for format)
  # Returns: covariance matrix between parameters
  
  tbl <- read_csv(filename)
  tbl <- tbl[, 2:ncol(tbl)] %>%
    as.matrix()
  
  tbl
}

read_sens_viralload <- function(filename = sens_viralload_filename) {
  # Reads csv file with data for test sensitivity logistic curve parameters (by viral load), and returns table of data in tidied format, indexed by day and test for fast lookup
  # Arguments: filename (string) - name of csv file containing test sensitivity parameters (see sensitivity_viralload.csv and sensitivity_viralload.txt for format)
  # Returns: data.table with one row per test type ("TestName", "b0", "b1"), with TestName as key
  
  tbl <- fread(filename) %>%
    setkey(TestName)
  
  tbl
}

read_test_sens <- function(filename=test_sens_filename) {
  # Reads csv file with test sensitivity data and returns table of data in tidied format, indexed by day and test for fast lookup
  # Arguments: filename (string) - name of csv file containing test sensitivities by day (see test_sensitivity.csv and test_sensitivity.txt for format)
  # Returns: data.table with one row per condition ("day" x "Test" x "Sensitivity"), with day and Test as keys

  tbl <- read_csv(filename)
  test_names <- names(tbl)[!(names(tbl) %in% "day")]
  tbl <- tbl %>%
    pivot_longer(all_of(test_names), names_to = "Test", values_to = "Sensitivity") %>%
    setDT(keep.rownames = FALSE, key = c("day", "Test"))

  tbl
}

read_policy <- function(option=policy_option, filename, policy_lst) {
  # Reads named list or csv file (depending on option selected) of quarantine and testing policy specifications for travel from Country B to Country A and returns data.table
  # of data, indexed by strategy number for fast lookup
  # Arguments: option (int: 1, 2, 3 or 4) - number denoting the method of input of quarantine and testing policy as described under 'Quarantine and testing parameters';
  # filename (string) - name of csv file containing policy description (see policy.csv and policy.txt or policy_reg.csv and policy_reg.txt for format);
  # policy_lst (named list) - named list of tests and days of quarantine (see Option 2 or 4 under 'Quarantine and testing parameters' for format)
  # Returns: data.table with one row per quarantine and testing policy, columns "StrategyNo", "StrategyDescription", "PreTest", "PreTestDelay", "EntryTest", "EntryTestDelay",
  # "ExitTest", "ExitTestDelay", "QDaysMin", "QDaysMax", "QTest_check_key_XX", "QTest_wait_key_XX" and if get_output_day_diagnosed = TRUE, also "QTest_testday_key_XX" 
  # (where XX is each test type); with StrategyNo as key

  if (option == 1) {                            # Make a data.table from csv
    tbl <- read_csv(filename)
    max_q_days <- max(tbl[["QDaysMax"]])
    tbl_names_char <- c("StrategyDescription", "PreTest", "EntryTest", "ExitTest")
    if (max_q_days > 0) {
      tbl_names_char <- c(tbl_names_char, gsub(" ", "", paste("QTestName", 1:max_q_days, sep = "")))
    }
    extra_tbl_names_char <- setdiff(tbl_names_char, colnames(tbl))
    tbl_names_num <- c("StrategyNo", "PreTestDelay", "EntryTestDelay", "ExitTestDelay", "QDaysMin", "QDaysMax")
    if (max_q_days > 0) {
      tbl_names_num <- c(tbl_names_num, gsub(" ", "", paste("QTestDelay", 1:max_q_days, sep = "")))
    }
    extra_tbl_names_num <- setdiff(tbl_names_num, colnames(tbl))
    tbl[, c(extra_tbl_names_char, extra_tbl_names_num)] <- NA
    tbl <- tbl %>%
      mutate(across(all_of(tbl_names_num), as.integer)) %>%                      # Set columns to integer types
      mutate(across(all_of(tbl_names_char), as.character))                       # Set columns to character types

    qtest_names <- as.character()
    for (col in (tbl %>% select(starts_with("QTestName")))) {
      qtest_names <- unique(c(qtest_names, col))                                # Compile set of unique test names
    }
    qtest_names <- qtest_names[!is.na(qtest_names)]
    for (qtest_name in qtest_names) {
      check_key <- paste("QTest_check_key_", qtest_name, sep = "")
      tbl[[check_key]] <- with(tbl, str_pad("", QDaysMax+1, "left", "0"))               # QTest_check_key_XX: for each policy, make a string of 0s and
                                                                                        # 1s representing days on which tests are performed (1 if performed
                                                                                        # on a given day) - starts the day before quarantine day 1
      wait_key_colname <- paste("QTest_wait_key_", qtest_name, sep = "")
      wait_key <- with(tbl, lapply(QDaysMax + 1, function(x) rep(-1, x)))               # QTest_wait_key_XX: for each policy, make a vector of ints
                                                                                        # representing min quarantine length needed to get results, for test type XX,
                                                                                        # when performed on any given day (starting the day before quarantine day 1)
      for (day in 2:(max_q_days+1)) {
        day_name <- tbl[[paste("QTestName", day-1, sep = "")]]
        day_delay <- replace_na(tbl[[paste("QTestDelay", day-1, sep = "")]], 0)
        actual_test_day <- ifelse(day_name %in% qtest_name & day - day_delay >= 1, day - day_delay, max_q_days + 2)
        substr(tbl[[check_key]], actual_test_day, actual_test_day) <- '1'
        wait_key_indices <- ifelse(day_name %in% qtest_name & day - day_delay >= 1, day - day_delay, 0)
        for (i in 1:length(wait_key_indices)) {
          idx <- wait_key_indices[i]
          if (idx > 0 & idx <= tbl[[i, "QDaysMax"]] + 1) {
            if (wait_key[[i]][[idx]] == -1) {
              wait_key[[i]][[idx]] <- day - 1
            } else {
              wait_key[[i]][[idx]] <- min(wait_key[[i]][[idx]], day - 1)
            }
          }
        }
      }
      wait_key <- sapply(wait_key, function(x) paste(as.character(x), collapse = " "))
      tbl[[wait_key_colname]] <- wait_key
    }

    tbl[c(gsub(" ", "", paste("QTestName", 1:max_q_days, sep = "")), gsub(" ", "", paste("QTestDelay", 1:max_q_days, sep = "")))] <- NULL
                                                                                        # Drop unnecessary columns

    tbl <- tbl %>% setDT(keep.rownames = FALSE, key = "StrategyNo")

    tbl
  }
  else if (option == 2) {                       # Make a single-row data.table with the necessary columns
    tbl_names_char <- c("StrategyDescription", "PreTest", "EntryTest", "ExitTest")
    tbl_names_num <- c("StrategyNo", "PreTestDelay", "EntryTestDelay", "ExitTestDelay", "QDaysMin", "QDaysMax")
    max_q_days <- policy_lst[["QDaysMax"]]
    if (max_q_days > 0) {
      tbl_names_char <- c(tbl_names_char, gsub(" ", "", paste("QTestName", 1:max_q_days, sep = "")))       # Make variables of char type including all days up to q_length
      tbl_names_num <- c(tbl_names_num, gsub(" ", "", paste("QTestDelay", 1:max_q_days, sep = "")))        # Make variables of num type
    }
    vc <- setNames(rep("", length(tbl_names_char)), tbl_names_char)
    vc[tbl_names_char] <- policy_lst[tbl_names_char]
    vc[sapply(vc, is.null)] <- NA                                  # Assign NA to char variables not provided by user
    vn <- setNames(rep(0, length(tbl_names_num)), tbl_names_num)
    vn[tbl_names_num] <- policy_lst[tbl_names_num]
    vn[sapply(vn, is.null)] <- NA                                  # Assign NA to num variables not provided by user
    v <- c(vc, vn)
    tbl <- bind_rows(v)                                                  # Make single row tibble from all variables
    tbl <- tbl %>%
      mutate(across(all_of(tbl_names_num), as.integer)) %>%                      # Set columns to integer types
      mutate(across(all_of(tbl_names_char), as.character))                       # Set columns to character types

    qtest_names <- as.character()
    for (col in (tbl %>% select(starts_with("QTestName")))) {
      qtest_names <- unique(c(qtest_names, col))                                # Compile set of unique test names
    }
    qtest_names <- qtest_names[!is.na(qtest_names)]
    for (qtest_name in qtest_names) {
      check_key <- paste("QTest_check_key_", qtest_name, sep = "")
      tbl[[check_key]] <- with(tbl, str_pad("", QDaysMax+1, "left", "0"))               # QTest_check_key_XX: for each policy, make a string of 0s and
                                                                                        # 1s representing days on which tests are performed (1 if performed
                                                                                        # on a given day) - starts the day before quarantine day 1
      wait_key_colname <- paste("QTest_wait_key_", qtest_name, sep = "")
      wait_key <- with(tbl, lapply(QDaysMax + 1, function(x) rep(-1, x)))               # QTest_wait_key_XX: for each policy, make a vector of ints
                                                                                        # representing min quarantine length needed to get results, for test type XX,
                                                                                        # when performed on any given day (starting the day before quarantine day 1)
      for (day in 2:(max_q_days+1)) {
        day_name <- tbl[[paste("QTestName", day-1, sep = "")]]
        day_delay <- replace_na(tbl[[paste("QTestDelay", day-1, sep = "")]], 0)
        actual_test_day <- ifelse(day_name %in% qtest_name & day - day_delay >= 1, day - day_delay, max_q_days + 2)
        substr(tbl[[check_key]], actual_test_day, actual_test_day) <- '1'
        wait_key_indices <- ifelse(day_name %in% qtest_name & day - day_delay >= 1, day - day_delay, 0)
        for (i in 1:length(wait_key_indices)) {
          idx <- wait_key_indices[i]
          if (idx > 0 & idx <= tbl[[i, "QDaysMax"]] + 1) {
            if (wait_key[[i]][[idx]] == -1) {
              wait_key[[i]][[idx]] <- day - 1
            } else {
              wait_key[[i]][[idx]] <- min(wait_key[[i]][[idx]], day - 1)
            }
          }
        }
      }
      wait_key <- sapply(wait_key, function(x) paste(as.character(x), collapse = " "))
      tbl[[wait_key_colname]] <- wait_key
    }

    tbl[c(gsub(" ", "", paste("QTestName", 1:max_q_days, sep = "")), gsub(" ", "", paste("QTestDelay", 1:max_q_days, sep = "")))] <- NULL
                                                                                        # Drop unnecessary columns

    tbl <- tbl %>% setDT(keep.rownames = FALSE, key = "StrategyNo")

    tbl
  }
  else if (option == 3) {                       # Make a data.table from csv with regular patterns
    tbl <- read_csv(filename)
    max_q_days <- max(tbl[["QDaysMax"]])
    qtestname_cols <- tbl %>%                                                             # Prepare column names for QTest variables
      select(starts_with("QTestName")) %>%
      names()
    qtestdelay_cols <- gsub("Name", "Delay", qtestname_cols)
    qteststart_cols <- gsub("Name", "Start", qtestname_cols)
    qtestend_cols <- gsub("Name", "End", qtestname_cols)
    qtestfreq_cols <- gsub("Name", "Freq", qtestname_cols)
    tbl_names_char <- c("StrategyDescription", "PreTest", "EntryTest", "ExitTest", qtestname_cols)
    extra_tbl_names_char <- setdiff(tbl_names_char, colnames(tbl))
    tbl_names_num <- c("StrategyNo", "PreTestDelay", "EntryTestDelay", "ExitTestDelay", "QDaysMin", "QDaysMax",
                       qtestdelay_cols, qteststart_cols, qtestend_cols, qtestfreq_cols)
    extra_tbl_names_num <- setdiff(tbl_names_num, colnames(tbl))
    tbl[, c(extra_tbl_names_char, extra_tbl_names_num)] <- NA
    tbl <- tbl %>%
      mutate(across(all_of(tbl_names_num), as.integer)) %>%                      # Set columns to integer types
      mutate(across(all_of(tbl_names_char), as.character))                       # Set columns to character types

    qtest_names <- as.character()
    for (col in (tbl %>% select(starts_with("QTestName")))) {
      qtest_names <- unique(c(qtest_names, col))                                # Compile set of unique test names
    }
    qtest_names <- qtest_names[!is.na(qtest_names)]
    for (qtest_name in qtest_names) {
      check_key <- paste("QTest_check_key_", qtest_name, sep = "")
      tbl[[check_key]] <- with(tbl, str_pad("", QDaysMax+1, "left", "0"))               # QTest_check_key_XX: for each policy, make a string of 0s and
                                                                                        # 1s representing days on which tests are performed (1 if performed
                                                                                        # on a given day) - starts the day before quarantine day 1
      wait_key_colname <- paste("QTest_wait_key_", qtest_name, sep = "")
      wait_key <- with(tbl, lapply(QDaysMax + 1, function(x) rep(-1, x)))               # QTest_wait_key_XX: for each policy, make a vector of ints
                                                                                        # representing min quarantine length needed to get results, for test type XX,
                                                                                        # when performed on any given day (starting the day before quarantine day 1)
      for (qtestname_col in qtestname_cols) {
        patternname <- tbl[[qtestname_col]]
        patterndelay <- replace_na(tbl[[gsub("Name", "Delay", qtestname_col)]], 0)
        patternstart <- replace_na(tbl[[gsub("Name", "Start", qtestname_col)]], 0)
        patternend <- replace_na(tbl[[gsub("Name", "End", qtestname_col)]], max_q_days + 2)
        patternfreq <- replace_na(tbl[[gsub("Name", "Freq", qtestname_col)]], 1)
        for (day in 2:(max_q_days+1)) {
          actual_test_day <- ifelse(patternname %in% qtest_name & day - 1 >= patternstart & (day - 1 - patternstart) %% patternfreq == 0 & day - patterndelay >= 1
                                    & day - 1 <= patternend, day - patterndelay, max_q_days + 2)
          substr(tbl[[check_key]], actual_test_day, actual_test_day) <- '1'
          wait_key_indices <- ifelse(patternname %in% qtest_name & day - 1 >= patternstart & (day - 1 - patternstart) %% patternfreq == 0 & day - patterndelay >= 1
                                     & day - 1 <= patternend, day - patterndelay, 0)
          for (i in 1:length(wait_key_indices)) {
            idx <- wait_key_indices[i]
            if (idx > 0 & idx <= tbl[[i, "QDaysMax"]] + 1) {
              if (wait_key[[i]][[idx]] == -1) {
                wait_key[[i]][[idx]] <- day - 1
              } else {
                wait_key[[i]][[idx]] <- min(wait_key[[i]][[idx]], day - 1)
              }
            }
          }
        }
        wait_key <- sapply(wait_key, function(x) paste(as.character(x), collapse = " "))
        tbl[[wait_key_colname]] <- wait_key
      }
    }

    tbl[c(qtestname_cols, qtestdelay_cols, qteststart_cols, qtestend_cols, qtestfreq_cols)] <- NULL
                                                                                        # Drop unnecessary columns

    tbl <- tbl %>% setDT(keep.rownames = FALSE, key = "StrategyNo")

    tbl
  }
  else if (option == 4) {                       # Make a single-row data.table with regular patterns

    qtestname_cols <- names(policy_lst)[startsWith(names(policy_lst), "QTestName")]               # Prepare column names for QTest variables
    qtestdelay_cols <- gsub("Name", "Delay", qtestname_cols)
    qteststart_cols <- gsub("Name", "Start", qtestname_cols)
    qtestend_cols <- gsub("Name", "End", qtestname_cols)
    qtestfreq_cols <- gsub("Name", "Freq", qtestname_cols)
    tbl_names_char <- c("StrategyDescription", "PreTest", "EntryTest", "ExitTest", qtestname_cols)                # Make variables of char type including all days up to q_length
    tbl_names_num <- c("StrategyNo", "PreTestDelay", "EntryTestDelay", "ExitTestDelay", "QDaysMin", "QDaysMax",   # Make variables of num type
                       qtestdelay_cols, qteststart_cols, qtestend_cols, qtestfreq_cols)
    vc <- setNames(rep("", length(tbl_names_char)), tbl_names_char)
    vc[tbl_names_char] <- policy_lst[tbl_names_char]
    vc[sapply(vc, is.null)] <- NA                                  # Assign NA to char variables not provided by user
    vn <- setNames(rep(0, length(tbl_names_num)), tbl_names_num)
    vn[tbl_names_num] <- policy_lst[tbl_names_num]
    vn[sapply(vn, is.null)] <- NA                                  # Assign NA to num variables not provided by user
    v <- c(vc, vn)
    tbl <- bind_rows(v)                                                  # Make single row tibble from all variables
    tbl <- tbl %>%
      mutate(across(all_of(tbl_names_num), as.integer)) %>%                      # Set columns to integer types
      mutate(across(all_of(tbl_names_char), as.character))                       # Set columns to character types

    qtest_names <- as.character()
    for (col in (tbl %>% select(starts_with("QTestName")))) {
      qtest_names <- unique(c(qtest_names, col))                                # Compile set of unique test names
    }
    qtest_names <- qtest_names[!is.na(qtest_names)]
    for (qtest_name in qtest_names) {
      check_key <- paste("QTest_check_key_", qtest_name, sep = "")
      tbl[[check_key]] <- with(tbl, str_pad("", QDaysMax+1, "left", "0"))               # QTest_check_key_XX: for each policy, make a string of 0s and
                                                                                        # 1s representing days on which tests are performed (1 if performed
                                                                                        # on a given day) - starts the day before quarantine day 1
      wait_key_colname <- paste("QTest_wait_key_", qtest_name, sep = "")
      wait_key <- with(tbl, lapply(QDaysMax + 1, function(x) rep(-1, x)))               # QTest_wait_key_XX: for each policy, make a vector of ints
                                                                                        # representing min quarantine length needed to get results, for test type XX,
                                                                                        # when performed on any given day (starting the day before quarantine day 1)
      for (qtestname_col in qtestname_cols) {
        patternname <- tbl[[qtestname_col]]
        patterndelay <- replace_na(tbl[[gsub("Name", "Delay", qtestname_col)]], 0)
        patternstart <- replace_na(tbl[[gsub("Name", "Start", qtestname_col)]], 0)
        patternend <- replace_na(tbl[[gsub("Name", "End", qtestname_col)]], max_q_days + 2)
        patternfreq <- replace_na(tbl[[gsub("Name", "Freq", qtestname_col)]], 1)
        for (day in 2:(max_q_days+1)) {
          actual_test_day <- ifelse(patternname %in% qtest_name & day - 1 >= patternstart & (day - 1 - patternstart) %% patternfreq == 0 & day - patterndelay >= 1
                                    & day - 1 <= patternend, day - patterndelay, max_q_days + 2)
          substr(tbl[[check_key]], actual_test_day, actual_test_day) <- '1'
          wait_key_indices <- ifelse(patternname %in% qtest_name & day - 1 >= patternstart & (day - 1 - patternstart) %% patternfreq == 0 & day - patterndelay >= 1
                                     & day - 1 <= patternend, day - patterndelay, 0)
          for (i in 1:length(wait_key_indices)) {
            idx <- wait_key_indices[i]
            if (idx > 0 & idx <= tbl[[i, "QDaysMax"]] + 1) {
              if (wait_key[[i]][[idx]] == -1) {
                wait_key[[i]][[idx]] <- day - 1
              } else {
                wait_key[[i]][[idx]] <- min(wait_key[[i]][[idx]], day - 1)
              }
            }
          }
        }
        wait_key <- sapply(wait_key, function(x) paste(as.character(x), collapse = " "))
        tbl[[wait_key_colname]] <- wait_key
      }
    }

    tbl[c(qtestname_cols, qtestdelay_cols, qteststart_cols, qtestend_cols, qtestfreq_cols)] <- NULL
                                                                                        # Drop unnecessary columns

    tbl <- tbl %>% setDT(keep.rownames = FALSE, key = "StrategyNo")

    tbl
  }
  else {
    stop("Invalid policy input option: please set policy_option = 1, 2, 3 or 4 under 'Quarantine and testing parameters'")
  }
}

get_test_types <- function(policy_tbl) {
  # Gets all test types from policy table, returning them in a vector
  # Arguments: policy_tbl (data.table) - table of quarantine and testing policies followed by all travellers going from Country B to Country A
  # Returns: vector containing the names of all test types in policy_tbl
  
  test_names <- as.character()
  
  for (test_name in policy_tbl[["PreTest"]]) {
    if (!is.na(test_name)) {
      test_names <- unique(c(test_names, test_name))
    }
  }
  for (test_name in policy_tbl[["EntryTest"]]) {
    if (!is.na(test_name)) {
      test_names <- unique(c(test_names, test_name))
    }
  }
  for (test_name in policy_tbl[["ExitTest"]]) {
    if (!is.na(test_name)) {
      test_names <- unique(c(test_names, test_name))
    }
  }
  qtest_colnames <- names(policy_tbl)[startsWith(names(policy_tbl), "QTest_check_key_")]
  test_names <- unique(c(test_names, gsub("QTest_check_key_", "", qtest_colnames)))
  
  test_names
}

read_transrisk_vl <- function(filename = transrisk_vl_filename) {
  # Reads csv file with data for transmission risk link functions (by viral load), and returns table of data, indexed by link function for fast lookup
  # Arguments: filename (string) - name of csv file containing link function parameters (see transrisk_vl_covid.csv and transrisk_vl.txt for format)
  # Returns: data.table with one row per test type ("Link", "b0", "b1"), with Link as key
  
  tbl <- fread(filename) %>%
    setkey(Link)
  
  tbl
}

### Simulation functions ###

init_travellers <- function(policy_no, policy_tbl_BA, viralload_time_params_sympt_tbl, viralload_time_params_asympt_tbl, viralload_time_covar_sympt_mat, 
                            viralload_time_covar_asympt_mat) {
  # For a given policy in policy_tbl, initialize table of infected travellers' infection / vaccine characteristics and policy information
  # Arguments: policy_no (int) - StrategyNo to initialize table for
  # policy_tbl_BA (data.table) - table of quarantine and testing policies followed by all travellers going from Country B to Country A;
  # viralload_time_params_sympt_tbl (data.table) - table of parameter means and distributions for viral load in symptomatic travellers by day of infection; 
  # viralload_time_params_asympt_tbl (data.table) - table of parameter means and distributions for viral load in asymptomatic travellers by day of infection; 
  # viralload_time_covar_sympt_mat (matrix) - covariance matrix for viral load parameters in symptomatic travellers by day of infection; 
  # viralload_time_covar_asympt_mat (matrix) - covariance matrix for viral load parameters in asymptomatic travellers by day of infection; 
  # Returns: tibble with one row per traveller not filtered out and columns "StrategyNo", "incubation_d", "infectious_d",
  # "total_inf_d", "day_of_inf", "remaining_d", "sympt", "Spread_factor", "QDaysMin_BA", "QDaysMax_BA", "PreTest_BA", "PreTestDelay_BA", "EntryTest_BA", "EntryTestDelay_BA",
  # "ExitTest_BA", "ExitTestDelay_BA", "QTest_check_key_BA_XX", "QTest_wait_key_BA_XX"
  # (where XX is each test type)

  set.seed(init_travellers_seed)
  
  N <- simulants_per_tbl                                                                                    # Total no. of rows in tibble

  travellers_tbl <- tibble(StrategyNo = rep(policy_no, N))                                                  # Create tibble with (max_inf_day) * no_infected rows

  day_of_inf <- sample(max_inf_day, size = N, replace = TRUE, prob = rep(1/max_inf_day, max_inf_day))       # randomly distribute day of infection; uniform distribution
  
  viralload_initial_sympt_mean <- viralload_time_params_sympt_tbl[["VL_Initial_Mean"]]
  viralload_growth_sympt_mean <- viralload_time_params_sympt_tbl[["VL_Growth_Mean"]]
  viralload_peak_sympt_mean <- viralload_time_params_sympt_tbl[["VL_Peak_Mean"]]
  viralload_decline_sympt_mean <- viralload_time_params_sympt_tbl[["VL_Decline_Mean"]]
  viralload_peaktosympt_sympt_mean <- viralload_time_params_sympt_tbl[["Delay_PeakSympt_Mean"]]
  p_symptomatic <- viralload_time_params_sympt_tbl[["Sympt_Prob"]]
  viralload_initial_asympt_mean <- viralload_time_params_asympt_tbl[["VL_Initial_Mean"]]
  viralload_growth_asympt_mean <- viralload_time_params_asympt_tbl[["VL_Growth_Mean"]]
  viralload_peak_asympt_mean <- viralload_time_params_asympt_tbl[["VL_Peak_Mean"]]
  viralload_decline_asympt_mean <- viralload_time_params_asympt_tbl[["VL_Decline_Mean"]]
  viralload_peaktosympt_asympt_mean <- viralload_time_params_asympt_tbl[["Delay_PeakSympt_Mean"]]
  
  means_sympt <- c(viralload_initial_sympt_mean, viralload_growth_sympt_mean, viralload_peak_sympt_mean, viralload_decline_sympt_mean, viralload_peaktosympt_sympt_mean)
  draws_sympt <- rmvnorm(N, mean = means_sympt, sigma = viralload_time_covar_sympt_mat)
  viralload_initial_sympt_is_lognorm <- viralload_time_params_sympt_tbl[["VL_Initial_Dist"]] == "Lognormal"
  if (viralload_initial_sympt_is_lognorm) {
    viralload_initial_sympt <- pmax(exp(draws_sympt[,1]), viralload_initial_min)
  } else {
    viralload_initial_sympt <- pmax(draws_sympt[,1], viralload_initial_min)
  }
  viralload_growth_sympt_is_lognorm <- viralload_time_params_sympt_tbl[["VL_Growth_Dist"]] == "Lognormal"
  if (viralload_growth_sympt_is_lognorm) {
    viralload_growth_sympt <- pmax(exp(draws_sympt[,2]), viralload_growth_min)
  } else {
    viralload_growth_sympt <- pmax(draws_sympt[,2], viralload_growth_min)
  }
  viralload_peak_sympt_is_lognorm <- viralload_time_params_sympt_tbl[["VL_Peak_Dist"]] == "Lognormal"
  if (viralload_peak_sympt_is_lognorm) {
    viralload_peak_sympt <- pmax(exp(draws_sympt[,3]), viralload_peak_min)
  } else {
    viralload_peak_sympt <- pmax(draws_sympt[,3], viralload_peak_min)
  }
  viralload_decline_sympt_is_lognorm <- viralload_time_params_sympt_tbl[["VL_Decline_Dist"]] == "Lognormal"
  if (viralload_decline_sympt_is_lognorm) {
    viralload_decline_sympt <- pmax(exp(draws_sympt[,4]), viralload_decline_min)
  } else {
    viralload_decline_sympt <- pmax(draws_sympt[,4], viralload_decline_min)
  }
  viralload_peaktosympt_sympt_is_lognorm <- viralload_time_params_sympt_tbl[["Delay_PeakSympt_Dist"]] == "Lognormal"
  if (viralload_peaktosympt_sympt_is_lognorm) {
    viralload_peaktosympt_sympt <- round(exp(draws_sympt[,5]))
  } else {
    viralload_peaktosympt_sympt <- round(draws_sympt[,5])
  }
  
  means_asympt <- c(viralload_initial_asympt_mean, viralload_growth_asympt_mean, viralload_peak_asympt_mean, viralload_decline_asympt_mean, viralload_peaktosympt_asympt_mean)
  draws_asympt <- rmvnorm(N, mean = means_asympt, sigma = viralload_time_covar_asympt_mat)
  viralload_initial_asympt_is_lognorm <- viralload_time_params_asympt_tbl[["VL_Initial_Dist"]] == "Lognormal"
  if (viralload_initial_asympt_is_lognorm) {
    viralload_initial_asympt <- pmax(exp(draws_asympt[,1]), viralload_initial_min)
  } else {
    viralload_initial_asympt <- pmax(draws_asympt[,1], viralload_initial_min)
  }
  viralload_growth_asympt_is_lognorm <- viralload_time_params_asympt_tbl[["VL_Growth_Dist"]] == "Lognormal"
  if (viralload_growth_asympt_is_lognorm) {
    viralload_growth_asympt <- pmax(exp(draws_asympt[,2]), viralload_growth_min)
  } else {
    viralload_growth_asympt <- pmax(draws_asympt[,2], viralload_growth_min)
  }
  viralload_peak_asympt_is_lognorm <- viralload_time_params_asympt_tbl[["VL_Peak_Dist"]] == "Lognormal"
  if (viralload_peak_asympt_is_lognorm) {
    viralload_peak_asympt <- pmax(exp(draws_asympt[,3]), viralload_peak_min)
  } else {
    viralload_peak_asympt <- pmax(draws_asympt[,3], viralload_peak_min)
  }
  viralload_decline_asympt_is_lognorm <- viralload_time_params_asympt_tbl[["VL_Decline_Dist"]] == "Lognormal"
  if (viralload_decline_asympt_is_lognorm) {
    viralload_decline_asympt <- pmax(exp(draws_asympt[,4]), viralload_decline_min)
  } else {
    viralload_decline_asympt <- pmax(draws_asympt[,4], viralload_decline_min)
  }
  viralload_peaktosympt_asympt_is_lognorm <- viralload_time_params_asympt_tbl[["Delay_PeakSympt_Dist"]] == "Lognormal"
  if (viralload_peaktosympt_asympt_is_lognorm) {
    viralload_peaktosympt_asympt <- round(exp(draws_asympt[,5]))
  } else {
    viralload_peaktosympt_asympt <- round(draws_asympt[,5])
  }
  
  days_to_peak_sympt <- ceiling((viralload_peak_sympt - viralload_initial_sympt)/viralload_growth_sympt)
  viralload_peak_day_sympt <- 1 + days_to_peak_sympt
  days_to_end_sympt <- ceiling((viralload_peak_sympt - viralload_initial_sympt)/viralload_decline_sympt)
  viralload_inf_end_day_sympt <- viralload_peak_day_sympt + days_to_end_sympt
  
  days_to_peak_asympt <- ceiling((viralload_peak_asympt - viralload_initial_asympt)/viralload_growth_asympt)
  viralload_peak_day_asympt <- 1 + days_to_peak_asympt
  days_to_end_asympt <- ceiling((viralload_peak_asympt - viralload_initial_asympt)/viralload_decline_asympt)
  viralload_inf_end_day_asympt <- viralload_peak_day_asympt + days_to_end_asympt

  travellers_tbl <- travellers_tbl %>%
    mutate(sympt = rbinom(N, 1, p_symptomatic)) %>%                                                         # sympt (values 0, 1): 0 if asymptomatic, 1 if symptomatic
    mutate(viralload_growth = ifelse(sympt, viralload_growth_sympt, viralload_growth_asympt)) %>%
    mutate(viralload_peak = ifelse(sympt, viralload_peak_sympt, viralload_peak_asympt)) %>%
    mutate(viralload_decline = ifelse(sympt, viralload_decline_sympt, viralload_decline_asympt)) %>%
    mutate(viralload_peak_day = ifelse(sympt, viralload_peak_day_sympt, viralload_peak_day_asympt)) %>%
    mutate(viralload_peaktosympt = ifelse(sympt, viralload_peaktosympt_sympt, viralload_peaktosympt_asympt)) %>%
    mutate(viralload_inf_end_day = ifelse(sympt, viralload_inf_end_day_sympt, viralload_inf_end_day_asympt)) %>%
    mutate(incubation_d = pmin(pmax(viralload_peak_day + viralload_peaktosympt - 1, 0), 
                               viralload_inf_end_day)) %>%                                                  # incubation_d (log-normally distributed): incubation period = Inc
                                                                                                            # days, from days 1 to Inc
    mutate(infectious_d = pmax(viralload_inf_end_day - incubation_d, 0)) %>%                                # infectious_d (normally distributed): infectious period = Inf
                                                                                                            # days, from days Inc + 1 to Inc + Inf
    mutate(total_inf_d = incubation_d + infectious_d) %>%                                                   # total_inf_d (incubation_d + infectious_d): full infection period
    mutate(day_of_inf = day_of_inf) %>%                                                                     # day_of_inf (values from 1 to max_inf_day):
                                                                                                            # day of infection; day_of_inf = n means a traveller was infected
                                                                                                            # n-1 days ago (on day of travel)
    mutate(remaining_d = total_inf_d - day_of_inf)                                                          # remaining_d (total_inf_d - day_of_inf): remaining days to
                                                                                                            # recovery


  travellers_tbl[["QDaysMin_BA"]] <- policy_tbl_BA[.(policy_no),][["QDaysMin"]]                             # Add BA policy information
  travellers_tbl[["QDaysMax_BA"]] <- policy_tbl_BA[.(policy_no),][["QDaysMax"]]
  travellers_tbl[["PreTest_BA"]] <- policy_tbl_BA[.(policy_no),][["PreTest"]]
  travellers_tbl[["PreTestDelay_BA"]] <- policy_tbl_BA[.(policy_no),][["PreTestDelay"]]
  travellers_tbl[["EntryTest_BA"]] <- policy_tbl_BA[.(policy_no),][["EntryTest"]]
  travellers_tbl[["EntryTestDelay_BA"]] <- policy_tbl_BA[.(policy_no),][["EntryTestDelay"]]
  travellers_tbl[["ExitTest_BA"]] <- policy_tbl_BA[.(policy_no),][["ExitTest"]]
  travellers_tbl[["ExitTestDelay_BA"]] <- policy_tbl_BA[.(policy_no),][["ExitTestDelay"]]

  qtest_checkkeys <- policy_tbl_BA %>% select(starts_with("QTest_check_key")) %>% names()                 # Compile set of distinct quarantine test names in policy table
  for (qtest_checkkey in qtest_checkkeys) {                                                               # "QTest_check_key_BA_XX": check key for BA quarantine test
    BA_checkkey <- gsub("QTest_check_key", "QTest_check_key_BA", qtest_checkkey)
    BA_waitkey <- gsub("QTest_check_key", "QTest_wait_key_BA", qtest_checkkey)
    qtest_waitkey <- gsub("QTest_check_key", "QTest_wait_key", qtest_checkkey)
    if (!is.null(policy_tbl_BA[.(policy_no),][[qtest_checkkey]])) {
      travellers_tbl[[BA_checkkey]] <- policy_tbl_BA[.(policy_no),][[qtest_checkkey]]
      travellers_tbl[[BA_waitkey]] <- policy_tbl_BA[.(policy_no),][[qtest_waitkey]]
    } else {
      travellers_tbl[[BA_checkkey]] <- str_pad("", travellers_tbl[["QDaysMax_BA"]] + 1, "left", "0")
      travellers_tbl[[BA_waitkey]] <- with(travellers_tbl, paste(rep("-1", QDaysMax_BA), collapse = " "))
    }
  }

  travellers_tbl
}

filter_sympt_infectious <- function(travellers_tbl) {
  # Filter out infected travellers that are symptomatic and in their infectious period on day of travel (current day_of_inf)
  # Arguments: travellers_tbl (tibble) - initialized table of travellers' infection / vaccine characteristics and policy information
  # Returns: travellers_tbl, with aforementioned travellers (rows) removed
  
  travellers_tbl <- travellers_tbl %>% filter(sympt == 0 | day_of_inf <= incubation_d)                    # Keep only travellers who are asymptomatic or have not left
                                                                                                          # incubation period
  
  travellers_tbl
}

filter_recovered <- function(travellers_tbl) {
  # Filter out infected travellers that have recovered by the day of travel (current day_of_inf)
  # Arguments: travellers_tbl (tibble) - initialized table of travellers' infection / vaccine characteristics and policy information
  # Returns: travellers_tbl, with aforementioned travellers (rows) removed
  
  travellers_tbl <- travellers_tbl %>% filter(day_of_inf <= total_inf_d)                                  # Keep only travellers who have not recovered by day of travel
  
  travellers_tbl
}

get_viral_loads <- function(viralload_growth, viralload_decline, viralload_peak, viralload_peak_day, day) {
  # Compute vector of viral loads from given viral load parameters and days of infection
  # Arguments: viralload_growth (vector of doubles) - viral load parameter viralload_growth;
  # viralload_decline (vector of doubles) - viral load parameter viralload_decline;
  # viralload_peak (vector of doubles) - viral load parameter viralload_peak;
  # viralload_peak_day (vector of ints) - viral load parameter viralload_peak_day;
  # day (vector of ints) - day of infection to receive viral load for
  # Returns: vector of doubles, representing viral loads of travellers on given days
  
  N <- length(viralload_growth)
  
  viral_loads <- ifelse(day <= viralload_peak_day,
                        ifelse(day < 1, 0, pmax(viralload_peak - viralload_growth * (viralload_peak_day - day), 0)),
                        pmax(viralload_peak - viralload_decline * (day - viralload_peak_day), 0))
  viral_loads
}

get_test_probability <- function(sens_viralload_tbl, viral_loads, test_type) {
  # Compute vector of test positive probabilities from test sensitivity table given viral loads and test type
  # Arguments: sens_viralload_tbl (data.table) - table of test sensitivities by viral load;
  # viral_loads (vector of doubles) - viral loads to obtain sensitivities for;
  # test_type (string) - name of test to get sensitivities for
  # Returns: vector of test sensitivities (probabilities of testing positive)
  
  b0 <- sens_viralload_tbl[.(test_type),][["b0"]]
  b1 <- sens_viralload_tbl[.(test_type),][["b1"]]
  viral_loads <- ifelse(viral_loads > 0, viral_loads, 0)
  y <- b0 + b1 * viral_loads
  probs <- exp(y) / (1 + exp(y))
  
  probs
}

draw_test_results <- function(travellers_tbl, sens_viralload_tbl) {
  # Pre-sample test results that would be received on each possible test day for each infected traveller in a given table of travellers
  # Arguments: travellers_tbl (tibble) - initialized table of travellers' infection / vaccine characteristics and policy information;
  # sens_viralload_tbl (data.table) - table of test sensitivities by viral load
  # Returns: travellers_tbl, with added columns "PreTest_results_BA", "EntryTest_results_BA", "ExitTest_results_BA", "QTest_results_BA_XX" (where XX is each test type)

  N <- nrow(travellers_tbl)
  day_of_inf <- travellers_tbl[["day_of_inf"]]
  total_inf_d <- travellers_tbl[["total_inf_d"]]
  incubation_d <- travellers_tbl[["incubation_d"]]
  viralload_growth <- travellers_tbl[["viralload_growth"]]
  viralload_decline <- travellers_tbl[["viralload_decline"]]
  viralload_peak <- travellers_tbl[["viralload_peak"]]
  viralload_peak_day <- travellers_tbl[["viralload_peak_day"]]
  qdaysmax_BA <- travellers_tbl[["QDaysMax_BA"]][1]

  pretest_BA <- travellers_tbl[["PreTest_BA"]][1]
  if (!is.na(pretest_BA)) {                                                                               # Add PreTest_BA_results only if there is a PreTest_BA
    pretest_BA_delay <- travellers_tbl[["PreTestDelay_BA"]][1]
    pretest_BA_day <- day_of_inf - pretest_BA_delay
    pretest_BA_viral_loads <- get_viral_loads(viralload_growth, viralload_decline, viralload_peak, viralload_peak_day, pretest_BA_day)
    pretest_BA_probabilities <- get_test_probability(sens_viralload_tbl, pretest_BA_viral_loads, pretest_BA)
    pretest_BA_probabilities <- ifelse(day_of_inf < pretest_BA_delay, 0, pretest_BA_probabilities)
    pretest_BA_results <- rbinom(rep(1, N), rep(1, N), pretest_BA_probabilities)
  } else {
    pretest_BA_results <- rep(0, N)
  }
  travellers_tbl[["PreTest_BA_results"]] <- pretest_BA_results                                            # Add PreTest_BA_results (tested on day of travel - PreTestDelay_BA)

  entrytest_BA <- travellers_tbl[["EntryTest_BA"]][1]
  if (!is.na(entrytest_BA)) {                                                                             # Add EntryTest_BA_results only if there is an EntryTest_BA
    entrytest_BA_day <- day_of_inf
    entrytest_BA_viral_loads <- get_viral_loads(viralload_growth, viralload_decline, viralload_peak, viralload_peak_day, entrytest_BA_day)
    entrytest_BA_probabilities <- get_test_probability(sens_viralload_tbl, entrytest_BA_viral_loads, entrytest_BA)
    entrytest_BA_results <- rbinom(rep(1, N), rep(1, N), entrytest_BA_probabilities)
    travellers_tbl[["EntryTest_BA_results"]] <- entrytest_BA_results                                      # Add EntryTest_BA_results (tested on day of travel - EntryTestDelay_BA)
  } else {
    entrytest_BA_results <- rep(0, N)
  }
  travellers_tbl[["EntryTest_BA_results"]] <- entrytest_BA_results                                        # Add EntryTest_BA_results (tested on day of travel - EntryTestDelay_BA)

  exittest_BA <- travellers_tbl[["ExitTest_BA"]][1]
  if (!is.na(exittest_BA)) {                                                                              # Add ExitTest_BA_results only if there is an ExitTest_BA
    exittest_BA_delay <- travellers_tbl[["ExitTestDelay_BA"]][1]
    exittest_BA_results <- rep("", N)
    for (day in 0:qdaysmax_BA) {
      exittest_BA_day <- day_of_inf + day - exittest_BA_delay
      exittest_BA_viral_loads <- get_viral_loads(viralload_growth, viralload_decline, viralload_peak, viralload_peak_day, exittest_BA_day)
      exittest_BA_probabilities <- get_test_probability(sens_viralload_tbl, exittest_BA_viral_loads, exittest_BA)
      exittest_BA_probabilities <- ifelse(day_of_inf + day - exittest_BA_delay > total_inf_d | day - exittest_BA_delay < 1, 0, exittest_BA_probabilities)
      exittest_BA_results <- paste(exittest_BA_results, as.character(rbinom(rep(1, N), rep(1, N), exittest_BA_probabilities)), sep = "")
    }
  } else {
    exittest_BA_results <- strrep("0", rep(qdaysmax_BA + 1, N))
  }
  travellers_tbl[["ExitTest_BA_results"]] <- exittest_BA_results                                          # Add ExitTest_BA_results (range of days starts on 
                                                                                                          # day of travel - ExitTestDelay_BA)

  qtest_checkkeys_BA <- travellers_tbl %>% select(starts_with("QTest_check_key_BA")) %>% names()
  qtest_types_BA <- gsub("QTest_check_key_BA_", "", qtest_checkkeys_BA)
  for (qtest_type in qtest_types_BA) {
    qtest_BA_results <- rep("", N)
    for (day in 0:qdaysmax_BA) {
      qtest_BA_day <- day_of_inf + day
      qtest_BA_viral_loads <- get_viral_loads(viralload_growth, viralload_decline, viralload_peak, viralload_peak_day, qtest_BA_day)
      qtest_BA_probabilities <- get_test_probability(sens_viralload_tbl, qtest_BA_viral_loads, qtest_type)
      qtest_BA_probabilities <- ifelse(day_of_inf + day > total_inf_d, 0, qtest_BA_probabilities)
      qtest_BA_results <- paste(qtest_BA_results, as.character(rbinom(rep(1, N), rep(1, N), qtest_BA_probabilities)), sep = "")
    }
    qtest_BA_colname <- paste("QTest_BA_results_", qtest_type, sep = "")
    travellers_tbl[[qtest_BA_colname]] <- qtest_BA_results                                                # Add QTest_BA_results_XX (starts on day of travel)
  }

  travellers_tbl
}

test_outcomes <- function(travellers_tbl, test_names_BA) {
  # Apply entry tests, quarantines and quarantine tests to table of travellers infected at B, by quarantine length combinations; 
  # track numbers of tests used (by test type)
  # Arguments: travellers_tbl (tibble) - table of pretested travellers' infection / vaccine characteristics and strategies being tested;
  # test_names_BA (character vector) - vector of test names for tests to be used on travel from B to A
  # Returns: named list with names "Sympt_filtered_BA_j", "PreTest_diagnosed_BA_j", "Diagnosed_BA_j", "Missed_BA_j", "TestsUsed_PreTest_BA_XX_j", 
  # "TestsUsed_EntryTest_BA_XX_j", "TestsUsed_QTest_BA_XX_j", "TestsUsed_ExitTest_BA_XX_j", for quarantine lengths j from QDaysMin_BA to QDaysMax_BA, test types XX
  
  N <- nrow(travellers_tbl)
  qdaysmin_BA <- travellers_tbl[["QDaysMin_BA"]][1]
  qdaysmax_BA <- travellers_tbl[["QDaysMax_BA"]][1]
  qtest_results_colnames <- travellers_tbl %>% select(starts_with("QTest_BA_results_")) %>% names()
  sympt <- travellers_tbl[["sympt"]]
  day_of_inf <- travellers_tbl[["day_of_inf"]]
  incubation_d <- travellers_tbl[["incubation_d"]]
  total_inf_d <- travellers_tbl [["total_inf_d"]]
  final_inf_d <- total_inf_d - day_of_inf
  
  outcomes <- list()
  
  for (j in qdaysmin_BA:qdaysmax_BA) {
    for (test_name in test_names_BA) {
      count_name <- paste("TestsUsed_PreTest_BA_", test_name, "_", j, sep = "")
      outcomes[[count_name]] <- 0L                                                                        # TestsUsed_PreTest_BA_XX_j: Number of tests of type XX used for 
                                                                                                          # BA pretests, for BA quarantine length j
      count_name <- paste("TestsUsed_EntryTest_BA_", test_name, "_", j, sep = "")
      outcomes[[count_name]] <- 0L                                                                        # TestsUsed_EntryTest_BA_XX_j: Number of tests of type XX used for 
                                                                                                          # BA entry tests, for BA quarantine length j
      count_name <- paste("TestsUsed_QTest_BA_", test_name, "_", j, sep = "")
      outcomes[[count_name]] <- 0L                                                                        # TestsUsed_QTest_BA_XX_j: Number of tests of type XX used for 
                                                                                                          # BA quarantine tests, for BA quarantine length j
      count_name <- paste("TestsUsed_ExitTest_BA_", test_name, "_", j, sep = "")
      outcomes[[count_name]] <- 0L                                                                        # TestsUsed_ExitTest_BA_XX_j: Number of tests of type XX used for 
                                                                                                          # BA exit tests, for BA quarantine length j
    }
    
    keep_checking_name <- paste("KeepChecking_", j, sep = "")
    outcomes[[keep_checking_name]] <- as.bit(rep(1, N))                                                   # KeepChecking_j: Auxiliary column indicating that traveller
                                                                                                          # has not been diagnosed or recovered yet
  }
  
  for (j in qdaysmin_BA:qdaysmax_BA) {
    sympt_filtered_name <- paste("Sympt_filtered_BA_", j, sep = "")
    keep_checking_name <- paste("KeepChecking_", j, sep = "")
    outcomes[[sympt_filtered_name]] <- as.bit(sympt & (day_of_inf > incubation_d))
    outcomes[[keep_checking_name]] <- outcomes[[keep_checking_name]] & !outcomes[[sympt_filtered_name]]
  }
  
  pretest_diagnosed <- travellers_tbl[["PreTest_BA_results"]] == "1"
  pretest_name <- travellers_tbl[["PreTest_BA"]][1]
  for (j in qdaysmin_BA:qdaysmax_BA) {
    pretest_diagnosed_name <- paste("PreTest_diagnosed_BA_", j, sep = "")
    keep_checking_name <- paste("KeepChecking_", j, sep = "")
    outcomes[[pretest_diagnosed_name]] <- as.bit(pretest_diagnosed) & outcomes[[keep_checking_name]]
                                                                                                          # PreTest_diagnosed_BA_j: Has been diagnosed by pretest on
                                                                                                          # travel from B to A, for quarantine length j
    if (!is.na(pretest_name)) {
      increment <- sum(as.integer(outcomes[[keep_checking_name]]))
      count_name <- paste("TestsUsed_PreTest_BA_", pretest_name, "_", j, sep = "")
      outcomes[[count_name]] <- outcomes[[count_name]] + increment                                        # increment count for pretest's test type
    }
    outcomes[[keep_checking_name]] <- outcomes[[keep_checking_name]] & !outcomes[[pretest_diagnosed_name]]
  }
  
  entrytest_diagnosed <- travellers_tbl[["EntryTest_BA_results"]] == "1"
  entrytest_name <- travellers_tbl[["EntryTest_BA"]][1]
  entrytest_delay <- travellers_tbl[["EntryTestDelay_BA"]][1]
  for (j in qdaysmin_BA:qdaysmax_BA) {
    diagnosed_name <- paste("Diagnosed_BA_", j, sep = "")
    keep_checking_name <- paste("KeepChecking_", j, sep = "")
    outcomes[[diagnosed_name]] <- as.bit(rep(0, N))
                                                                                                          # Diagnosed_BA_j: Has been diagnosed by tests or symptoms during 
                                                                                                          # quarantine on travel from B to A, for quarantine length j
    if (!is.na(entrytest_name)) {
      outcomes[[diagnosed_name]] <- as.bit(entrytest_diagnosed & j >= entrytest_delay) & outcomes[[keep_checking_name]]
      increment <- sum(as.integer(as.bit(j >= entrytest_delay) & outcomes[[keep_checking_name]]))
      count_name <- paste("TestsUsed_EntryTest_BA_", entrytest_name, "_", j, sep = "")
      outcomes[[count_name]] <- outcomes[[count_name]] + increment                                        # increment count for entry test's test type
    }
    outcomes[[keep_checking_name]] <- outcomes[[keep_checking_name]] & !outcomes[[diagnosed_name]]
  }
  
  earliest_day_diagnosed <- ifelse(as.integer(entrytest_diagnosed), entrytest_delay, qdaysmax_BA + 1)     # earliest_day_diagnosed: first day of quarantine on which 
                                                                                                          # travellers are diagnosed (inclusive of day 0 by entry test)
  
  for (day in 0:qdaysmax_BA) {                                                                            # compute earliest_day_diagnosed for each traveller
    for (colname in qtest_results_colnames) {
      qtest_results <- travellers_tbl[[colname]]
      checkkey_colname <- gsub("QTest_BA_results_", "QTest_check_key_BA_", colname)
      waitkey_colname <- gsub("QTest_BA_results_", "QTest_wait_key_BA_", colname)
      qtest_checkkey <- travellers_tbl[[checkkey_colname]][1]
      qtest_waitkey <- as.integer(unlist(strsplit(travellers_tbl[[waitkey_colname]][1], " ")))
      result_today <- substr(qtest_results, day + 1, day + 1) == "1"
      checkkey_today <- substr(qtest_checkkey, day + 1, day + 1) == "1"
      waitkey_today <- qtest_waitkey[day + 1]
      earliest_day_diagnosed <- ifelse(result_today & checkkey_today & waitkey_today >= 0, pmin(waitkey_today, earliest_day_diagnosed), earliest_day_diagnosed)
                                                                                                          # Earliest day diagnosed is no later than a day on which
                                                                                                          # positive test result is obtained
    }                                                                                                     
    earliest_day_diagnosed <- ifelse(sympt & (day_of_inf + day > incubation_d), pmin(day, earliest_day_diagnosed), earliest_day_diagnosed) 
  }                                                                                                       # Earliest day diagnosed can include day that symptoms emerge
  
  exittest_name <- travellers_tbl[["ExitTest_BA"]][1]
  exittest_delay <- travellers_tbl[["ExitTestDelay_BA"]][1]
  for (j in qdaysmin_BA:qdaysmax_BA) {
    keep_checking_name <- paste("KeepChecking_", j, sep = "")
    for (colname in qtest_results_colnames) {
      testname <- gsub("QTest_BA_results_", "", colname)
      checkkey_colname <- paste("QTest_check_key_BA_", testname, sep = "")
      waitkey_colname <- paste("QTest_wait_key_BA_", testname, sep = "")
      qtest_checkkey <- travellers_tbl[[checkkey_colname]][1]
      qtest_waitkey <- as.integer(unlist(strsplit(travellers_tbl[[waitkey_colname]][1], " ")))
      increment <- 0L
      for (today in 0:j) {
        checkkey_today <- substr(qtest_checkkey, today + 1, today + 1) == "1"
        if (!is.na(entrytest_name)) {
          if (testname == entrytest_name & today == 0 & entrytest_delay <= j) {
            checkkey_today <- FALSE
          }
        }
        if (!is.na(exittest_name)) {
          if (testname == exittest_name & today == j - exittest_delay & today >= 1) {
            checkkey_today <- FALSE
          }
        }
        waitkey_today <- qtest_waitkey[today + 1]
        increment <- increment + as.integer((today <= earliest_day_diagnosed) & checkkey_today & (waitkey_today >= 0 & waitkey_today <= j))
                                                                                                          # increment test count for tests performed each day; test results
                                                                                                          # must fall within quarantine length for test to be performed and  
                                                                                                          # tests cannot be past earliest_day_diagnosed
      }
      total_increment = sum(increment * as.integer(outcomes[[keep_checking_name]]))
      count_name <- paste("TestsUsed_QTest_BA_", testname, "_", j, sep = "")
      outcomes[[count_name]] <- outcomes[[count_name]] + total_increment
    }
    
    exittest_diagnosed <- (substr(travellers_tbl[["ExitTest_BA_results"]], j + 1, j + 1) == "1") & (j - exittest_delay >= 1)
    if (!is.na(exittest_name)) {
      increment <- sum(as.integer(outcomes[[keep_checking_name]] & as.bit(j - exittest_delay <= earliest_day_diagnosed) &
                                    as.bit(j - exittest_delay >= 1)))
      count_name <- paste("TestsUsed_ExitTest_BA_", exittest_name, "_", j, sep = "")
      outcomes[[count_name]] <- outcomes[[count_name]] + increment
    }
    diagnosed_name <- paste("Diagnosed_BA_", j, sep = "")
    outcomes[[diagnosed_name]] <- outcomes[[diagnosed_name]] | (outcomes[[keep_checking_name]] & (as.bit(j >= earliest_day_diagnosed) | exittest_diagnosed))
    outcomes[[keep_checking_name]] <- outcomes[[keep_checking_name]] & !outcomes[[diagnosed_name]] & as.bit(j < final_inf_d)
    missed_name <- paste("Missed_BA_", j, sep = "")
    outcomes[[missed_name]] <- outcomes[[keep_checking_name]]                                             # Missed_BA_j: Infected traveller was missed by quarantine and
                                                                                                          # tests on travel from B to A (and did not recover)
    recovered_name <- paste("Recovered_BA_", j, sep = "")
    outcomes[[recovered_name]] <- !(outcomes[[pretest_name]] | outcomes[[entrytest_name]] | outcomes[[diagnosed_name]] | outcomes[[missed_name]])
  }
  
  for (j in qdaysmin_BA:qdaysmax_BA) {
    keep_checking_name <- paste("KeepChecking_", j, sep = "")
    outcomes <- outcomes[names(outcomes) != keep_checking_name]
  }
  
  outcomes
}

logit_transrisk <- function(viral_loads, transrisk_vls) {
  # Calculates transmission risk on a given day from viral loads and table of link function parameters
  # Arguments: viral_loads (vector of doubles) - log viral loads on a given day;
  # transrisk_vls (data.table) - table of transmission risk link functions by viral loads
  # Returns: vector of estimated transmission risks for a given day
  
  b0 <- transrisk_vls[.("Logit"),][["b0"]]
  b1 <- transrisk_vls[.("Logit"),][["b1"]]
  transrisks <- 1/(1 + exp(-(b0 + b1 * viral_loads)))
  
  transrisks
}

get_total_tr <- function(travellers_tbl, day, transrisk_vls, link) {
  # Calculates vector of total transmission risk based on log viral loads for each traveller in a table from given days of infection
  # Arguments: travellers_tbl (tibble) - table of travellers' infection characteristics and strategies being tested; 
  # day (vector of ints) - day of infection to receive total transmission risk from; 
  # transrisk_vls (data.table) - table of transmission risk link functions by viral loads; 
  # link (function) - link function to calculate transmission risks from viral loads
  # Returns: vector of doubles, representing total transmission risk from given day onwards
  
  N <- nrow(travellers_tbl)
  viralload_growth <- travellers_tbl[["viralload_growth"]]
  viralload_peak <- travellers_tbl[["viralload_peak"]]
  viralload_decline <- travellers_tbl[["viralload_decline"]]
  viralload_peak_day <- travellers_tbl[["viralload_peak_day"]]
  total_inf_d <- travellers_tbl[["total_inf_d"]]
  
  remaining_days <- total_inf_d - day
  last_d <- max(remaining_days)
  total_tr <- rep(0, N)
  
  for (d in 0:last_d) {
    today <- day + d
    today_vl <- get_viral_loads(viralload_growth, viralload_decline, viralload_peak, viralload_peak_day, today)
    today_tr <- link(today_vl, transrisk_vls)
    total_tr <- total_tr + today_tr
  }
  
  total_tr
}

infdays_trs <- function(travellers_tbl, transrisk_vls) {
  # Creates table of days of infection remaining for travellers after quarantine (with symptoms and without), as well as total transmission risk over 
  # days of remaining infection after quarantine (with symptoms and without)
  # Arguments: travellers_tbl (tibble) - table of travellers' infection characteristics and strategies being tested; 
  # transrisk_vls (data.table) - table of transmission risk link functions by viral loads;
  # Returns: named list with "PostQ_InfDays_Sympt_j", "PostQ_InfDays_NoSympt_j", "PostQ_TR_Sympt_j" and "PostQ_TR_NoSympt_j", for BA quarantine lengths j from QDaysMin_BA to 
  # QDaysMax_BA
  
  N <- nrow(travellers_tbl)
  qdaysmin_BA <- travellers_tbl[["QDaysMin_BA"]][1]
  qdaysmax_BA <- travellers_tbl[["QDaysMax_BA"]][1]
  sympt <- travellers_tbl[["sympt"]]
  day_of_inf <- travellers_tbl[["day_of_inf"]]
  incubation_d <- travellers_tbl[["incubation_d"]]
  total_inf_d <- travellers_tbl [["total_inf_d"]]
  
  infdays_trs_lst <- list()
  for (j in qdaysmin_BA:qdaysmax_BA) {
    infdays_sympt <- ifelse(sympt, pmin(total_inf_d - incubation_d, pmax(total_inf_d - (day_of_inf + j), 0)), 0)
    infdays_sympt_name <- paste("PostQ_InfDays_Sympt_", j, sep = "")
    infdays_trs_lst[[infdays_sympt_name]] <- infdays_sympt
    infdays_nosympt <- ifelse(sympt, pmax(incubation_d - (day_of_inf + j), 0), pmax(total_inf_d - (day_of_inf + j), 0))
    infdays_nosympt_name <- paste("PostQ_InfDays_NoSympt_", j, sep = "")
    infdays_trs_lst[[infdays_nosympt_name]] <- infdays_nosympt
    tr_sympt_logit <- ifelse(sympt, get_total_tr(travellers_tbl, pmax(day_of_inf + j, incubation_d + 1), transrisk_vls, logit_transrisk), 0)
    tr_sympt_logit_name <- paste("PostQ_TR_Sympt_Logit_", j, sep = "")
    infdays_trs_lst[[tr_sympt_logit_name]] <- tr_sympt_logit
    tr_nosympt_logit <- get_total_tr(travellers_tbl, day_of_inf + j, transrisk_vls, logit_transrisk) - tr_sympt_logit
    tr_nosympt_logit_name <- paste("PostQ_TR_NoSympt_Logit_", j, sep = "")
    infdays_trs_lst[[tr_nosympt_logit_name]] <- tr_nosympt_logit
  }
  
  infdays_trs_lst
}



count_totals <- function(travellers_tbl) {
  # Count total number of simulated travellers in a table (for use after initialization)
  # Arguments: travellers_tbl (tibble) - just-initialized and pre-filtered table of travellers' infection characteristics and strategies being tested
  # Returns: tibble with columns "StrategyNo" and "total"
  
  tbl <- travellers_tbl %>%
    group_by(StrategyNo) %>%
    summarise(total = n())
  
  tbl
}

### Final processing functions ###

count_vals <- function(travellers_tbl, outcomes, infdays_trs_lst, test_names) {
  # Create table of numbers of diagnosed / recovered / missed travellers, total remaining person-days and estimated transmission risk over days spent outside of quarantine, 
  # and tests used for missed travellers, for a given policy and for each quarantine length
  # Arguments: travellers_tbl (tibble) - table of travellers' infection characteristics and strategies being tested;
  # outcomes (named list) - as output returned by test_outcomes;
  # infdays_trs_lst (named list) - as output returned by infdays_trs;
  # test_names (character vector) - vector of test names
  # Returns: tibble with columns "StrategyNo", "QDays", "Pretest_diagnosed_unrecovered", "Pretest_diagnosed_recovered", "Sympt_filtered_unrecovered", "Sympt_filtered_recovered", 
  # "Quarantine_diagnosed_unrecovered", "Quarantine_diagnosed_recovered", "Missed", "PostQ_InfDays_Sympt", "PostQ_InfDays_NoSympt", "PostQ_TR_Sympt_Logit", "PostQ_TR_NoSympt_Logit", 
  # "TestsUsed_XX", for each type of test XX in the policy table
  
  strategy_no <- travellers_tbl[["StrategyNo"]][1]
  qdaysmin_BA <- travellers_tbl[["QDaysMin_BA"]][1]
  qdaysmax_BA <- travellers_tbl[["QDaysMax_BA"]][1]
  
  summary_tbl <- tibble()
  N <- nrow(travellers_tbl)
  for (day in qdaysmin_BA:qdaysmax_BA) {
    infdays_sympt_name <- paste("PostQ_InfDays_Sympt_", day, sep = "")
    infdays_sympt <- infdays_trs_lst[[infdays_sympt_name]]
    infdays_nosympt_name <- paste("PostQ_InfDays_NoSympt_", day, sep = "")
    infdays_nosympt <- infdays_trs_lst[[infdays_nosympt_name]]
    
    tr_sympt_logit_name <- paste("PostQ_TR_Sympt_Logit_", day, sep = "")
    tr_sympt_logit <- infdays_trs_lst[[tr_sympt_logit_name]]
    tr_nosympt_logit_name <- paste("PostQ_TR_NoSympt_Logit_", day, sep = "")
    tr_nosympt_logit <- infdays_trs_lst[[tr_nosympt_logit_name]]
    
    pretest_diagnosed_colname <- paste("PreTest_diagnosed_BA_", day, sep = "")
    pretest_diagnosed <- outcomes[[pretest_diagnosed_colname]]
    sympt_filtered_colname <- paste("Sympt_filtered_BA_", day, sep = "")
    sympt_filtered <- outcomes[[sympt_filtered_colname]]
    diagnosed_colname <- paste("Diagnosed_BA_", day, sep = "")
    diagnosed <- outcomes[[diagnosed_colname]]
    
    recovered_in_quarantine <- ifelse(infdays_sympt + infdays_nosympt == 0, TRUE, FALSE)
    
    pretest_diagnosed_unrecovered <- mean(as.integer(pretest_diagnosed & !(recovered_in_quarantine)))
    pretest_diagnosed_recovered <- mean(as.integer(pretest_diagnosed & recovered_in_quarantine))
    sympt_filtered_unrecovered <- mean(as.integer(sympt_filtered & !(recovered_in_quarantine)))
    sympt_filtered_recovered <- mean(as.integer(sympt_filtered & recovered_in_quarantine))
    diagnosed_unrecovered <- mean(as.integer(diagnosed & !(recovered_in_quarantine)))
    diagnosed_recovered <- mean(as.integer(diagnosed & recovered_in_quarantine))
    
    missed_colname <- paste("Missed_BA_", day, sep = "")
    missed <- as.integer(outcomes[[missed_colname]])
    
    infdays_sympt <- mean(as.integer(infdays_sympt & missed))
    infdays_nosympt <- mean(as.integer(infdays_nosympt & missed))
    tr_sympt_logit_missed <- sum(tr_sympt_logit * missed)
    tr_nosympt_logit_missed <- sum(tr_nosympt_logit * missed)
    
    missed <- mean(missed)
    
    summary_rows <- tibble(StrategyNo = strategy_no, QDays = day, Pretest_diagnosed_unrecovered = pretest_diagnosed_unrecovered, 
                           Pretest_diagnosed_recovered = pretest_diagnosed_recovered, Sympt_filtered_unrecovered = sympt_filtered_unrecovered, 
                           Sympt_filtered_recovered = sympt_filtered_recovered, Quarantine_diagnosed_unrecovered = diagnosed_unrecovered, 
                           Quarantine_diagnosed_recovered = diagnosed_recovered, Missed = missed, PostQ_InfDays_Sympt = infdays_sympt, 
                           PostQ_InfDays_NoSympt = infdays_nosympt, PostQ_TR_Sympt_Logit = tr_sympt_logit_missed, PostQ_TR_NoSympt_Logit = tr_nosympt_logit_missed)
    
    for (test_name in test_names) {
      tests_used_name <- paste("TestsUsed_", test_name, sep = "")
      tests_used_total <- 0
      for (test_part in c("TestsUsed_PreTest_", "TestsUsed_EntryTest_", "TestsUsed_QTest_", "TestsUsed_ExitTest_")) {
        tests_used_colname <- paste(test_part, "BA_", test_name, "_", day, sep = "")
        tests_used <- outcomes[[tests_used_colname]]
        tests_used_total <- tests_used_total + tests_used
      }
      summary_rows <- summary_rows %>% mutate(!!tests_used_name := tests_used_total / N)
    }
    
    summary_tbl <- rbind(summary_tbl, summary_rows)
  }
  
  summary_tbl
}

calc_tests_uninfected <- function(policy_tbl_BA, test_names) {
  # Calculate numbers of tests used by uninfected travellers, for each policy, quarantine length and test type
  # Arguments: policy_tbl_BA (tibble) - table of quarantine and testing policies followed by all travellers going from Country B to Country A;
  # test_names (character vector) - vector of test names
  # Returns: data.table with columns "StrategyNo", "QDays" and "UninfTestsUsed_XX", for each test type XX in policy table, with StrategyNo and QDays as keys
  
  count_tbl <- tibble()
  strategy_nos <- policy_tbl_BA[["StrategyNo"]]
  for (strategy_no in strategy_nos) {
    qdaysmin_BA <- policy_tbl_BA[.(strategy_no),][["QDaysMin"]]
    qdaysmax_BA <- policy_tbl_BA[.(strategy_no),][["QDaysMax"]]
    
    for (day in qdaysmin_BA:qdaysmax_BA) {
      row_tbl <- tibble(StrategyNo = strategy_no, QDays = day)                                # Prepare single-row table for binding
      for (name in test_names) {
        count_colname <- paste("UninfTestsUsed_", name, sep = "")
        row_tbl[[count_colname]] <- 0
        
        if (!is.na(policy_tbl_BA[.(strategy_no),][["PreTest"]])) {                                             # Increment test count for pretest
          if (policy_tbl_BA[.(strategy_no),][["PreTest"]] == name) {
            row_tbl[[count_colname]] <- row_tbl[[count_colname]] + 1
          }
        }
        
        if (!is.na(policy_tbl_BA[.(strategy_no),][["EntryTest"]])) {
          if (policy_tbl_BA[.(strategy_no),][["EntryTest"]] == name) {                                         # Increment test count for entry test
            row_tbl[[count_colname]] <- row_tbl[[count_colname]] + 1
          }
        }
        
        checkkey_colname <- paste("QTest_check_key_", name, sep = "")
        waitkey_colname <- paste("QTest_wait_key_", name, sep = "")
        increment <- 0
        if (!is.null(policy_tbl_BA[.(strategy_no),][[checkkey_colname]])) {
          qtest_checkkey <- policy_tbl_BA[.(strategy_no),][.(strategy_no),][[checkkey_colname]]
          qtest_waitkey <- as.numeric(unlist(strsplit(policy_tbl_BA[.(strategy_no),][.(strategy_no),][[waitkey_colname]], " ")))
          for (today in 0:day) {
            checkkey_today <- as.numeric(substr(qtest_checkkey, today + 1, today + 1) == "1")
            waitkey_today <- qtest_waitkey[today + 1]
            increment <- increment + ifelse(checkkey_today & waitkey_today <= day & waitkey_today >= 0, 1, 0)      # Increment test count for quarantine tests
          }
        }
        row_tbl[[count_colname]] <- row_tbl[[count_colname]] + increment
        
        if (!is.na(policy_tbl_BA[.(strategy_no),][["ExitTest"]])) {
          if ((policy_tbl_BA[.(strategy_no),][["ExitTest"]] == name) & (day - policy_tbl_BA[.(strategy_no),][["ExitTestDelay"]] >= 1)) {
            row_tbl[[count_colname]] <- row_tbl[[count_colname]] + 1
          }
        }
      }
      count_tbl <- rbind(count_tbl, row_tbl)
    }
  }
  
  count_tbl
}

attach_descriptions <- function(summary_tbl, policy_tbl_BA) {
  # Attach descriptions of strategies to summary table from policy table
  # Arguments: summary_tbl (tibble) - summary_tbl from calc_props or calc_per_million for outcomes;
  # policy_tbl_BA (data.table) - table of quarantine and testing policies followed by all travellers going from Country B to Country A
  # Returns: summary_tbl with added column "StrategyDescription" of descriptions of policies according to StrategyNo
  
  policy_desc_tbl <- policy_tbl_BA %>%
    transmute(StrategyNo, StrategyDescription)
  
  summary_tbl <- policy_desc_tbl %>%
    left_join(summary_tbl, by = "StrategyNo")
  
  summary_tbl
}

### Main loop ###

if (!test_mode) {

sens_viralload_tbl <- read_sens_viralload(filename = sens_viralload_filename)                               # Load required tables
if (policy_option == 1 | policy_option == 2) {
  policy_tbl <- read_policy(option = policy_option, filename = policy_BA_filename, policy_lst = policy_BA)
} else {
  policy_tbl <- read_policy(option = policy_option, filename = policy_reg_BA_filename, policy_lst = policy_reg_BA)
}
viralload_time_params_sympt_tbl <- read_viralload_time_params(filename = viralload_time_params_sympt_filename)
viralload_time_params_asympt_tbl <- read_viralload_time_params(filename = viralload_time_params_asympt_filename)
viralload_time_covar_sympt_mat <- read_viralload_time_covar(filename = viralload_time_covar_sympt_filename)
viralload_time_covar_asympt_mat <- read_viralload_time_covar(filename = viralload_time_covar_asympt_filename)
transrisk_vls <- read_transrisk_vl(filename = transrisk_vl_filename)
test_names <- get_test_types(policy_tbl = policy_tbl)
tests_uninfected_tbl <- calc_tests_uninfected(policy_tbl_BA = policy_tbl, test_names = test_names)

policy_nos <- policy_tbl[["StrategyNo"]]
output_prop_tbl <- tibble()
for (policy_no in policy_nos) {
  travellers_tbl <- init_travellers(policy_no = policy_no, policy_tbl_BA = policy_tbl, viralload_time_params_sympt_tbl = viralload_time_params_sympt_tbl, 
                                    viralload_time_params_asympt_tbl = viralload_time_params_asympt_tbl, viralload_time_covar_sympt_mat = viralload_time_covar_sympt_mat, 
                                    viralload_time_covar_asympt_mat = viralload_time_covar_asympt_mat)
                                                                                                            # Simulate travellers to Country A
  travellers_tbl <- travellers_tbl %>%
    filter_recovered() %>%                                                                                  # Remove recovered travellers
    filter_sympt_infectious() %>%                                                                           # Now remove symptomatic infectious travellers, who wouldn't embark
    draw_test_results(sens_viralload_tbl)
  
  outcomes <- test_outcomes(travellers_tbl, test_names)
  
  infdays_trs_lst <- infdays_trs(travellers_tbl, transrisk_vls)
  
  output_subtbl <- count_vals(travellers_tbl = travellers_tbl, outcomes = outcomes, infdays_trs = infdays_trs_lst, test_names = test_names)
  
  output_prop_tbl <- rbind(output_prop_tbl, output_subtbl)
  
  print(paste(policy_no, " has been processed", sep = ""))
}

output_prop_tbl <- attach_descriptions(output_prop_tbl, policy_tbl)                                         # Attach policy descriptions to output tables
output_prop_tbl <- output_prop_tbl %>% 
  left_join(tests_uninfected_tbl)
if (save_output_prop) {
  write_csv(output_prop_tbl, output_prop_filename)                                                          # Save output table (in proportions)
}

}

### For quick testing ###

if (test_mode & quick_test) {


}