#' DO Analysis
#'
#' Assesses dissolved oxygen data against the relevant standard and calculates DO sat where applicable
#' @param df dataframe with DO data, including temperature and elevation columns
#' @param datetime_column POSIXCT column name containing sample datetimes
#' @param spawn_start_column mm/dd spawn start date (as characters)
#' @param spawn_end_column mm/dd spawn end date (as characters)
#' @param result_column numeric results column name
#' @param temp_column column with numeric temperature values in degrees celsius (NAs for no data)
#' @param elev_column column with numeric elevation values in feet
#' @return a dataframe with relevant Ecoli criteria and excursion variables added
#' @export
#' @examples function(df = your_ecoli_data, datetime_column = "sample_datetime")
#'
DO_assessment <- function(df, datetime_column = "sample_datetime", spawn_start_column = "spawn_start",
                          spawn_end_column = "spawn_end", result_column = "Result_cen",
                          temp_column = "temperature", elev_column = "ELEV_Ft"){

  library(lubridate)
  library(zoo)

  print("Preparing data...")

  # Year round --------------------------------------------------------------

  sample_datetime <- as.symbol(datetime_column)
  spawn_start <- as.symbol(spawn_start_column)
  spawn_end <- as.symbol(spawn_end_column)
  result <- as.symbol(result_column)
  temperature <- as.symbol(temp_column)
  elevation <- as.symbol(elev_column)

  df$DO_Class <- LU_DOCode[match(df$DO_code, LU_DOCode$DO_code), "DO_Class"]
  df <- df %>% mutate(sample_date = as.Date(sample_datetime, format = "%m/%d/%Y"))

  # add spawn start and end dates as dates, include indicator if actdate is within spawn
  # add critical period start and end dates, include indicator is actdate is within critperiod
  data <- df %>% filter(Char_Name == "Dissolved oxygen (DO)") %>%
    mutate(
      # Add columns for Critcal period start and end date
      critstart = mdy(paste0("7/1/",year(sample_datetime) )),
      critend = mdy(paste0("9/30/",year(sample_datetime) )),
      # Append spawn start and end dates with year
      Start_spawn = ifelse(!is.na(spawn_start), paste0(spawn_start,"/",year(sample_datetime)), NA ) ,
      End_spawn = ifelse(!is.na(spawn_end), paste0(spawn_end,"/",year(sample_datetime)), NA ),
      # Make spwnmn start and end date date format
      Start_spawn = mdy(Start_spawn),
      End_spawn = mdy(End_spawn),
      # If Spawn dates span a calendar year, account for year change in spawn end date
      End_spawn = if_else(End_spawn < Start_spawn & sample_datetime >= End_spawn, End_spawn + years(1), # add a year if in spawn period carrying to next year
                          End_spawn), # otherwise, keep End_spawn as current year
      Start_spawn = if_else(End_spawn < Start_spawn & sample_datetime <= End_spawn, Start_spawn - years(1), # subtract a year if in spawn period carrying from previous year
                            Start_spawn), # otherwise, keep Start_spawn as current year
      # Flag for results in spawning and/or critical period
      Spawn_type = if_else(sample_datetime >= Start_spawn & sample_datetime <= End_spawn & !is.na(Start_spawn), "Spawn", "Not_Spawn" ),
      is.crit = if_else(sample_datetime >= critstart & sample_datetime <= critend, 1, 0 ))

  data <- data %>% dplyr::mutate(yr_exc_30DADMean = as.numeric(NaN),
                                 yr_exc_7DADMin = as.numeric(NaN),
                                 yr_exc_inst = as.numeric(NaN),
                                 yr_exc_min = as.numeric(NaN),
                                 spwn_exc_inst = as.numeric(NaN),
                                 spwn_exc_7DADMean = as.numeric(NaN),
                                 spwn_exc_min = as.numeric(NaN),
                                 startdate7 = NA,
                                 startdate30 = NA)

  # data <- data %>%
  #   filter(Statistical_Base %in% c("30DADMean", "7DADMin", '7DADMean', "Minimum", NA)) %>%
  #   mutate(yr_excursion = if_else(is.na(Statistical_Base) & Result_cen < Do_crit_instant, 1,
  #                              if_else(Statistical_Base == "30DADMean" & Result_cen < Do_crit_30D, 1,
  #                                      if_else(Statistical_Base == "7DADMin" & Result_cen < Do_crit_7Mi, 1,
  #                                              if_else(Statistical_Base == "Minimum" & Result_cen < DO_crit_min, 1, 0)))))
  # data <- data %>%
  #   mutate(spawn_excursion = if_else(in_spawn == 1 & Statistical_Base %in% c("7DADMean", "Minimum", NA) & Result_cen < 11, 1, 0))

  # Store 7DADMin data
  DO_7DADMin <- data %>%
    filter(Statistical_Base == '7DADMin') %>%
    mutate(yr_exc_7DADMin = if_else(Result_cen < Do_crit_7Mi, 1, 0))

# Calculate 30DADMean DO sat values ---------------------------------------

  sat_data_30d <- data %>% filter(Statistical_Base %in% c("30DADMean")
                                  # ,yr_excursion == 1, DO_Class == "Cold Water"
                                  )
  if(nrow(sat_data_30d) > 0){
    print("Calculating 30DADMean DO sat values...")
    sat_data_30d <- merge(sat_data_30d, select(filter(df, Char_Name == "Dissolved oxygen saturation", Statistical_Base == "Mean"),
                                               MLocID, sample_date, DO_sat_mean = Result_cen),
                          by = c("MLocID", "sample_date"), all.x = TRUE, all.y = FALSE)
    sat_data_30d <- merge(sat_data_30d, select(filter(df, Char_Name == "Temperature, water", Statistical_Base == "Mean"),
                                               MLocID, sample_date, temp_mean = Result_Numeric),
                          by = c("MLocID", "sample_date"), all.x = TRUE, all.y = FALSE)
    sat_data_30d <- mutate(sat_data_30d,
                           DO_sat = if_else(is.na(DO_sat_mean),
                                            mapply(DO_sat_calc, Result_cen, temp_mean, ELEV_Ft, USE.NAMES = FALSE),
                                            DO_sat_mean),
                           DO_sat = if_else(DO_sat > 100, 100, DO_sat))

    sat_data_30d <- sat_data_30d %>%
      arrange(MLocID, sample_date) %>%
      group_by(MLocID) %>%
      mutate(startdate30 = sample_date - 30)

    sat_data_30d$DO_sat_30DADM <- pbmapply(DO_30DADMean_calc, sat_data_30d$MLocID,
                                           sat_data_30d$startdate30, sat_data_30d$sample_date,
                                           MoreArgs = list(df = sat_data_30d), USE.NAMES = FALSE, SIMPLIFY = TRUE)

    sat_data_30d <- sat_data_30d %>% mutate(yr_exc_30DADMean = if_else(Result_cen < Do_crit_30D,
                                                                  if_else(DO_Class == "Cold Water",
                                                                          if_else(is.na(DO_sat_30DADM) | DO_sat_30DADM < 90, 1, 0),
                                                                          1),
                                                                  0)
                                            )
  }

# Calculate 7DADMean DO sat values ----------------------------------------

  sat_data_7d <- data %>% filter(Statistical_Base %in% c("7DADMean")
                                 # , spawn_excursion == 1
  )

  if(nrow(sat_data_7d) > 0){
    print("Calculating 7DADMean DO sat values...")

    sat_data_7d <- merge(sat_data_7d, select(filter(df, Char_Name == "Dissolved oxygen saturation", Statistical_Base == "Mean"),
                                             MLocID, sample_date, DO_sat_mean = Result_cen),
                         by = c("MLocID", "sample_date"), all.x = TRUE, all.y = FALSE)
    sat_data_7d <- merge(sat_data_7d, select(filter(df, Char_Name == "Temperature, water", Statistical_Base == "Mean"),
                                             MLocID, sample_date, temp_mean = Result_Numeric),
                         by = c("MLocID", "sample_date"), all.x = TRUE, all.y = FALSE)
    sat_data_7d <- mutate(sat_data_7d,
                          DO_sat = if_else(is.na(DO_sat_mean),
                                           pbmapply(DO_sat_calc, Result_cen, temp_mean, ELEV_Ft, USE.NAMES = FALSE),
                                           DO_sat_mean),
                          DO_sat = if_else(DO_sat > 100, 100, DO_sat))
    sat_data_7d <- sat_data_7d %>%
      arrange(MLocID, sample_date) %>%
      group_by(MLocID) %>%
      mutate(startdate7 = lag(sample_date, 6, order_by = sample_date),
             # flag out which result gets a moving average calculated
             calc7ma = ifelse(startdate7 == (sample_date - 6), 1, 0 ),
             DO_sat_7DADMean = ifelse(calc7ma == 1, round(rollmean(x = DO_sat, 7, align = "right", fill = NA),1) , NA ))

    sat_data_7d <- sat_data_7d %>% mutate(spwn_exc_7DADMean = if_else(Result_cen < 11,
                                                                      if_else(is.na(DO_sat_7DADMean) | DO_sat_7DADMean < 95,
                                                                              1, 0),
                                                                      0)
    )
  }

# Calculate daily minimum DO sat values -----------------------------------


  sat_data_min <- data %>% filter(Statistical_Base %in% c("Minimum")
                                  # , spawn_excursion == 1
                                  )
  if(nrow(sat_data_min) > 0){
    print("Calculating minimum DO sat values...")
    sat_data_min <- merge(sat_data_min, select(filter(df, Char_Name == "Dissolved oxygen saturation", Statistical_Base == "Minimum"),
                                               MLocID, sample_date, DO_sat_min = Result_cen),
                          by = c("MLocID", "sample_date"), all.x = TRUE, all.y = FALSE)
    sat_data_min <- merge(sat_data_min, select(filter(df, Char_Name == "Temperature, water", Statistical_Base == "Minimum"),
                                               MLocID, sample_date, temp_min = Result_Numeric),
                          by = c("MLocID", "sample_date"), all.x = TRUE, all.y = FALSE)
    sat_data_min <- mutate(sat_data_min,
                           DO_sat_min = if_else(is.na(DO_sat_min),
                                            pbmapply(DO_sat_calc, Result_cen, temp_min, ELEV_Ft, USE.NAMES = FALSE),
                                            DO_sat_min),
                           DO_sat_min = if_else(DO_sat_min > 100, 100, DO_sat_min))

    sat_data_min <- sat_data_min %>% mutate(yr_exc_min = if_else(Result_cen < DO_crit_min, 1, 0),
                                            spwn_exc_min = if_else(Result_cen < 11,
                                                                   if_else(is.na(DO_sat_min) | DO_sat_min < 95,
                                                                           1, 0),
                                                                   0)
                                            )
  }

# Calculate instantaneous DO sat values -----------------------------------

  sat_data_inst <- data %>% filter(is.na(Statistical_Base)
                                      # , yr_excursion == 1
                                      )
  if(nrow(sat_data_inst) > 0){
    print("Calculating instantaneous DO sat values...")
    sat_data_inst <- merge(sat_data_inst, select(filter(df, Char_Name == "Dissolved oxygen saturation", is.na(Statistical_Base)),
                                                 MLocID, sample_date, DO_sat = Result_cen),
                           by = c("MLocID", "sample_date"), all.x = TRUE, all.y = FALSE)
    sat_data_inst <- merge(sat_data_inst, select(filter(df, Char_Name == "Temperature, water", is.na(Statistical_Base)),
                                                 MLocID, sample_date, temperature = Result_Numeric),
                           by = c("MLocID", "sample_date"), all.x = TRUE, all.y = FALSE)
    sat_data_inst <- mutate(sat_data_inst,
                            DO_sat = if_else(is.na(DO_sat),
                                             pbmapply(DO_sat_calc, Result_cen, temperature, ELEV_Ft, USE.NAMES = FALSE),
                                             DO_sat),
                            DO_sat = if_else(DO_sat > 100, 100, DO_sat))

    sat_data_inst <- sat_data_inst %>% mutate(yr_exc_inst = if_else(Result_cen < Do_crit_instant, 1, 0),
                                              spwn_exc_inst = if_else(Result_cen < 11,
                                                                     if_else(is.na(DO_sat) | DO_sat < 95,
                                                                             1, 0),
                                                                     0)
    )
  }

  print("Evaluating excursions...")
  data <- bind_rows(sat_data_30d, sat_data_7d, sat_data_inst, sat_data_min, DO_7DADMin)

  data <- mutate(data,
                 yr_excursion = if_else(1 %in% c(yr_exc_30DADMean, yr_exc_7DADMin, yr_exc_inst, yr_exc_min), 1, 0),
                 spawn_excursion = if_else((Spawn_type == "Spawn") & (1 %in% c(spwn_exc_inst, spwn_exc_7DADMean, spwn_exc_min)), 1, 0),
                 excursion_cen = if_else((Spawn_type == "Spawn") & (spawn_excursion == 1),
                                         1,
                                         if_else((Spawn_type == "Spawn") & (spawn_excursion == 0),
                                                 0,
                                                 if_else(yr_excursion == 1,
                                                         1,
                                                         0
                                                 )
                                         )
                 )
  ) %>%
    select(-startdate30, -startdate7) %>%
    ungroup()

  # data <- data %>% mutate(excursion_cen = if_else(Statistical_Base == "30DADMean" & yr_excursion == 1 & DO_Class == "Cold Water",
  #                                                 if_else(is.na(DO_sat_30DADM) | DO_sat_30DADM < 90, 1, 0),
  #                                                 if_else(Statistical_Base %in% c("7DADMin", "Minimum", NA) & yr_excursion == 1,
  #                                                         1,
  #                                                         if_else(Statistical_Base %in% c("7DADMean", "Minimum", NA) & spawn_excursion == 1,
  #                                                                 if_else(is.na()))))
  # data$excursion_cen <- if_else(data$yr_excursion == 1 | data$spawn_excursion == 1, 1, 0)


  print("Analysis complete, returning data...")
  return(data)
}


#' DO saturation percent calculator
#'
#' This function will calculate DO saturation percentage
#' based on DO mg/L values, temperature in C, and elevcation in ft
#' This function is based on the equation found in The Dissolved
#' Ocygen Water Quality Standard Implementatiion Guidence.
#' This function differs from the oxySol function in the wql package
#' because it calcultaes the percentage dirrectly and incorporates elevation,
#' as opposed to pressure
#'
#' @param DO DO value in mg/L
#' @param TempC Temperature value in C
#' @param elevft Elevation value in feet
#' @export
#' @examples
#' DO_sat_calc()
#'

DO_sat_calc <- function(DO, TempC, elevft) {
  DO / (exp(-139.34411 + ((1.575701*10^5)/(TempC+273.15)) -
              ((6.642308 * 10^7)/((TempC+273.15)^2)) +
              ((1.243800 * 10^10)/((TempC+273.15)^3)) -
              ((8.621949 * 10^11)/((TempC+273.15)^4))) *
          (1 - (0.0001148 * elevft/3.281 ))) * 100
}

#' 30 Day Average Daily Mean DO saturation percent calculator
#'
#' This function will calculate 30DADMean DO saturation percentage
#' based on daily mean DO sat values. It takes a station, start date,
#' end date (totaling 30 days), and a dataframe with mean DO sat values.
#'
#' @param station The station to calculate
#' @param start_date The first day in the 30 day period
#' @param end The the last day in the 30 day period
#' @param df The dataframe containing DO sat values for subsetting
#' @param n The number of days required in the 30 day period
#' @export
#' @examples
#' DO_30DADMean_calc(station, start_date, end_date, df, n = 29)
#'

DO_30DADMean_calc <- function(station, start_date, end_date, df, n = 29){
  data <- df %>%
    filter(between(sample_date, start_date, end_date), MLocID == station) %>%
    group_by(sample_date) %>%
    summarise(count = if_else(any(!is.na(DO_sat)), 1, 0),
              DO_sat_mean = mean(DO_sat, na.rm = TRUE)
    )

  DO_sat_30DADMean <- if_else(sum(data$count, na.rm = TRUE) >= n, mean(data$DO_sat_mean, na.rm = TRUE), NA_real_)
  return(DO_sat_30DADMean)
}
