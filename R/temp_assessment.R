#' Temperature Assessment
#' 
#' Creates three new variables for each sample summarizing the data within the last 3 years at that station; the number of
#' excursions, the number of samples during the critical period, and the number of samples during the spawning period.
#' @param df dataframe with temperature data
#' @param datetime_column POSIXCT column name containing sample datetimes
#' @param spawn_start mm/dd spawn start date (as characters)
#' @param spawn_end mm/dd spawn end date (as characters)
#' @param results numeric results column name
#' @param criteria name of column containing temperature criteria values
#' @return a dataframe with relevant temperature spawning, criteria, and excursion variables added
#' @export
#' @examples
#' temp_assessment(df = temperature-data, datetime_column = "sample_datetime",
#' spawn_start = "spawn_start", spawn_end = "spawn_end",
#' results = "Results_cen", criteria = "temp_crit")

temp_assessment <- function(df, datetime_column = "sample_datetime", spawn_start_column = "spawn_start", spawn_end_column = "spawn_end",
                            result_column = "Result_cen", criteria_column = "temp_crit"){
  
  sample_datetime <- as.symbol(datetime_column)
  spawn_start <- as.symbol(spawn_start_column)
  spawn_end <- as.symbol(spawn_end_column)
  result <- as.symbol(result_column)
  criteria <- as.symbol(criteria_column)
  
  temp_analysis <- df %>%
    dplyr::mutate(
      # Add columns for Critcal period start and end date
      Crit_period_start = lubridate::mdy(paste0("7/1/",lubridate::year(sample_datetime))),
      Crit_period_end = lubridate::mdy(paste0("9/30/",lubridate::year(sample_datetime))),
      # Append spawn start and end dates with year
      Start_spawn = ifelse(!is.na(spawn_start), paste0(spawn_start,"/",lubridate::year(sample_datetime)), NA ) ,
      End_spawn = ifelse(!is.na(spawn_end), paste0(spawn_end,"/",lubridate::year(sample_datetime)), NA ),
      # Make spwnmn start and end date date format
      Start_spawn = lubridate::mdy(Start_spawn),
      End_spawn = lubridate::mdy(End_spawn),
      # If Spawn dates span a calendar year, account for year change in spawn end date
      End_spawn = dplyr::if_else(End_spawn < Start_spawn & sample_datetime >= End_spawn, End_spawn + lubridate::years(1), # add a year if in spawn period carrying to next year
                                 End_spawn), # otherwise, keep End_spawn as current year
      Start_spawn = dplyr::if_else(End_spawn < Start_spawn & sample_datetime <= End_spawn, Start_spawn - lubridate::years(1), # subtract a year if in spawn period carrying from previous year
                                   Start_spawn), # otherwise, keep Start_spawn as current year
      # Flag for results in critical period
      In_crit_period = ifelse(sample_datetime >= Crit_period_start & sample_datetime <= Crit_period_end, 1, 0 ),
      # Print if result is in spawn or out of spawn
      Spawn_type = ifelse((sample_datetime >= Start_spawn & sample_datetime <= End_spawn & !is.na(Start_spawn)),  "Spawn", "Not_Spawn")
    ) %>%
    dplyr::arrange(sample_datetime)
  
  # Flag if result does not attain standard,  use 13 for during spawn dates, else use criteria
  temp_analysis$excursion = dplyr::if_else((temp_analysis$Spawn_type == "Spawn" & temp_analysis[,result_column] > 13), 1,
                                           dplyr::if_else(temp_analysis$Spawn_type == "Not_Spawn" & temp_analysis[,result_column] > temp_analysis[, criteria_column], 1, 0)
  )
  # Flag for if excursion was in spawn period
  temp_analysis$Spawn_excursion = dplyr::if_else(temp_analysis$Spawn_type == "Spawn" & temp_analysis$excursion == 1, 1, 0 )
  
  # For each observation, determine the number of excursions within the last 3 years, the number of samples
  # within the critical period, and the number of samples in the spawning period, all for the same station.
  cat(paste("Checking within a 3 year period from sample date for the number of previous excursions",
            "and the number of samples taken during spawning and critical periods...\n"))
  temp_analysis <- dplyr::bind_rows(
    pbapply::pblapply(unique(temp_analysis$MLocID), function(x, df){
      cat(paste("\nStation: ", x, "(", which(unique(df$MLocID) == x), "out of", length(unique(df$MLocID)), ")\n"))
      subData <- dplyr::filter(df, MLocID == x)
      subData <- dplyr::bind_cols(
        subData,
        dplyr::bind_rows(pbapply::pblapply(subData$sample_datetime, 
                                           function(x, df){
                                             start3yr <- x - lubridate::years(3)
                                             end3yr <- x
                                             data <- dplyr::filter(df[,c(datetime_column, "excursion", "In_crit_period", "Spawn_type")], 
                                                                   dplyr::between(sample_datetime, start3yr, end3yr))
                                             data <- dplyr::summarise(data, 
                                                                      n_excursions_3yr = sum(data$excursion),
                                                                      samples_in_crit_period = sum(data$In_crit_period),
                                                                      samples_in_spawn_period = sum(data$Spawn_type == "Spawn"))
                                             return(data[, c("n_excursions_3yr", "samples_in_crit_period", "samples_in_spawn_period")])
                                           }, 
                                           df = subData[,c(datetime_column, "excursion", "In_crit_period", "Spawn_type")])
        )
      )
    }, df = temp_analysis)
  )
  
  temp_analysis$excursion_cen <- ifelse(temp_analysis$excursion == 1 & temp_analysis$n_excursions_3yr >= 3, 1, 0)
  temp_analysis$spawning_crit <- 13
  
  return(temp_analysis)
}
