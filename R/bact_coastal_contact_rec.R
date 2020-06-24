#' Bacteria Coastal Contact Recreation Analysis
#'
#' Assesses Ecoli data against the standard
#' @param df dataframe with Ecoli data
#' @param datetime_column POSIXCT column name containing sample datetimes
#' @return a dataframe with relevant Ecoli criteria and excursion variables added
#' @export
#' @examples function(df = your_ecoli_data, datetime_column = "sample_datetime")

Coastal_Contact_rec <- function(df, datetime_column = "sample_datetime"){

  print("Begin coastal contact rec analysis")

  #create lists to get data out of for loops
  geomeanlist = list()

  SampleStartDate <- as.symbol(datetime_column)

  # Water Contact Recreation - Coastal -----------------------------------

  # Get filter down to data needed only for coastal contact rec data
  # Bacteria code #2 and Entero
  Coastal <- df %>%
    dplyr::filter(BacteriaCode %in%  c(1, 3),
                  Char_Name == "Enterococcus") %>%
    #add blank columns to be filled in during analysis phase
    dplyr::mutate(geomean = as.numeric(NaN),
                  count_period = as.numeric(NaN),
                  n_above_crit = as.numeric(NaN),
                  perc_above_crit_10 = as.numeric(NaN),
                  n_samples_greater_perc_crit = as.numeric(NaN),
                  less_5 = as.numeric(NaN),
                  Max_value = as.numeric(NaN),
                  SS_Crit = NaN,
                  Geomean_Crit = 35,
                  Perc_Crit = 130)


  if(length(unique(Coastal$MLocID)) == 0) {
    # stop("No Enterococcus Data")
    return(Coastal %>% dplyr::mutate(excursion_cen = NA))
  }

  # Geometric mean calculations --------------------------------------------


  # Process the geometirc means
  # These for loops first filter data down to individual monitoring stations
  # and sets a variable for each sampling date that indicates the start of a 90 day geomean window.
  # The second for loop loops through each activity date and creates a table of all activity dates in that
  # 90 day window and calculates the geomettric mean. It then assigns the geomeans into the single location table
  # created in the first loop, if there are more than 5 sampling dates in that window.
  # The end of the first loop puts the single location table into a list which is used to bring
  # the data out of the for loop by binding it together after the loop into table "Coastal_singlestation"

  print("Begin analysis")

  pb <- txtProgressBar(0, length(unique(Coastal$MLocID)), style = 3)

  for(i in 1:length(unique(Coastal$MLocID))){
    setTxtProgressBar(pb, i)
    station <- unique(Coastal$MLocID)[i]

    # Filter table down to single station
    Coastal_singlestation <- Coastal %>%
      dplyr::filter(MLocID == station) %>%
      dplyr::mutate(geomean_start_date = as.Date(!!SampleStartDate)-90)

    for(j in 1:nrow(Coastal_singlestation)){

      #start of 90 day window
      geomean_date <- Coastal_singlestation$geomean_start_date[j]
      # end of 90 day window
      enddate <- Coastal_singlestation[[SampleStartDate]][j]

      #create table for only samples in that window
      geomean_period <- Coastal_singlestation %>%
        dplyr::filter(!!SampleStartDate <= enddate & !!SampleStartDate >= geomean_date )

      count_period = nrow(geomean_period)

      #get geomeans if number of samples in that window is 5 or greater
      Coastal_singlestation[j,"geomean"] <- ifelse(nrow(geomean_period) >= 5, geo_mean(geomean_period$Result_cen), NA)
      #get count of 90 day period
      Coastal_singlestation[j,"count_period"] <- count_period
      # get number that are above 130 criterion
      Coastal_singlestation[j,"n_above_crit"] <- sum(geomean_period$Result_cen > 130)
      # get percent that are above criteria if more than 10 samples in 90 day period
      Coastal_singlestation[j,"perc_above_crit_10"] <- ifelse(count_period >= 10, sum(geomean_period$Result_cen > 130) /count_period, NA)
      # get lowest value in 90 day window if 5-9 samples in 90 day window
      Coastal_singlestation[j,"n_samples_greater_perc_crit"]  <- ifelse(count_period < 10 & count_period >= 5, sum(geomean_period$Result_cen > 130), NA )
      # flag if less than 5 in 90 day window
      Coastal_singlestation[j,"less_5"] <- ifelse(nrow(geomean_period) < 5, 1, 0)
      #Max Value
      Coastal_singlestation[j,"Max_value"] <- max(geomean_period$Result_cen)


    }

    geomeanlist[[i]] <- Coastal_singlestation

  }

  close(pb)

  #Bind list to dataframe and ensure numbers are numeric
  Coastal_analysis <- dplyr::bind_rows(geomeanlist) %>%
    dplyr::mutate(geomean = as.numeric(geomean),
                  count_period = as.numeric(count_period),
                  n_above_crit = as.numeric(n_above_crit),
                  perc_above_crit_10 = as.numeric(perc_above_crit_10),
                  n_samples_greater_perc_crit = as.numeric(n_samples_greater_perc_crit ),
                  less_5 = as.numeric(less_5),
                  geomean_excursion = ifelse(geomean > 35, 1, 0),
                  perc_excursion = ifelse(perc_above_crit_10 > 0.1, 1, 0),
                  ss_excursion = ifelse(count_period >= 5 & count_period <= 9 & n_above_crit > 1, 1, 0),
                  excursion_cen = ifelse(less_5 == 1,
                                         2,
                                         ifelse(geomean_excursion == 1 | perc_excursion == 1 | ss_excursion == 1,
                                                1,
                                                0))
    )

  print("Finish coastal contact rec analysis")
  return(Coastal_analysis)

}
