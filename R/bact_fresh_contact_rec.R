#' Bacteria Freshwater Contact Recreation Analysis
#'
#' Assesses Ecoli data against the standard
#' @param df dataframe with Ecoli data
#' @param datetime_column POSIXCT column name containing sample datetimes
#' @return a dataframe with relevant Ecoli criteria and excursion variables added
#' @export
#' @examples function(df = your_ecoli_data, datetime_column = "sample_datetime")

Fresh_Contact_rec <- function(df, datetime_column = "sample_datetime"){
  print("Begin fresh contact rec analysis")

  #create lists to get data out of for loops
  geomeanlist = list()

  SampleStartDate <- as.symbol(datetime_column)

  fresh_contact <- df %>%
    dplyr::filter(BacteriaCode == 2,
                  Char_Name == "Escherichia coli") %>%
    #add blank columns to be filled in during analysis phase
    dplyr::mutate(geomean = "",
                  count_period = "")


  if(length(unique(fresh_contact$MLocID)) == 0) {
    print("No E coli Data")
    return(fresh_contact)
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

  pb <- txtProgressBar(0, length(unique(fresh_contact$MLocID)), style = 3)

  for(i in 1:length(unique(fresh_contact$MLocID))){
    utils::setTxtProgressBar(pb, i)
    station <- unique(fresh_contact$MLocID)[i]

    # Filter table down to single station
    fresh_singlestation <- fresh_contact %>%
      dplyr::filter(MLocID == station) %>%
      dplyr::mutate(geomean_start_date = as.Date(!!SampleStartDate)-90)
    #print(paste("i = ", i))

    for(j in 1:nrow(fresh_singlestation)){

      # print(paste("j = ", j))
      #start of 90 day window
      geomean_date <- fresh_singlestation$geomean_start_date[j]
      # end of 90 day window
      enddate <- as.Date(fresh_singlestation[,as.character(SampleStartDate)][j])

      #create table for only samples in that window
      geomean_period <- fresh_singlestation %>%
        dplyr::filter(!!SampleStartDate <= enddate & !!SampleStartDate >= geomean_date )

      count_period = nrow(geomean_period)

      #get geomeans if number of samples in that window is 5 or greater
      fresh_singlestation[j,"geomean"] <- ifelse(nrow(geomean_period) >= 5, geo_mean(geomean_period$Result_cen), NA)
      #get count of 90 day period
      fresh_singlestation[j,"count_period"] <- count_period
      # get number of samples that are above criteria if more than 5 samples in 90 day period
      fresh_singlestation[j,"n_above_crit"] <- ifelse(count_period >= 5, sum(geomean_period$Result_cen > 406), NA)

    }

    geomeanlist[[i]] <- fresh_singlestation

  }

  close(pb)

  #Bind list to dataframe and ensure numbers are numeric
  fresh_analysis <- dplyr::bind_rows(geomeanlist) %>%
    dplyr::mutate(geomean = as.numeric(geomean),
                  count_period = as.numeric(count_period),
                  critical_excursions = odeqassessment::excursions_req(count_period),
                  ss_excursion = if_else(Result_cen > bact_crit_ss, 1, 0),
                  geomean_excursion = if_else(is.na(geomean), 0,
                                              if_else(geomean > bact_crit_geomean, 1, 0)),
                  perc_exceed = if_else(count_period >= 5 & n_above_crit > critical_excursions, 1, 0),
                  excursion_cen = if_else(count_period < 5,
                                          NaN,
                                          if_else(ss_excursion == 1 | geomean_excursion == 1 | perc_exceed == 1,
                                                  1,
                                                  0)
                                          )
    )

  print("Finish fresh contact rec analysis")
  return(fresh_analysis)

}
