#' Assess TP data
#'
#' Asssesses TP data against the appropriate criteria and produces columns for excursions
#' @param TP_data dataframe with TP data and criteria values
#' @return dataframe with the following columns added; 'excursion_cen'
#' @export
#' @examples
#' TP_assessment(TP_data = "your-TP-data")
#'

TP_assessment <- function(TP_data) {
  options(scipen = 999)
  # tmdl_lookup <- read.csv("//deqhq1/wqnps/Status_and_Trend_Reports/Lookups_Statewide/TMDL_targets.csv")
  
  tp_sum_medians <- data.frame()
  if("median" %in% TP_data$stat.base){
    tp_sum_medians <- TP_data %>%
      dplyr::mutate(tp_year = lubridate::year(sample_datetime),
                    tp_month = dplyr::if_else(tp_summer == 1,
                                              lubridate::month("2019-07-01 00:00:00", label = TRUE, abbr = TRUE),
                                              lubridate::month("2019-01-01 00:00:00", label = TRUE, abbr = TRUE)),
                    sample_datetime = dplyr::if_else(tp_summer == 1,
                                                     as.POSIXct(paste0(tp_year, "-07-01 00:00:00")),
                                                     as.POSIXct(paste0(tp_year, "-01-01 00:00:00")))
      ) %>%
      dplyr::filter(stat.base == "median") %>% dplyr::group_by(MLocID, Char_Name, StationDes, Statistical_Base, tp_year, tp_month, tp_summer, sample_datetime) %>%
      dplyr::summarise(tp_median = stats::median(Result_cen, na.rm = TRUE),
                       Result_cen = tp_median,
                       TP_crit = unique(TP_crit),
                       excursion_cen = dplyr::if_else(tp_median > TP_crit, 1, 0)
      ) %>% ungroup()
  }
  
  TP_summary <- TP_data %>% dplyr::filter(!MLocID %in% unique(tp_sum_medians$MLocID)) %>%
    dplyr::mutate(excursion_cen = dplyr::if_else(!is.na(TP_crit),
                                                 dplyr::if_else(Result_cen > TP_crit, 1, 0),
                                                 NaN)
    )
  
  if(nrow(tp_sum_medians) > 0){
    TP_data_combined <- dplyr::bind_rows(TP_summary, tp_sum_medians)
  } else {TP_data_combined <- TP_summary}
  
  # TP_data_combined[is.na(TP_data_combined$excursion_cen),"excursion_cen"] <- 0
  
  return(TP_data_combined)
}
