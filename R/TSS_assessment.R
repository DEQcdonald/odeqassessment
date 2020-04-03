#' Assess TSS data
#'
#' Asssesses TSS data against the appropriate criteria and produces columns for excursions
#' @param TSS_data dataframe with TSS data and criteria values
#' @return dataframe with the following columns added; 'excursion_cen'
#' @export
#' @examples
#' TP_assessment(TSS_data = "your-TSS-data")
#'

TSS_assessment <- function(TSS_data) {
  TSS_summary <- TSS_data %>%
    dplyr::mutate(excursion_cen = ifelse(!is.na(TSS_crit),
                                         dplyr::if_else(Result_cen > TSS_crit, 1, 0),
                                         NaN)
    )
  
  # TSS_summary[is.na(TSS_summary$excursion_cen),"excursion_cen"] <- 0
  
  return(TSS_summary)
}
