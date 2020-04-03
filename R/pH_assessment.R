#' Assess pH data
#' 
#' Asssesses pH data against the appropriate criteria and produces columns for excursions, including whether they were high or low
#' @param pH_data dataframe with pH data and criteria values
#' @return dataframe with the following columns added; 'pH_violation', 'pH_violation_high', 'pH_violation_low.'
#' @export
#' @examples
#' pH_assessment(pH_data = "your-pH-data")
#' 

pH_assessment <- function(pH_data) {
  pH_summary <- pH_data %>%
    dplyr::mutate(pH_excursion = dplyr::if_else(Result_cen < pH_Min | Result_cen > pH_Max, 1, 0 ),
                  pH_excursion_high = dplyr::if_else(Result_cen > pH_Max, 1, 0 ),
                  pH_excursion_low = dplyr::if_else(Result_cen < pH_Min, 1, 0 )
    )
  
  pH_summary$excursion_cen <- pH_summary$pH_excursion
  
  return(pH_summary)
}
