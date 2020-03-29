#' Censored Data
#' 
#' This function deals with censored data as described in the white paper "Use of Censored Data" from 5/11/2018
#' When the QL > criteria value, 1/2 the value of the criteria is substituted for any sample reported as censored
#' When the QL < Criteria value, 1/2 the QL value will be substituted
#' Samples greater than the Max QL, use QL value
#' Input is a dataframe
#' Adds a column named Result_cen for the value to use
#' Column names should be strings.
#' @param df dataframe to modify. 
#' @param results name of column in dataframe with the result to be compared against. Defaults to 'Result_Numeric'.
#' @param results_qualifier name of column in dataframe for the qualifier. Defaults to 'Result_Operator'.
#' @param criteria name of column in dataframe with applicable criteria value. In the case of parameters with 
#' multiple applicable criteria (i.e. pH) the criteria should be lowest criteria. Defaults to 'criteria'.
#' @return dataframe with additional "Results_cen" column
#' @export
#' 

Censored_data <- function(df, results = "Result_Numeric", results_qualifier = "Result_Operator", criteria = 'criteria') {
  results = as.symbol(results)
  results_qualifier = as.symbol(results_qualifier)
  criteria  = as.symbol(criteria)
  
  
  # Perform censored data modifications
  Results_censored <- df %>%
    # Get lowest criteria value to set censored results
    dplyr::mutate(Result_cen = ifelse(!!results_qualifier == "=", as.numeric(!!results),
                                      ifelse(!!results_qualifier == ">", as.numeric(!!results), 
                                             ifelse(!!results_qualifier == "<", ifelse(!!results > as.numeric(!!criteria), 0.5 * as.numeric(!!criteria) , 0.5 * as.numeric(!!results) ), "ER" )))) %>%
    dplyr::mutate(Result_cen = as.numeric(Result_cen))
  
  return(Results_censored)
  
}
