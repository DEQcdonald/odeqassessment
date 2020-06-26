#' excursions_req
#'
#' This function returns the number of excursions required to list as impaired for conventional pollutants,
#' @param n number of results
#' @export
#' @examples
#' excursions_req()


excursions_req <- function(n){

  x = ifelse(n <= 11, 2, stats::qbinom(0.90, n, 0.10, lower.tail = TRUE)+1 )
  return(x)
}
