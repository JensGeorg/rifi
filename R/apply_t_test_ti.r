#' apply_t_test_ti: compares the mean of two neighboring TI fragments
#' within the same TU.
#' apply_t_test_ti: this function uses the statistical t_test to check
#' if two neighboring TI fragments are significant.
#' @param data dataframe: the probe based data frame.
#' 
#' @return the probe data frame with the columns regarding statistics:
#' \describe{
#'   \item{ID:}{The bin/probe specific ID}
#'   \item{position:}{The bin/probe specific position}
#'   \item{strand:}{The bin/probe specific strand}
#'   \item{flag:}{Information on which fitting model is applied}
#'   \item{position_segment:}{The position based segment}
#'   \item{TI_termination_factor:}{The termination factor of the bin/probe}
#'   \item{TU:}{The overarching transcription unit}
#'   \item{TI_termination_fragment:}{The TI fragment the bin belongs to}
#'   \item{TI_mean_termination_factor:}{The mean termination factor of the
#'   respective TI fragment}
#'   \item{pausing_site:}{}
#'   \item{iTSS_I:}{}
#'   \item{ps_ts_fragment:}{}
#'   \item{event_ps_itss_p_value_Ttest:}{}
#'   \item{p_value_slope:}{}
#'   \item{delay_frg_slope:}{}
#'   \item{velocity_ratio:}{}
#'   \item{event_duration:}{}
#'   \item{event_position:}{}
#'   \item{FC_HL:}{}
#'   \item{FC_fragment_HL:}{}
#'   \item{p_value_HL:}{}
#'   \item{FC_intensity:}{}
#'   \item{FC_fragment_intensity:}{}
#'   \item{p_value_intensity:}{}
#'   \item{FC_HL_intensity:}{}
#'   \item{FC_HL_intensity_fragment:}{}
#'   \item{FC_HL_adapted:}{}
#'   \item{synthesis_ratio:}{}
#'   \item{synthesis_ratio_event:}{}
#'   \item{p_value_Manova:}{}
#'   \item{p_value_TI:}{}
#'   \item{TI_fragments_p_value:}{}
#' }
#' 
#' @examples
#' data(stats_minimal)
#' apply_t_test_ti(data = stats_minimal)
#' 
#' @export
#' 
apply_t_test_ti <- function(data) {
  #new column is added
  data$p_value_TI <- NA
  data$TI_fragments_p_value <- NA
  #excluding outliers
  data_1 <-
    data[!grepl("_T|_O|_NA", data$TI_termination_fragment), ]
  #select unique TUs
  unique_TU <- unique(data$TU)
  #exclude TU terminales, outliers and NAs
  unique_TU <- na.omit(unique_TU[!grepl("_T|_O|_NA", unique_TU)])
  for (i in seq_along(unique_TU)) {
    #grep only TI fragments
    tu <- data_1[which(unique_TU[i] == data_1$TU),
                 c(
                   "ID",
                   "position",
                   "flag",
                   "TI_termination_fragment",
                   "TI_termination_factor",
                   "intensity"
                 )]
    tu <- tu[!is.na(tu$TI_termination_fragment), ]
    ti <- tu[grep("_TI_", tu$flag), ]
    #adjust the fragments for t-test
    if (nrow(ti) == 0) {
      next ()
    } else {
      ti_frag <- unique(tu$TI_termination_fragment)
      if (length(ti_frag) > 1) {
        for (k in seq_len(length(ti_frag) - 1)) {
          seg1 <- tu[which(ti_frag[k] == tu$TI_termination_fragment),
                     "TI_termination_factor"]
          seg2 <-
            tu[which(ti_frag[k + 1] == tu$TI_termination_fragment),
               "TI_termination_factor"]
          if (length(seg1) < 2 | length(seg2) < 2) {
            next ()
          }
          tryCatch({
            #t-test
            ti_test <- t.test(seg1,
                              seg2,
                              alternative = "two.sided",
                              var.equal = FALSE)
            #add 2 columns, fragments column and p_value from t-test
            p_value_tiTest <- ti_test[[3]]
            data[which(ti_frag[k] == data$TI_termination_fragment),
                 "TI_fragments_p_value"] <-
              paste0(ti_frag[k], ":", ti_frag[k + 1])
            data[which(ti_frag[k] == data$TI_termination_fragment),
                 "p_value_TI"] <- p_value_tiTest
          }, error = function(e) {
          })
        }
      }
    }
  }
  return(data)
}