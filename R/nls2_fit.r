#' =========================================================================
#' nls2_fit 
#' -------------------------------------------------------------------------
#'nls2_fit estimates decay for each probe or bin
#'
#' nls2_fit uses nls2 function to fit a probe or bin using intensities of the
#' time series data from different time point. nls2 uses different starting 
#' values through expand grid and selects the best fit. Different filters could 
#' be applied prior fitting to the model.
#' 
#' To apply nls2_fit function, prior filtration could applied.

#' 1. generic_filter_BG: filter probes with intensities below background using
#' threshold. Those probes are filtered.

#' 2. filtration_below_backg: additional functions exclusive to microarrays
#' could be applied. Its very strict to the background (not recommended in
#' usual case).

#' 3. filtration_above_backg: selects probes with a very high intensity and
#' above the background (recommended for special transcripts). Probes are
#' flagged with "_ABG_".

#' Those transcripts are usually related to a specific function in bacteria.
#' This filter selects all probes with the same ID, the mean is applied, 
#' the last time point is selected and compared to the threshold.
#'
#' The model used estimates the delay, decay, intensity of the first time
#' point (synthesis rate/decay) and the background.

#' The coefficients are gathered in vectors with the corresponding IDs.
#' Absence of the fit or a very bad fit are assigned with NA.

#' In case of probes with very high intensities and above the background,
#' the model used makes abstinence of background coefficient.

#' The output of all coefficients is saved in the metadata.

#' The fits are plotted using the function_plot_fit.r through rifi_fit.
#'
#' @param inp SummarizedExperiment: the input with correct format.
#' @param cores integer: the number of assigned cores for the task.
#' @param decay numeric vector: A sequence of starting values for the decay.
#' Default is seq(.08, 0.11, by=.02)
#' @param delay numeric vector: A sequence of starting values for the delay.
#' Default is seq(0,10, by=.1)
#' @param k numeric vector: A sequence of starting values for the synthesis
#' rate. Default is seq(0.1,1,0.2)
#' @param bg numeric vector: A sequence of starting values. Default is 0.2.
#'
#' @return the SummarizedExperiment object: with delay and decay added to the
#' rowRanges. The full fit data is saved in the metadata as "fit_STD".
#'   \describe{
#'   \item{delay:}{Integer, the delay value of the bin/probe}
#'   \item{half_life:}{Integer, the half-life of the bin/probe}
#'   }
#' 
#' @examples
#' data(preprocess_minimal)
#' nls2_fit(inp = preprocess_minimal, cores = 2)
#' 
#' @export

nls2_fit <-
  function(inp,
           cores = 1,
           decay = seq(.01, .11, by = .02),
           delay = seq(0, 10, by = 0.1),
           k = seq(0.1, 1, 0.2),
           bg = 0.2) {
    #order the input
    inp <- inp_order(inp)
    if(!"delay" %in% names(mcols(rowRanges(inp)))){
      rowRanges(inp)$delay <- as.numeric(NA)
    }
    if(!"half_life" %in% names(mcols(rowRanges(inp)))){
      rowRanges(inp)$half_life <- as.numeric(NA)
    }
    FLT_inp <- inp
    assay(FLT_inp)[decode_FLT(FLT_inp)] <- NA
    #normalize
    row_max <- apply(assay(FLT_inp), 1, max, na.rm = TRUE)
    assay(FLT_inp) <- assay(FLT_inp)/row_max
    #make the tmp_df
    tmp_df <- inp_df(FLT_inp, "ID", "position", "flag")
    #only STD
    tmp_df <- tmp_df[!grepl("_TI_", tmp_df$flag), ]
    #reset the values
    rowRanges(inp)$delay[rowRanges(inp)$ID %in% tmp_df$ID] <- NA
    rowRanges(inp)$half_life[rowRanges(inp)$ID %in% tmp_df$ID] <- NA
    #IDs
    ids_ABG<-tmp_df$ID[grepl("ABG",tmp_df$flag)]
    #time points
    time <- metadata(FLT_inp)$timepoints
    #start values
    st_STD <- expand.grid(decay = decay, delay = delay, k = k, bg = bg)
    st_ABG <- expand.grid(decay = decay, delay = delay, k = k)
    #boarders
    upper_STD <- list(decay = log(2)/(1/60), delay = max(time),
                       k = 1/(log(2)/(60)))
    lower_STD <- lower_STD <- list(decay = log(2)/(60), delay = 0.001, 
                                   k = 0.01, bg = 0.2)
    upper_ABG <- list(decay = log(2)/(1/60), delay = max(time),
                       k = 1/(log(2)/(60)))
    lower_ABG <- list(decay = log(2)/(60), delay = 0.001, k= 0.01)
    #models
    model_STD <- inty ~ I(time < delay) * I(k / decay + bg) + 
      (time >= delay) * I(bg + (k / decay) * (exp(-decay * (time - delay))))

    model_ABG <- inty ~ I(time < delay) * I(k / decay) + 
      (time >= delay) * I(k / decay * (exp(-decay * (time - delay))))
    
    n_fit <- mclapply(seq_len(nrow(tmp_df)), function(i) {
      #get the Data
      tmp_Data <- assay(FLT_inp)[rowRanges(FLT_inp)$ID %in% tmp_df$ID[i],]
      Data_fit <- data.frame(time = time, inty = as.numeric(tmp_Data))
      Data_fit <- na.omit(Data_fit)
      # probes with flag different from "_ABG_" are selected for the model with
      # background coefficient,
      # otherwise the model without background coefficient is applied.
      if (tmp_df$ID[i] %in% ids_ABG) {
        cc <- capture.output(type="message",
                             halfLE2 <- tryCatch({
                               halfLE2 <- nls2(
                                 model_ABG,
                                 data = Data_fit,
                                 algorithm = "port",
                                 control = list(warnOnly = TRUE),
                                 start = st_ABG,
                                 lower = lower_ABG
                                 #upper = upper_ABG
                               )},
                               error = function(e) {
                                 return(NULL)
                               }
                             ))
      } else {
        cc <- capture.output(type="message",
                             halfLE2 <- tryCatch({
                               halfLE2 <- nls2(
                                 model_STD,
                                 data = Data_fit,
                                 algorithm = "port",
                                 control = list(warnOnly = TRUE),
                                 start = st_STD,
                                 lower = lower_STD
                                 #upper = upper_STD
                               )},
                               error = function(e) {
                                 return(NULL)
                               }
                             ))
      }
      tryCatch({
        if (is.null(halfLE2)[1] | is.na(halfLE2)[1]) {
          decay_v <- NA
          delay_v <- NA
          bg_v <- 0
          k_v <- NA
        } else {
          decay_v <- coef(halfLE2)[1]
          delay_v <- coef(halfLE2)[2]
          k_v <- coef(halfLE2)[3]
          bg_v <- 0
          if (length(coef(halfLE2)) == 4) {
            bg_v <- coef(halfLE2)[4]
          }
        }
      },
      warning = function(war) {
        print(paste("my warning in processing HalfLE2:", i, war))
      },
      error = function(err) {
        print(paste("my error in processing HalfLE2:", i, err))
      }
      )
      data_c <- data.frame(tmp_df$ID[i], tmp_df$position[i], delay_v, decay_v,
                           k_v, bg_v)
      colnames(data_c) <-
        c("ID", "position", "delay", "decay", "k", "bg")
      return(data_c)
    }, mc.preschedule = FALSE, mc.cores = cores)
    fit_nls2 <- as.data.frame(do.call(rbind, n_fit))
    if (length(n_fit) == 0) {
      fit_nls2 <- data.frame(matrix(nrow = 0, ncol = 6))
      colnames(fit_nls2) <-
        c("ID", "position", "delay", "decay", "k", "bg")
    }
    inp <- inp[order(rowRanges(inp)$ID), ]
    fit_nls2 <- fit_nls2[order(fit_nls2$ID), ]
    
    metadata(inp)$fit_STD <- fit_nls2
    
    rowRanges(inp)$delay[rowRanges(inp)$ID %in% tmp_df$ID] <- fit_nls2$delay
    rowData(inp)$delay[!is.finite(rowData(inp)$delay)]<-NA
    rowRanges(inp)$half_life[rowRanges(inp)$ID %in% tmp_df$ID] <-
      log(2) / fit_nls2$decay
    rowData(inp)$half_life[!is.finite(rowData(inp)$half_life)]<-NA
    
    inp <- inp_order(inp)
    
    inp
  }
