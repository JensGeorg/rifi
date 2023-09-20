#' =========================================================================
#' TI_fit               
#' -------------------------------------------------------------------------
#' TI_fit estimates transcription interference and termination factor using
#' nls function for probe or bin flagged as "TI".
#' 
#' TI_fit uses nls2 function to fit the flagged probes or bins with "TI" found
#' using finding_TI.r.
#' It estimates the transcription interference level (referred later to TI) as
#' well as the transcription factor fitting the probes/bins with nls function
#' looping into several starting values.
#' 
#' 
#' To determine TI and termination factor, TI_fit function is applied to the
#' flagged probes and to the probes localized 1000 nucleotides upstream.
#' Before applying TI_fit function, some probes/bins are filtered out if they
#' are below the background using generic_filter_BG.
#' The model loops into a dataframe containing sequences of starting values and
#' the coefficients are extracted from the fit with the lowest
#' residuals. When many residuals are equal to 0, the lowest residual can not
#' be determined and the coefficients extracted could be wrong.
#' Therefore, a second filter was developed. First we loop into all starting
#' values, we collect nls objects and the corresponding residuals. They are
#' sorted and residuals non equal to 0 are collected in a vector. If the first
#' residuals are not equal to 0, 20 % of the best residuals are collected in
#' tmp_r_min vector and the minimum termination factor is selected. In case the
#' first residuals are equal to 0 then values between 0 to 20% of the values
#' collected in tmp_r_min vector are gathered. The minimum termination factor
#' coefficient is determined and saved. The coefficients are gathered in res
#' vector and saved as an object.
#'
#' @param inp SummarizedExperiment: the input with correct format.
#' @param cores integer: the number of assigned cores for the task.
#' @param restr numeric: a parameter that restricts the freedom of the fit
#' to avoid wrong TI-term_factors, ranges from 0 to 0.2.
#' @param k numeric vector: A sequence of starting values for the synthesis
#' rate. Default is seq(0, 1, by = 0.5).
#' @param decay numeric vector: A sequence of starting values for the decay
#' Default is c(0.05, 0.1, 0.2, 0.5, 0.6).
#' @param ti numeric vector: A sequence of starting values for the delay.
#' Default is  seq(0, 1, by = 0.5).
#' @param ti_delay numeric vector: A sequence of starting values for the
#' delay.
#' Default is seq(0, 2, by = 0.5).
#' @param rest_delay numeric vector: A sequence of starting values. Default
#' is seq(0, 2, by = 0.5).
#' @param bg numeric vector: A sequence of starting values. Default is 0.
#'
#' @return the SummarizedExperiment object: with delay, decay  and
#' TI_termination_factor added to the rowRanges. The full fit data is saved in
#' the metadata as "fit_TI".
#'
#' @examples 
#' data(preprocess_minimal)
#' TI_fit(inp = preprocess_minimal, cores=2, restr=0.01)
#'
#' @export

TI_fit <-
  function(inp,
           cores = 1,
           restr = 0.2,
           k = seq(0, 1, by = 0.5),
           decay = c(0.05, 0.1, 0.2, 0.5, 0.6),
           ti = seq(0, 1, by = 0.5),
           ti_delay = seq(0, 2, by = 0.5),
           rest_delay = seq(0, 2, by = 0.5),
           bg = 0) {
    inp <- inp_order(inp)
    if(!"delay" %in% names(mcols(rowRanges(inp)))){
      rowRanges(inp)$delay <- as.numeric(NA)
    }
    if(!"half_life" %in% names(mcols(rowRanges(inp)))){
      rowRanges(inp)$half_life <- as.numeric(NA)
    }
    if(!"TI_termination_factor" %in% names(mcols(rowRanges(inp)))){
      rowRanges(inp)$TI_termination_factor <- as.numeric(NA)
    }
    FLT_inp <- inp
    assay(FLT_inp)[decode_FLT(FLT_inp)] <- NA
    #normalize
    row_max <- apply(assay(FLT_inp), 1, max, na.rm = TRUE)
    assay(FLT_inp) <- assay(FLT_inp)/row_max
    #make the tmp_df
    tmp_df <- inp_df(FLT_inp, "ID", "position", "flag")
    #only STD
    tmp_df <- tmp_df[grepl("_TI_", tmp_df$flag), ]
    #reset the values
    rowRanges(inp)$delay[rowRanges(inp)$ID %in% tmp_df$ID] <- NA
    rowRanges(inp)$half_life[rowRanges(inp)$ID %in% tmp_df$ID] <- NA
    rowRanges(inp)$TI_termination_factor[rowRanges(inp)$ID %in% tmp_df$ID] <- NA
    #IDs
    ids_ABG <- tmp_df$ID[grepl("ABG",tmp_df$flag)]
    #time points
    time <- metadata(FLT_inp)$timepoints
    #start values
    st_STD <- expand.grid(decay = decay, ti_delay = ti_delay, k = k,
                          rest_delay = rest_delay, ti = ti, bg = bg)
    st_ABG <- expand.grid(decay = decay, ti_delay = ti_delay, k = k,
                          rest_delay = rest_delay, ti = ti)
    
    ###################################################################
    #boarders
    # upper_STD <- list(decay = log(2)/(1/60), ti_delay = max(time),
    #                   k = 1/(log(2)/(60)), rest_delay = max(time),
    #                   ti = 1/(log(2)/(60)))
    lower_STD <- list(decay = log(2)/(60), ti_delay = 0, k = log(2)/(60),
                      rest_delay = 0, ti = 0, bg = 0)
    # upper_ABG <- list(decay = log(2)/(1/60), ti_delay = max(time),
    #                   k = 1/(log(2)/(60)), rest_delay = max(time),
    #                   ti = 1/(log(2)/(60)))
    lower_ABG <- list(decay = log(2)/(60), ti_delay = 0, k = log(2)/(60),
                      rest_delay = 0, ti = 0)
    #models
    model_STD <- inty ~ I(time < ti_delay) * I(k / decay - ti / decay + bg) +
      I(time < ti_delay + rest_delay & time >= ti_delay) *
      I(k / decay - ti / decay * exp(-decay * (time - ti_delay)) + bg) +
      I(time >= ti_delay + rest_delay) *
      I((k / decay - ti / decay * exp(-decay * rest_delay)) *
          exp(-decay * (time - (ti_delay + rest_delay))) + bg)
    model_ABG <- inty ~ I(time < ti_delay) * I(k / decay - ti / decay) +
      I(time < ti_delay + rest_delay & time >= ti_delay) *
      I(k / decay - ti / decay * exp(-decay * (time - ti_delay))) +
      I(time >= ti_delay + rest_delay) *
      I((k / decay - ti / decay * exp(-decay * rest_delay)) *
          exp(-decay * (time - (ti_delay + rest_delay))))

    n_fit <- mclapply(seq_len(nrow(tmp_df)), function(i) {
      #get the Data
      tmp_Data <- assay(FLT_inp)[rowRanges(FLT_inp)$ID %in% tmp_df$ID[i],]
      Data_fit <- data.frame(time = time, inty = as.numeric(tmp_Data))
      Data_fit <- na.omit(Data_fit)
      # probes with flag different from "_" are selected for the model with
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
                                 lower = lower_STD,
                                 #upper = upper_STD,
                                 all = TRUE
                               )},
                               error = function(e) {
                                 return(list(NULL))
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
                                 lower = lower_STD,
                                 # upper = upper_STD,
                                 all = TRUE
                               )},
                               error = function(e) {
                                 return(list(NULL))
                               }
                             ))
      }
     
      tryCatch({
        if (is.null(halfLE2)[1] | is.na(halfLE2)[1]) {
          decay_v <- NA
          ti_delay_v <- NA
          k_v <- NA
          rest_delay_v <- NA
          ti_v <- NA
          bg_v <- 0
        } else {
          #get the one with minimum ti in range of restr
          rss <- lapply(halfLE2, deviance)
          no_rss <- unlist(lapply(rss, is.null))
          halfLE2 <- halfLE2[!no_rss]
          min_rss <- min(unlist(rss)[unlist(rss) != 0])
          in_range <- which(unlist(rss) <= min_rss * (1 + restr))
          halfLE2 <- halfLE2[in_range]
          co <- lapply(halfLE2, function(x) coef(x)[5])
          min_co <- which.min(unlist(co))[1]
          halfLE2 <- halfLE2[[min_co]]
          decay_v <- coef(halfLE2)[1]
          ti_delay_v <- coef(halfLE2)[2]
          k_v <- coef(halfLE2)[3]
          rest_delay_v <- coef(halfLE2)[4]
          ti_v <- coef(halfLE2)[5]
          bg_v <- 0
          if (length(coef(halfLE2)) == 6) {
            bg_v <- coef(halfLE2)[6]
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
      data_c <- data.frame(tmp_df$ID[i], tmp_df$position[i], ti_delay_v,
                           rest_delay_v, decay_v, k_v, ti_v, bg_v)
      colnames(data_c) <-
        c("ID", "position", "ti_delay", "rest_delay", "decay", "k", "ti", "bg")
      return(data_c)
    }, mc.preschedule = FALSE, mc.cores = cores)
    fit_nls2 <- as.data.frame(do.call(rbind, n_fit))
    if (length(n_fit) == 0) {
      fit_nls2 <- data.frame(matrix(nrow = 0, ncol = 8))
      colnames(fit_nls2) <-
        c("ID", "position", "ti_delay", "rest_delay", "decay", "k", "ti", "bg")
    }
    inp <- inp[order(rowRanges(inp)$ID), ]
    fit_nls2 <- fit_nls2[order(fit_nls2$ID), ]
    
    metadata(inp)$fit_TI <- fit_nls2
    
    rowRanges(inp)$delay[rowRanges(inp)$ID %in% tmp_df$ID] <-
      fit_nls2$ti_delay + fit_nls2$rest_delay
    rowData(inp)$delay[!is.finite(rowData(inp)$delay)]<-NA
    rowRanges(inp)$half_life[rowRanges(inp)$ID %in% tmp_df$ID] <-
      log(2) / fit_nls2$decay
    rowData(inp)$half_life[!is.finite(rowData(inp)$half_life)]<-NA
    rowRanges(inp)$TI_termination_factor[rowRanges(inp)$ID %in% tmp_df$ID] <-
      fit_nls2$ti / fit_nls2$k
    rowData(inp)$TI_termination_factor[
      !is.finite(rowData(inp)$TI_termination_factor)]<-NA
    
    inp <- inp_order(inp)
    
    inp
  }
