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

nls2_fit<-
  function(inp,
           cores = 1,
           decay = seq(.01, .11, by = .02),
           delay = seq(0, 10, by = 0.1),
           k = seq(0.1, 1, 0.2),
           bg = 0.2) {
    inp <- inp_order(inp)
    
    # Safe checks natively against rowData instead of mcols
    if(!"delay" %in% names(rowData(inp))){
      rowData(inp)$delay <- as.numeric(NA)
    }
    if(!"half_life" %in% names(rowData(inp))){
      rowData(inp)$half_life <- as.numeric(NA)
    }
    FLT_inp <- inp
    assay(FLT_inp)[decode_FLT(FLT_inp)] <- NA
    row_max <- apply(assay(FLT_inp), 1, max, na.rm = TRUE)
    assay(FLT_inp) <- assay(FLT_inp)/row_max
    tmp_df <- inp_df(FLT_inp, "ID", "position", "flag")
    tmp_df <- tmp_df[!grepl("_TI_", tmp_df$flag), ]
    
    rowData(inp)$delay[rowData(inp)$ID %in% tmp_df$ID] <- NA
    rowData(inp)$half_life[rowData(inp)$ID %in% tmp_df$ID] <- NA
    ids_ABG <- tmp_df$ID[grepl("ABG", tmp_df$flag)]
    time <- metadata(FLT_inp)$timepoints
    
    st_STD <- expand.grid(decay = decay, delay = delay, k = k, bg = bg)
    st_ABG <- expand.grid(decay = decay, delay = delay, k = k)
    upper_STD <- list(decay = log(2)/(1/60), delay = max(time),
                      k = 1/(log(2)/(60)))
    lower_STD <- list(decay = log(2)/(60), delay = 0.001,
                      k = 0.01, bg = 0.2)
    upper_ABG <- list(decay = log(2)/(1/60), delay = max(time),
                      k = 1/(log(2)/(60)))
    lower_ABG <- list(decay = log(2)/(60), delay = 0.001, k = 0.01)
    model_STD <- inty ~ I(time < delay) * I(k / decay + bg) +
      (time >= delay) * I(bg + (k / decay) * (exp(-decay * (time - delay))))
    model_ABG <- inty ~ I(time < delay) * I(k / decay) +
      (time >= delay) * I(k / decay * (exp(-decay * (time - delay))))

    extract_se <- function(fit, param, is_ABG = FALSE) {
      null_result <- list(
        se_decay = NA, t_decay = NA, p_decay = NA,
        se_delay = NA, t_delay = NA, p_delay = NA,
        se_k     = NA, t_k     = NA, p_k     = NA,
        se_bg    = NA, t_bg    = NA, p_bg    = NA
      )
      if (is.null(fit) || any(is.na(fit))) return(null_result)
      tryCatch({
        coef_tab <- summary(fit)$coefficients
        get_val <- function(par, col) {
          if (par %in% rownames(coef_tab)) coef_tab[par, col] else NA
        }
        list(
          se_decay = get_val("decay", "Std. Error"),
          t_decay  = get_val("decay", "t value"),
          p_decay  = get_val("decay", "Pr(>|t|)"),
          se_delay = get_val("delay", "Std. Error"),
          t_delay  = get_val("delay", "t value"),
          p_delay  = get_val("delay", "Pr(>|t|)"),
          se_k     = get_val("k",     "Std. Error"),
          t_k      = get_val("k",     "t value"),
          p_k      = get_val("k",     "Pr(>|t|)"),
          se_bg    = if (is_ABG) NA else get_val("bg", "Std. Error"),
          t_bg     = if (is_ABG) NA else get_val("bg", "t value"),
          p_bg     = if (is_ABG) NA else get_val("bg", "Pr(>|t|)")
        )
      },
      error   = function(e) {
        message("summary() failed for ID: ", param, " — ", conditionMessage(e))
        null_result
      },
      warning = function(w) {
        message("summary() warning for ID: ", param, " — ", conditionMessage(w))
        null_result
      })
    }

    n_fit <- mclapply(seq_len(nrow(tmp_df)), function(i) {
      tmp_Data <- assay(FLT_inp)[rowData(FLT_inp)$ID %in% tmp_df$ID[i],]
      Data_fit <- data.frame(time = time, inty = as.numeric(tmp_Data))
      Data_fit <- na.omit(Data_fit)
      is_ABG <- tmp_df$ID[i] %in% ids_ABG

      if (is_ABG) {
        cc <- capture.output(type = "message",
                             halfLE2 <- tryCatch({
                               nls2(model_ABG, data = Data_fit,
                                    algorithm = "port",
                                    control = list(warnOnly = TRUE),
                                    start = st_ABG, lower = lower_ABG)
                             }, error = function(e) NULL))
      } else {
        cc <- capture.output(type = "message",
                             halfLE2 <- tryCatch({
                               nls2(model_STD, data = Data_fit,
                                    algorithm = "port",
                                    control = list(warnOnly = TRUE),
                                    start = st_STD, lower = lower_STD)
                             }, error = function(e) NULL))
      }

      tryCatch({
        if (is.null(halfLE2)[1] | is.na(halfLE2)[1]) {
          decay_v <- NA; delay_v <- NA; bg_v <- 0; k_v <- NA
        } else {
          decay_v <- coef(halfLE2)[1]
          delay_v <- coef(halfLE2)[2]
          k_v     <- coef(halfLE2)[3]
          bg_v    <- 0
          if (length(coef(halfLE2)) == 4) bg_v <- coef(halfLE2)[4]
        }
      },
      warning = function(war) {
        print(paste("my warning in processing HalfLE2:", i, war))
      },
      error = function(err) {
        print(paste("my error in processing HalfLE2:", i, err))
      })

      se_vals <- extract_se(halfLE2, param = tmp_df$ID[i], is_ABG = is_ABG)

      data_c <- data.frame(
        ID       = tmp_df$ID[i],
        position = tmp_df$position[i],
        delay    = delay_v,
        decay    = decay_v,
        k        = k_v,
        bg       = bg_v,
        se_decay = se_vals$se_decay,
        t_decay  = se_vals$t_decay,
        p_decay  = se_vals$p_decay,
        se_delay = se_vals$se_delay,
        t_delay  = se_vals$t_delay,
        p_delay  = se_vals$p_delay,
        se_k     = se_vals$se_k,
        t_k      = se_vals$t_k,
        p_k      = se_vals$p_k,
        se_bg    = se_vals$se_bg,
        t_bg     = se_vals$t_bg,
        p_bg     = se_vals$p_bg
      )
      return(data_c)
    }, mc.preschedule = FALSE, mc.cores = cores)

    all_cols <- c("ID", "position", "delay", "decay", "k", "bg",
                  "se_decay", "t_decay", "p_decay",
                  "se_delay", "t_delay", "p_delay",
                  "se_k",     "t_k",     "p_k",
                  "se_bg",    "t_bg",    "p_bg")
    numeric_cols <- setdiff(all_cols, c("ID", "position"))

    if (length(n_fit) == 0) {
      fit_nls2 <- as.data.frame(
        matrix(nrow = 0, ncol = length(all_cols)),
        stringsAsFactors = FALSE)
      colnames(fit_nls2) <- all_cols
    } else {
      fit_nls2 <- as.data.frame(
        setNames(
          lapply(all_cols, function(col) {
            vals <- sapply(n_fit, function(row) {
              v <- row[[col]]
              if (is.null(v) || length(v) == 0) return(NA)
              v[[1]]
            })
            if (col %in% numeric_cols) as.numeric(vals) else as.character(vals)
          }),
          all_cols),
        stringsAsFactors = FALSE)
    }

    inp <- inp[order(rowData(inp)$ID), ]
    fit_nls2 <- fit_nls2[order(fit_nls2$ID), ]

    metadata(inp)$fit_STD <- fit_nls2

    # Standardize column access exclusively to rowData
    rowData(inp)$delay[rowData(inp)$ID %in% tmp_df$ID]     <- fit_nls2$delay
    rowData(inp)$delay[!is.finite(rowData(inp)$delay)]         <- NA
    rowData(inp)$half_life[rowData(inp)$ID %in% tmp_df$ID] <- log(2) / fit_nls2$decay
    rowData(inp)$half_life[!is.finite(rowData(inp)$half_life)] <- NA

    se_cols <- c("se_decay", "t_decay", "p_decay",
                 "se_delay", "t_delay", "p_delay",
                 "se_k",     "t_k",     "p_k",
                 "se_bg",    "t_bg",    "p_bg")

    for (col in se_cols) {
      if (!col %in% names(rowData(inp))) {
        rowData(inp)[[col]] <- rep(as.numeric(NA), nrow(inp))
      }
    }

    idx <- rowData(inp)$ID %in% tmp_df$ID
    
    for (col in se_cols) {
      vals <- as.numeric(fit_nls2[[col]])
      vals[!is.finite(vals)] <- NA
      rowData(inp)[[col]][idx] <- vals
    }

    inp <- inp_order(inp)
    inp
  }
