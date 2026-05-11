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
    
    # Use rowData instead of mcols(rowRanges(...))
    if(!"delay" %in% names(rowData(inp))){
      rowData(inp)$delay <- as.numeric(NA)
    }
    if(!"half_life" %in% names(rowData(inp))){
      rowData(inp)$half_life <- as.numeric(NA)
    }
    if(!"TI_termination_factor" %in% names(rowData(inp))){
      rowData(inp)$TI_termination_factor <- as.numeric(NA)
    }
    FLT_inp <- inp
    assay(FLT_inp)[decode_FLT(FLT_inp)] <- NA
    row_max <- apply(assay(FLT_inp), 1, max, na.rm = TRUE)
    assay(FLT_inp) <- assay(FLT_inp)/row_max
    tmp_df <- inp_df(FLT_inp, "ID", "position", "flag")
    tmp_df <- tmp_df[grepl("_TI_", tmp_df$flag), ]
    
    # Safely assign to rowData
    rowData(inp)$delay[rowData(inp)$ID %in% tmp_df$ID] <- NA
    rowData(inp)$half_life[rowData(inp)$ID %in% tmp_df$ID] <- NA
    rowData(inp)$TI_termination_factor[rowData(inp)$ID %in% tmp_df$ID] <- NA
    ids_ABG  <- tmp_df$ID[grepl("ABG", tmp_df$flag)]
    time     <- metadata(FLT_inp)$timepoints
    st_STD   <- expand.grid(decay = decay, ti_delay = ti_delay, k = k,
                            rest_delay = rest_delay, ti = ti, bg = bg)
    st_ABG   <- expand.grid(decay = decay, ti_delay = ti_delay, k = k,
                            rest_delay = rest_delay, ti = ti)
    lower_STD <- list(decay = log(2)/(60), ti_delay = 0, k = log(2)/(60),
                      rest_delay = 0, ti = 0, bg = 0)
    lower_ABG <- list(decay = log(2)/(60), ti_delay = 0, k = log(2)/(60),
                      rest_delay = 0, ti = 0)
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

    # helper: extract SE and t/p-values from a single TI fit object safely
    extract_se_TI <- function(fit, id, is_ABG = FALSE) {
      null_result <- list(
        se_decay      = NA, t_decay      = NA, p_decay      = NA,
        se_ti_delay   = NA, t_ti_delay   = NA, p_ti_delay   = NA,
        se_k          = NA, t_k          = NA, p_k          = NA,
        se_rest_delay = NA, t_rest_delay = NA, p_rest_delay = NA,
        se_ti         = NA, t_ti         = NA, p_ti         = NA,
        se_bg         = NA, t_bg         = NA, p_bg         = NA
      )
      if (is.null(fit) || any(is.na(fit))) return(null_result)
      tryCatch({
        coef_tab <- summary(fit)$coefficients
        get_val <- function(par, col) {
          if (par %in% rownames(coef_tab)) coef_tab[par, col] else NA
        }
        list(
          se_decay      = get_val("decay",      "Std. Error"),
          t_decay       = get_val("decay",      "t value"),
          p_decay       = get_val("decay",      "Pr(>|t|)"),
          se_ti_delay   = get_val("ti_delay",   "Std. Error"),
          t_ti_delay    = get_val("ti_delay",   "t value"),
          p_ti_delay    = get_val("ti_delay",   "Pr(>|t|)"),
          se_k          = get_val("k",          "Std. Error"),
          t_k           = get_val("k",          "t value"),
          p_k           = get_val("k",          "Pr(>|t|)"),
          se_rest_delay = get_val("rest_delay", "Std. Error"),
          t_rest_delay  = get_val("rest_delay", "t value"),
          p_rest_delay  = get_val("rest_delay", "Pr(>|t|)"),
          se_ti         = get_val("ti",         "Std. Error"),
          t_ti          = get_val("ti",         "t value"),
          p_ti          = get_val("ti",         "Pr(>|t|)"),
          se_bg         = if (is_ABG) NA else get_val("bg", "Std. Error"),
          t_bg          = if (is_ABG) NA else get_val("bg", "t value"),
          p_bg          = if (is_ABG) NA else get_val("bg", "Pr(>|t|)")
        )
      },
      error = function(e) {
        message("summary() failed for ID: ", id, " — ", conditionMessage(e))
        null_result
      },
      warning = function(w) {
        message("summary() warning for ID: ", id, " — ", conditionMessage(w))
        null_result
      })
    }

    n_fit <- mclapply(seq_len(nrow(tmp_df)), function(i) {
      tmp_Data <- assay(FLT_inp)[rowData(FLT_inp)$ID %in% tmp_df$ID[i],]
      Data_fit <- data.frame(time = time, inty = as.numeric(tmp_Data))
      Data_fit <- na.omit(Data_fit)
      is_ABG   <- tmp_df$ID[i] %in% ids_ABG

      if (is_ABG) {
        cc <- capture.output(type = "message",
                             halfLE2 <- tryCatch({
                               nls2(model_ABG, data = Data_fit,
                                    algorithm = "port",
                                    control = list(warnOnly = TRUE),
                                    start = st_ABG, lower = lower_ABG,
                                    all = TRUE)
                             }, error = function(e) return(list(NULL))))
      } else {
        cc <- capture.output(type = "message",
                             halfLE2 <- tryCatch({
                               nls2(model_STD, data = Data_fit,
                                    algorithm = "port",
                                    control = list(warnOnly = TRUE),
                                    start = st_STD, lower = lower_STD,
                                    all = TRUE)
                             }, error = function(e) return(list(NULL))))
      }

      best_fit <- NULL
      tryCatch({
        if (is.null(halfLE2)[1] | is.na(halfLE2)[1]) {
          decay_v <- NA; ti_delay_v <- NA; k_v <- NA
          rest_delay_v <- NA; ti_v <- NA; bg_v <- 0
        } else {
          rss      <- lapply(halfLE2, deviance)
          no_rss   <- unlist(lapply(rss, is.null))
          halfLE2  <- halfLE2[!no_rss]
          min_rss  <- min(unlist(rss)[unlist(rss) != 0])
          in_range <- which(unlist(rss) <= min_rss * (1 + restr))
          halfLE2  <- halfLE2[in_range]
          co       <- lapply(halfLE2, function(x) coef(x)[5])
          min_co   <- which.min(unlist(co))[1]
          best_fit <- halfLE2[[min_co]]     
          decay_v      <- coef(best_fit)[1]
          ti_delay_v   <- coef(best_fit)[2]
          k_v          <- coef(best_fit)[3]
          rest_delay_v <- coef(best_fit)[4]
          ti_v         <- coef(best_fit)[5]
          bg_v         <- 0
          if (length(coef(best_fit)) == 6) bg_v <- coef(best_fit)[6]
        }
      },
      warning = function(war) {
        print(paste("my warning in processing HalfLE2:", i, war))
      },
      error = function(err) {
        print(paste("my error in processing HalfLE2:", i, err))
      })

      se_vals <- extract_se_TI(best_fit, id = tmp_df$ID[i], is_ABG = is_ABG)

      data_c <- data.frame(
        ID           = tmp_df$ID[i],
        position     = tmp_df$position[i],
        ti_delay     = ti_delay_v,
        rest_delay   = rest_delay_v,
        decay        = decay_v,
        k            = k_v,
        ti           = ti_v,
        bg           = bg_v,
        se_decay     = se_vals$se_decay,
        t_decay      = se_vals$t_decay,
        p_decay      = se_vals$p_decay,
        se_ti_delay  = se_vals$se_ti_delay,
        t_ti_delay   = se_vals$t_ti_delay,
        p_ti_delay   = se_vals$p_ti_delay,
        se_k         = se_vals$se_k,
        t_k          = se_vals$t_k,
        p_k          = se_vals$p_k,
        se_rest_delay = se_vals$se_rest_delay,
        t_rest_delay  = se_vals$t_rest_delay,
        p_rest_delay  = se_vals$p_rest_delay,
        se_ti        = se_vals$se_ti,
        t_ti         = se_vals$t_ti,
        p_ti         = se_vals$p_ti,
        se_bg        = se_vals$se_bg,
        t_bg         = se_vals$t_bg,
        p_bg         = se_vals$p_bg
      )
      return(data_c)
    }, mc.preschedule = FALSE, mc.cores = cores)

    all_cols <- c(
      "ID", "position", "ti_delay", "rest_delay", "decay", "k", "ti", "bg",
      "se_decay",      "t_decay",      "p_decay",
      "se_ti_delay",   "t_ti_delay",   "p_ti_delay",
      "se_k",          "t_k",          "p_k",
      "se_rest_delay", "t_rest_delay", "p_rest_delay",
      "se_ti",         "t_ti",         "p_ti",
      "se_bg",         "t_bg",         "p_bg")
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

    # Order correctly referencing rowData
    inp     <- inp[order(rowData(inp)$ID), ]
    fit_nls2 <- fit_nls2[order(fit_nls2$ID), ]

    metadata(inp)$fit_TI <- fit_nls2

    # ALL metadata assignment strictly defaults to rowData
    rowData(inp)$delay[rowData(inp)$ID %in% tmp_df$ID] <-
      fit_nls2$ti_delay + fit_nls2$rest_delay
    rowData(inp)$delay[!is.finite(rowData(inp)$delay)] <- NA
    
    rowData(inp)$half_life[rowData(inp)$ID %in% tmp_df$ID] <-
      log(2) / fit_nls2$decay
    rowData(inp)$half_life[!is.finite(rowData(inp)$half_life)] <- NA
    
    rowData(inp)$TI_termination_factor[rowData(inp)$ID %in% tmp_df$ID] <-
      fit_nls2$ti / fit_nls2$k
    rowData(inp)$TI_termination_factor[
      !is.finite(rowData(inp)$TI_termination_factor)] <- NA

    # Initialise SE/t/p columns safely in rowData DataFrame directly
    se_cols <- c("se_decay",      "t_decay",      "p_decay",
                 "se_ti_delay",   "t_ti_delay",   "p_ti_delay",
                 "se_k",          "t_k",          "p_k",
                 "se_rest_delay", "t_rest_delay", "p_rest_delay",
                 "se_ti",         "t_ti",         "p_ti",
                 "se_bg",         "t_bg",         "p_bg")
                 
    for (col in se_cols) {
      if (!col %in% names(rowData(inp))) {
        rowData(inp)[[col]] <- rep(as.numeric(NA), nrow(inp))
      }
    }

    # Assign SE/t/p directly into the rowData
    idx <- rowData(inp)$ID %in% tmp_df$ID
    for (col in se_cols) {
      vals <- as.numeric(fit_nls2[[col]])
      vals[!is.finite(vals)] <- NA
      rowData(inp)[[col]][idx] <- vals
    }

    inp <- inp_order(inp)
    inp
  }
