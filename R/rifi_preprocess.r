#' =========================================================================
#' rifi_preprocess
#' -------------------------------------------------------------------------
#' rifi_preprocess wraps conveniently all pre-processing steps
#'
#' rifi_preprocess wraps the functions:

#' 1. check_input

#' 2. make_df

#' 3. function_seg

#' 4. finding_PDD

#' 5. finding_TI
#' 

#' rifi_preprocess allows for the optional integration of filter functions.

#' Filter functions mark replicates with TRUE. Those are then not considered
#' in the fit!

#' FUN_filter is a general filter usually to exclude probes with low
#' expression or "bad" patterns.
#'
#' @param inp SummarizedExperiment: the input.
#' @param cores integer: the number of assigned cores for the task.
#' @param FUN_filter function: A function of x, returning a logical.
#' x is the numeric vector of the intensity from all time points for a specific
#' replicate.
#' @param bg numeric: threshold over which the last time point has to be to be
#' fitted with the above background mode.
#' @param rm_FLT logical: remove IDs where all replicates are marked as filtered
#' by the background check. Default is FALSE.
#' @param thrsh_check numeric: the minimal allowed intensity for time point "0".
#' Advised to be kept at 0! Default is 0.
#' @param dista integer: the amount of nucleotides defining the gap. Default is
#' 300.
#' @param run_PDD logical: running the PDD flag function
#' @param pen_PDD numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer fragments. Advised to be kept at 2. Default
#' is 2.
#' @param pen_out_PDD numeric: an internal parameter for the dynamic
#' programming. Higher values result in fewer possible outliers. Advised to
#' be kept at 1. Default is 1.
#' @param thrsh_PDD numeric: an internal parameter that allows fragments with
#' slopes steeper than the threshold to be flagged with "_PDD_". Higher values
#' result in fewer candidates . Advised to be kept at 0.001. Default is 0.001.
#' @param pen_TI numeric: an internal parameter for the dynamic programming.
#' Higher values result in fewer fragments. Advised to be kept at 10. Default
#' is 10.
#' @param thrsh_TI numeric: an internal parameter that allows fragments with a
#' certain amount of IDs with higher relative intensities at time points later
#' than "0" to be flagged as "_TI_". Higher values result in fewer candidates.
#' -0.5 is 25 %, 0 is 50%, 0.5 is 75%. Advised to be kept at 0.5. Default
#' is 0.5.
#' @param add integer: range of nucleotides before a potential TI event where
#' in IDs are fitted with the TI fit.
#' 
#' @return The SummarizedExperiment object: checked, and with position, ID,
#' intensity, probe_TI, position_segment, flag and filtration added to the
#' rowRanges.
#'
#' @seealso `check_input`
#' @seealso `make_df`
#' @seealso `segment_pos`
#' @seealso `finding_PDD`
#' @seealso `finding_TI`
#'
#' @examples
#' data(example_input_minimal)
#' rifi_preprocess(
#'   inp = example_input_minimal, cores = 2, bg = 100, rm_FLT = FALSE,
#'   thrsh_check = 0, dista = 300, run_PDD = FALSE
#'   )
#'      
#' @export

rifi_preprocess <-
  function(inp,
           cores,
           FUN_filter = function(x) {
             FALSE
           },
           bg = 0,
           rm_FLT = FALSE,
           thrsh_check = 0,
           dista = 300,
           run_PDD = FALSE,
           pen_PDD = 2,
           pen_out_PDD = 1,
           thrsh_PDD = 0.001,
           pen_TI = 10,
           thrsh_TI = 0.5,
           add = 1000) {

    message("running check_input...")
    inp <- check_input(inp = inp, thrsh = thrsh_check)
    inp_save <- inp
    inp <- tryCatch({
        for(i in seq_along(metadata(inp)$replicate)){
          tmp_inp <- inp[,colData(inp)$replicate == i]
          logi <- apply(assay(tmp_inp), 1, FUN = FUN_filter)
          rows <- which(logi)
          inp <- encode_FLT(obj = inp, rows = rows, rep = i)
        }
        inp
      },
      error = function(e) {
        writeLines(
          paste(
            "An unknown error has appeared!\n",
            e,
            "An emergency output was returned!\n The given filtration
            malfunctioned!"
          )
        )
        inp <- inp_save
        return(inp)
      }
    )

    message("running make_df...")
    inp <- make_df(
      inp = inp,
      cores = cores,
      bg = bg,
      rm_FLT = rm_FLT
    )
    message("running segment_pos...")
    inp <- segment_pos(inp = inp, dista = dista)
    
    if (run_PDD == TRUE) {
      message("running finding_PDD...")
      tryCatch({
          inp <-
            finding_PDD(
              inp = inp,
              pen = pen_PDD,
              pen_out = pen_out_PDD,
              thrsh = thrsh_PDD,
              cores = cores
            )
        },
        error = function(e) {
          writeLines(
            paste(
              "An unknown error has appeared!\n",
              e,
              "An emergency output was returned!\n Please rerun finding_PDD
              manually!"
            )
          )
        }
      )
    }

    message("running finding_TI...")
    tryCatch({
        inp <-
          finding_TI(inp,
            pen = pen_TI,
            thrsh = thrsh_TI,
            cores,
            add = add
          )
      },
      error = function(e) {
        writeLines(
          paste(
            "An unknown error has appeared!\n",
            e,
            "An emergency output was returned!\n Please rerun finding_TI
            manually!"
          )
        )
      }
    )
    metadata(inp)$sessioninfo<-sessionInfo()
    metadata(inp)$date_of_analysis<-date()
    res <- inp
    res
  }
