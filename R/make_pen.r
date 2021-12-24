#' make_pen: automatically assigns a penalties.
#' 'make_pen' calls one of four available penalty functions to automatically
#' assign penalties for the dynamic programming.
#' The four functions to be called are:
#' 1. fragment_delay_pen
#' 2. fragment_HL_pen
#' 3. fragment_inty_pen
#' 4. fragment_TI_pen.
#' These functions return the amount of statistically correct and statistically
#' wrong splits at a specific pair of penalties.
#' 'make_pen' iterates over many penalty pairs and picks the most suitable pair
#' based on the difference between wrong and correct splits.
#' The sample size, penalty range and resolution as well as the number of cycles
#' can be customized.
#' The primary start parameters create a matrix with n = rez_pen rows and
#' n = rez_pen_out columns with values between sta_pen/sta_pen_out
#' and end_pen/end_pen_out. The best penalty pair is picked. If dept is bigger
#' than 1 the same process is repeated with a new matrix of the same size
#' based on the result of the previous cycle. Only position segments with length
#' within the sample size range are considered for the penalties to increase
#' run time.
#' Returns a penalty object (list of 4 objects) the first being the logbook.
#'
#' @param probe data frame: the probe based data frame.
#' @param FUN function: one of the four bottom level functions (see details)
#' @param cores integer: the number of assigned cores for the task
#' @param logs numeric vector: the logbook vector.
#' @param dpt integer: the number of times a full iteration cycle is repeated
#' with a more narrow range based on the previous cycle. Default is 2.
#' @param smpl_min integer: the smaller end of the sampling size. Default is 10.
#' @param smpl_max integer: the larger end of the sampling size. Default is 100.
#' @param sta_pen numeric: the lower starting penalty. Default is 0.5.
#' @param end_pen numeric: the higher starting penalty. Default is 4.5.
#' @param rez_pen numeric: the number of penalties iterated within the penalty
#' range. Default is 9.
#' @param sta_pen_out numeric: the lower starting outlier penalty.
#' Default is 0.5.
#' @param end_pen_out numeric: the higher starting outlier penalty.
#' Default is 3.5.
#' @param rez_pen_out numeric: the number of outlier penalties iterated within
#' the outlier penalty range. Default is 7.
#'
#' @return A list with 4 items:
#'     \describe{
#'       \item{logbook:}{The logbook vector containing all penalty information}
#'       \item{penalties:}{a vetor with the respective penalty and outlier
#'       penalty}
#'       \item{correct:}{a matrix of the correct splits}
#'       \item{wrong:}{a matrix of the incorrect splits}
#'     }
#'     
#' @examples
#' data(fit_minimal)
#' make_pen(
#'   probe = fit_minimal, FUN = rifi:::fragment_HL_pen, cores = 2,
#'   logs = as.numeric(rep(NA, 8)), dpt = 1, smpl_min = 10, smpl_max = 50,
#'   sta_pen = 0.5, end_pen = 4.5, rez_pen = 9, sta_pen_out = 0.5,
#'   end_pen_out = 3.5, rez_pen_out = 7
#' )
#' 
#' @export

make_pen <- function(probe,
                     FUN,
                     cores = 1,
                     logs,
                     dpt = 1,
                     smpl_min = 10,
                     smpl_max = 100,
                     sta_pen = 0.5,
                     end_pen = 4.5,
                     rez_pen = 9,
                     sta_pen_out = 0.5,
                     end_pen_out = 3.5,
                     rez_pen_out = 7) {
  num_args <-
    list(
      cores,
      dpt,
      smpl_min,
      smpl_max,
      sta_pen,
      end_pen,
      rez_pen,
      sta_pen_out,
      end_pen_out,
      rez_pen_out
    )
  names(num_args) <-
    c(
      "cores",
      "dpt",
      "smpl_min",
      "smpl_max",
      "sta_pen",
      "end_pen",
      "rez_pen",
      "sta_pen_out",
      "end_pen_out",
      "rez_pen_out"
    )
  assert(all(unlist(lapply(
    num_args,
    FUN = function(x) {
      (is.numeric(x) &
         length(x) == 1)
    }
  ))),
  paste0("'", names(which(
    unlist(lapply(
      num_args,
      FUN = function(x) {
        (is.numeric(x) &
           length(x) == 1)
      }
    )) == FALSE
  ))[1], "' must be numeric of length one"))
  assert(cores > 0, "'cores' must be a positive integer")
  assert(
    is.numeric(logs) &
      is.vector(logs) &
      length(logs) == 8,
    "'logs' must be a numeric vector of length 8"
  )
  req_cols_probe <-
    c("ID", "position", "strand", "position_segment")
  assert(
    all(req_cols_probe %in% colnames(probe)),
    paste0("'", req_cols_probe[which(!req_cols_probe %in% colnames(probe))],
           "' must be a column in 'probe'!")
  )
  assert(
    is.function(FUN),
    "'FUN' must be a function from fragment_delay_pen, fragment_HL_pen,
    fragment_inty_pen, or fragment_TI_pen"
  )
  res2 <- vector("list", dpt)
  res3 <- vector("list", dpt)
  step_pen <- (end_pen - sta_pen) / (rez_pen - 1)
  step_pen_out <- (end_pen_out - sta_pen_out) / (rez_pen_out - 1)
  pen <- seq(sta_pen, end_pen, step_pen)
  pen_out <- seq(sta_pen_out, end_pen_out, step_pen_out)
  for (i in seq_len(dpt)) {
    correct <- matrix(, length(pen), length(pen_out))
    wrong <- matrix(, length(pen), length(pen_out))
    rownames(correct) <- pen
    colnames(correct) <- pen_out
    rownames(wrong) <- pen
    colnames(wrong) <- pen_out
    for (j in seq_along(pen)) {
      tmp_pen_out <- pen_out[pen_out >= 0.4 * pen[j]]
      for (k in seq_along(tmp_pen_out)) {
        tmp <-
          FUN(probe, pen[j], tmp_pen_out[k], smpl_min, smpl_max, cores = cores)
        correct[j, as.character(tmp_pen_out[k])] <- tmp[[1]]
        wrong[j, as.character(tmp_pen_out[k])] <- tmp[[2]]
      }
    }
    res2[[i]] <- correct
    res3[[i]] <- wrong
    dif <- correct - wrong
    if (any(dif >= wrong, na.rm = TRUE)) {
      ind <-
        which(dif == max(dif[dif >= wrong], na.rm = TRUE), arr.ind = TRUE)[1, ]
    } else {
      ind <- which(dif == max(dif, na.rm = TRUE), arr.ind = TRUE)[1, ]
      warning("no set of penalties suits the criteria, continued with next best
              option")
    }
    sta_pen <- as.numeric(rownames(dif)[ind[1]]) - step_pen
    end_pen <- as.numeric(rownames(dif)[ind[1]]) + step_pen
    sta_pen_out <- as.numeric(colnames(dif)[ind[2]]) - step_pen_out
    end_pen_out <- as.numeric(colnames(dif)[ind[2]]) + step_pen_out
    step_pen <- (end_pen - sta_pen) / (rez_pen - 1)
    step_pen_out <- (end_pen_out - sta_pen_out) / (rez_pen_out - 1)
    pen <- seq(sta_pen, end_pen, step_pen)
    pen_out <- seq(sta_pen_out, end_pen_out, step_pen_out)
    res1 <-
      c(as.numeric(rownames(dif)[ind[1]]), as.numeric(colnames(dif)[ind[2]]))
  }
  logs[c(
    paste0(names(tmp)[1], "_penalty"),
    paste0(names(tmp)[1], "_outlier_penalty")
  )] <- c(res1[1], res1[2])
  names(res1) <- names(tmp)
  res <- list(logs, res1, res2, res3)
  res
}