fragment_function <- function(seg) {
  segs <- c()
  unique_Input <- unique(seg)
  for (j in seq_along(unique_Input)) {
    ind_input <- which(unique_Input[j] == seg)
    if (length(ind_input) == 1) {
      next ()
    } else{
      segs <- c(segs, unique_Input[j])
    }
  }
  return(segs)
}

t_test_function <-
  function(data, seg, param, o, tu, threshold) {
    if (length(seg) > 1) {
      for (j in seq_along(seg) - 1) {
        frag <- paste0(o, "_fragment")
        frg.1.h <- tu[which(tu[, frag] == seg[j]),]
        frg.2.h <- tu[which(tu[, frag] == seg[j + 1]),]
        #if the dataframe has less than 3 probes/bins, will not be subjected to
        # t-test
        if (nrow(frg.1.h) < 2 | nrow(frg.2.h) < 2) {
          next ()
        } else if(o == "HL" & nrow(frg.1.h) >= 3 & nrow(frg.2.h) == 0){
          quot.hl <-1
          frag_hl <- paste0(seg[j], ":", seg[j])
          #assign NA to a p_value
          t_h <- NA
        } else if (abs(last(frg.1.h$position) - frg.2.h$position[1]) <=
                   threshold) {
          frg.1.h <- tu[which(tu[, frag] == seg[j]), param]
          frg.2.h <- tu[which(tu[, frag] == seg[j + 1]), param]
          t_h <-
            t.test(frg.1.h,
                   frg.2.h,
                   alternative = "two.sided",
                   var.equal = FALSE)
          #write fragments name separated by ":" subjected to t-test
          frag_hl <- paste0(seg[j], ":", seg[j + 1])
          #extract the p_value from t-test
          t_h <- t_h[[3]]
          #Fold change of the neighboring fragments, the mean of the
          #fragment2/the mean of the fragment1
          quot.hl <- log2(mean(frg.2.h) / mean(frg.1.h))
          #assign the fragments name to the corresponding column
          Column <- as.data.frame(rowRanges(data)[,
                 frag])
				  Column<-Column[,ncol(Column)]
          rows <- which(Column %in% seg[j])
          if (o == "HL") {
            rowRanges(data)$FC_fragment_HL[rows] <-
              rep(frag_hl, times = length(which(Column %in% seg[j])))
            #assign the ratio of FC to the corresponding column
            rowRanges(data)$FC_HL[
              which(rowRanges(data)$HL_fragment == seg[j])] <-
              rep(quot.hl, times = length(which(
                rowRanges(data)$HL_fragment == seg[j]
              )))
            #assign the p_value to the corresponding column
            rowRanges(data)$p_value_HL[
              which(rowRanges(data)$HL_fragment == seg[j])] <-
              rep(t_h, times = length(which(
                rowRanges(data)$HL_fragment == seg[j]
              )))
          } else if (o == "intensity") {
            rowRanges(data)$FC_fragment_intensity[rows] <-
              rep(frag_hl, times = length(which(Column %in% seg[j])))
            #assign the ratio of FC to the corresponding column
            rowRanges(data)$FC_intensity[
              which(rowRanges(data)$intensity_fragment == seg[j])] <-
              rep(quot.hl, times = length(which(
                rowRanges(data)$intensity_fragment == seg[j]
              )))
            #assign the p_value to the corresponding column
            rowRanges(data)$p_value_intensity[
              which(rowRanges(data)$intensity_fragment == seg[j])] <-
              rep(t_h, times = length(which(
                rowRanges(data)$intensity_fragment == seg[j]
              )))
            }
          } 
        }
      } else {
        next ()
      }
    return(data)
  }

manova_function <- function(x, y, z, data) {
  res.man <- manova(cbind(x, y) ~ z, data = data)
  p_value_Manova <- summary(res.man)[[4]]
  p_value_Manova <- p_value_Manova[1, 6]
  return(p_value_Manova)
}
