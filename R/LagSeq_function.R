## TODO
## - not converting codes to numbers
## - handle the case that not all codes appearing in a sequence


#' A utility function for merging same codes appearing together
#' 
#' This function reads a vector of data and then merges same codes that appear together (i.e., without being interrupted by other codes)
#' @param vec A vector of data representing a sequence
#' @return The resulting vector after merging
#' @examples
#' vec = c(1, 2, 2, 3, 4)
#' MergeSameCodes(vec)
MergeSameCodes <- function(vec) {
  
  if(length(vec) == 1) return(vec)
  
  vec <- as.vector(vec)
  v <- c(vec[1])
  for(i in 2:length(vec)) {
    if(vec[i] != v[length(v)]) {
      v <- c(v, vec[i])
    }
  }
  return(v)
}

#' A basic function that converts data to a transitional frequency matrix
#' 
#' This function reads a vector of data and then computes and returns a square output matrix (the cells of which correspond to the cells of the transitional frequency matrix) for each ofthe following: transitional frequencies with row and column totals, ex- pected frequencies, transitional probabilities, adjusted residuals and significance levels, Yule's Q values, transformed kappas (Wampold , 1989, 1992, 1995), z values for the kappas, and significance levels.
#' @param vec A vector of data representing a sequence
#' @param ncodes Integer. The number of codes (or types of events/behaviors) expected in the vector. Optional.
#' @param lag Integer. The lag to be applied to the conversion. Default: 1, denoting converting based on immediate adjacency.
#' @param merge Boolean. Whether to merge the same codes appearing without interruption. Default: FALSE.
#' @return A list of a bunch of matrices representing different statistics of transitional frequency
#' @export
#' @examples
#' vec = c(1, 2, 3, 4, 3, 2, 1, 2)
#' LagSeq(vec, ncodes = 4)
LagSeq <- function(vec, ncodes=0, lag=1, merge=FALSE) {
  
  tmp = as.factor(vec)
  codes_levels = levels(tmp)
  vec <- as.integer(tmp) # TODO
  if(merge==TRUE) vec <- MergeSameCodes(vec)
  
  if(is.null(ncodes) || is.na(ncodes) || ncodes == 0)
    ncodes = length(unique(vec))

  ## Transitional frequency matrix
  freqs = matrix(0, ncodes, ncodes)
  for(c in 1:length(vec)) {
    if(c+lag <= length(vec))
      freqs[vec[c], vec[c+lag]] <- freqs[vec[c], vec[c+lag]] + 1
  }
  
  ## Expected frequency matrix and Adjusted residuals
  rowtots  = rowSums(freqs) # sums of rows
  coltots  = rowSums(freqs) # sums of cols
  ntrans   = sum(rowtots) # total transitions
  prows    = rowtots / ntrans # probability for each row
  pcols    = coltots / ntrans # probability for each column
  expfreq  = matrix(-1.001, ncodes,ncodes) # Expected Values/Frequencies
  zadjres  = matrix(-1.001, ncodes,ncodes) # Adjusted Residuals
  
  for(i in 1:ncodes) {
    for(j in 1:ncodes) {
      if (!merge) {
        expfreq[i,j] = rowtots[i] * coltots[j] / ntrans 
      }
      if (merge && (ntrans - rowtots[j]) > 0 ) {
        expfreq[i,j] = (rowtots[i] * coltots[j]) / (ntrans - rowtots[j])
      }
      
      if ( (expfreq[i,j]*(1-pcols[j])*(1-prows[i])) > 0) {
        zadjres[i,j]=(freqs[i,j]-expfreq[i,j]) / sqrt( expfreq[i,j]*(1-pcols[j]) * (1-prows[i]) )
      }
    }
  }
  
  ## Yule's Q values
  yulesq   = matrix(-1.001, ncodes,ncodes) # Yule's Q Values
  for(i in 1:ncodes) {
    for(j in 1:ncodes) {
      a = freqs[i,j]
      b = rowtots[i] - freqs[i,j]
      c = coltots[j] - freqs[i,j]
      d = ntrans - rowtots[i] - coltots[j] + freqs[i,j]
      if ( (a*d + b*c) > 0 ) {
        yulesq[i,j] = ( a*d - b*c ) / (a*d + b*c)
      }
    }
  }
  
  return(list(freq = freqs, expfreq = expfreq, adjres = zadjres, yulesq = yulesq))
}

#' Compare transitional patterns between two groups
#' 
#' Given two groups, each of which contains multiple sequences of codes (or types of events/behaviors), compare these two groups whether there is any significant difference for each pair of codes for given trasational relationship measure(s).
#' @param df A data frame containing required data. Data should be strict tabular format, with at least the following columns---group membership, sequence membership, codes.
#' @param group Index of the column representing group membership.
#' @param seq Index of the column representing sequence membership.
#' @param codes Index of the column representing codes.
#' @param measure A vector containing measures to compare between two groups. Possible values include `freq`, `adr`, and `yule`.
#' @param ncodes Integer. The number of codes (or types of events/behaviors) expected in the vector. Optional.
#' @param lag Integer. The lag to be applied to the conversion. Default: 1, denoting converting based on immediate adjacency.
#' @param merge Boolean. Whether to merge the same codes appearing without interruption. Default: FALSE.
#' @export
#' @return Nothing. But results of comparisons will be printed.
#' @examples
#' load("lagseq_example_data.Rdata")
#' Lag_Seq_Groups(df, group=6, seq=1, codes=5)
LagSeq_Groups <- function(df, 
                           group, seq, codes,
                           measure="freq", 
                           ncodes=0, lag=1, merge=FALSE) {
  
  options(stringsAsFactors=FALSE)
  
  if(is.null(ncodes) || is.na(ncodes) || ncodes == 0)
    ncodes = length(unique(df[, codes]))
  
  ## convert codes to integers
  tmp = as.factor(df[, codes])
  (codes_levels = levels(tmp))
  df[, codes] = as.integer(tmp)
  ## make sequences unique by pasting seq with group
  df[, seq] <- paste(df[, group], df[, seq], sep="_")
  
  lag_measures <- data.frame(rep(c(), 1 + ncodes ** 2))
  seqs <- levels(factor(df[, seq]))
  groups = rep(NA, length(seqs))
  for(i in seq(1, length(seqs))) {
    s = seqs[i]
    df_sub <- df[df[, seq] == s, ]
    if(nrow(df_sub) <= 1)
      lag_measures <- rbind(lag_measures, c(1, rep(NA, ncodes ** 2)))
    else {
      matrices = LagSeq(df_sub[, codes], ncodes, lag, merge) # ncodes problematic here
      v_m = matrices$freq
      count <- sum(rowSums(v_m))
      if(measure == "adr")
        v_m <- matrices$adjres
      if(measure == "yule")
        v_m <- matrices$yulesq
      
      lag_measures <- rbind(lag_measures, c(count, as.vector(t(v_m))))
    }
    
    groups[i] = unique(df_sub[, group])
  }
  names(lag_measures) <- c("count", sapply(1:ncodes, function(x) paste(1:ncodes, x, sep="_")))
  lag_measures$seq = seqs
  lag_measures$group = groups
  
  # describe first
  require(psych)
  print(describeBy(lag_measures[, 1:(ncol(lag_measures)-2)], lag_measures$group))
  
  # t-tests
  (groups_u = unique(groups))
  for(c in 1:(ncol(lag_measures)-2)) {
    tryCatch({
      t = t.test(lag_measures[lag_measures$group == groups_u[1], c],
                 lag_measures[lag_measures$group == groups_u[2], c])
      if(!is.na(t) && t$p.value < 0.1) {
        cat("t-test for", names(lag_measures)[c], ":\n")
        print(t)
      }
    }, error=function(cond) {
      # message("Here's the original error message:")
      # message(cond)
    })
  }
}
