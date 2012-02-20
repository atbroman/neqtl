
findpeaks <-
function(results, lodcolumn=7, window=5, n.peaks=20)
{
  results$chr <- as.character(results$chr)
  for(i in 1:n.peaks) {
    temp <- max(results, lodcolumn=lodcolumn)
    if(i==1) 
      output <- temp
    else
      output <- rbind(output, temp)
    
    chr <- as.character(temp[[1]])
    pos <- temp[[2]]
    results[results$chr==chr & results$pos >= pos-window &
            results$pos <= pos+window,lodcolumn+2] <- NA
  }
  output
}
      
