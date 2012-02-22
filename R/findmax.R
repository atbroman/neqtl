
findmax <-
function(neqtl.out, lodcolumn=1, window=5, n.max=10)
{
  neqtl.out$chr <- as.character(neqtl.out$chr)
  for(i in 1:n.max) {
    temp <- max(neqtl.out, lodcolumn=lodcolumn)
    if(length(temp[[3]])==0) break else{
    if(i==1) 
      output <- temp
    else
      output <- rbind(output, temp)
    
    chr <- as.character(temp[[1]])
    pos <- temp[[2]]
    neqtl.out[neqtl.out$chr==chr & neqtl.out$pos >= pos-window &
            neqtl.out$pos <= pos+window,lodcolumn+2] <- NA
  }}
  output
}
      
