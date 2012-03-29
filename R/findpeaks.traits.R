## Find transcripts that map to a hotspot
## chr & pos are vectors of the same length
## chr must be as.character()
findpeaks.traits <- function(chr,pos,sigpos.out,
                         maxlod.out,
                         peaks.out,win=5){
  chr <- as.character(chr)
  results.out <- NULL
  for(i in 1:length(chr)){
        print(i)
        temp <- sigpos.out[[chr[i]]]
        wh <- temp >= pos[i]-(win/2) & temp <= pos[i]+(win/2)
        if(sum(wh)>0){
          results <- data.frame(index=i,chr=rep(chr[i], sum(wh)),
               pos=rep(pos[i], sum(wh)),
               id=names(temp)[wh],peak.pos=temp[wh],
               stringsAsFactors=FALSE)
          if(!missing(maxlod.out))
            results$lod <- maxlod.out$maxlod[chr[i],names(temp)[wh]]
          if(!missing(peaks.out))
            results$cis <- ifelse(is.na(match(names(temp)[wh],
                      names(peaks.out$maxpos.cis[[chr[i]]]))),
                      ifelse(is.na(match(names(temp)[wh],
                      names(peaks.out$maxpos.trans[[chr[i]]]))),NA,0),1)
          if(is.null(results.out)) results.out <- results else
               results.out <- rbind(results.out,results)
        } else next }
  rownames(results.out) <- NULL
  results.out
}   
