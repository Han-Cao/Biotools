library(plotrix)

extract_output <- function(outtab) {
  dim(outtab) <- c(1, length(outtab))
  ncol = dim(outtab)[2]
  nsnp = dim(outtab)[1]
  mat <- as.matrix(outtab[,17:ncol]) # p-m data
  if (nsnp == 1) {
    dim(mat) <- c(1, length(mat))
  }
  ncol <- ncol - 16
  p <- as.numeric(mat[,1:(ncol/2)])
  h <- as.numeric(mat[,(ncol/2+1):ncol])
  FE_P <- as.numeric(outtab[, 3])
  FE_P <- formatC(FE_P,format = "E", digits = 2)
  RE2_P <- as.numeric(outtab[, 9])
  RE2_P <- formatC(RE2_P,format = "E", digits = 1)
  FE_beta <- as.numeric(outtab[,4]) # fixed effects model summary
  FE_se <- as.numeric(outtab[5])
  RE_beta <- as.numeric(outtab[,7]) # random effects model summary
  RE_se <- as.numeric(outtab[,8])
  I_square <- as.numeric(outtab[,13])
  I_square <- paste0(round(I_square,1), "%")
  Q_pval <- as.numeric(outtab[,15])
  Q_pval <- formatC(Q_pval, format = "E", digits = 1)
  return(list(p=p, h=h, RE2_P=RE2_P, FE_P=FE_P,
              FE_beta=FE_beta, FE_se=FE_se,
              RE_beta=RE_beta, RE_se=RE_se,
              I_square=I_square, Q_pval=Q_pval))
}

extract_input <- function(intab) {
  dim(intab) <- c(1, length(intab))
  # read data for forest plot
  #data <- intab[sel, ]
  ncol = dim(intab)[2]
  nsnp = dim(intab)[1]
  mat <- as.matrix(intab[,2:ncol])
  if (nsnp) {
    dim(mat) <- c(1, length(mat))
    effects <- mat[,c(TRUE,FALSE)]
    stderrs <- mat[,c(FALSE,TRUE)]
    dim(effects) <- c(1, length(effects))
    dim(stderrs) <- c(1, length(stderrs))
  }
  else {
    effects <- mat[,c(TRUE,FALSE)]
    stderrs <- mat[,c(FALSE,TRUE)]
  }
  return(list(effects=effects, stderrs=stderrs))
}

meta.colors<-function(all.elements,box="black",lines="gray",summary="black",zero="lightgray",
                      mirror="lightblue",text="black", axes="black",background="white"){
  
  if (missing(all.elements)){
    return(list(box=box, lines=lines, summary=summary,
                zero=zero, mirror=mirror, text=text,
                axes=axes, background=background))
  }
  
  if (is.null(all.elements))
    all.elements<-par("fg")
  
  return(list(box=all.elements, lines=all.elements,
              summary=all.elements, zero=all.elements,
              mirror=all.elements, text=all.elements,
              axes=all.elements, background=NA))
  
}



metaplot <- function( mn, se, nn=NULL, labels=NULL, conf.level = .95,
                      xlab = "Odds ratio", ylab = NA, study_group = NULL,
                      xlim = NULL, summn = NULL, pooled = FALSE,
                      sumse = NULL, sumnn = NULL,
                      summlabel = "Summary", heterolabel = NULL,
                      logeffect = FALSE, lwd = 2, boxsize = 1,
                      zero = as.numeric(logeffect),
                      colors=meta.colors(), xaxt="s", logticks=TRUE,
                      ... ) {
  nth<-function(x,i){
    x[ (i-1) %% length(x) +1]
  }
  ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
  ok <- is.finite( mn + se )
  if ( is.null( xlim ) )
    xlim <- c( min( mn[ok] - ci.value * se[ok], na.rm = TRUE )-0.00,
               max( mn[ok] + ci.value * se[ok], na.rm = TRUE ) )
  ##par( pty="s" )
  n <- length( mn )
  if ( logeffect ) {
    xlog <- "x"
    nxlim <- exp( xlim )
  }
  else {
    xlog <- ""
    nxlim <- xlim
  }
  
  leftedge<-nxlim[1]
  lalim = nxlim
  
  if ( !is.null( labels ) ) {
    if ( logeffect )
      nxlim[1] <- nxlim[1] / sqrt( nxlim[2] / nxlim[1] )
    else
      #nxlim[1] <- nxlim[1] - 0.5 * ( nxlim[2] - nxlim[1] )
      #lalim[1] <- nxlim[1] - 0.6 * ( nxlim[2] - nxlim[1] )
      lalim[1] <- nxlim[1] - 0.6 * ( nxlim[2] - nxlim[1] )
    nxlim[1] <- nxlim[1] - 0.7 * ( nxlim[2] - nxlim[1] )
    
    labels<-as.character(labels)
    
  }
  par( xaxt = "n",yaxt = "n", bg=colors$background )
  plot( nxlim,c( 1,-n-2-3 * !is.null( summn ) ),
        type = "n", bty = "n", xaxt = "n", yaxt = "n",
        log = xlog, xlab=xlab, ylab=ylab,..., col.lab=colors$axes )
  
  par( xaxt = "s" )
  if (xaxt=="s"){
    if (logeffect) {
      if (logticks){
        ats<-round( 10 ^ pretty( log( exp( xlim ),10), 8,min.n=6  ), 2 )
        ats<-ats[ats> exp(xlim[1]) & ats< 10^(par("usr")[2])]
        axis( 1, at = ats, col= colors$axes, col.axis= colors$axes)
      } else {
        ats<-pretty(exp(xlim),8, min.n=6)
        ats<-ats[ats> exp(xlim[1]) & ats <10^(par("usr")[2])]
        axis( 1, at=ats, col= colors$axes, col.axis= colors$axes)
      }
    }  else {
      ats<-pretty(xlim, 6)
      ##ats<-ats[ats> xlim[1] & ats <xlim[2]]
      axis( 1, at=ats, col= colors$axes, col.axis= colors$axes)
    }
  }
  
  if ( !is.null( zero )&& zero>leftedge )
    abline( v = zero, lty = 2, lwd = 2 ,col=colors$zero)
  
  ci.value <- -qnorm( ( 1 - conf.level ) / 2 )
  lower <- mn - ci.value * se
  upper <- mn + ci.value * se
  if ( logeffect ){
    lower <- exp( lower )
    upper <- exp( upper )
  }
  for ( i in 1:n ){
    if ( is.na( lower[i]+upper[i] ) )
      next
    lines( c( lower[i], upper[i] ), c( -i, -i ), lwd = lwd, col=nth(colors$lines,i),... )
  }
  
  if ( !is.null( labels ) )
    #text( rep( nxlim[1], n ), -( 1:n ), labels,..., col=rep(colors$text,length.out=n),adj=0 )
    #text( rep( lalim[1], n ), -( 1:n ), labels,..., col=rep(colors$text,length.out=n),adj=-0.4 )
    #text( rep( lalim[1]+0.4, n ), -( 1:n ), labels,..., col=rep(colors$text,length.out=n),adj=0 )
    ## added by yurang.park
    studies = labels
  if(max(nchar(studies))>5){
    mxstrlenlabel <- max(strwidth(labels))-0.4    		    
  }else{
    mxstrlenlabel <- 0         
  }        
  #text( rep( lalim[1]-mxstrlenlabel, n ), -( 1:n ), labels,..., col=rep(colors$text,length.out=n),adj=0)
  #Edit by Hanc            
  text( rep( lalim[1], n ), -( 1:n ), labels,..., col=rep(colors$text,length.out=n),adj=0)
  
  
  if ( is.null( nn ) )
    nn <- se ^ -2
  yscale <- 0.3 * boxsize / max( sqrt( nn ), na.rm = TRUE )
  
  if ( logeffect ) {
    scale <- ( nxlim[2] / nxlim[1] ) ^ ( yscale / ( 4 + n ) )
    xl <- exp( mn ) * ( scale ^ -sqrt( nn ) )
    xr <- exp( mn ) * ( scale ^ sqrt( nn ) )
  }
  else {
    scale <- yscale * ( nxlim[2] - nxlim[1] ) / ( 4 + n )
    xl <- mn - scale * sqrt( nn )
    xr <- mn + scale * sqrt( nn )
  }
  yb <- ( 1:n ) - yscale * sqrt( nn )
  yt <- ( 1:n ) + yscale * sqrt( nn )
  for ( i in 1:n ) {
    if ( !is.finite( mn[i] ) )
      next
    rect( xl[i], -yb[i], xr[i], -yt[i], col = nth(colors$box,i),border=nth(colors$box,i))
  }
  if ( !is.null( summn ) ) {
    # meta summary
    for (i in 1:length(summn)) {
      if ( logeffect ) {
        x0 <- exp( summn[i] )
        xl <- exp( summn[i] - ci.value * sumse[i] )
        xr <- exp( summn[i] + ci.value * sumse[i] )
      }
      else{
        x0 <- summn[i]
        xl <- summn[i] - ci.value * sumse[i]
        xr <- summn[i] + ci.value * sumse[i]
      }
      if(pooled){
        n <- n + 1
        text( lalim[1], -n, study_group,..., col=rep(colors$text,length.out=n),adj=0)
      }
      else
        text( lalim[1], 0, study_group,..., col=rep(colors$text,length.out=n),adj=0)
      
      
      y0 <- n + i + 1
      yb <- y0 - sqrt( sumnn[i] ) * yscale
      yt <- y0 + sqrt( sumnn[i] ) * yscale
      polygon( c( xl, x0, xr, x0 ), -c( y0, yt, y0, yb ),
               col = colors$summary, border = colors$summary )
      text( lalim[1], -y0, labels = summlabel[[i]], adj = 0,col=colors$text )
    }
    # hetero label
    text( lalim[1], -n - 1, labels = heterolabel, adj = 0,col=colors$text )
  }
}

forestplot <- function(intab, outtab, labels, xlim, study_group, RE = FALSE, pooled= FALSE, P_list=NULL, beta_list=NULL, se_list=NULL){
  beta <- as.numeric(extract_input(intab)$effects)
  stderr <- as.numeric(extract_input(intab)$stderrs)
  FE_beta <- as.numeric(extract_output(outtab)$FE_beta)
  FE_se <- as.numeric(extract_output(outtab)$FE_se)
  RE_beta <- as.numeric(extract_output(outtab)$RE_beta)
  RE_se <- as.numeric(extract_output(outtab)$RE_se)
  FE_P <- extract_output(outtab)$FE_P
  RE2_p <- extract_output(outtab)$RE2_P
  I_square <- extract_output(outtab)$I_square
  Q_pval <- extract_output(outtab)$Q_pval
  if (is.null(P_list)){
    if (RE) {
      summary_beta = RE_beta
      summary_se = RE_se
      summary_P = RE2_P
    }
    else{
      summary_beta = FE_beta
      summary_se = FE_se
      summary_P = FE_P
    }
  }
  else {
    summary_beta = beta_list
    summary_se = se_list
    summary_P = P_list
  }

  summlabel <- lapply(summary_P, function(x) paste0("P = ", x))
  heterolabel <- paste0("Hetero: I2 = ", I_square, " P = ", Q_pval)
  metaplot(beta, stderr, labels = labels, xlim = xlim, 
           summn = summary_beta, sumse = summary_se, 
           sumnn=1/(as.numeric(summary_se)^2),
           summlabel = summlabel, heterolabel = heterolabel,
           study_group = study_group, pooled = pooled)
}
