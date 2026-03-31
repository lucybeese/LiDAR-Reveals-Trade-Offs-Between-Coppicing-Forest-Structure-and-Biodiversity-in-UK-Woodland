library(Hmisc)

wvioplot = function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
                     horizontal = FALSE, col = "grey", border = "black", lty = 1, 
                     lwd = 1, rectCol = "black", colMed = "white", pchMed = 19, 
                     at, add = FALSE, wex = 1, drawRect = TRUE, centpt = "median", 
                     weights = NULL, adjust = 3, clip = TRUE) 
{
  
  datas <- list(x, ...)
  n <- length(datas)
  if (missing(at)) 
    at <- 1:n
  upper <- vector(mode = "numeric", length = n)
  lower <- vector(mode = "numeric", length = n)
  q1 <- vector(mode = "numeric", length = n)
  q3 <- vector(mode = "numeric", length = n)
  med <- vector(mode = "numeric", length = n)
  base <- vector(mode = "list", length = n)
  height <- vector(mode = "list", length = n)
  baserange <- c(Inf, -Inf)
  
  for (i in 1:n) {
    if(is.list(weights) & !is.null(weights)){
      wgts = weights[[i]]
    }else{
      wgts = weights}
    if(is.null(weights)){
      wgts = rep(1, length(datas[[i]]))
    }		
    data <- datas[[i]]
    wgts = wgts[!is.na(data)]
    data = data[!is.na(data)]
    data.min <- min(data, na.rm=T)
    data.max <- max(data, na.rm=T)
    q1[i] <- wtd.quantile(data, weights = wgts, probs = 0.25, na.rm=T)
    
    q3[i] <- wtd.quantile(data, weights = wgts, probs = 0.75, na.rm=T)
    
    if(centpt=="median"){
      med[i] <- wtd.quantile(data, weights = wgts, probs = .5, na.rm=T)
    }
    if(centpt=="mean"){
      med[i] <- wtd.mean(data, weights = wgts, na.rm=T)
    }
    
    iqd <- q3[i] - q1[i]
    upper[i] <- min(q3[i] + range * iqd, data.max, na.rm=T)
    lower[i] <- max(q1[i] - range * iqd, data.min, na.rm=T)
    est.xlim <- c(min(lower[i], data.min, na.rm=T), max(upper[i], 
                                                        data.max, na.rm=T))
    if(clip){
      args <- list(from = est.xlim[1], to= est.xlim[2])
    }else{args = list()}
    #		print(wgts)
    smout <- do.call("density", c(list(data, 
                                       adjust = adjust, 
                                       weights = (wgts/sum(wgts)),
                                       na.rm=T),args ) )
    hscale <- 0.4/max(smout$y) * wex
    base[[i]] <- smout$x
    height[[i]] <- smout$y * hscale
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1], t[1])
    baserange[2] <- max(baserange[2], t[2])
  }
  if (!add) {
    xlim <- if (n == 1) 
      at + c(-0.5, 0.5)
    else range(at) + min(diff(at))/2 * c(-1, 1)
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  }else {
    label <- names
  }
  boxwidth <- 0.05 * wex
  if (!add) 
    plot.new()
  if (!horizontal) {
    if (!add) {
      plot.window(xlim = xlim, ylim = ylim)
      axis(2)
      axis(1, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
              c(base[[i]], rev(base[[i]])), col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
              lty = lty)
        rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
             q3[i], col = rectCol)
        points(at[i], med[i], pch = pchMed, col = colMed)
      }
    }
  } else {
    if (!add) {
      plot.window(xlim = ylim, ylim = xlim)
      axis(1)
      axis(2, at = at, label = label)
    }
    box()
    for (i in 1:n) {
      polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
                                              rev(at[i] + height[[i]])), col = col, border = border, 
              lty = lty, lwd = lwd)
      if (drawRect) {
        lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
              lty = lty)
        rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
               boxwidth/2, col = rectCol)
        points(med[i], at[i], pch = pchMed, col = colMed)
      }
    }
  }
  invisible(list(upper = upper, lower = lower, median = med, 
                 q1 = q1, q3 = q3))
}