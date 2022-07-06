# Calculate denominators
denom <- function(x, y, popfrac) {
  
  cont.val <- ifelse(class(x) == "dcm", "nruns", "nsims")
  if (popfrac == TRUE) {
    for (i in 1:length(y)) {
      dname <- paste(strsplit(y[i], "[.]")[[1]][-1], collapse = ".")
      x$epi[[y[i]]] <- x$epi[[y[i]]] / x$epi[[dname]]
    }
  }
  if (popfrac == FALSE && x$control[[cont.val]] == 1) {
    for (j in 1:length(y)) {
      x$epi[[y[j]]] <- data.frame(x$epi[[y[j]]])
    }
  }
  
  return(x)
}

addmediansci_model1<-function(mod){
  mod$epi$s.num$median<-apply(mod$epi$s.num, 1, median)
  mod$epi$s.num$low<-apply(mod$epi$s.num, 1, quantile, probs=0.25)
  mod$epi$s.num$high<-apply(mod$epi$s.num, 1, quantile, probs=0.75)
  mod$epi$s.num$min<-apply(mod$epi$s.num, 1, min)
  mod$epi$s.num$max<-apply(mod$epi$s.num, 1, max)
  
  mod$epi$e.num$median<-apply(mod$epi$e.num, 1, median)
  mod$epi$e.num$low<-apply(mod$epi$e.num, 1, quantile, probs=0.25)
  mod$epi$e.num$high<-apply(mod$epi$e.num, 1, quantile, probs=0.75)
  mod$epi$e.num$min<-apply(mod$epi$e.num, 1, min)
  mod$epi$e.num$max<-apply(mod$epi$e.num, 1, max)
  
  mod$epi$i.num$median<-apply(mod$epi$i.num, 1, median)
  mod$epi$i.num$low<-apply(mod$epi$i.num, 1, quantile, probs=0.25)
  mod$epi$i.num$high<-apply(mod$epi$i.num, 1, quantile, probs=0.75)
  
  mod$epi$rplus.num$median<-apply(mod$epi$rplus.num, 1, median)
  mod$epi$rplus.num$low<-apply(mod$epi$rplus.num, 1, quantile, probs=0.25)
  mod$epi$rplus.num$high<-apply(mod$epi$rplus.num, 1, quantile, probs=0.75)
  
  mod$epi$f.num$median<-apply(mod$epi$f.num, 1, median)
  mod$epi$f.num$low<-apply(mod$epi$f.num, 1, quantile, probs=0.25)
  mod$epi$f.num$high<-apply(mod$epi$f.num, 1, quantile, probs=0.75)
  
  mod$epi$se.flow$median<-apply(mod$epi$se.flow, 1, median)
  mod$epi$se.flow$low<-apply(mod$epi$se.flow, 1, quantile, probs=0.25)
  mod$epi$se.flow$high<-apply(mod$epi$se.flow, 1, quantile, probs=0.75)
  
  mod$epi$ei.flow$median<-apply(mod$epi$ei.flow, 1, median)
  mod$epi$ei.flow$low<-apply(mod$epi$ei.flow, 1, quantile, probs=0.25)
  mod$epi$ei.flow$high<-apply(mod$epi$ei.flow, 1, quantile, probs=0.75)
  
  mod$epi$irplus.flow$median<-apply(mod$epi$irplus.flow, 1, median)
  mod$epi$irplus.flow$low<-apply(mod$epi$irplus.flow, 1, quantile, probs=0.25)
  mod$epi$irplus.flow$high<-apply(mod$epi$irplus.flow, 1, quantile, probs=0.75)
  
  mod$epi$isminus.flow$median<-apply(mod$epi$isminus.flow, 1, median)
  mod$epi$isminus.flow$low<-apply(mod$epi$isminus.flow, 1, quantile, probs=0.25)
  mod$epi$isminus.flow$high<-apply(mod$epi$isminus.flow, 1, quantile, probs=0.75)
  
  mod$epi$if.flow$median<-apply(mod$epi$if.flow, 1, median)
  mod$epi$if.flow$low<-apply(mod$epi$if.flow, 1, quantile, probs=0.25)
  mod$epi$if.flow$high<-apply(mod$epi$if.flow, 1, quantile, probs=0.75)
  
  mod$epi$rpluss.flow$median<-apply(mod$epi$rpluss.flow, 1, median)
  mod$epi$rpluss.flow$low<-apply(mod$epi$rpluss.flow, 1, quantile, probs=0.25)
  mod$epi$rpluss.flow$high<-apply(mod$epi$rpluss.flow, 1, quantile, probs=0.75)
  
  return(mod)
}


addmediansci_model2<-function(mod){
  mod$epi$s.num$median<-apply(mod$epi$s.num, 1, median)
  mod$epi$s.num$low<-apply(mod$epi$s.num, 1, quantile, probs=0.25)
  mod$epi$s.num$high<-apply(mod$epi$s.num, 1, quantile, probs=0.75)
  
  mod$epi$e.num$median<-apply(mod$epi$e.num, 1, median)
  mod$epi$e.num$low<-apply(mod$epi$e.num, 1, quantile, probs=0.25)
  mod$epi$e.num$high<-apply(mod$epi$e.num, 1, quantile, probs=0.75)
  
  mod$epi$i.num$median<-apply(mod$epi$i.num, 1, median)
  mod$epi$i.num$low<-apply(mod$epi$i.num, 1, quantile, probs=0.25)
  mod$epi$i.num$high<-apply(mod$epi$i.num, 1, quantile, probs=0.75)
  
  mod$epi$rplus.num$median<-apply(mod$epi$rplus.num, 1, median)
  mod$epi$rplus.num$low<-apply(mod$epi$rplus.num, 1, quantile, probs=0.25)
  mod$epi$rplus.num$high<-apply(mod$epi$rplus.num, 1, quantile, probs=0.75)
  
  mod$epi$f.num$median<-apply(mod$epi$f.num, 1, median)
  mod$epi$f.num$low<-apply(mod$epi$f.num, 1, quantile, probs=0.25)
  mod$epi$f.num$high<-apply(mod$epi$f.num, 1, quantile, probs=0.75)
  
  mod$epi$v1.num$median<-apply(mod$epi$v1.num, 1, median)
  mod$epi$v1.num$low<-apply(mod$epi$v1.num, 1, quantile, probs=0.25)
  mod$epi$v1.num$high<-apply(mod$epi$v1.num, 1, quantile, probs=0.75)
  
  mod$epi$v2.num$median<-apply(mod$epi$v2.num, 1, median)
  mod$epi$v2.num$low<-apply(mod$epi$v2.num, 1, quantile, probs=0.25)
  mod$epi$v2.num$high<-apply(mod$epi$v2.num, 1, quantile, probs=0.75)
  
  mod$epi$se.flow$median<-apply(mod$epi$se.flow, 1, median)
  mod$epi$se.flow$low<-apply(mod$epi$se.flow, 1, quantile, probs=0.25)
  mod$epi$se.flow$high<-apply(mod$epi$se.flow, 1, quantile, probs=0.75)
  
  mod$epi$ei.flow$median<-apply(mod$epi$ei.flow, 1, median)
  mod$epi$ei.flow$low<-apply(mod$epi$ei.flow, 1, quantile, probs=0.25)
  mod$epi$ei.flow$high<-apply(mod$epi$ei.flow, 1, quantile, probs=0.75)
  
  mod$epi$irplus.flow$median<-apply(mod$epi$irplus.flow, 1, median)
  mod$epi$irplus.flow$low<-apply(mod$epi$irplus.flow, 1, quantile, probs=0.25)
  mod$epi$irplus.flow$high<-apply(mod$epi$irplus.flow, 1, quantile, probs=0.75)
  
  mod$epi$isminus.flow$median<-apply(mod$epi$isminus.flow, 1, median)
  mod$epi$isminus.flow$low<-apply(mod$epi$isminus.flow, 1, quantile, probs=0.25)
  mod$epi$isminus.flow$high<-apply(mod$epi$isminus.flow, 1, quantile, probs=0.75)
  
  mod$epi$if.flow$median<-apply(mod$epi$if.flow, 1, median)
  mod$epi$if.flow$low<-apply(mod$epi$if.flow, 1, quantile, probs=0.25)
  mod$epi$if.flow$high<-apply(mod$epi$if.flow, 1, quantile, probs=0.75)
  
  mod$epi$rpluss.flow$median<-apply(mod$epi$rpluss.flow, 1, median)
  mod$epi$rpluss.flow$low<-apply(mod$epi$rpluss.flow, 1, quantile, probs=0.25)
  mod$epi$rpluss.flow$high<-apply(mod$epi$rpluss.flow, 1, quantile, probs=0.75)
  
  mod$epi$v1v2.flow$median<-apply(mod$epi$v1v2.flow, 1, median)
  mod$epi$v1v2.flow$low<-apply(mod$epi$v1v2.flow, 1, quantile, probs=0.25)
  mod$epi$v1v2.flow$high<-apply(mod$epi$v1v2.flow, 1, quantile, probs=0.75)
  
  mod$epi$sv1.flow$median<-apply(mod$epi$sv1.flow, 1, median)
  mod$epi$sv1.flow$low<-apply(mod$epi$sv1.flow, 1, quantile, probs=0.25)
  mod$epi$sv1.flow$high<-apply(mod$epi$sv1.flow, 1, quantile, probs=0.75)
  
  mod$epi$v2e.flow$median<-apply(mod$epi$v2e.flow, 1, median)
  mod$epi$v2e.flow$low<-apply(mod$epi$v2e.flow, 1, quantile, probs=0.25)
  mod$epi$v2e.flow$high<-apply(mod$epi$v2e.flow, 1, quantile, probs=0.75)
  
  mod$epi$v1e.flow$median<-apply(mod$epi$v1e.flow, 1, median)
  mod$epi$v1e.flow$low<-apply(mod$epi$v1e.flow, 1, quantile, probs=0.25)
  mod$epi$v1e.flow$high<-apply(mod$epi$v1e.flow, 1, quantile, probs=0.75)
  
  return(mod)
}

plot.dcm2 <- function(x, y, popfrac = FALSE, run, col, lwd, lty, alpha = 0.9,
                      legend, leg.name, leg.cex = 0.8, axs = "r", grid = FALSE,
                      add = FALSE, y_names=y_names, ...) {
  
  ## Set missing flags
  noy <- ifelse(missing(y), TRUE, FALSE)
  norun <- ifelse(missing(run), TRUE, FALSE)
  nocol <- ifelse(missing(col), TRUE, FALSE)
  nolwd <- ifelse(missing(lwd), TRUE, FALSE)
  nolty <- ifelse(missing(lty), TRUE, FALSE)
  noleg <- ifelse(missing(legend), TRUE, FALSE)
  
  
  ## Dot args
  da <- list(...)
  
  
  ## Model dimensions
  nsteps <- x$control$nsteps
  nruns <- x$control$nruns
  print(nruns)
  groups <- x$param$groups
  dis.type <- x$control$type
  
  ## Main title default
  if (is.null(da$main)) {
    main <- "today"
  } else {
    main <- da$main
  }
  
  
  ## Defaults for missing y
  if (noy == TRUE && nruns == 1) {
    y <- grep(".num", names(x$epi), value = TRUE)
  }
  if (noy == TRUE && nruns > 1) {
    y <- grep("i.num", names(x$epi), value = TRUE)
  }
  if (all(y %in% names(x$epi)) == FALSE) {
    stop("Specified y is unavailable", call. = FALSE)
  }
  lcomp <- length(y)
  print(lcomp)
  
  ## Prevalence calculations
  x <- denom(x, y, popfrac)
  
  
  ## Compartment ymax calculations
  if (popfrac == FALSE) {
    allmax <- sapply(1:lcomp, function(i) max(x$epi[[y[i]]], na.rm = TRUE))
    ymax <- ceiling(max(allmax))
  } else {
    ymax <- 1
  }
  
  
  ## Defaults for ylim, xlim
  if (is.null(da$ylim)) {
    ylim <- c(0, ymax)
  } else {
    ylim <- da$ylim
  }
  if (is.null(da$xlim)) {
    xlim <- c(0, nsteps)
  } else {
    xlim <- da$xlim
  }
  
  
  ## Defaults for lwd
  if (nolwd == FALSE && lcomp > 1 && length(lwd) < lcomp) {
    lwd <- rep(lwd, lcomp)
  }
  if (nolwd == FALSE && lcomp == 1 && length(lwd) < nruns) {
    lwd <- rep(lwd, nruns)
  }
  if (nolwd == TRUE) {
    lwd <- rep(2.5, lcomp * nruns)
  }
  
  
  ## Defaults for lty
  if (nolty == FALSE && lcomp > 1 && length(lty) < lcomp) {
    lty <- rep(lty, lcomp)
  }
  if (nolty == FALSE && lcomp == 1 && length(lty) < nruns) {
    lty <- rep(lty, nruns)
  }
  if (nolty == TRUE) {
    lty <- rep(1, lcomp * nruns)
    if (groups == 2 && noy == TRUE) {
      lty <- rep(1:2, each = lcomp / 2)
    }
  }
  
  ## Defaults for xlab and ylab
  if (is.null(da$xlab)) {
    xlab <- "Time"
  } else {
    xlab <- da$xlab
  }
  
  if (is.null(da$ylab)) {
    if (popfrac == FALSE) {
      ylab <- "Number"
    } else {
      ylab <- "Prevalence"
    }
  } else {
    ylab <- da$ylab
  }
  
  
  ## Main plot window
  if (add == FALSE) {
    plot(1, 1, type = "n", bty = "n",
         xaxs = axs, yaxs = axs, xlim = xlim, ylim = ylim,
         xlab = 'Days', ylab = 'Individuals', main = main)
  }
  
  
  ## Default line colors
  pal <- NULL
  # Missing col
  if (nocol == TRUE) {
    if (lcomp == 1) {
      if (nruns == 1) {
        col <- "black"
      }
      if (nruns > 1) {
        col <- "Set1"
      }
      if (nruns > 5) {
        col <- "Spectral"
      }
      if (norun == FALSE && length(run) == 1) {
        col <- "black"
      }
    }
    if (lcomp > 1) {
      col <- "Set1"
    }
  }
  
  
  # Test if using a RColorBrewer palette
  if (length(col) == 1 && col %in% row.names(brewer.pal.info)) {
    use.brewer <- TRUE
  } else {
    use.brewer <- FALSE
  }
  
  # Set color palette
  if (is.null(pal)) {
    if (lcomp == 1) {
      if (use.brewer == TRUE) {
        if (nruns < 6) {
          pal <- adjustcolor(brewer.pal(5, col)[1:nruns], alpha)
        } else {
          pal <- adjustcolor(brewer_ramp(nruns, col), alpha)
        }
      }
      if (use.brewer == FALSE) {
        pal <- adjustcolor(rep(col, nruns), alpha)
      }
    }
    if (lcomp > 1) {
      if (use.brewer == TRUE) {
        if (lcomp > 4) {
          pal <- adjustcolor(brewer_ramp(lcomp, col), alpha)
        } else {
          pal <- adjustcolor(brewer.pal(max(c(lcomp, 4)), col), alpha)
          fixpal <- pal
          fixpal[1] <- pal[2]; fixpal[2] <- pal[1]
          pal <- fixpal
        }
        if (groups == 2 && noy == TRUE) {
          pal <- adjustcolor(brewer.pal(3, col), alpha)
          fixpal <- pal
          fixpal[1] <- pal[2]; fixpal[2] <- pal[1]
          pal <- fixpal
          if (dis.type != "SIR") {
            pal <- pal[1:2]
          }
          pal <- rep(pal, times = lcomp / 2)
        }
      }
      if (use.brewer == FALSE) {
        pal <- adjustcolor(rep(col, lcomp), alpha)
        if (groups == 2 && noy == TRUE) {
          pal <- adjustcolor(rep(col, times = 2), alpha)
        }
      }
    }
  }
  
  
  ## Plot lines
  
  if (lcomp > 1) {
    if (nruns == 1) {
      for (i in 1:lcomp) {
        lines(x$control$timesteps, x$epi[[y[i]]][, nruns+1],
              lwd = lwd, lty = lty[i], col = pal[i])
      }
    }
    if (nruns > 1) {
      if (norun == TRUE) {
        for (i in 1:lcomp) {
          run <- 1
          lines(x$control$timesteps, x$epi[[y[i]]][, nruns+1],
                lwd = lwd[i], lty = lty[i], col = pal[i])
          
        }
      }
      if (norun == FALSE) {
        if (length(run) > 1) {
          stop("Plotting multiple runs of multiple y is not supported",
               call. = FALSE)
        }
        for (i in 1:lcomp) {
          lines(x$control$timesteps, x$epi[[y[i]]][, nruns+1],
                lwd = lwd[i], lty = lty[i], col = pal[i])
          polygon(c(x$control$timesteps, rev(x$control$timesteps)), c(x$epi[[y[i]]][, nruns+3], rev(x$epi[[y[i]]][, nruns+2])),
                  col = alpha(pal[i],0.2), border=FALSE)
          #lines(x$control$timesteps, x$epi[[y[i]]][, nruns+2],
          #     lwd = lwd[i], lty = lty[i], col = alpha(pal[i],0.2))
          #lines(x$control$timesteps, x$epi[[y[i]]][, nruns+3],
          #     lwd = lwd[i], lty = lty[i], col = alpha(pal[i],0.2))
        }
      }
    }
  }
  
  ## Grid
  if (grid == TRUE) {
    grid()
  }
  
  ## Legend
  
  # Default legend type
  if (noleg == TRUE) {
    legend <- "n"
    if (lcomp == 1 & nruns < 3) {
      legend <- "full"
    }
    if (lcomp == 1 & nruns >= 3) {
      legend <- "lim"
    }
    if (lcomp > 1) {
      legend <- "full"
    }
    if (noy == FALSE) {
      legend <- "n"
    }
  } else {
    if (legend == "lim" & nruns < 3) {
      legend <- "full"
    }
    if (legend == "lim" & lcomp == 2) {
      legend <- "full"
    }
  }
  
  # Default legend names
  if (missing(leg.name)) {
    if (nruns == 1) {
      leg.names <- y
    }
    if (nruns > 1) {
      if (norun == TRUE & lcomp == 1) {
        leg.names <- names(x$epi[[y[1]]])
      }
      if (norun == FALSE & lcomp == 1) {
        if (length(run) == 1) {
          leg.names <- y_names
        }
        if (length(run) > 1) {
          leg.names <- names(x$epi[[y[1]]][run])
        }
      }
      if (lcomp > 1) {
        leg.names <- y_names
        print(leg.names)
      }
    }
  } else {
    if (lcomp == 1) {
      leg.names <- paste(leg.name, 1:nruns)
    }
    if (lcomp > 1) {
      leg.names <- y
      warning("Legend names ignored for multiple y plots of multiple run
              models", call. = FALSE)
    }
  }
  
  # Legend
  if (norun == TRUE) {
    if (legend == "full") {
      legend("topright", legend = leg.names,
             bg = "white", lty = lty, lwd = lwd,
             col = pal, cex = leg.cex)
    }
    if (legend == "lim") {
      legend("topright",
             legend = c(leg.names[1], "...", leg.names[nruns]),
             bg = "white",
             lty = c(lty[1], 1, lty[nruns]), lwd = lwd + 1,
             col = c(pal[1], "white", pal[nruns]), cex = leg.cex)
    }
  }
  if (norun == FALSE & legend != "n") {
    if (lcomp == 1) {
      legend("topright", legend = leg.names,
             bg = "white", lty = lty[1:length(run)],
             lwd = lwd[1:length(run)],
             col = pal[1:length(run)], cex = leg.cex)
    }
    if (lcomp > 1) {
      legend("topright", legend = leg.names,
             bg = "white", lty = lty, lwd = lwd,
             col = pal, cex = leg.cex)
      print('legend')
    }
  }
  
}



