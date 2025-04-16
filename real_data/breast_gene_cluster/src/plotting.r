# ################################################################
# ################################################################
## PLOTTING FUNCTIONS:

cbPalette <- c("#999999", "#E69F00", "#56B4E9",
               "#009E73", "#F0E442", "#0072B2",
               "#D55E00", "#CC79A7")

###################################
## Let us create a function that takes a
##    number of plots desired,
##    an ordering vector that gives a score to each variable,
## and generates a multiple plot that shows a sample of the
## variables in the data, ordering the sample by the given
## scores.
plot.mardist <- function(data = NULL,
                        n, p,
                        scores = NULL,
                        random = TRUE,
                        minmax = TRUE,
                        nprow, npcol,
                        lwd = 1, lcol = "black",
                        pwd = 0.5, pcol = "blue",
                        main = NULL, cex.omain = 3, # nolint
                        cex.main = 1) { # nolint
  .nplot <- nprow * npcol

  #### Step 1: Select which variables to plot.
  ####
  ## Order variables by their score.
  .orderscores <- order(scores,
                       decreasing = TRUE)
  .indexes <- rep(NA, .nplot)

  par(mfrow = c(1, 1))
  if (minmax && random) { ## If we choose randomly, including min/max variables.
    .indexes[1]      <- .orderscores[1]
    .indexes[.nplot] <- .orderscores[p]

    ## Choosing which positions to select
    .ordind <- sample(x = 2:(p - 1),
                        size = .nplot - 2,
                        replace = FALSE)
    .ordind <- sort(.ordind,
                   decreasing = FALSE)
    .indexes[2:(.nplot - 1)] <- .orderscores[.ordind]
  } else if (!minmax && random) { ## If we choose variables fully randomly:
    .ordind <- sample(x = 1:p,
                     size = .nplot,
                     replace = FALSE)
    .ordind <- sort(.ordind,
                   decreasing = FALSE)
    .indexes <- .orderscores[.ordind]
    } else if (!random) { ## If we select a deterministic grid of values:
    a <- floor(p / (.nplot - 1))
    r <- p - a * .nplot

    .ordind <- a * (1:(.nplot - 1)) + r
    .ordind <- c(1, .ordind)

    .indexes <- .orderscores[.ordind]
  }

  ## Once we determined the set of variables
  ## to plot, we plot them.
  par(mfrow = c(nprow, npcol),
      oma = c(0, 0, 4, 0))

  for (.var in .indexes){
    .missing <- is.na(data[, .var])
    .nvarcol <- sum(!.missing)
    .varcol <- data[!.missing, .var]
    .pmain <- paste(wrap_sentence(string = colnames(data)[.var],
                                 width = 30), "\n",
                   round(scores[.var], digits = 2))

    .densityplot <- density(.varcol)
    plot(.densityplot,
         main = .pmain,
         cex.main = cex.main,
         col = lcol,
         lwd = lwd)

    abline(v = 0)

    ## Create cloud plot:
    .ymax <- max(.densityplot$y)
    points(x = .varcol,
           y = runif(.nvarcol,
                     min = .ymax / 3,
                     max = 2 * .ymax / 3),
           pch = 19,
           col = pcol,
           cex = pwd)
  }
  # title
  mtext(text = main,
        side = 3,
        outer = TRUE,
        cex = cex.omain)

  ## restore original settings
}

wrap_sentence <- function(string, width) {
  words <- unlist(strsplit(string, "_"))
  fullsentence <- ""
  checklen <- ""
  for (i in 1:length(words)) { # nolint
    checklen <- paste(checklen, words[i])
    if (nchar(checklen) > (width + 1)) {
      fullsentence <- paste0(fullsentence, "\n")
      checklen <- ""
    }
    fullsentence <- paste(fullsentence, words[i])
  }
  fullsentence <- sub("^\\s", "", fullsentence)
  fullsentence <- gsub("\n ", "\n", fullsentence)
  return(fullsentence)
}



###################################
plot.pdm <- function( # nolint
  data, main, breaks = 21, diagonal = FALSE,
  cex = 0.05, cex.main = 1.5, entryrange = NULL) {
  .d <- dim(data)[1]

  if (is.null(entryrange)) {
    entryrange <- 0
    if (diagonal) {
      entryrange <- max(abs(data))
    } else {
      for (.i in 1:.d) {
        for (.j in (1:.d)[-.i]) {
          entryrange <- max(entryrange, abs(data[.i, .j]))
        }
      }
    }
  }

  rb_pal <- colorRampPalette(c("red", "gray95", "blue"))
  color <- rb_pal(breaks)
  plot(1, 1,
       col = "white",
       xlim = c(1, .d),
       ylim = c(1, .d),
       xlab = "Row Index",
       ylab = "Column Index",
       xaxt = "none",
       yaxt = "none",
       main = main,
       cex.main = cex.main)
  axis(2, at = c(1, seq(from = 0, to = .d, length.out = 6)[-1]),
       labels = c(seq(from = .d, to = 0, length.out = 6)[-6], 1))
  axis(1, at = c(1, seq(from = 0, to = .d, length.out = 6)[-1]),
       labels = c(1, seq(from = 0, to = .d, length.out = 6)[-1]))
  
  for (.i in 1:.d) {
    if (diagonal) {
      for (.j in 1:.d) {
        colorindex <- floor( 
          breaks * (data[.i, .j] + entryrange) / (2 * entryrange) + 0.99
        )
        points(x = c(.j), y = c(.d - .i + 1), col = color[colorindex],
               pch = 19, cex = cex)
      }
    } else {
      for (.j in (1:.d)[-.i]) {
        colorindex <- floor( 
          breaks * (data[.i, .j] + entryrange) / (2 * entryrange) + 0.99
        )
        points(x = c(.j), y = c(.d - .i + 1), col = color[colorindex],
               pch = 19, cex = cex)
      }
    }
  }
}