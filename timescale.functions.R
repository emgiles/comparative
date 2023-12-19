# Two timescale plotting functions
# carl simpson, 4-16-2015

# download these files so that the functions know interval boundaries
periods <- data.frame(tpx=c("Cm","O","S","D","C","P","Tr","J","K","Pg","Ng","R"),
                      periods=c(541,485,444,419,359,299,252,201,145,66,23,0))
epochs <- data.frame(epoch.age=c(541.00,521.00,509.00,497.00,485.00,470.00,458.00,444.00,433.00,
            427.00,423.00,419.00,
            393.00,383.00,359.00,347.00,331.00,323.00,315.00,307.00,299.00,272.00,260.00,252.00,
            247.00,237.00,201.00,174.00,164.00,145.00,100.00,66.00,56.00,33.90,23.00,5.30,2.60,
            0.01,0.00))
# Funciton for period scale on ribbon and backdrop
tscales.period <- function(top, bot, s.bot, col, lwd, col.border,...){
  bases <- periods
  cc <- rep(c("white","white "),length(bases))

  rect(xleft = bases[-12, 2], ybottom = rep(bot,11), xright = bases[-1, 2],
       ytop = rep(top, 11), col = cc, border=NA)

  rect(xleft = bases[-12, 2], ybottom = rep(bot,11), xright = bases[-1, 2],
       ytop = rep(s.bot, 11), border = col.border, col=col, lwd = lwd)

  bt <- (bot+s.bot)/2
  tpl <- bases[,2]+c(diff(bases[,2])/2,0)

  text(x=tpl[-12], y=bt[-12], labels = bases[1:11, 1])
}

# Funciton for period scale on ribbon and epoch scale on backdrop
tscales.epoch <- function(top, bot, s.bot){
  labels <- periods
  bases <- epochs
  end <- dim(bases)[1]
  cc <- rep(c("white","white"),length(bases))

  rect(xleft = bases[-end,], ybottom = rep(bot, end), xright = bases[-1,],
       ytop = rep(top, end), col = cc, border=NA)

  rect(xleft = labels[-12, 2], ybottom = rep(bot,11), xright = labels[-1, 2],
       ytop = rep(s.bot, 11), border = 'grey90')

  bt <- (bot+s.bot)/2
  tpl <- labels[,2]+c(diff(labels[,2])/2,0)

  text(x=tpl[-12], y=bt[-12], labels = labels[1:11, 1])
}
