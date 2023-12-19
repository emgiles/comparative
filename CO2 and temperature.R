rm(list=ls())
source("timescale.functions.R")
co2 <- read.table("CO2_Foster2017.txt", h=T, sep = "\t")
co2 <- co2[co2$age<=300,]
co2$age2 <- co2$age*-1
sealevel <- read.table("SeaLevel_Hannisdal2011.txt", h=T, sep = "\t")
sealevel <- sealevel[sealevel$age<=300,]
sealevel$age2 <- sealevel$age*-1

#pdf("pCO2 and temperature.pdf", width = 10, height = 8)
par(mfrow=c(1,1), mar=c(6,6,2,5.8)) 
plot(1, type="n", xlab="", ylab="", cex.axis=1.3, yaxt="n", bty="n",
     xlim=rev(c(0, 300)), ylim=c(-700, 4000), las=1, xaxt="n", col.axis="blue")

tscales.period(4000, -600, -900, lwd=1.5, col.border="gray75",
                col = "gray90")                       
#col = c("palegreen4", "dodgerblue3", "green","olivedrab", 
#"white","orange", "brown"))

with(co2,polygon(c(age,rev(age)),c(up95,rev(lw95)), 
                                         col="grey80", border=NA))
with(co2,polygon(c(age,rev(age)),c(up68,rev(lw68)), 
                                         col="gray60", border=NA))
lines(pCO2 ~ age, data = co2, type="l", las=1, col="blue", lwd=3)
#box(lwd=2.5)
mtext("pCO2",side = 2, cex=2.5, font=4, line = 3.5, col = "blue")
axis(side = 2, at = seq(500,2500, 500), tck=-0.01, lwd = 3, las=1, cex.axis=0.85,
     labels = c("500","","1500","","2500"),
     col.ticks = "blue", col.axis="blue")
axis(side = 2, at = seq(0,3000, 1000), labels = seq(0,3000, 1000), lwd = 3, las=1,
     col="blue", col.axis="blue")

par(new = TRUE)
plot(1, type="n", xlab="", ylab="", cex.axis=1.3, yaxt="n", bty="n",
     xlim=rev(c(0, 300)), ylim=c(-30, 50), las=1, xaxt="n", col.axis="blue")
lines(temp_Epstein1953 ~ age, data = sealevel, type="l", las=1, col="red", lwd=3)
axis(side = 4, at = seq(15,45, 10), tck=-0.01, lwd = 3, las=1, cex.axis=0.7, labels = F,
     col.ticks = "red", col.axis="red", col="red")
axis(side = 4, at = seq(10,50, 10), lwd = 3, las=1,
     col="red", col.axis="red")
mtext("Temperature (°C)",side = 4, cex=1.5, font=4, line = 3.5, col = "red", adj = 0.9)

axis(1, lwd=2)
mtext("Time before present (Ma)",side = 1, cex=1.5, font=4, line = 3, col = "black", adj = 0.5)
#dev.off()

