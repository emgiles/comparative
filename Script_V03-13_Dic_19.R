# checking packages ---------------------------------------------------------------
source("setup_packages.R")

# loading dependences -------------------------------------------------------------
library(ape)
library(phytools)
library(phangorn)
library(geiger)
library(l1ou)
library(ips)

# reading and ordering comparative datasets ---------------------------------------
amphibia <- read.table("data_amphibia.txt", h=T, sep="\t")
str(amphibia)

reptiles <- read.table("data_reptilia.txt", h=T, sep="\t")
str(reptiles)

phylo.Hedges <- read.tree("~/Dropbox/Filogenias/Hedges et al 2015/time-tree.txt")

tr.amphibia <- treedata(phylo.Hedges, 
           setNames(amphibia$sp, amphibia$sp))[[1]]
tr.reptiles <- treedata(phylo.Hedges, 
           setNames(reptiles$sp, reptiles$sp))[[1]]

# l1ou in amphibia ---------------------------------------------------------------
vars.amphibia <- amphibia
rownames(vars.amphibia) <- as.character(vars.amphibia$sp)
nums <- unlist(lapply(vars.amphibia, is.numeric)) # which columns are numeric? 
vars.amphibia <- vars.amphibia[,nums] # selecting only those columns in vars.amphibia that are numeric
vars.amphibia <- vars.amphibia[,-1] # omit ID column or another variables that we should omit.
vars.amphibia <- vars.amphibia[match(tr.amphibia$tip.label, rownames(vars.amphibia)),] # reordering dataset 
str(vars.amphibia)
rownames(vars.amphibia)
colnames(vars.amphibia)
all.equal(rownames(vars.amphibia), tr.amphibia$tip.label)

# creating several empty list for each model that we will test.
  amphibia.l1ou <- adjust_data(tr.amphibia, vars.amphibia[,"CT_max"], normalize = FALSE)
  eModel.full <- estimate_shift_configuration(amphibia.l1ou$tree, 
                        amphibia.l1ou$Y, criterion="BIC", 
                        max.nShifts=45, nCores = 3)
  plot(eModel.full)

# 
sp.opt <- tr.amphibia$tip.label[unlist(Descendants(tr.amphibia, node = 126, type="tips"))]
amphibia[amphibia$sp %in% sp.opt,]
mat.amphibia <- as.matrix(amphibia[,c("Latitude", "Longitude")])
rownames(mat.amphibia) <- amphibia$sp
colnames(mat.amphibia) <- c("Lat", "Long")
obj <- phylo.to.map(tr.amphibia, mat.amphibia, plot=FALSE)
#pdf(file="AmphibiosMaps.pdf",width=22,height=19)
plot(obj,direction = "downwards", colors="gray60",cex.points=c(0,2), lwd=c(4,2),ftype="i", fsize=1.3)
#dev.off()

plot(tr.amphibia)
# columns = models; rows = variables/traits
K.amphibia <- data.frame(matrix(ncol = 1, nrow = 6))
colnames(K.amphibia) <- c("CT_max")
rownames(K.amphibia) <- c("Nshift", "alpha", "sigma", "BIC", "Tetha_1", "Tetha_2")

K.amphibia$CT_max <- round(c(eModel.full$nShifts, 
        eModel.full$alpha,
        eModel.full$sigma2,
        eModel.full$score, 
  as.numeric(unique(as.character(round(eModel.full$optima, 2)))[1]), 
  as.numeric(unique(as.character(round(eModel.full$optima, 2)))[2])), 2)

# AMPHIBIA: PLOT PHYLOGENY - CTMAX - OPTIMA --------------------------------------
amp.vector <- setNames(vars.amphibia[,1], rownames(vars.amphibia))

plot(eModel.full, asterisk = F, edge.shift.ann = F)
nodelabels(frame = "n") ## The shift is located on the node 126

dec.amphibia <- tr.amphibia$tip.label[Descendants(tr.amphibia, node = 126)[[1]]]
cols.amphibia <- setNames(rep("gray", length(tr.amphibia$tip.label)), tr.amphibia$tip.label)
cols.amphibia[dec.amphibia] <- "red"

tr.amphibia2 <- paintSubTree(tr.amphibia, node = 126, state = "Shift_1", stem = T)
all.equal(tr.amphibia2$tip.label, names(amp.vector))
cols <- setNames(c("gray", "red"), c("1","Shift_1"))

#pdf("Fig shifts Amphibia CT_Max.pdf", width = 8, height = 9)
ss<-getStates(tr.amphibia2,"tips")
barcols<-setNames(sapply(ss,function(x,y) y[which(names(y)==x)],
                         y=cols),names(ss))

plotTree.barplot(tr.amphibia2,amp.vector[tr.amphibia2$tip.label],
                 args.plotTree=list(ftype="i", fsize=0.7),
                 args.barplot=list(col=barcols,border=barcols,
                                   xlab="CT_Max (°C)",
                                   xlim=c(20,60)),
                 args.axis=list(at=seq(20,50,by=10)))
abline(v=31.24506, col="red", lty=2, lwd=1.5)
abline(v=36.60905, col="grey40", lty=2, lwd=1.5)

## this is how I go back to my tree frame
par(mfg=c(1,1))
plot(tr.amphibia2,cols, ftype='i', mar=c(5.1,1.1,2.1,0), fsize=0.7, offset=-50)
#dev.off()

# l1ou in reptiles ---------------------------------------------------------------
vars.reptiles <- reptiles
rownames(vars.reptiles) <- as.character(vars.reptiles$sp)
nums <- unlist(lapply(vars.reptiles, is.numeric)) # which columns are numeric? 
vars.reptiles <- vars.reptiles[,nums] # selecting only those columns in vars.reptiles that are numeric
vars.reptiles <- vars.reptiles[,-1] # omit ID column or another variables that we should omit.
vars.reptiles <- vars.reptiles[match(tr.reptiles$tip.label, rownames(vars.reptiles)),] # reordering dataset 
str(vars.reptiles)
rownames(vars.reptiles)
colnames(vars.reptiles)
all.equal(rownames(vars.reptiles), tr.reptiles$tip.label)

# creating several empty list for each model that we will test.
reptiles.l1ou <- adjust_data(tr.reptiles, vars.reptiles[,"CT_max"], normalize = FALSE)
eModel.full <- estimate_shift_configuration(reptiles.l1ou$tree, 
                                            reptiles.l1ou$Y, criterion="BIC", 
                                            max.nShifts=45, nCores = 3)
eModel.full <- estimate_convergent_regimes(eModel.full, criterion="BIC")
plot(eModel.full, edge.shift.ann = F,asterisk = F, plot.bar = F)
nodelabels(cex=0.5, frame="none")
# to map
sp.opt.131 <- tr.reptiles$tip.label[unlist(Descendants(tr.reptiles, node = 131, type="tips"))]
sp.opt.147 <- tr.reptiles$tip.label[unlist(Descendants(tr.reptiles, node = 147, type="tips"))]
sp.opt.133 <- tr.reptiles$tip.label[unlist(Descendants(tr.reptiles, node = 133, type="tips"))]
sp.opt.131 <- sp.opt.131[!sp.opt.131 %in% c(sp.opt.147, sp.opt.133)]
plot(eModel.full, edge.shift.ann = F,asterisk = F, plot.bar = T)
reptiles[reptiles$sp %in% c(sp.opt.147, sp.opt.133),] # grupo con valores bajos

mat.reptiles <- as.matrix(reptiles[,c("Latitude", "Longitude")])
rownames(mat.reptiles) <- reptiles$sp
colnames(mat.reptiles) <- c("Lat", "Long")
obj <- phylo.to.map(tr.reptiles, mat.reptiles, plot=FALSE)
#pdf(file="ReptilesMaps.pdf",width=22,height=19)
plot(obj,direction = "downwards", colors="gray60",cex.points=c(0,2), lwd=c(4,2),ftype="i", fsize=1.3)
#dev.off()

plot(tr.reptiles)
# columns = models; rows = variables/traits
K.reptiles <- data.frame(matrix(ncol = 1, nrow = 6))
colnames(K.reptiles) <- c("CT_max")
rownames(K.reptiles) <- c("Nshift", "alpha", "sigma", "BIC", "Theta_1", "Theta_2")

K.reptiles$CT_max <- round(c(eModel.full$nShifts, 
                             eModel.full$alpha,
                             eModel.full$sigma2,
                             eModel.full$score, 
                             as.numeric(unique(as.character(round(eModel.full$optima, 2)))[1]), 
                             as.numeric(unique(as.character(round(eModel.full$optima, 2)))[2])), 2)

# REPTILES: PLOT PHYLOGENY - CTMAX - OPTIMA --------------------------------------
rep.vector <- setNames(vars.reptiles[,1], rownames(vars.reptiles))

plot(eModel.full, asterisk = F, edge.shift.ann = F, plot.bar = T, 
     show.tip.label=F)
nodelabels(frame = "n") 
## The shifts are located on the node 131 (blue), 147 and 133 (red, stem=T)

dec.amphibia <- tr.amphibia$tip.label[Descendants(tr.amphibia, node = 126)[[1]]]
cols.amphibia <- setNames(rep("gray", length(tr.amphibia$tip.label)), tr.amphibia$tip.label)
cols.amphibia[dec.amphibia] <- "red"

tr.reptiles2 <- paintSubTree(tr.reptiles, node = 131, state = "Shift_1", stem = T)
tr.reptiles2 <- paintSubTree(tr.reptiles2, node = 147, state = "Shift_2",
                             stem = T)
tr.reptiles2 <- paintSubTree(tr.reptiles2, node = 133, state = "Shift_2",
                             stem = T)
all.equal(tr.amphibia2$tip.label, names(amp.vector))
cols <- setNames(c("gray", "red", "blue"), c("1","Shift_1", "Shift_2"))
plotSimmap(tr.reptiles2, colors = cols)

#pdf("Fig shifts Reptiles CT_Max.pdf", width = 8, height = 10)
ss<-getStates(tr.reptiles2,"tips")
barcols<-setNames(sapply(ss,function(x,y) y[which(names(y)==x)],
                         y=cols),names(ss))

plotTree.barplot(tr.reptiles2,rep.vector[tr.reptiles2$tip.label],
                 args.plotTree=list(ftype="i", fsize=0.5),
                 args.barplot=list(col=barcols,border=barcols,
                                   xlab="CT_Max (°C)",
                                   xlim=c(20,60)),
                 args.axis=list(at=seq(20,50,by=10)))
abline(v=mean(c(39.53708, 38.07875)), col="darkblue", lty=2, lwd=1.5)
abline(v=43.71827, col="darkred", lty=2, lwd=1.5)
abline(v=41.52152, col="grey40", lty=2, lwd=1.5)

## this is how I go back to my tree frame
par(mfg=c(1,1))
plot(tr.reptiles2,cols, ftype='i', mar=c(5.1,1.1,2.1,0), fsize=0.5, offset=-50)
#dev.off()
