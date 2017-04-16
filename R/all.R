rm(list=ls(all=T))
cat("\014") 

require(fmsb)
require(car)
require(boot)
require(coefplot)
require(MASS)

# Set working directory to script location!
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

source('coefplot.r')

save2disk = T
options(scipen=0, digits=4)

if (save2disk) sink("../Results/weil-obwohl.txt")

################## DATA IMPORT + EDIT ##################

weilobwohl <- read.delim("../Data/Sample2015.csv")
weilobwohl$Senttype <- factor(weilobwohl$Senttype, levels=c("VL", "V2"))
weilobwohl <- weilobwohl[-which(is.na(weilobwohl$Senttype)),]
weilobwohl$Hypo <- factor(ifelse(weilobwohl$Hypo==1, "1", "0"), levels=c("0", "1"))
weilobwohl$Majus <- factor(weilobwohl$Majus, levels=c("0", "1"))
weilobwohl <- weilobwohl[-which(is.na(weilobwohl$Majus)),]
weilobwohl$Left <- factor(weilobwohl$Left, levels=c("Comma", "Word", "End", "Paro", "Three", "Hyphen", "Emo", "Colon"))
weilobwohl <- weilobwohl[-which(is.na(weilobwohl$Left)),]
weilobwohl$Right[which(weilobwohl$Right=="Emo" | weilobwohl$Right=="End" | weilobwohl$Right=="Other" | weilobwohl$Right=="Parc" | weilobwohl$Right=="Paro" | weilobwohl$Right=="Quot")] <- "Other"
weilobwohl$Right <- factor(weilobwohl$Right, levels=c("Word", "Comma", "Three", "Hyphen", "Colon"))
weilobwohl <- weilobwohl[-which(is.na(weilobwohl$Right)),]
weilobwohl$ModBin <- factor(ifelse(weilobwohl$ModBin==1, "1", "0"), levels=c("0", "1"))
weilobwohl$Mod <- factor(weilobwohl$Mod, levels=c("0", "Nur", "Und", "Unddas", "Other", "Nicht", "Sondern", "Aber", "Eben", "Gerade", "Vorallem", "Auch", "Denn", "Vielleicht"))
weilobwohl$Mod[which(is.na(weilobwohl$Mod))] <- "Other"
weilobwohl$Isolated <- factor(weilobwohl$Isolated)

# Create factor encoding simply whether the clause is independent or not.
weilobwohl$Independent <- factor(ifelse(weilobwohl$Isolated=="1" | (weilobwohl$Senttype == "V2" & weilobwohl$Left %in% c("End","Emo","Paro")), 1, 0))

weil <- weilobwohl[which(weilobwohl$Target=="Weil"),]
weil$Left <- factor(weil$Left, levels=c("Comma", "Word", "End", "Three", "Paro", "Hyphen", "Emo", "Colon"))
weil$Right <- factor(weil$Right, levels=c("Word", "Hyphen", "Colon", "Three", "Comma"))
weil$Mod <- factor(weil$Mod, levels=c("0", "Nur", "Und", "Other", "Nicht", "Sondern", "Eben", "Gerade", "Vorallem", "Auch", "Vielleicht", "Aber"))

weil.vl <- weil[which(weil$Senttype=="VL"),]
weil.vl$Right <- factor(weil.vl$Right, levels=c("Word", "Hyphen"))
weil.vl <- weil.vl[-which(is.na(weil.vl$Right)),]

weil.v2 <- weil[which(weil$Senttype=="V2"),]
weil.v2$Left <- factor(weil.v2$Left, levels=c("Comma", "End", "Word", "Emo", "Three", "Hyphen", "Paro", "Colon"))
weil.v2$Mod <- factor(weil.v2$Mod, levels=c("0", "Nicht", "Auch", "Nur"))

obwohl <- weilobwohl[which(weilobwohl$Target=="Obwohl"),]
obwohl$Left <- factor(obwohl$Left, levels=c("Comma", "End", "Word", "Paro", "Three", "Hyphen", "Emo", "Colon"))
obwohl$Right <- factor(obwohl$Right, levels=c("Word", "Comma", "Three", "Colon", "Hyphen"))
obwohl$Mod <- factor(obwohl$Mod, levels=c("0", "Unddas", "Und", "Aber", "Denn", "Other"))
obwohl$Mod[which(is.na(obwohl$Mod))] <- "Other"

obwohl.vl <- obwohl[which(obwohl$Senttype=="VL"),]
obwohl$Left <- factor(obwohl$Left, levels=c("Comma", "Word", "End", "Paro", "Three", "Hyphen", "Emo", "Colon"))

obwohl.v2 <- obwohl[which(obwohl$Senttype=="V2"),]
obwohl.v2$Left <- factor(obwohl.v2$Left, levels=c("End", "Comma", "Three", "Emo", "Word", "Paro", "Hyphen", "Colon"))

v2 <- rbind(weil.v2, obwohl.v2)
weilobwohl.vl <- rbind(weil.vl, obwohl.vl)


# =================================================
# ================== DECRIPTIVE ===================

cat("\n\n--------- DESCRIPTIVE STATISTICS TABLE ---------\n\n")

percentize <- function(v) round(v/sum(v)*100, 2)

cat("\n Left\n")
LeftTab <- t(table(weilobwohl$Target,weilobwohl$Left))
LeftTab <- cbind(LeftTab[,1], percentize(LeftTab[,1]), LeftTab[,2], percentize(LeftTab[,2]))
LeftTab <- rbind(LeftTab, apply(LeftTab, 2, sum))
colnames(LeftTab) <- c("Obw", "Obw %", "Weil", "Weil %")
print(LeftTab)

cat("\n Right\n")
RightTab <- t(table(weilobwohl$Target,weilobwohl$Right))
RightTab <- cbind(RightTab[,1], percentize(RightTab[,1]), RightTab[,2], percentize(RightTab[,2]))
RightTab <- rbind(RightTab, apply(RightTab, 2, sum))
colnames(RightTab) <- c("Obw", "Obw %", "Weil", "Weil %")
print(RightTab)

cat("\n Mod\n")
ModTab <- t(table(weilobwohl$Target,weilobwohl$ModBin))
ModTab <- cbind(ModTab[,1], percentize(ModTab[,1]), ModTab[,2], percentize(ModTab[,2]))
ModTab <- rbind(ModTab, apply(ModTab, 2, sum))
colnames(ModTab) <- c("Obw", "Obw %", "Weil", "Weil %")
print(ModTab)

cat("\n Hypo\n")
HypoTab <- t(table(weilobwohl$Target,weilobwohl$Hypo))
HypoTab <- cbind(HypoTab[,1], percentize(HypoTab[,1]), HypoTab[,2], percentize(HypoTab[,2]))
HypoTab <- rbind(HypoTab, apply(HypoTab, 2, sum))
colnames(HypoTab) <- c("Obw", "Obw %", "Weil", "Weil %")
print(HypoTab)

cat("\n Senttype\n")
SenttypeTab <- t(table(weilobwohl$Target,weilobwohl$Senttype))
SenttypeTab <- cbind(SenttypeTab[,1], percentize(SenttypeTab[,1]), SenttypeTab[,2], percentize(SenttypeTab[,2]))
SenttypeTab <- rbind(SenttypeTab, apply(SenttypeTab, 2, sum))
colnames(SenttypeTab) <- c("Obw", "Obw %", "Weil", "Weil %")
print(SenttypeTab)

# =================================================
# ================= OVERVIEW PLOT =================

# plt the extra bars for section 4
rnams <- c("obwohl", "weil")
colos <- c("#555555", "#CCCCCC")
cnams.l <- c("Comma", "End", "Word", "Emo", "Three", "Hyphen", "Paro", "Colon")
cnams.l.pretty <- c("Comma", "End", "Word", "Emo", "Ellipsis", "Dash", "Paro", "Colon")
cnams.r <- c("Word", "Hyphen", "Colon", "Three", "Comma")
cnams.r.pretty <- c("Word", "Dash", "Colon", "Ellipsis", "Comma")
cexo <- 0.8
cexa <- 1.1
offso.l <- 5
offso.r <- 10
offsi <- 0

# Plot LEFT

adhoc.l.w <- function(intp) { length(which(v2$Left==intp & v2$Target=="Weil")) }
adhoc.l.o <- function(intp) { length(which(v2$Left==intp & v2$Target=="Obwohl")) }
leftObwCount <- as.numeric(lapply(cnams.l,adhoc.l.o))
leftObw <- round(leftObwCount/sum(leftObwCount)*100, 2)
leftObwLab <- paste(leftObw, "% (", leftObwCount, ")", sep="")
leftWeiCount <- as.numeric(lapply(cnams.l,adhoc.l.w))
leftWei <- round(leftWeiCount/sum(leftWeiCount)*100, 2)
leftWeiLab <- paste(leftWei, "% (", leftWeiCount, ")", sep="")
v2.left <- rbind(leftObw, leftWei)
rownames(v2.left) <- rnams
colnames(v2.left) <- cnams.l.pretty

if (save2disk) svg("../Results/v2.left.svg")
v2.l.bp <- barplot(v2.left, beside=TRUE, col=colos, ylim=c(0,59), cex.axis = cexa, cex.names = cexa)
legend(15 , 45, rnams, fill=colos)
text(x = v2.l.bp[1,]-offsi, y = v2.left[1,]+offso.l, labels = leftObwLab, cex=cexo, srt=90)
text(x = v2.l.bp[2,]-offsi, y = v2.left[2,]+offso.l, labels = leftWeiLab, cex=cexo, srt=90)
if (save2disk) dev.off()

# Plot RIGHT

adhoc.r.w <- function(intp) { length(which(v2$Right==intp & v2$Target=="Weil")) }
adhoc.r.o <- function(intp) { length(which(v2$Right==intp & v2$Target=="Obwohl")) }

#v2.right <- rbind(as.numeric(lapply(cnams.r,adhoc.r.o)), as.numeric(lapply(cnams.r,adhoc.r.w)))
#xv2.right <- rbind(as.numeric(lapply(cnams.r,adhoc.r.o)), as.numeric(lapply(cnams.r,adhoc.r.w)))

rightObwCount <- as.numeric(lapply(cnams.r,adhoc.r.o))
rightObw <- round(rightObwCount/sum(rightObwCount)*100, 2)
rightObwLab <- paste(rightObw, "% (", rightObwCount, ")", sep="")
rightWeiCount <- as.numeric(lapply(cnams.r,adhoc.r.w))
rightWei <- round(rightWeiCount/sum(rightWeiCount)*100, 2)
rightWeiLab <- paste(rightWei, "% (", rightWeiCount, ")", sep="")
v2.right <- rbind(rightObw, rightWei)

rownames(v2.right) <- rnams
colnames(v2.right) <- cnams.r.pretty

if (save2disk) svg("../Results/v2.right.svg")
v2.r.bp <- barplot(v2.right, beside=TRUE, col=colos, ylim=c(0,109), cex.axis = cexa, cex.names = cexa)
legend(10, 85, rnams, fill=colos)
text(x = v2.r.bp[1,]-offsi, y = v2.right[1,]+offso.r, labels = rightObwLab, cex=cexo, srt=90)
text(x = v2.r.bp[2,]-offsi, y = v2.right[2,]+offso.r, labels = rightWeiLab, cex=cexo, srt=90)
if (save2disk) dev.off()

# =================================================



# Delete one data point which messes up the GLM 10-CV.
v2 <- v2[-which(v2$Left=="Colon"),]
v2$Left <- droplevels(v2$Left)


#### AUX STUDY Fn. 30 ######

cat("\n\n---- Distribution of Hypo ----- \n")

obw.hyp.tab <- table(obwohl$Hypo, obwohl$Senttype)
obw.hyp.csum <- apply(obw.hyp.tab, 2, sum)
obw.hyp.tab.perc <- rbind(round(obw.hyp.tab[1,]/obw.hyp.csum*100, 2), round(obw.hyp.tab[2,]/obw.hyp.csum*100, 2))
obw.hyp.ft <- fisher.test(obw.hyp.tab)

cat("\n obwohl \n")
print(obw.hyp.tab)
print(obw.hyp.tab.perc)
print(obw.hyp.ft)

wei.hyp.tab <- table(weil$Hypo, weil$Senttype)
wei.hyp.csum <- apply(wei.hyp.tab, 2, sum)
wei.hyp.tab.perc <- rbind(round(wei.hyp.tab[1,]/wei.hyp.csum*100, 2), round(wei.hyp.tab[2,]/wei.hyp.csum*100, 2))
wei.hyp.ft <- fisher.test(wei.hyp.tab)

cat("\n weil \n")
print(wei.hyp.tab)
print(wei.hyp.tab.perc)
print(wei.hyp.ft)

################## GLMs ##################

cat("\n\n--------- GLMs ---------\n\n")




# ================= OBWOHL =================

# Remove cases with parentheticals after the particle.
obwohl <- obwohl[which(obwohl$Hypo=="0"),]

obwohl.glm.tmp <- glm(Senttype ~ Right + Independent * ModBin, family=binomial, data=obwohl)
obwohl.glm <- stepAIC(obwohl.glm.tmp, trace=F)

obwohl.glm.0 <- glm(Senttype~1, family=binomial, data=obwohl)
obwohl.lr.test <- lr.test(obwohl.glm, obwohl.glm.0)
obwohl.glm.r2 <- NagelkerkeR2(obwohl.glm)$R2

cat("\n\n--------- OBWOHL ---------\n\n Response: +1 : V2, -1 : VL\n\n")
cat("n=", nrow(obwohl))
cat("\n")
print(summary(obwohl.glm))
cat("\n\nOdds ratios:\n")
print(cbind(exp(coef(obwohl.glm))))

cat("\n\nNagelkerke R² = ", round(obwohl.glm.r2, 3), "\n")
cat("Dispersion  φ = ", round(phi.glm(obwohl.glm), 3), "\n")
cat("LR Test     LR = ", obwohl.lr.test$lr, "   df = ",  obwohl.lr.test$df, "   p = ",  obwohl.lr.test$p, "\n")
cat("10-CV    Delta = ", cv.glm(obwohl, obwohl.glm, K=10)$delta[2], "\n")
cat("\nVariance inflation / colinearity:\n")
print(vif(obwohl.glm))

if (save2disk) svg("../Results/obwohl.svg")
plot.coef(obwohl.glm, main="GLM: Obwohl\n(positive coefficients: V2)")
if (save2disk) dev.off()


if (save2disk) svg("../Results/obwohl-boot.svg")
replicates <- 1e4
cat("\nBootstrap with ", replicates, " replicates\n")

load("../data/obwohl.boot.dat")
# obwohl.boot <- Boot(obwohl.glm, R=replicates)
# save(obwohl.boot, file="../data/obwohl.boot.dat")

print(summary(obwohl.boot))

load("../data/obwohl.boot.ci.dat")
# obwohl.boot.ci <- confint(obwohl.boot)
# save(obwohl.boot.ci, file="../data/obwohl.boot.ci.dat")

obwohl.bo <- order(names(obwohl.boot$t0), decreasing = T)

obw.boot.labs <- rownames(obwohl.boot.ci)[obwohl.bo]
obw.boot.labs <- gsub("ModBin", obw.boot.labs, replacement="Mod")
obw.boot.labs <- gsub("Three", obw.boot.labs, replacement="Ellipsis")
obw.boot.labs <- gsub("Hyphen", obw.boot.labs, replacement="Dash")

dotchart(unlist(summary(obwohl.boot)["bootMed"])[obwohl.bo],
         xlim=c(min(unlist(obwohl.boot.ci[,1])), max(obwohl.boot.ci[,2])),
         pch=18, labels=obw.boot.labs, cex=1)
lapply(1:nrow(obwohl.boot.ci), function(x) lines(c(obwohl.boot.ci[obwohl.bo[x],1], obwohl.boot.ci[obwohl.bo[x],2]), c(x,x), lwd=2, col="black") )
lines(c(0,0), c(0,nrow(obwohl.boot.ci)+1), lwd=2, lty=2, col="black")
points(unlist(summary(obwohl.boot)["bootMed"])[obwohl.bo], 1:nrow(obwohl.boot.ci), pch=20, cex=2)
if (save2disk) dev.off()



# ================= WEIL =================

# Remove cases with parentheticals after the particle.
weil <- weil[which(weil$Hypo=="0"),]

weil.glm.tmp <- glm(Senttype ~ Right + Independent * ModBin, family=binomial, data=weil)
weil.glm <- stepAIC(weil.glm.tmp, trace=F)
weil.glm.0 <- glm(Senttype~1, family=binomial, data=weil)

weil.lr.test <- lr.test(weil.glm, weil.glm.0)
weil.glm.r2 <- NagelkerkeR2(weil.glm)$R2

cat("\n\n--------- WEIL ---------\n\n Response: + = V2, - = VL\n\n")
cat("n=", nrow(weil))
cat("\n")
print(summary(weil.glm))
cat("\n\nOdds ratios:\n")
print(cbind(exp(coef(weil.glm))))

cat("\n\nNagelkerke R²= ", round(weil.glm.r2, 3), "\n")
cat("Dispersion  φ = ", round(phi.glm(weil.glm), 3), "\n")
cat("LR Test     LR = ", weil.lr.test$lr, "   df = ",  weil.lr.test$df, "   p = ",  weil.lr.test$p, "\n")
cat("Variance inflation / colinearity:\n")
print(vif(weil.glm))

# ================= V2 =================

# Remove cases with parentheticals after the particle.
v2 <- v2[which(v2$Hypo=="0"),]

v2.glm.tmp <- glm(Target ~ Right + Left * ModBin, family=binomial, data=v2)
v2.glm <- stepAIC(v2.glm.tmp, trace=F)
v2.glm.0 <- glm(Target~1, family=binomial, data=v2)

v2.lr.test <- lr.test(v2.glm, v2.glm.0)
v2.glm.r2 <- NagelkerkeR2(v2.glm)$R2

cat("\n\n--------- V2 ---------\n\n Response: +1 : Weil, -1 : Obwohl\n\n")
cat("n=", nrow(v2))
cat("\n")
print(summary(v2.glm))
cat("\n\nOdds ratios:\n")
print(cbind(exp(coef(v2.glm))))

cat("\n\nNagelkerke R² = ", round(v2.glm.r2, 3), "\n")
cat("Dispersion  φ = ", round(phi.glm(v2.glm), 3), "\n")
cat("LR Test     LR = ", v2.lr.test$lr, "   df = ",  v2.lr.test$df, "   p = ",  v2.lr.test$p, "\n")

cat("10-CV    Delta = ", cv.glm(v2, v2.glm, K=10)$delta[2], "\n")
cat("Variance inflation / colinearity:\n")
print(vif(v2.glm))

if (save2disk) svg("../Results/weil-obwohl-v2.svg")
plot.coef(v2.glm, main="GLM: V2\n(positive coefficients: weil)")
if (save2disk) dev.off()

if (save2disk) svg("../Results/v2-boot.svg")
replicates <- 1e4
cat("\nBootstrap with ", format(replicates, big.mark=",", scientific=FALSE), " replicates\n")

load("../data/v2.boot.dat")
# v2.boot <- Boot(v2.glm, R=replicates)
# save(v2.boot, file="v2.boot.dat")

print(summary(v2.boot))

load("../data/v2.boot.ci.dat")
# v2.boot.ci <- confint(v2.boot)
# save(v2.boot.ci, file="v2.boot.ci.dat")

v2.bo <- order(names(v2.boot$t0), decreasing = T)

v2.boot.labs <- rownames(v2.boot.ci)[v2.bo]
v2.boot.labs <- gsub("ModBin", v2.boot.labs, replacement="Mod")
v2.boot.labs <- gsub("Three", v2.boot.labs, replacement="Ellipsis")
v2.boot.labs <- gsub("Hyphen", v2.boot.labs, replacement="Dash")

dotchart(unlist(summary(v2.boot)["bootMed"])[v2.bo], xlim=c(min(unlist(v2.boot.ci[,1])), max(v2.boot.ci[,2])), pch=18, labels=v2.boot.labs, cex=1)
lapply(1:nrow(v2.boot.ci), function(x) lines(c(v2.boot.ci[v2.bo[x],1], v2.boot.ci[v2.bo[x], 2]), c(x,x), lwd=2, col="black") )
lines(c(0,0), c(0,nrow(v2.boot.ci)+1), lwd=2, lty=2, col="black")
points(unlist(summary(v2.boot)["bootMed"])[v2.bo], 1:nrow(v2.boot.ci), pch=20, cex=2)
if (save2disk) dev.off()



# ================= VL =================

# Remove cases with parentheticals after the particle.
weilobwohl.vl <- weilobwohl.vl[which(weilobwohl.vl$Hypo=="0"),]

weilobwohl.vl.glm.tmp <- glm(Target ~ Right + Independent * ModBin, family=binomial, data=weilobwohl.vl)
weilobwohl.vl.glm <- stepAIC(weilobwohl.vl.glm.tmp, trace=F)
weilobwohl.vl.glm.0 <- glm(Target~1, family=binomial, data=weilobwohl.vl)

weilobwohl.vl.lr.test <- lr.test(weilobwohl.vl.glm, weilobwohl.vl.glm.0)
weilobwohl.vl.glm.r2 <- NagelkerkeR2(weilobwohl.vl.glm)$R2

cat("\n\n--------- VL ---------\n\n Response: +1 : Weil, -1 : Obwohl\n\n")
cat("n=", nrow(weilobwohl.vl))
cat("\n")
print(summary(weilobwohl.vl.glm))
cat("\n\nOdds ratios:\n")
print(cbind(exp(coef(weilobwohl.vl.glm))))

cat("\n\nNagelkerke R² = ", round(weilobwohl.vl.glm.r2, 3), "\n")
cat("Dispersion  φ = ", round(phi.glm(weilobwohl.vl.glm), 3), "\n")
cat("LR Test     LR = ", weilobwohl.vl.lr.test$lr, "   df = ",  weilobwohl.vl.lr.test$df, "   p = ",  weilobwohl.vl.lr.test$p, "\n")

cat("\nVariance inflation / colinearity:\n")
print(vif(weilobwohl.vl.glm))





# ========= Aux study (Sec 5) ==========

cat("\n\n===============================================\n")
cat("Auxiliary study: total isolation of VL clauses.\n\n")

weil.initial <- weilobwohl.vl[which(weilobwohl.vl$Target=="Weil" & weilobwohl.vl$Majus==1 & (weilobwohl.vl$Left=="End" | weilobwohl.vl$Left=="Emo") & weilobwohl.vl$Isolated!="2" ),]
weil.initial$Isolated <- droplevels(weil.initial$Isolated)
obwohl.initial <- weilobwohl.vl[which(weilobwohl.vl$Target=="Obwohl" & weilobwohl.vl$Majus==1 & (weilobwohl.vl$Left=="End" | weilobwohl.vl$Left=="Emo") & weilobwohl.vl$Isolated!="2" ),]
obwohl.initial$Isolated <- droplevels(obwohl.initial$Isolated)
initial <- rbind(weil.initial,obwohl.initial)
initial.tab <- table(initial$Target, initial$Isolated)
colnames(initial.tab) <- c("Matrix", "Nomatrix")
print(initial.tab)
cat("\n\n Percentages\n")
initial.tab.perc <- rbind(percentize(initial.tab[1,]), percentize(initial.tab[2,]))
rownames(initial.tab.perc) <- rownames(initial.tab)
print(initial.tab.perc)
print(fisher.test(initial.tab))

cat("\n\nSample sizes in aux study:\n\n")

weil.n <- length(which(weilobwohl.vl$Target=="Obwohl"))
obwohl.n <- length(which(weilobwohl.vl$Target=="Weil"))

cat("Obwohl:", nrow(obwohl.initial), "of", obwohl.n, "(=", round(nrow(obwohl.initial)/obwohl.n*100, 2), "%)\n")
cat("weil:", nrow(weil.initial), "of", weil.n, "(=", round(nrow(weil.initial)/weil.n*100, 2), "%)\n")

cat("\n\n===============================================\n")
cat("Remark in Sec. 5: VL with Right=Comma\n\n")

cat("obwohl VL with Right=Comma:", length(which(weilobwohl.vl$Target=="Obwohl" & weilobwohl.vl$Right=="Comma" )), "\n\n")
cat("weil VL with Right=Comma:", length(which(weilobwohl.vl$Target=="Weil" & weilobwohl.vl$Right=="Comma" )), "\n\n")

cat("\n\n===============================================\n")
cat("Aux study in Section 5.2 (top sentence-inital words with PM)\n\n")

inits.colon <- read.delim("../data/12qs.init.colon.csv", header=FALSE, stringsAsFactors=FALSE, quote="")
inits.comma <- read.delim("../data/12qs.init.comma.csv", header=FALSE, stringsAsFactors=FALSE, quote="")
inits.dash <- read.delim("../data/12qs.init.dash.csv", header=FALSE, stringsAsFactors=FALSE, quote="")
inits.ellipsis <- read.delim("../data/12qs.init.ellipsis.csv", header=FALSE, stringsAsFactors=FALSE, quote="")

summary(inits.ellipsis)

inits.colon.n <- sum(inits.colon[,1])
inits.comma.n <- sum(inits.comma[,1])
inits.dash.n <- sum(inits.dash[,1])
inits.ellipsis.n <- sum(inits.ellipsis[,1])
inits.colon <- cbind(inits.colon[,2], round(inits.colon[,1]/inits.colon.n*100, 2), inits.colon[,1])
inits.comma <- cbind(inits.comma[,2], round(inits.comma[,1]/inits.comma.n*100, 2), inits.comma[,1])
inits.dash <- cbind(inits.dash[,2], round(inits.dash[,1]/inits.dash.n*100, 2), inits.dash[,1])
inits.ellipsis <- cbind(inits.ellipsis[,2], round(inits.ellipsis[,1]/inits.ellipsis.n*100, 2), inits.ellipsis[,1])
colnames(inits.colon) <- c("Word", "Percent", "Count")
colnames(inits.comma) <- c("Word", "Percent", "Count")
colnames(inits.dash) <- c("Word", "Percent", "Count")
colnames(inits.ellipsis) <- c("Word", "Percent", "Count")

cat(" === COLON === \n\n")
inits.colon <- inits.colon[which(inits.colon[,1] > 10),]
print(inits.colon[1:20,])
cat("\n Total: ", inits.colon.n, "\n\n")

cat(" === COMMA === \n\n")
inits.comma <- inits.comma[which(inits.comma[,1] > 10),]
print(inits.comma[1:20,])
cat("\n\n Total: ", inits.comma.n, "\n\n")

cat(" === DASH === \n\n")
inits.dash <- inits.dash[which(inits.dash[,1] > 10),]
print(inits.dash[1:20,])
cat("\n\n Total: ", inits.dash.n, "\n\n")

cat(" === ELLIPSIS === \n\n")
inits.ellipsis <- inits.ellipsis[which(inits.ellipsis[,1] > 10),]
print(inits.ellipsis[1:20,])
cat("\n\n Total: ", inits.ellipsis.n, "\n\n")



cat("===============================================\n")
cat("See you then. Be good! Bye-bye!\n")
cat("===============================================\n\n\n")

if (save2disk) sink()
