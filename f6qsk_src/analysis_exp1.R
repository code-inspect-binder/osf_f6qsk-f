# --------------------------------------------------------------------- #
# "Perceiving location of moving objects across eye blinks"
# Maus, Goh, & Lisi
# analysis of Experiment 1
# --------------------------------------------------------------------- #
# load libraries

rm(list=ls())
setwd("~/git_local/motion-blink/")
library(ggplot2)
library(mlisi) # available here: https://github.com/mattelisi/mlisi

library(ggplot2)
# nicer white theme
nice_theme <- theme_bw()+theme(text=element_text(family="Helvetica",size=9),panel.border=element_blank(),strip.background = element_rect(fill="white",color="white",size=0),strip.text=element_text(size=rel(0.8)),panel.grid.major.x=element_blank(),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),axis.line.x=element_line(size=.4),axis.line.y=element_line(size=.4),axis.text.x=element_text(size=7,color="black"),axis.text.y=element_text(size=7,color="black"),axis.line=element_line(size=.4), axis.ticks=element_line(color="black"))

# --------------------------------------------------------------------- #
# Load dataset and prepare for analysis

# load
d <- read.table("./data/Experiment1.csv",header=T,sep=",")

# eclude trials with 'inappropriate' blinking
# i.e. trials with blink in no-blink conditions, and trials
# without blinks in blink conditions
d$condition <- ifelse(d$Block==1|d$Block==3,"blink","no-blink")
d$include <- ifelse((d$condition=="no-blink" & d$Blink==1)|(d$condition=="blink" & d$Blink==0), 0, 1)
(1-mean(d$include))*100 # 4.082031
d <- d[d$include==1,]

# trials with blink dur = 0 in blink blocks
(sum(d$BlinkDuration==0 & d$condition=="blink") / nrow(d)) * 100 # 9.570352
d <- d[d$condition=="no-blink" | d$BlinkDuration>0,]

# custom functions to do cut-off based on >3SD from mean (response error)
outfilter_blink <- function (rtv, sv, bv, nsd = 3) {
  if (length(rtv) == length(sv) & length(sv) == length(bv)) {
    index <- {}
    subjects <- unique(sv)
    for (i in 1:length(subjects)) {
      aRT <- rtv[which(sv == subjects[i])]
      c_i <- bv[which(sv == subjects[i])]
      aM <- mean(aRT[c_i=="blink"])
      aSD <- sd(aRT[c_i=="blink"])
      index_i <- ifelse(c_i=="no-blink",0, aRT > (aM + nsd * aSD) | aRT < (aM - nsd * aSD))
      index <- c(index, index_i)
    }
    invisible(index)
  }
  else {
    stop("Error: the two vectors have different length!")
  }
}
out_1 <- outfilter_blink(d$BlinkDuration, d$Subject, d$condition, 3)
sum(out_1)/nrow(d)*100 # 0.2702094
d <- d[-out_1,]

# filter outliers in response error
out_2 <- outfilter(d$ResponseError, d$Subject, 3)
length(out_2)/nrow(d) * 100
d <- d[-out_2,]

# there are some unusually short blink durations; we removed them
d <- d[which((d$BlinkDuration>15 & d$condition=="blink") | d$condition=="no-blink"),]

# store table with mean blink durations (to plot on real scale later on)
tbdur <- with(d[d$condition=="blink",],tapply(BlinkDuration,Subject,mean))

# center blink distributions so that their mean is zero (for modelling below)
for(i in unique(d$Subject)){
  d$BlinkDuration[d$Subject==i & d$condition=="blink"] <- scale(d$BlinkDuration[d$Subject==i & d$condition=="blink"],scale=F,center=T)  
}
# sd_centered_bdur <- sd(d$BlinkDuration[d$condition=="blink"])

# --------------------------------------------------------------------- #
### Make a nice plot of response error (FIG. 1A)

d$blink_label <- ifelse(d$condition=="blink","blink trials","no-blink trials")
d$vel_label <- ifelse(d$Velocity==240,"240 deg/sec","180 deg/sec")
dag1 <- aggregate(ResponseError ~ blink_label + Subject, d, mean)
dag2 <- aggregate(ResponseError ~ blink_label, dag1, mean)
dag2$se <- aggregate(ResponseError ~ blink_label, dag1, bootMeanSE)$ResponseError

dag1 <- aggregate(ResponseError ~ blink_label + Subject + vel_label, d, mean)
dag1$se <- aggregate(ResponseError ~ blink_label + Subject + vel_label, d, bootMeanSE)$ResponseError
dag2 <- aggregate(ResponseError ~ blink_label + vel_label, dag1, mean)
dag2$se <- aggregate(ResponseError ~ blink_label + vel_label, dag1, bootMeanSE)$ResponseError

# simple barplot
ggplot(dag2, aes(x=vel_label, y=ResponseError, fill=vel_label, color=vel_label))+geom_col(alpha=0.5,width=0.4)+geom_errorbar(aes(ymin=ResponseError-se, ymax=ResponseError+se),width=0.2)+facet_grid(.~blink_label)+nice_theme+geom_hline(yintercept=0,size=0.4,lty=2)+scale_color_manual(values=c("black","dark grey"),name="")+scale_fill_manual(values=c("black", "dark grey"),name="")+ theme(panel.spacing = unit(2, "lines"))

# make a nicer scatterplot with individual values
library(tidyr)
dag_x1 <- pivot_wider(data=dag1, names_from=blink_label, values_from=c(ResponseError, se))
dag_x2 <- pivot_wider(data=dag2, names_from=blink_label, values_from=c(ResponseError, se))
dag_x1$x_lb <- dag_x1$"ResponseError_blink trials" - dag_x1$"se_blink trials" 
dag_x1$x_ub <- dag_x1$"ResponseError_blink trials" + dag_x1$"se_blink trials" 
dag_x1$y_lb <- dag_x1$"ResponseError_no-blink trials" - dag_x1$"se_no-blink trials" 
dag_x1$y_ub <- dag_x1$"ResponseError_no-blink trials" + dag_x1$"se_no-blink trials" 
dag_x2$x_lb <- dag_x2$"ResponseError_blink trials" - dag_x2$"se_blink trials" 
dag_x2$x_ub <- dag_x2$"ResponseError_blink trials" + dag_x2$"se_blink trials" 
dag_x2$y_lb <- dag_x2$"ResponseError_no-blink trials" - dag_x2$"se_no-blink trials" 
dag_x2$y_ub <- dag_x2$"ResponseError_no-blink trials" + dag_x2$"se_no-blink trials" 

p1 <- ggplot(dag_x1, aes(x=get("ResponseError_blink trials"), y=get("ResponseError_no-blink trials"), color=vel_label, fill= vel_label, xmin=x_lb, xmax=x_ub, ymin=y_lb, ymax=y_ub))+facet_grid(.~vel_label)+coord_equal(xlim=c(-7.5,20),ylim=c(-7.5,20))+nice_theme+scale_color_manual(values=c("dark grey","black"),name="", guide=F)+scale_fill_manual(values=c("dark grey", "black"),name="", guide=F)+geom_abline(slope=1,intercept=0,size=0.3,lty=2)+geom_vline(xintercept=0,size=0.25,lty=2)+geom_hline(yintercept=0,size=0.3,lty=2)+geom_errorbar(width=0,size=0.2)+geom_errorbarh(height=0,size=0.2)+geom_point(pch=21,fill=NA,stroke=0.3)+ theme(panel.spacing = unit(2, "lines"))+geom_errorbar(data=dag_x2,size=1,width=0)+geom_errorbarh(data=dag_x2,size=0.6,height=0)+labs(x="response error in blink trials [deg]", y="response error in no-blink trials [deg]")+geom_point(data=dag_x2,size=3,pch=19)

pdf("exp1_i.pdf",width=5,height=3)
p1
dev.off()


# --------------------------------------------------------------------- #
# fit within subject model
library(lme4)
library(mlisi)
library(BayesFactor)

#
d$bdur <- d$BlinkDuration/1000 # transform in seconds
d$vel1 <- ifelse(d$Velocity==240,1,0) # dummy coding of velocity
d$cond1 <- ifelse(d$condition=="blink",1,0) # dummy coding of condition

# fit model
m0 <- lmList(ResponseError~cond1*vel1 + cond1:bdur| Subject, d)
coef(m0)

# test individual coefficients
t.test(coef(m0)[,2])
t.test(coef(m0)[,3])
t.test(coef(m0)[,4])
t.test(coef(m0)[,5])

# bayes factor for trial-by-trial influence of blink duration
bf = ttestBF(x = coef(m0)[,5], rscale=1)
1/bf

# --------------------------------------------------------------------- #
# This make a table of coefficients with standardized measures of effect size

ktab <- coef(m0)
colnames(ktab) <- c("beta_0","beta_1","beta_2","beta_3","beta_4")

# formula from Hedges and Olkin (1985)
# Statistical Methods in Meta-Analysis
# January 1985Journal of Educational Statistics 20(1)
# DOI: 10.2307/1164953
cohen.d.onesample <- function(x){
  mean(x)/sd(x)
}
cohen.d.ci <- function(x, alpha=0.05){
  d <- cohen.d.onesample(x)
  n <- length(x)
  d.se <- sqrt(1/n + (d^2)/(2*n))
  d.ci <- qnorm(1-alpha/2) * d.se
  return(list(estimate=d,conf.int =c(d-d.ci, d+d.ci), error=d.se))
}

tab2 <- {}
for(i in 1:ncol(ktab)){
  t_ <- t.test(ktab[,i])
  efsz <- cohen.d.ci(ktab[,i])
  bf_ = ttestBF(x = ktab[,i], rscale=1)
  line_ <- data.frame(colnames(ktab)[i], mean(ktab[,i]), sd(ktab[,i]), 
                      t_$statistic, t_$parameter, t_$p.value,
                      extractBF(bf_)$bf, efsz$estimate, efsz$conf.int[1], efsz$conf.int[2])
  
  colnames(line_) <-  c("parameter","mean","sd","t","df","p-value","BF_10","Cohen's d", "Cohen's d lower bound" , "Cohen's d upper bound")
  tab2 <- rbind(tab2,line_)
}

# print the table
print(tab2, digits=2)

# --------------------------------------------------------------------- #
# This additional analysis test the influence of trial-by-trial blink duration
# on response error by computing individual Bayes factors
#
# the Bayes factor is calculated by examining the residuals of a null model
# and how they correlate with blink-durations
# the procedure use a 'default' prior for calculating bayes factors for correlations
# code here:
source("BF_correlations.R")

bf_null <- rep(NA, length(unique(d$Subject)))
dur_slope <- rep(NA, length(unique(d$Subject)))
for(i in 1:length(bf_null)){
  
  d_i <- d[d$Subject==unique(d$Subject)[i] & d$cond1==1,]
  m0 <- lm(ResponseError~vel1, d_i)
  m1 <- lm(ResponseError~vel1 + cond1:bdur, d_i)
  dur_slope[i] <- coef(m1)[3]
  
  bf_null[i] <- 1/bf10JeffreysIntegrate(n=nrow(d_i), r=cor(residuals(m0), d_i$bdur))  
}

# mean BF_01
mean(bf_null)
bootMeanCI(bf_null)

sum(bf_null > 10^(1/2)) # n subject for showing strong support for the null
length(bf_null) # out of 

# median and range of BF supporting the null
round(median(bf_null[bf_null > 10^(1/2)]),digits=2)
round(range(bf_null[bf_null > 10^(1/2)]),digits=2)

# median and range of BF not-supporting the null (inconclusive)
round(median(bf_null[bf_null < 10^(1/2)]),digits=2)
round(range(bf_null[bf_null < 10^(1/2)]),digits=2)

# plot of individual BF (not included in paper)
plot( dur_slope, bf_null, log="y", ylim=c(1/(10^(1/1.8)), 13), cex=1.2,col="blue",xlab="slope [deg/sec]", ylab=expression("BF"["01"]), main="Exp. 1") 
abline(h=c(10^(1/2), 1/(10^(1/2))),lty=2); 
abline(h=1,lty=1)

# ------------------------------------------------------------------------- #  
# additional analysis with Direction (CW vs CCW) as predictor
d$dirCW <- ifelse(d$Direction==1,1,0)
m0 <- lmList(ResponseError~cond1*vel1+dirCW + cond1:bdur | Subject, d)
coef(m0)
t.test(coef(m0)[,4])
t.test(coef(m0)[,6])

# ------------------------------------------------------------------------- #  
# plot split by quartiles of blink duration (FIG. 1B)

# calculate individual quantiles
d$blinkBin <- NA
for(i in unique(d$Subject)){
  d$blinkBin[d$Subject==i& d$condition=="blink"] <- cut(d$BlinkDuration[d$Subject==i& d$condition=="blink"],breaks=quantile(d$BlinkDuration[d$Subject==i & d$condition=="blink"],probs = seq(0, 1, 0.25)),labels=1:4)
}
d$blinkBin <- ifelse(d$condition=="no-blink", "no-blink", d$blinkBin)

# make aggregated dataset with pooled standard errors
dag <- aggregate(cbind(BlinkDuration,ResponseError)~condition+Velocity+blinkBin+Subject, d, mean)
dag$n <- aggregate(ResponseError~condition+Velocity+blinkBin+Subject, d, length)$ResponseError

# given that blink duration was centered, add again the individual mean to get the correct duration
dag$bd <- NA
for(i in 1:nrow(dag)){
  dag$bd[i] <- ifelse(dag$condition[i]=="blink", dag$BlinkDuration[i] + tbdur[names(tbdur)==dag$Subject[i]], 0)
}
dag2 <- aggregate(cbind(BlinkDuration,ResponseError,bd)~condition+Velocity+blinkBin, dag, mean)

# To account for some variability in the number of good trials between subjects
# I am using pooled mean and pooled standard error (e.g. weighted by the relative number
# of trials available for each subject in each bin). This should represent more sensibly the trend in
# the data because subject with less trials would have more variable measures are weighted slightly less.
# Note however that the resuls is basically identical with what obtained from standard mean/se,
# which can be run by running this instead:
# dag2$r_se <- aggregate(ResponseError~condition+Velocity+blinkBin, dag, mean)$ResponseError

# this function calculated a bootstrapped standard error for the pooled mean
bootPoolMN <- function (M, N, nsim = 1000) {
  d <- data.frame(M, N)
  poolFoo <- function(d){
    M_pool <- sum(d$M * (d$N)) / sum(d$N)
    return(M_pool)
  }
  bootFoo <- function(d, i) poolFoo(d[i,])
  bootRes <- boot::boot(d, bootFoo, nsim)
  return(sd(bootRes$t, na.rm = T))
}

# calculated pooled mean and SEM
dag2$b_se <- NA
dag2$r_se <- NA
for(i in 1:nrow(dag2)){
  index <- dag$condition==dag2$condition[i] & dag$Velocity==dag2$Velocity[i] & dag$blinkBin==dag2$blinkBin[i]

  M <- dag$ResponseError[index]
  N <- dag$n[index]
  dag2$ResponseError[i] <- sum(M * N) / sum(N)
  dag2$r_se[i] <- bootPoolMN(M,N)

  M <- dag$bd[index]
  dag2$bd[i] <- sum(M * N) / sum(N)
  dag2$bd_se[i] <- bootPoolMN(M,N)
}

# adjust labels for the plot
dag2$Velocity <- paste(dag2$Velocity,"deg/sec")
dag2$group2 <- paste(dag2$condition,dag2$Velocity,sep="_")

dag$Velocity <- paste(dag$Velocity,"deg/sec")
dag$group2 <- paste(dag$condition,dag$Velocity,sep="_")
dag$group3 <- paste(dag$Subject, dag$group2, sep="_")

# do plot
p2 <-   ggplot(dag2, aes(x=bd,y=ResponseError,color=Velocity,group=group2,shape=condition))+geom_hline(yintercept=0,lty=2,size=0.4)+geom_errorbar(aes(ymin=ResponseError-r_se, ymax=ResponseError+r_se),width=0)+geom_point(size=3)+geom_line()+nice_theme+facet_grid(.~Velocity)+scale_color_manual(values=c("dark grey", "black"),guide=F)+labs(y="response error [deg]", x="blink duration [ms]")+scale_x_continuous(limits=c(-20,260))+scale_y_continuous(limits=c(-5,15))+scale_shape_manual(values=c(19,19),guide=F)+ theme(panel.spacing = unit(2, "lines"))

# save in  pdf
pdf("exp1",width=5,height=3)
p2
dev.off()

# make one sigle big plot as in the paper
library(ggpubr)
pdf("exp1_all",width=5,height=5.5)
ggarrange(p1, p2, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
dev.off()

