# --------------------------------------------------------------------- #
# "Perceiving location of moving objects across eye blinks"
# Maus, Goh, & Lisi
# analysis of Experiment 2
# --------------------------------------------------------------------- #
# load libraries

rm(list=ls())
setwd("~/git_local/motion-blink/") # set to local working directory
library(ggplot2)
library(mlisi) # available here: https://github.com/mattelisi/mlisi

library(ggplot2)
# nicer white theme
nice_theme <- theme_bw()+theme(text=element_text(family="Helvetica",size=9),panel.border=element_blank(),strip.background = element_rect(fill="white",color="white",size=0),strip.text=element_text(size=rel(0.8)),panel.grid.major.x=element_blank(),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),axis.line.x=element_line(size=.4),axis.line.y=element_line(size=.4),axis.text.x=element_text(size=7,color="black"),axis.text.y=element_text(size=7,color="black"),axis.line=element_line(size=.4), axis.ticks=element_line(color="black"))

# others
library(lme4)
library(BayesFactor)
library(sjstats)

# --------------------------------------------------------------------- #
# Load dataset and prepare for analysis

# load data
d <- read.table("./data/Experiment2.csv",header=T,sep=",")

# eclude trials with inappropriate blinking
d$condition <- ifelse(d$Block==1|d$Block==3,"blink","no-blink")
d$include <- ifelse((d$condition=="no-blink" & d$Blink==1)|(d$condition=="blink" & d$Blink==0), 0, 1)
round((1-mean(d$include))*100, digits=2) # 2.310268
d <- d[d$include==1,]

# cut-offs as Gerrit's analysis
outfilter_blink <- function (rtv, sv, bv, nsd = 2) {
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
round(sum(out_1)/nrow(d)*100,digits=2)
d <- d[-out_1,]

# additional filter 
round(sum(d$BlinkDuration<15 & d$condition=="blink")/nrow(d) * 100,digits=2) # 0.1485205
d <- d[which((d$BlinkDuration>15 & d$condition=="blink") | d$condition=="no-blink"),]

# center and normalize blink distributions
for(i in unique(d$Subject)){
  d$BlinkDuration[d$Subject==i & d$condition=="blink"] <- scale(d$BlinkDuration[d$Subject==i & d$condition=="blink"],scale=F,center=T)  
}

# prepare variables for modelling
d$bdur <- d$BlinkDuration/1000 # transform in seconds
d$Velocity <- factor(d$Velocity)
d$condition <- factor(d$condition)
contrasts(d$Velocity)
contrasts(d$condition)
d$cond1 <- ifelse(d$condition=="blink",1,0)

# --------------------------------------------------------------------- #
# evaluate goodness of fit of psychometric function using deviance test
d_fit <-{}
for(i in unique(d$Subject)){
  m_1 <- glm(Response ~ JumpSize * Velocity +BlinkDuration , d[d$Subject==i & d$cond1==1,], family=binomial(probit))
  m_0<- glm(Response ~ Velocity +BlinkDuration , d[d$Subject==i & d$cond1==1,], family=binomial(probit))
  LRT <- anova(m_0, m_1,test="LRT")
  d_fit <- rbind(d_fit, data.frame(id=i, D=LRT$Deviance[2], df=LRT$Df[2], p=LRT$`Pr(>Chi)`[2]))
}
print(d_fit,digits=2,row.names=F)
round(d_fit$p,digits=10)
# this was added to address a reviewer comments. It confirms that the 
# proportion of "forward" responses was modulated by the jump size for '
# all observers
d_fit$p < 10^(-9)

# --------------------------------------------------------------------- #
# use lmList to fit individual GLM models in one go
mo <- lmList(Response ~ cond1 * JumpSize * Velocity +bdur:cond1 | Subject, d, family=binomial(probit))

# store coefficients
B <- coef(mo)

# --------------------------------------------------------------------- #
# test the dependence in blink duration at group level
# slope_dur is the expecte rate of change in PSE due to 
# the deviation of blink duration from the individual mean
slope_bdur <- B[,8] / (((B[,3]+B[,5]) + (B[,3]+B[,7]+B[,5]+B[,9]))/2)
t.test(slope_bdur)
1/ttestBF(x = slope_bdur, rscale=1)

# --------------------------------------------------------------------- #
# transform the beta coefficients in PSE, then do ANOVA
# -beta_0/beta_1
nobli180 <- -B[,1]/B[,3]
nobli240 <- -(B[,1]+B[,4])/(B[,3]+B[,7])
blink180 <- -(B[,1]+B[,2])/(B[,3]+B[,5])
blink240 <- -(B[,1]+B[,4]+B[,2]+B[,6])/(B[,3]+B[,7]+B[,5]+B[,9])
dPSE <- data.frame(PSE=c(nobli180,nobli240,blink180,blink240), condition=c(rep("no-blink",16),rep("blink",16)),Velocity=as.factor(rep(c(rep(180,8),rep(240,8)),2)), id=as.factor(rep(1:8,4)))
m0 <- aov(PSE~condition*Velocity + Error(id/(condition*Velocity)),dPSE)

## ANOVA results
summary(m0)

##  ANOVA results with standardized effect size (partial.etasq)
anova_stats(m0)[,c(1:3,6:7,9)]

# --------------------------------------------------------------------- #
# calculate individual BF for the effect of duration
# since there isn't an established way to calculate "default"
# bayes factor, we fit a fully Bayesian model in Stan, then
# use Savage-Dickey density ratio

# load custom functions
source("BayesianGLM_functions.R")

# This code runs it iteratively for all participants
# see main text for explanation of how we set the prior for the 
# trial-by-trial effect of blink duration
# out_all <- {}
# for(i in rownames(coef(mo))){
#   m_sig_i <- glm(Response ~ JumpSize, d[d$Subject==i & d$condition=="blink",], family=binomial(probit))
#   #prior_sd <- unname(mean(c(180,240)) / (1/coef(m_sig_i)[2]))
#   prior_sd <- unname(240 / (1/coef(m_sig_i)[2]))
#   out_i <- bayesian_GLM_varSD(d[d$Subject==i,], prior_sd )
#   out_i$prior_sd <- prior_sd 
#   out_all <-rbind(out_all, out_i)
# }
# saveRDS(out_all,file="out_all_exp2_varSD.RDS")

# this load the results (if already run)
outB <- readRDS("out_all_exp2_varSD.RDS")

# visualize Bayes factors and parameter estimates
par(mfrow=c(1,2))

slope_bdur <- B[,8] / (((B[,3]+B[,5]) + (B[,3]+B[,7]+B[,5]+B[,9]))/2)
slope_bayes <- outB$beta_5 / (((outB$beta_2) + (outB$beta_2+outB$beta_4))/2)

# sanity check to ensure that we get same estimates from frequentist
# or Bayesian analysis
plot(slope_bdur, slope_bayes, ylab="slope (Bayesian)",xlab="slope (frequentist)", pch=19)
abline(a=0,b=1,lty=2)
cor.test(slope_bdur, slope_bayes)

# plot BF_01 as a function of parameter estimate
plot( slope_bayes, outB$BF01,  cex=1.2, col="blue",xlab="slope [deg/sec]", ylab=expression("BF"["01"])) 
abline(h=c(10^(1/2), 1/(10^(1/2))),lty=2); 
abline(h=1,lty=1)


## summaries

# median and range of BF supporting the null
length(outB$BF01[outB$BF01 > 10^(1/2)])
round(median(outB$BF01[outB$BF01 > 10^(1/2)]),digits=2)
round(range(outB$BF01[outB$BF01> 10^(1/2)]),digits=2)

# median and range of BF not-supporting the null (inconclusive)
length(outB$BF01[1/(10^(1/2)) < outB$BF01  & outB$BF01 < 10^(1/2)])
round(median(outB$BF01[1/(10^(1/2)) < outB$BF01  & outB$BF01 < 10^(1/2)]),digits=2)
round(range(outB$BF01[1/(10^(1/2)) < outB$BF01 & outB$BF01 < 10^(1/2)]),digits=2)

# median and range of BF not-supporting the alternative
length(outB$BF01[outB$BF01 < 1/(10^(1/2))])
outB$BF01[outB$BF01 < 1/(10^(1/2))]

# --------------------------------------------------------------------- #
# ANOVA on psychometric slopes
nobli180 <- 1/B[,3]
nobli240 <- 1/(B[,3]+B[,7])
blink180 <- 1/(B[,3]+B[,5])
blink240 <- 1/(B[,3]+B[,7]+B[,5]+B[,9])
dSIGMA <- data.frame(sigma=c(nobli180,nobli240,blink180,blink240), condition=c(rep("no-blink",16),rep("blink",16)),Velocity=as.factor(rep(c(rep(180,8),rep(240,8)),2)), id=as.factor(rep(1:8,4)))
m0 <- aov(sigma~condition*Velocity + Error(id/(condition*Velocity)),dSIGMA)

# results
summary(m0)

# with standardized effect sizes
anova_stats(m0)[,c(1:3,6:7,9)]

# --------------------------------------------------------------------- #
# plot of psychometric functions
mo <- lmList(Response ~ cond1 * JumpSize * Velocity +BlinkDuration:cond1 | Subject, d, family=binomial(probit))
nd <- expand.grid(JumpSize=seq(-65, 65, length.out=100), Velocity=unique(d$Velocity), cond1=unique(d$cond1), Subject=unique(d$Subject), BlinkDuration=0)
nd$Response <- pnorm(predict(mo, newdata=nd, type="response"))

dag0 <- aggregate(Response ~ JumpSize + Velocity + cond1 + Subject, d, mean)
dag1 <- aggregate(Response ~ JumpSize + Velocity + cond1, dag0, mean)
dag1$se <- aggregate(Response ~ JumpSize + Velocity + cond1, dag0, bootMeanSE)$Response

nd$plot_g <- paste(nd$Subject,nd$cond1,sep="_")
dag1$plot_g <- NA
nd$cond1 <- ifelse(nd$cond1==1,"blink","no-blink")
dag1$cond1 <- ifelse(dag1$cond1==1,"blink","no-blink")

# make plot
ggplot(nd, aes(x=JumpSize, y=Response, color=cond1, group=plot_g))+geom_line(alpha=0.5)+facet_grid(cond1~Velocity) + scale_color_manual(values=c("black","dark grey"),guide=F)+nice_theme+geom_errorbar(data=dag1, aes(ymin=Response-se, ymax=Response+se),width=3,size=1)+geom_point(data=dag1, size=3)+labs(x="jump size [deg]",y="p(forward)")



# --------------------------------------------------------------------- #
# plot of individual PSE values
nobli180 <- -B[,1]/B[,3]
nobli240 <- -(B[,1]+B[,4])/(B[,3]+B[,7])
blink180 <- -(B[,1]+B[,2])/(B[,3]+B[,5])
blink240 <- -(B[,1]+B[,4]+B[,2]+B[,6])/(B[,3]+B[,7]+B[,5]+B[,9])
dPSE <- data.frame(PSE=c(nobli180,nobli240,blink180,blink240), condition=c(rep("no-blink",16),rep("blink",16)),Velocity=as.factor(rep(c(rep(180,8),rep(240,8)),2)), id=as.factor(rep(1:8,4)))
str(dPSE)

# add the standard errors
Bse <- summary(mo)$coefficients[,2,]
se_nobli180 <- sqrt((Bse[,1]/B[,1])^2 + (Bse[,3]/B[,3])^2)
se_nobli240 <- sqrt((sqrt(Bse[,1]^2 + Bse[,4]^2)/(B[,1]+B[,4]))^2 +(sqrt(Bse[,3]^2 + Bse[,7]^2)/(B[,3]+B[,7]))^2)
se_blink180 <- sqrt((sqrt(Bse[,1]^2 + Bse[,2]^2)/(B[,1]+B[,2]))^2 +(sqrt(Bse[,3]^2 + Bse[,5]^2)/(B[,3]+B[,5]))^2)
se_blink240 <- sqrt((sqrt(Bse[,1]^2 + Bse[,2]^2 + Bse[,4]^2 + Bse[,6]^2)/(B[,1]+B[,2]+B[,4]+B[,6]))^2 +(sqrt(Bse[,3]^2 +Bse[,5]^2+Bse[,7]^2+Bse[,9]^2)/(B[,3]+B[,5]+B[,7]+B[,9]))^2)
dPSE$se <- c(se_nobli180,se_nobli240,se_blink180,se_blink240)
dPSE$blink_label <- ifelse(dPSE$condition=="blink","blink trials","no-blink trials")
dPSE$vel_label <- ifelse(dPSE$Velocity==240,"240 deg/sec","180 deg/sec")

# mean PSE
dPSE_ag <- aggregate(PSE ~ condition + Velocity + blink_label + vel_label, dPSE, mean)
dPSE_ag$se <- aggregate(PSE ~ condition + Velocity + blink_label + vel_label, dPSE, bootMeanSE)$PSE

# set the data right for plotting
library(tidyr)
dPSE_ag$condition <- {}
dPSE$condition <- {}
dag_x1 <- pivot_wider(data=dPSE, names_from=blink_label, values_from=c(PSE, se))
dag_x2 <- pivot_wider(data=dPSE_ag, names_from=blink_label, values_from=c(PSE, se))
dag_x1$x_lb <- dag_x1$"PSE_blink trials" - dag_x1$"se_blink trials" 
dag_x1$x_ub <- dag_x1$"PSE_blink trials" + dag_x1$"se_blink trials" 
dag_x1$y_lb <- dag_x1$"PSE_no-blink trials" - dag_x1$"se_no-blink trials" 
dag_x1$y_ub <- dag_x1$"PSE_no-blink trials" + dag_x1$"se_no-blink trials" 
dag_x2$x_lb <- dag_x2$"PSE_blink trials" - dag_x2$"se_blink trials" 
dag_x2$x_ub <- dag_x2$"PSE_blink trials" + dag_x2$"se_blink trials" 
dag_x2$y_lb <- dag_x2$"PSE_no-blink trials" - dag_x2$"se_no-blink trials" 
dag_x2$y_ub <- dag_x2$"PSE_no-blink trials" + dag_x2$"se_no-blink trials" 

# OK finally here's the plot
(ggplot(dag_x1, aes(x=get("PSE_blink trials"), y=get("PSE_no-blink trials"), color=vel_label, fill= vel_label, xmin=x_lb, xmax=x_ub, ymin=y_lb, ymax=y_ub))
  +facet_grid(.~vel_label)
  +nice_theme
  +scale_color_manual(values=c("dark grey","black"),name="", guide=F)
  +scale_fill_manual(values=c("dark grey", "black"),name="", guide=F)
  +geom_abline(slope=1,intercept=0,size=0.3,lty=2)
  +geom_vline(xintercept=0,size=0.25,lty=2)
  +geom_hline(yintercept=0,size=0.3,lty=2)
  +geom_errorbar(width=0,size=0.2)
  +geom_errorbarh(height=0,size=0.2)
  +geom_point(pch=21,fill=NA,stroke=0.3)
  +theme(panel.spacing = unit(2, "lines"))
  +geom_errorbar(data=dag_x2,size=1,width=0)
  +geom_errorbarh(data=dag_x2,size=0.6,height=0)
  +labs(x="response error in blink trials [deg]", y="response error in no-blink trials [deg]")
  +geom_point(data=dag_x2,size=3,pch=19)
  +coord_equal(xlim=c(-65,10),ylim=c(-65,10)))



