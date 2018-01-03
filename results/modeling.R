library(tidyverse)
library(BayesFactor)
library(cowplot)
setwd("~/github/physNumComparisonTask/results/")
source("waldFunctions.R")

rawdata <- read_csv("processed.csv")

data = rawdata %>%
  filter(correct==1) %>%
  mutate(rt=response_time/1000) %>%
  filter(rt>median(rt)-3*mad(rt) & rt < median(rt)+6*mad(rt))

data %>%
  ggplot(aes(x=rt, group=congruity)) +
  geom_density(aes(fill=congruity), alpha=0.5) +
  theme_classic(14)

# size-congruity effect on RTs
detach(package:plyr)
agg = data %>%
  group_by(subject_nr,congruity) %>%
  summarize(medRT = median(rt), meanRT=mean(rt), sdRT = sd(rt))



t.test(agg$medRT[agg$congruity=="incongruent"], agg$medRT[agg$congruity=="congruent"], paired=TRUE)
t.test(agg$meanRT[agg$congruity=="incongruent"], agg$meanRT[agg$congruity=="congruent"], paired=TRUE)
t.test(agg$sdRT[agg$congruity=="incongruent"], agg$sdRT[agg$congruity=="congruent"], paired=TRUE)




# fit Wald model for overall RT
fit = swapply(dat=data, obsvar="rt", facs = c("congruity", "subject_nr")) # fit the data

# plot diagnostics and parameters
library(plyr)
swplotws(fit$vars,c("subject_nr"))
swchecks(fit$vars, data, xqnts = fit$xqnts, pxqnts = fit$pxqnts)

params = fit$vars
params$subject_nr <- as.factor(params$subject_nr)
params$congruity=factor(params$congruity, labels=c("congruent","incongruent"))

detach(package:plyr)
gammaPlot = params %>%
  group_by(congruity) %>%
  summarize(mean = mean(gamma), se = sd(gamma)/sqrt(23)) %>%
  ggplot(aes(x=congruity, y=mean)) +
  geom_bar(stat="identity", fill="gray", colour="black", width=0.5) +
  geom_errorbar(width=0.1,aes(ymin=mean-se,ymax=mean+se)) +
  theme_classic(14) +
  ylab(expression(paste("Drift rate ",gamma)))
  
  
alphaPlot = params %>%
  group_by(congruity) %>%
  summarize(mean = mean(alpha), se = sd(alpha)/sqrt(23)) %>%
  ggplot(aes(x=congruity, y=mean)) +
  geom_bar(stat="identity", fill="gray", colour="black", width=0.5) +
  geom_errorbar(width=0.1,aes(ymin=mean-se,ymax=mean+se)) +
  theme_classic(14) +
  ylab(expression(paste("Response threshold ",alpha)))
  

thetaPlot = params %>%
  group_by(congruity) %>%
  summarize(mean = mean(theta), se = sd(theta)/sqrt(23)) %>%
  ggplot(aes(x=congruity, y=mean)) +
  geom_bar(stat="identity", fill="gray", colour="black", width=0.5) +
  geom_errorbar(width=0.1,aes(ymin=mean-se,ymax=mean+se)) +
  theme_classic(14) +
  ylab(expression(paste("Nondecision time ",theta)))
  

plot_grid(gammaPlot,alphaPlot,thetaPlot, nrow=1, ncol=3, labels="AUTO")

gammaCong = params$gamma[params$congruity=="congruent"]
gammaIncong = params$gamma[params$congruity=="incongruent"]
t.test(gammaCong,gammaIncong, paired=TRUE)
gammaBF=ttestBF(gammaCong,gammaIncong, paired=TRUE)
gammaBF

alphaCong = params$alpha[params$congruity=="congruent"]
alphaIncong = params$alpha[params$congruity=="incongruent"]
t.test(alphaCong,alphaIncong, paired=TRUE)
alphaBF=ttestBF(alphaCong,alphaIncong, paired=TRUE)
alphaBF

thetaCong = params$theta[params$congruity=="congruent"]
thetaIncong = params$theta[params$congruity=="incongruent"]
t.test(thetaCong,thetaIncong, paired=TRUE)
thetaBF=ttestBF(thetaCong,thetaIncong, paired=TRUE)
thetaBF
1/thetaBF
