library(MCMCpack)
library(BayesFactor)
library(msm)

setwd("~/Dropbox/papers/2022/lognormalConstraint/sims")
source("fitFunctions.R")

# parent distribution means and sds (from F, Vick, Bowman, 2018)
# note -- RT scale is seconds, not milliseconds
G = 3.91 # gamma
G.sd = 0.70
A = 0.92 # alpha
A.sd = 0.17
H = 0.32 # theta
H.sd = 0.05

runSim = function(nSub, nTrials, correctModel="unconstrained", numSims=10){
  
  # build empty vector to store BFs
  bf = numeric(1)
  bfLog = numeric(1)
  
  if(correctModel=="positive"){
    filename = sprintf("P-%d-%d.csv", nSub, nTrials)
  }
  else{
    filename = sprintf("U-%d-%d.csv", nSub, nTrials)
  }
  if(file.exists(filename)) dat = read.csv(filename) else dat = data.frame()
  
  for (j in 1:numSims){
    # each simulation run starts HERE
    # build simulated RT matrix
    # rows = subjects, columns = 2*trials
    rts = matrix(0, nrow=nSub, ncol=2*nTrials)
    
    # build matrix to store "target" parameters
    # rows = subjects, columns = gamma, alpha, theta
    targets1 = matrix(0, nrow=nSub, ncol=3)
    targets2 = matrix(0, nrow=nSub, ncol=3)
    
    # Step 1 - randomly draw target values from parent distribution for
    # congruent trials
    # (one unique parameter value for each subject)
    for (i in 1:nSub){
      targets1[i,1] = rnorm(1, mean=G, sd=G.sd) # draw random gamma
      targets1[i,2] = rnorm(1, mean=A, sd=A.sd) # draw random alpha
      targets1[i,3] = rnorm(1, mean=H, sd=H.sd) # draw random theta
    }  
     
    # from these target SW values, generate the distribution of observed RTs for each subject
    # for congruent trials
    for (i in 1:nSub){
      rts[i,1:nTrials] = rwald(nTrials, 
                      gamma = targets1[i,1], 
                      alpha = targets1[i,2], 
                      theta = targets1[i,3])
    }
    
    # Step 2 -- instantiate an effect in theta, where effect 
    # drawn from unconstrained normal, using settings from
    # Rouder, Kumar, & Haaf (2019)
    
    delta = numeric(nSub)
    
    for (i in 1:nSub){
      if(correctModel=="positive"){
        mu.delta = rnorm(1, mean=0.1, sd=0.05)
        s2.delta = rinvgamma(1, 2, 0.03^2)
        delta[i] = rtnorm(1, mean=mu.delta, sd=sqrt(s2.delta), lower=0)
      }
      else{
        mu.delta = rnorm(1, mean=0.1, sd=0.10)
        s2.delta = rinvgamma(1, 2, 0.03^2)
        delta[i] = rnorm(1, mean=mu.delta, sd=sqrt(s2.delta))
      }
    }
      
    
    # generate new targets with effect added to theta parameter
    for (i in 1:nSub){
      targets2[i,1] = targets1[i,1]
      targets2[i,2] = targets1[i,2]
      targets2[i,3] = targets1[i,3] + delta[i] 
    }  
    
    # generate new RTs for incongruent trials
    for (i in 1:nSub){
      rts[i,(nTrials+1):(2*nTrials)] = rwald(nTrials, 
                               gamma = targets2[i,1], 
                               alpha = targets2[i,2], 
                               theta = targets2[i,3])
    }
    
    # build data frame from simulated RTs
    Y = c(rts) # converts matrix to single vector
    sub = rep(1:nSub,2*nTrials)
    cond = c(rep(1,nSub*nTrials), rep(2,nSub*nTrials))
    
    prep = prep.models(sub,cond)
    BF = makeBF(y=Y, meanScale=1/3, effectScale=1/5, prep)
    BFlog = makeBF(y=log(Y-0.95*min(Y)), meanScale=1/3, effectScale=1/5, prep)
  
    newrun = data.frame(bf,bfLog)
    newrun$bf = log(BF$bf.pu)
    newrun$bfLog = log(BFlog$bf.pu)
    dat = rbind(dat, newrun)
    write.csv(dat, file = filename, row.names=FALSE)
    cat(sprintf("Finished run %d of %d\n", j, numSims))
  }
}

runSim(nSub=20, nTrials=50, correctModel="unconstrained", numSims=200)
runSim(nSub=20, nTrials=100, correctModel="unconstrained", numSims=200)
runSim(nSub=80, nTrials=50, correctModel="unconstrained", numSims=200)
runSim(nSub=80, nTrials=100, correctModel="unconstrained", numSims=200)

runSim(nSub=20, nTrials=50, correctModel="positive", numSims=200)
runSim(nSub=20, nTrials=100, correctModel="positive", numSims=200)
runSim(nSub=80, nTrials=50, correctModel="positive", numSims=200)
runSim(nSub=80, nTrials=100, correctModel="positive", numSims=200)


# load data

P.20.50 =  read.csv("P-20-50.csv")
P.20.100 = read.csv("P-20-100.csv")
P.20.200 = read.csv("P-20-200.csv")

U.20.50 =  read.csv("U-20-50.csv")
U.20.100 = read.csv("U-20-100.csv")
U.20.200 = read.csv("U-20-200.csv")

P.80.50 =  read.csv("P-80-50.csv")
P.80.100 = read.csv("P-80-100.csv")
P.80.200 = read.csv("P-80-200.csv")

U.80.50 =  read.csv("U-80-50.csv")
U.80.100 = read.csv("U-80-100.csv")
U.80.200 = read.csv("U-80-200.csv")



## check accuracy
sum(P.20.50$bf>0)/200
sum(P.20.50$bfLog>0)/200

sum(P.20.100$bf>0)/200
sum(P.20.100$bfLog>0)/200

sum(P.20.200$bf>0)/200
sum(P.20.200$bfLog>0)/200

sum(P.80.50$bf>0)/200
sum(P.80.50$bfLog>0)/200

sum(P.80.100$bf>0)/200
sum(P.80.100$bfLog>0)/200

sum(P.80.200$bf>0)/200
sum(P.80.200$bfLog>0)/200

sum(U.20.50$bf<0)/200
sum(U.20.50$bfLog<0)/200

sum(U.20.100$bf<0)/200
sum(U.20.100$bfLog<0)/200

sum(U.20.200$bf<0)/200
sum(U.20.200$bfLog<0)/200


sum(U.80.50$bf<0)/200
sum(U.80.50$bfLog<0)/200

sum(U.80.100$bf<0)/200
sum(U.80.100$bfLog<0)/200

sum(U.80.200$bf<0)/200
sum(U.80.200$bfLog<0)/200


# check consistency
sum(P.20.50$bf*P.20.50$bfLog>0)/200
sum(P.20.100$bf*P.20.100$bfLog>0)/200
sum(P.20.200$bf*P.20.200$bfLog>0)/200

sum(P.80.50$bf*P.80.50$bfLog>0)/200
sum(P.80.100$bf*P.80.100$bfLog>0)/200
sum(P.80.200$bf*P.80.200$bfLog>0)/200

sum(U.20.50$bf*U.20.50$bfLog>0)/200
sum(U.20.100$bf*U.20.100$bfLog>0)/200
sum(U.20.200$bf*U.20.200$bfLog>0)/200

sum(U.80.50$bf*U.80.50$bfLog>0)/200
sum(U.80.100$bf*U.80.100$bfLog>0)/200
sum(U.80.200$bf*U.80.200$bfLog>0)/200
