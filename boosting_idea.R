library(aciccomp2016)

input_2016 = aciccomp2016::input_2016
parameterNum = 1:77
simulationNum = 1:100
dat = dgp_2016(input_2016,parameters =12,random.seed = 20)

## true value
tau.true =  dat$y.1 - dat$y.0

## MSE
mse = function(a,hat_a){
  mean((a-hat_a)^2)
}

library(gbm)
library(tidyverse)

for(i in 1:ncol(input_2016)){
  if(class(input_2016[[i]]) == "integer"){
    if(length(unique(input_2016[[i]])) > 5)
      input_2016[[i]] = as.numeric(input_2016[[i]])
    else
      input_2016[[i]] = as.factor(input_2016[[i]])
  }
}
sapply(input_2016, class)

## T learner
T_boosting = function(X, z, y, seed = 1, ntrees = 100, ...){
  set.seed(seed)
  dat = data.frame(X,z,y)
  
  X1 = dat %>% filter(z==1) %>% dplyr::select(-z, -y)
  X0 = dat %>% filter(z==0) %>% dplyr::select(-z, -y)
  
  boost_mu1 = gbm(Y~., data = data.frame(X1, Y=y[z==1]), 
                  n.trees = ntrees, ...)
  boost_mu0 = gbm(Y~., data = data.frame(X0, Y=y[z==0]), 
                  n.trees = ntrees, ...)
  
  y1.pred = predict(boost_mu1, newdata = X, n.trees = ntrees)
  y0.pred = predict(boost_mu0, newdata = X, n.trees = ntrees)
  
  tau.pred = y1.pred - y0.pred
  
  ans = list(tau.pred = tau.pred)
  return(ans)
}

tau.boostT = T_boosting(input_2016,dat$z,dat$y, 
                        interaction.depth = 3, shrinkage = 0.1)
mse.boostT = mse(tau.true, tau.boostT$tau.pred)

## S learner
S_boosting = function(X, z, y, seed = 1, ntrees = 100, ...){
  set.seed(seed)
  dat = data.frame(X,Z=z,Y=y)
  
  boost_mu = gbm(Y~., data = dat, n.trees = ntrees, ...)
                         
  dat1 = dat; dat1$Z = 1;
  dat0 = dat; dat0$Z = 0;
  
  y1.pred = predict(boost_mu, newdata = dat1, n.trees = ntrees)
  y0.pred = predict(boost_mu, newdata = dat0, n.trees = ntrees)
  
  tau.pred = y1.pred - y0.pred
  
  ans = list(tau.pred = tau.pred)
  return(ans)
}

tau.boostS = S_boosting(input_2016,dat$z,dat$y, 
                        interaction.depth = 3, shrinkage = 0.1)
mse.boostS = mse(tau.true, tau.boostS$tau.pred)

## X learner
X_boosting = function(X, z, y, seed = 1, ntrees = 100, ...){
  set.seed(seed)
  # ps.pred = dat$e
  dat = data.frame(X,z,y)
  
  # propensity score
  ps_boost = gbm(Z~., distribution = "bernoulli", 
                 data=data.frame(X,Z=z), n.trees = ntrees,...)
  ps.pred = plogis(ps_boost$fit)
  
  X1 = dat %>% filter(z==1) %>% dplyr::select(-z, -y)
  X0 = dat %>% filter(z==0) %>% dplyr::select(-z, -y)
  
  ## response surface
  boost_mu1 = gbm(Y~., data = data.frame(X1, Y=y[z==1]), 
                  n.trees = ntrees, ...)
  boost_mu0 = gbm(Y~., data = data.frame(X0, Y=y[z==0]), 
                  n.trees = ntrees, ...)
  
  y1.pred.ctrl = predict(boost_mu1, newdata = X0, n.trees = ntrees)
  y0.pred.trt = predict(boost_mu0, newdata = X1, n.trees = ntrees)
  
  tau.impute.trt = y[z==1] - y0.pred.trt
  tau.impute.ctrl = y1.pred.ctrl - y[z==0]
  
  impute_trt_boost = gbm(Y~., data = data.frame(X1, Y=tau.impute.trt), 
                         n.trees = ntrees, ...)
  tau.impute.trt.pred = predict(impute_trt_boost, newdata = X, n.trees = ntrees)
  
  impute_ctrl_boost = gbm(Y~., data = data.frame(X0, Y=tau.impute.ctrl), 
                         n.trees = ntrees, ...)
  tau.impute.ctrl.pred = predict(impute_ctrl_boost, newdata = X, n.trees = ntrees)
  
  tau.pred = ps.pred * tau.impute.ctrl.pred + (1-ps.pred) * tau.impute.trt.pred
  
  ans = list(tau.pred = tau.pred)
  return(ans)
}

tau.boostX = X_boosting(input_2016,dat$z,dat$y, 
                        interaction.depth = 3, shrinkage = 0.1)
mse.boostX = mse(tau.true, tau.boostX$tau.pred)


## R_learner
R_boosting = function(X, z, y, seed = 1, fold = 10, ntrees = 100, ...){
  set.seed(seed)
  
  ## propensity score and mean (10 fold)
  ps.pred = rep(NA, nrow(X)); m.pred = rep(NA, nrow(X))
  if(fold > 0){
    randnum = sample(nrow(X)) / nrow(X)
    for(i in 1:fold){
      ind = (randnum <= i/fold) & (randnum > (i-1)/fold)
      ps_boost = gbm(Z ~ ., distribution = "bernoulli",
                     data=data.frame(X[!ind,],Z=z[!ind]),n.trees = ntrees,...)
      m_boost = gbm(Y ~ ., data=data.frame(X[!ind,],Y=y[!ind]),
                    n.trees = ntrees,...)
      ps.pred[ind] = predict(ps_boost, X[ind,], type = "response", n.trees = ntrees)
      m.pred[ind] = predict(m_boost, X[ind,], n.trees = ntrees)
    }  
  }
  else{
    ps_boost = gbm(Z ~ ., distribution = "bernoulli",
                   data=data.frame(X,Z=z),n.trees = ntrees,...)
    m_boost = gbm(Y ~ ., data=data.frame(X,Y=y), n.trees = ntrees,...)
    ps.pred = plogis(ps_boost$fit)
    m.pred = m_boost$fit
  }
  
  transform_y = (y-m.pred)/(z-ps.pred)
  oracle_boost = gbm(Y~., data=data.frame(X,Y=transform_y), 
                     weights = (z-ps.pred)^2, n.trees = ntrees, ...)
  tau.pred = oracle_boost$fit
  
  ans = list(tau.pred = tau.pred)
  return(ans)
}

tau.boostR = R_boosting(input_2016,dat$z,dat$y, fold = 0,
                        interaction.depth = 3, shrinkage = 0.1)
mse.boostR = mse(tau.true, tau.boostR$tau.pred)

## DR_learner: augmented IPTW
DR_boosting = function(X, z, y, seed = 1, ntrees = 100, ...){
  
  set.seed(seed)
  dat = data.frame(X,z,y)
  
  X1 = dat %>% filter(z==1) %>% dplyr::select(-z, -y)
  X0 = dat %>% filter(z==0) %>% dplyr::select(-z, -y)
  
  ## propensity score
  ps_boost = gbm(Z~., distribution = "bernoulli", 
                 data=data.frame(X,Z=z), n.trees = ntrees,...)
  ps.pred = plogis(ps_boost$fit)
  
  ## mu1 and mu0 
  boost_mu1 = gbm(Y~., data = data.frame(X1, Y=y[z==1]), 
                  n.trees = ntrees, ...)
  boost_mu0 = gbm(Y~., data = data.frame(X0, Y=y[z==0]), 
                  n.trees = ntrees, ...)
  y1.pred = predict(boost_mu1, newdata = X, n.trees = ntrees)
  y0.pred = predict(boost_mu0, newdata = X, n.trees = ntrees)
  
  ## Augmented IPTW
  tau.pseudo = (z- ps.pred)/(ps.pred * (1-ps.pred)) * (y- (z*y1.pred +(1-z) * y0.pred)) + 
    y1.pred - y0.pred
  
  boost_tau = gbm(Y~., data=data.frame(X, Y=tau.pseudo), n.trees = ntrees, ...)
  tau.pred =  boost_tau$fit
  
  ans = list(tau.pred = tau.pred)
  return(ans)
}

tau.boostDR = DR_boosting(input_2016,dat$z,dat$y,
                          interaction.depth = 3, shrinkage = 0.1)
mse.boostDR = mse(tau.true, tau.boostDR$tau.pred)

## double boosting
R_Double_Boosting = function(X, z, y, seed = 1, ntrees = 100, ninit = 10, ...){
  
  set.seed(seed)
  
  ## propensity score
  ps_boost = gbm(Z ~ ., distribution = "bernoulli",
                 data=data.frame(X,Z=z),n.trees = ninit,...)
  ps.pred = plogis(ps_boost$fit)
  
  ## initial
  m_trees = gbm(Y ~ ., data = data.frame(X, Y=y), n.trees = ninit, ...)
  y_m = (y - m_trees$fit) / (z - ps.pred)
  tao_trees = gbm(Y ~ ., data = data.frame(X, Y=y_m), 
                  weights = (z-ps.pred)^2, n.trees = ninit, ...)
  
  for(t in 1:ntrees){
    y_tao = y - (z-ps.pred) * tao_trees$fit
    y_m = (y - m_trees$fit) / (z - ps.pred)
    m_trees = gbm.more(m_trees, n.new.trees = 1, data = data.frame(X, Y = y_tao))
    tao_trees = gbm.more(tao_trees, n.new.trees = 1, 
                         data = data.frame(X, Y = y_m),
                         weights = (z-ps.pred)^2)  
  }
  tau.pred = tao_trees$fit
  
  ans = list(tau.pred = tau.pred, model = list(tao=tao_trees,m=m_trees,ps=ps_boost))
  return(ans)
}

tau.dR = R_Double_Boosting(input_2016,dat$z,dat$y, interaction.depth = 3, shrinkage = 0.1,
                           ntrees = 10, ninit = 100)
mse.dR = mse(tau.true, tau.dR$tau.pred)

## plot
result = data.frame(tau.true,
                    tau.T = tau.boostT$tau.pred,
                    tau.S = tau.boostS$tau.pred,
                    tau.X = tau.boostX$tau.pred,
                    tau.R = tau.boostR$tau.pred,
                    tau.DR = tau.boostDR$tau.pred,
                    tau.dR = tau.dR$tau.pred)
result.gather = result %>%
  gather(algorithm,tau.pred,-tau.true)
result.gather %>% 
  ggplot(aes(x = tau.true,y = tau.pred))+
  geom_point(alpha= 0.2) + facet_grid(cols = vars(algorithm)) +
  geom_abline(slope = 1,intercept = 0,col="red") +
  theme_bw()


