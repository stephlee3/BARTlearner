## simulation data
# devtools::install_github("vdorie/aciccomp/2016")
library(aciccomp2016)
parameterNum = 1:77
simulationNum = 1:100
dat = dgp_2016(input_2016,parameters =12,random.seed = 20)

## BART package
# devtools::install_github("vdorie/bartCause")
library(bartCause)
library(tidyverse)

## BART learners for estimating HTE
## T learner
T_learner = function(X, z, y, ndpost = 200, seed = 1){
  set.seed(seed)
  dat = data.frame(X,z,y)
  
  X1 = dat %>% filter(z==1) %>% dplyr::select(-z, -y)
  X0 = dat %>% filter(z==0) %>% dplyr::select(-z, -y)
  
  bart_mu1 = dbarts::bart(x.train = X1, y.train = y[which(z==1)],
                          ndpost= ndpost,keeptrees = T)
  bart_mu0 = dbarts::bart(x.train = X0, y.train = y[which(z==0)],
                          ndpost= ndpost,keeptrees = T)
  y1.pred = predict(bart_mu1, newdata = dat)
  y0.pred = predict(bart_mu0, newdata = dat)
  
  
  
  tau.pred.mat = y1.pred - y0.pred
  tau.pred = apply(tau.pred.mat,2,mean)
  
  ans = list(tau.pred = tau.pred, tau.pred.mat = tau.pred.mat)
  return(ans)
}


## S learner
S_learner = function(X, z, y, ndpost = 200, seed = 1){
  set.seed(seed)
  dat = data.frame(X,z,y)
  
  bart_mu = dbarts::bart(x.train = dat %>% dplyr::select(-y) , y.train = y,
                         ndpost= ndpost,keeptrees = T)
  
  dat1 = dat; dat1$z = 1;
  dat0 = dat; dat0$z = 0;
  
  y1.pred = predict(bart_mu,newdata = dat1)
  y0.pred = predict(bart_mu,newdata = dat0)
  
  tau.pred.mat = y1.pred - y0.pred
  tau.pred = apply(tau.pred.mat,2,mean)
  
  ans = list(tau.pred = tau.pred, tau.pred.mat = tau.pred.mat)
  return(ans)
}


## X learner
X_learner = function(X, z, y, ndpost = 500, seed = 1){
  set.seed(seed)
  # ps.pred = dat$e
  dat = data.frame(X,z,y)
  
  # propensity score
  ps_bart = dbarts::bart(X,z,ndpost = ndpost,keeptrees=T)
  ps.pred = pnorm(as.vector((ps_bart$yhat.train)[ndpost,]))
  
  X1 = dat %>% filter(z==1) %>% dplyr::select(-z, -y)
  X0 = dat %>% filter(z==0) %>% dplyr::select(-z, -y)
  
  ## response surface
  bart_mu1 = dbarts::bart(x.train = X1, y.train = y[which(z==1)],
                          ndpost= ndpost,keeptrees = T)
  bart_mu0 = dbarts::bart(x.train = X0, y.train = y[which(z==0)],
                          ndpost= ndpost,keeptrees = T)
  y1.pred.ctrl = predict(bart_mu1, newdata = X0)
  y0.pred.trt = predict(bart_mu0, newdata = X1)
  
  tau.impute.trt = y[which(z==1)] - apply(y0.pred.trt,2,mean)
  tau.impute.ctrl = apply(y1.pred.ctrl,2, mean) - y[which(z==0)]
  
  impute_trt_bart = dbarts::bart(x.train = X1,y.train =tau.impute.trt,
                                 ndpost=ndpost,keeptrees = T)
  tau.impute.trt.pred.mat = predict(impute_trt_bart, newdata = X)
  tau.impute.trt.pred = apply(tau.impute.trt.pred.mat,2,mean)
  
  impute_ctrl_bart = dbarts::bart(x.train = X0,y.train =tau.impute.ctrl,
                                  ndpost=ndpost,keeptrees = T)
  tau.impute.ctrl.pred.mat = predict(impute_ctrl_bart, newdata = X)
  tau.impute.ctrl.pred = apply(tau.impute.ctrl.pred.mat,2,mean)
  
  tau.pred.mat = ps.pred * tau.impute.ctrl.pred.mat + (1-ps.pred) * tau.impute.trt.pred.mat
  
  tau.pred = ps.pred * tau.impute.ctrl.pred + (1-ps.pred) * tau.impute.trt.pred
  
  ans = list(tau.pred = tau.pred, tau.pred.mat = tau.pred.mat)
  return(ans)
  
}


## R_learner
R_learner = function(X, z, y, ndpost = 500, seed = 1, fold = 10){
  
  set.seed(seed)
  
  ## propensity score and mean (10 fold)
  ps.pred = rep(NA, nrow(X)); m.pred = rep(NA, nrow(X))
  if(fold > 0){
    randnum = sample(nrow(X)) / nrow(X)
    for(i in 1:fold){
      ind = (randnum <= i/fold) & (randnum > (i-1)/fold)
      ps_bart = dbarts::bart(X[!ind,],z[!ind],ndpost = ndpost,keeptrees=T)
      m_bart = dbarts::bart(X[!ind,],y[!ind],ndpost = ndpost,keeptrees = T)
      ps.pred[ind] = predict(ps_bart, X[ind,], type = "ev")[ndpost,]
      m.pred[ind] = predict(m_bart, X[ind,], type = "ev")[ndpost,]
    }  
  }
  else{
    ps_bart = dbarts::bart(X,z,ndpost = ndpost,keeptrees=T)
    ps.pred = pnorm(as.vector((ps_bart$yhat.train)[ndpost,]))
    
    m_bart = dbarts::bart(X,y,ndpost = ndpost,keeptrees = T)
    m.pred = m_bart$yhat.train.mean
  }
    
  transform_y = (y-m.pred)/(z-ps.pred)
  oracle_bart = dbarts::bart(X,transform_y, weights = (z-ps.pred)^2,
                             ndpost = ndpost, keeptrees = T)
  tau.pred.mat = oracle_bart$yhat.train
  tau.pred = oracle_bart$yhat.train.mean
  
  ans = list(tau.pred = tau.pred, tau.pred.mat = tau.pred.mat)
  return(ans)
}


## DR_learner: augmented IPTW
DR_learner = function(X, z, y, ndpost = 500, seed = 1){
  
  set.seed(seed)
  dat = data.frame(X,z,y)
  
  X1 = dat %>% filter(z==1) %>% dplyr::select(-z, -y)
  X0 = dat %>% filter(z==0) %>% dplyr::select(-z, -y)
  
  ## propensity score
  ps_bart = dbarts::bart(X,z,ndpost = ndpost,keeptrees=T)
  ps.pred = pnorm(as.vector((ps_bart$yhat.train)[ndpost,]))
  
  ## mu1 and mu0 
  bart_mu1 = dbarts::bart(x.train = X1, y.train = y[which(z==1)],
                          ndpost= ndpost,keeptrees = T)
  bart_mu0 = dbarts::bart(x.train = X0, y.train = y[which(z==0)],
                          ndpost= ndpost,keeptrees = T)
  y1.pred = apply(predict(bart_mu1, newdata = X),2,mean)
  y0.pred = apply(predict(bart_mu0, newdata = X),2,mean)
  
  ## Augmented IPTW
  tau.pseudo = (z- ps.pred)/(ps.pred * (1-ps.pred)) * (y- (z*y1.pred +(1-z) * y0.pred)) + 
    y1.pred - y0.pred
  
  bart_tau = dbarts::bart(x.train = X, y.train = tau.pseudo,
                          ndpost= ndpost,keeptrees = T)
  tau.pred.mat = bart_tau$yhat.train
  tau.pred =  bart_tau$yhat.train.mean
  
  ans = list(tau.pred = tau.pred, tau.pred.mat = tau.pred.mat)
  return(ans)
}

## Run the examples
tau.T = T_learner(input_2016,dat$z,dat$y)
tau.S = S_learner(input_2016,dat$z,dat$y)
tau.X = X_learner(input_2016,dat$z,dat$y)
tau.R = R_learner(input_2016,dat$z,dat$y, fold = 0)
tau.R.fold = R_learner(input_2016,dat$z,dat$y, fold = 10)
tau.DR = DR_learner(input_2016,dat$z,dat$y)

## Evaluation
## true value
tau.true =  dat$y.1 - dat$y.0

## MSE
mse = function(a,hat_a){
  mean((a-hat_a)^2)
}

## 95% CI and cover 
CI_eval = function(post.mat, method="norm", true.value = tau.true){ # A is the posterior sampling matrix
  post.mean = apply(post.mat,2,mean)
  post.sd = apply(post.mat,2,sd)
  if (method=="norm"){
    post.CI= cbind(post.mean-1.96 * post.sd, post.mean + 1.96 * post.sd)
    post.len = mean(post.CI)
    cover.ind = ifelse(post.CI[,1]<= true.value & post.CI[,2]>= true.value,1,0)
    cover.mean = mean(cover.ind)
  }
  
  if (method=="quantile"){
    post.CI = t(apply(A,2,quantile,probs= c(0.025,0.975)))
    post.len = mean(post.CI)
  }
  
  ans = list(len = post.len, cover = cover.mean)
  return(ans)
}

## Results
## MSE and CI
mse.T = mse(tau.true, tau.T$tau.pred)
mse.S = mse(tau.true, tau.S$tau.pred)
mse.X = mse(tau.true, tau.X$tau.pred)
mse.R = mse(tau.true, tau.R$tau.pred)
mse.R.fold = mse(tau.true, tau.R.fold$tau.pred)
mse.DR = mse(tau.true, tau.DR$tau.pred)

CI.T = CI_eval(tau.T$tau.pred.mat)
CI.S = CI_eval(tau.S$tau.pred.mat)
CI.X = CI_eval(tau.X$tau.pred.mat)
CI.R = CI_eval(tau.R$tau.pred.mat)
CI.R.fold = CI_eval(tau.R.fold$tau.pred.mat)
CI.DR = CI_eval(tau.DR$tau.pred.mat)

model_eval = data.frame(
  MSE = c(mse.T, mse.S, mse.X,mse.R,mse.DR,mse.R.fold),
  interval_length = c(CI.T$len,CI.S$len,CI.X$len,CI.R$len,CI.DR$len,CI.R.fold$len),
  cover_indicator = c(CI.T$cover,CI.S$cover,CI.X$cover,CI.R$cover,CI.DR$cover,CI.R.fold$cover)
)
rownames(model_eval)= c("T","S","X","R","DR","R.fold")
print(model_eval)

## plot
result = data.frame(tau.true,tau.T = tau.T$tau.pred,tau.S = tau.S$tau.pred,
                    tau.X = tau.X$tau.pred,tau.R = tau.R$tau.pred,
                    tau.DR = tau.DR$tau.pred,tau.R.fold = tau.R.fold$tau.pred)
result.gather = result %>%
  gather(algorithm,tau.pred,-tau.true)
result.gather %>% 
  ggplot(aes(x = tau.true,y = tau.pred))+
  geom_point(alpha= 0.2) + facet_grid(cols = vars(algorithm)) +
  geom_abline(slope = 1,intercept = 0,col="red") +
  theme_bw()



