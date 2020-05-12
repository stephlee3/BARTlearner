## iterative R_learner
iter_R_learner = function(X, z, y, ndpost = 500, seed = 10, niter = 25){
  
  set.seed(seed)
  
  ## propensity score 
  ps_bart = dbarts::bart(X,z,ndpost = ndpost,keeptrees=T)
  ps.pred = pnorm(as.vector((ps_bart$yhat.train)[ndpost,]))
  
  m_bart = dbarts::bart(X,y,ndpost = ndpost,keeptrees = T)
  m.pred = m_bart$yhat.train.mean
  
  transform_y = (y-m.pred)/(z-ps.pred)
  oracle_bart = dbarts::bart(X,transform_y, weights = (z-ps.pred)^2,
                             ndpost = ndpost, keeptrees = T)
  tau.pred.mat = oracle_bart$yhat.train
  tau.pred = oracle_bart$yhat.train.mean
  
  ## iteratively update m and tau
  mse = rep(0,niter)
  max.abs.err = rep(0,niter)
  
  for(i in 1:niter){
    print(paste0("iteration:",i))
    
    ## m 
    m.impute = y - (z-ps.pred)* tau.pred
    m_bart_iter = dbarts::bart(X,m.impute,ndpost = ndpost,keeptrees = T)
    m.pred.iter = m_bart_iter$yhat.train.mean
    
    ## tau
    y.impute = (y-m.pred.iter)/(z-ps.pred)
    tau_bart_iter = dbarts::bart(X,y.impute, weights = (z-ps.pred)^2,
                                 ndpost = ndpost, keeptrees = T)
    max.abs.err[i] = max(abs(tau.pred-tau_bart_iter$yhat.train.mean))
    mse[i] = mean((tau.pred - tau_bart_iter$yhat.train.mean)^2)
    
    tau.pred.mat = tau_bart_iter$yhat.train
    tau.pred = tau_bart_iter$yhat.train.mean
    
    
    print("--------------------------")
  }
  print(max.abs.err)
  print(mse)
  ans = list(tau.pred = tau.pred, tau.pred.mat = tau.pred.mat)
  return(ans)
}


tau.R.iter = iter_R_learner(input_2016,dat$z,dat$y)
mse.R.iter = mse(tau.true, tau.R.iter$tau.pred)
CI.R.iter = CI_eval(tau.R.iter$tau.pred.mat)


CI.R.subset = CI_eval(tau.R$tau.pred.mat[,tau.true>20], true.value = tau.true[tau.true>20])
CI.R.iter.subset = CI_eval(tau.R.iter$tau.pred.mat[,tau.true>20], true.value = tau.true[tau.true>20])


print(mse.R); print(mse.R.iter)
print(CI.R); print(CI.R.iter)
print(CI.R.subset); print(CI.R.iter.subset)

## R learner vs iter R learner
result = data.frame(tau.true,
                    tau.R = tau.R$tau.pred,
                    tau.R.iter = tau.R.iter$tau.pred)
result.gather = result %>%
  # filter(tau.true>20) %>%
  gather(algorithm,tau.pred,-tau.true)
result.gather %>% 
  ggplot(aes(x = tau.true,y = tau.pred,col = algorithm))+
  geom_point(alpha= 0.5)+
  geom_abline(slope = 1,intercept = 0,col="red")
