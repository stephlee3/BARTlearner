
# BARTlearner

<!-- badges: start -->
<!-- badges: end -->

The goal of BARTlearner is to unify meta-algorithms for estimating HTE via BART. We implemented S,T,X,R,DR learner and used the simulation data from 2016 Atlantic Causal Inference Competition to evaluate the model performance

## Simulation Data
You need to install the following R package to obtain the data
```
devtools::install_github("vdorie/aciccomp/2016")
```
We set the `parameters = 12` and `random.seed = 20` in our example. 

## BART
You need to install the following R package to fit BART models.
```
devtools::install_github("vdorie/bartCause")
```
Some other packages for BART: `bayestree`,`BART`,`dbarts`.

## Meta Algorithms
We wrapped up several meta algorithms into some built-in functions:  `T_learner`, `S_learner`, `X_learner`, `R_learner`,`DR_learner`. For example,

```
T_learner = function(X, z, y, ndpost = 200, seed = 1)
```
Your input is `X` matrix of baseline covariates , `z` vector of treatment assignment, `y` vector of outcome, `ndpost` the number of posterior draws after burn in, `seed` the random seed. Your output is `tau.pred.mat` the matrix of posterior sampling for treatment effect estimation (each row is a posterior draw, each column is an observation) and `tau.pred` the vector of means of posterior sampling matrix.


## Evaluation
We evaluate the model performance by the following metrics: 

- `MSE`: mean squared error
- `interval length`: the intervel length of 95% CI of the treatment effect estimation
- `cover indicator`: whether the 95% CI covers the true value.



