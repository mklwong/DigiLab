function result = optimFuncTest(obj,bnd,opts)

result = MCMC(obj,[],bnd,opts);
opts = MCMCOptimset(opts,'T',1,'prior',result);
result = MCMC(obj,[],bnd,opts);