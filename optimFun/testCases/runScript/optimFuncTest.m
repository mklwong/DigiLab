function pts = optimFuncTest(obj,bnd,opts)

[pts,logP,ptUnq] = MCMC(obj,[],bnd,opts);
result.logP  = logP;	result.pts = pts;	result.ptUn = ptUnq;
opts = MCMCOptimset(opts,'T',1,'prior',result);
[pts,~,~] = MCMC(obj,[],bnd,opts);