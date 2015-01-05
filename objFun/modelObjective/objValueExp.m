function resid = objValueExp(YSim,YExpMean,YExpStd,normNow)

% Get number of x steps
[nx,~] = size(YExpMean);

% Rescale YExp to Simulation units
YExpMean = YExpMean.*(ones(nx,1)*normNow);
YExpStd = YExpStd.*(ones(nx,1)*normNow);
YExpStd(YExpStd==0)=nan;

% Find Residual
resid = nansum(((YExpMean-YSim)./YExpStd).^2)/2; %root mean square
resid = nansum(resid);                           %sum over states