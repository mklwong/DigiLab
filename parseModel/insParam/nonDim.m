function [matNonDim,matRamp] = nonDim(inMat,tspan,nDimInd)

nonDimFac = tspan(2)-tspan(1);

matRamp = inMat;
matRamp(:,nDimInd) = 0;

matNonDim = inMat;
matNonDim(:,nDimInd) = inMat(:,nDimInd)*nonDimFac;
