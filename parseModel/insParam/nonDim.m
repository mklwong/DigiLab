function matNonDim = nonDim(inMat,tspan,nDimInd)

nonDimFac = tspan(end)-tspan(1);

matNonDim = inMat;
matNonDim(:,nDimInd) = inMat(:,nDimInd)*nonDimFac;
