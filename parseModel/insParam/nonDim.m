function matNonDim = nonDim(inMat,tspan,nDimInd)

nonDimFac = tspan(2)-tspan(1);

matNonDim = inMat;
matNonDim(:,nDimInd) = inMat(:,nDimInd)*nonDimFac;
