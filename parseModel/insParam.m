function [modelOut,modelRamp] = insParam(model,p,tspan,normInp)

% Putting concentration into tensors
freeInd = model.modSpc.pInd>0;
concVals = model.modSpc.matVal;
if ~isempty(freeInd)
	concVals(freeInd,end) = model.modSpc.matVal(freeInd,end).*p(model.modSpc.pInd(freeInd));
end
modelOut.modSpc = concVals;

% Putting compartments into tensors
freeInd = model.modComp.pInd>0;
compVals = model.modComp.matVal;
if ~isempty(freeInd)
	compVals(freeInd,end) = model.modComp.matVal(freeInd,end).*p(model.modComp.pInd(freeInd));
	compVals = compVals(model.modSpc.comp);
end
modelOut.comp = compVals;

for ii = 1:length(model.param)
	freeInd = model.param(ii).pInd>0;
	%Substitute in parameter
	model.param(ii).matVal(freeInd) = model.param(ii).matVal(freeInd).*p(model.param(ii).pInd(freeInd));
    modelOut.param(ii).matVal = model.param(ii).matVal;
	modelOut.param(ii).name   = model.param(ii).name;
end

modelRamp = model.rxnRules('compile',modelOut,[0 0]);
modelOut = model.rxnRules('compile',modelOut,tspan);

[~,ii] = intersect({modelOut.param.name},'k0');
modelOut.sigma = @(t) normInp(t) + modelOut.param(ii).matVal*ones(1,length(t)); 