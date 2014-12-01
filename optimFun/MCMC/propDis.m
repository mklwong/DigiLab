function [pt1,pdfBias,logScl] = propDis(pt0,bnd,opts)

% pdfBias is q(x,y)/q(y,x)

% Make p a column
if isrow(pt0)
	pt0 = pt0';
end

if isrow(opts.step)
	opts.step = opts.step';
end

%initialise stuff
pt1 = 0*pt0;

% Determine log scale
logScl = ((log10(bnd(:,2))-log10(bnd(:,1)))>1) &(bnd(:,1)>=0);

reRoll = true(size(pt0));
n = 7; %scaling factor to turn range of randn to 1.

% Generate random variable
while sum(reRoll)
	prop = randn(size(pt0))/n.*opts.step;
	pt1(~logScl&reRoll) = pt0(~logScl&reRoll)+prop(~logScl&reRoll);
	pt1(logScl&reRoll)  = pt0(logScl&reRoll).*10.^(prop(logScl&reRoll));
	reRoll = pt1<bnd(:,1)|pt1>bnd(:,2);
end

% Calculate biasing by truncating the cumulative distribution function.
% This is for the hastings ratio
pdfBias = biasCalc(pt0,pt1,bnd,n,logScl,opts.step);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% End Function %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pdfBias = biasCalc(pt0,pt1,bnd,n,logScl,step)
xyBias = 0*pt0;
yxBias = 0*pt0;

xyLin = biasLin(bnd(:,1)-pt0,bnd(:,2)-pt0,n./step);
yxLin = biasLin(bnd(:,1)-pt1,bnd(:,2)-pt1,n./step);

try
	xyLog = biasLin(log10(bnd(logScl,1)./pt0(logScl)),log10(bnd(logScl,2)./pt0(logScl)),n./step(logScl));
	yxLog = biasLin(log10(bnd(logScl,1)./pt1(logScl)),log10(bnd(logScl,2)./pt1(logScl)),n./step(logScl));
catch
	keyboard
end
xyBias(~logScl) = xyLin(~logScl);
xyBias(logScl) = xyLog;
yxBias(~logScl) = yxLin(~logScl); 
yxBias(logScl) = yxLog;

pdfBias = prod(xyBias)/prod(yxBias);
end

%%%%%%%%%%%%

function rat = biasLin(x1,x2,n)

rat = 1./((erf(x2/sqrt(2).*n)+1)/2-(erf(x1/sqrt(2).*n)+1)/2);
    %Cumulative function to ub take away cumulative function to lb
end