function [pt1,pdfBias] = propDis(pt0,bnd,opts)

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

reRoll = true(size(pt0));
n = 7; %scaling factor to turn range of randn to 1.

bndRng = (bnd(:,2)-bnd(:,1));
bndRng(isinf(bndRng)) = 1;

%Determine scale
logTest = log(bnd(:,2))-log(bnd(:,1));
logScale = logTest>1&imag(logTest)==0;

%%
% Generate random variable. Undirected
while sum(reRoll)
	prop = randn(size(pt0))/n.*opts.step.*bndRng;
	pt1(reRoll&~logScale) = pt0(reRoll&~logScale)+prop(reRoll&~logScale);
	pt1(reRoll&logScale) = pt0(reRoll&logScale).*exp(prop(reRoll&logScale));
	reRoll = pt1<bnd(:,1)|pt1>bnd(:,2);
end

% Generate random variable. Directed
% nPt = length(pt0);
% curBasis = opts.basis;
% try
% 	basis = [curBasis null(curBasis')];
% catch
% 	keyboard
% end
% pow = opts.basisSkew;
% alpha = 0.75;
% 
% while sum(reRoll)
% 	pt1(reRoll) = pt0(reRoll) + propFunc(randn(sum(reRoll),1),skew)/n.*opts.step;
% 	reRoll = pt1<bnd(:,1)|pt1>bnd(:,2);
% end

%%
% Calculate biasing by truncating the cumulative distribution function.
% This is for the hastings ratio

%Create mixed log and lin scale points
pt0_M = pt0;
pt0_M(logScale) = log(pt0(logScale));
pt1_M = pt1;
pt1_M(logScale) = log(pt1(logScale));
bnd_M = bnd;
bnd_M(logScale,:) = log(bnd(logScale,:));

pdfBias = biasCalc(pt0_M,pt1_M,bnd_M,n,opts.step);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% End Function %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pdfBias = biasCalc(pt0,pt1,bnd,n,step)
	if ~isempty(bnd)
		xy = biasLin(bnd(:,1)-pt0,bnd(:,2)-pt0,n./step);
		yx = biasLin(bnd(:,1)-pt1,bnd(:,2)-pt1,n./step);
		
		pdfBias = prod(xy)/prod(yx);
	else
		pdfBias = 1;   % Unbounded MCMC run is by definition not biased in it's proposal distribution.
	end
end

%%%%%%%%%%%%

function rat = biasLin(x1,x2,n)

rat = 1./((erf(x2/sqrt(2).*n)+1)/2-(erf(x1/sqrt(2).*n)+1)/2);
    %Cumulative function to ub take away cumulative function to lb
end

function pt = propFunc(x,skew)
n = 1;
skewCoeff = 1;

skew = (skew/skewCoeff).^n;

pt = 2*normpdf(x,0,1).*normcdf(3.1*skew./(1+abs(skew)).*x,0,1);

%the skew = (skew/skewCoeff).^n makes the skew scale like a hill function.
%the larger the skew, the more biased in that direction, up to the point
%where only 10% of the sampled points will be in the non-skewed direction.

end