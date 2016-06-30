function [pt1,pdfBias] = propDis(runVar)

% pdfBias is q(x,y)/q(y,x)

%initialise stuff
pt1 = 0*runVar.pt;

reRoll = true(size(runVar.pt));
n = 6; %scaling factor to turn range of randn to the width of the boundary.

if ~isempty(runVar.bnd)
	% Determine whether to use logarithmic scale to sample boundary or to use
	% linear
	bndRng = (runVar.bnd(:,2)-runVar.bnd(:,1));
	bndRng(isinf(bndRng)) = 1; %Unbounded ones are set to 1
else
	bndRng = ones(size(pt1));
end

%%
% Generate random variable. Undirected
while any(reRoll)
	prop = randn(size(runVar.pt))/n.*runVar.step.*bndRng;
	pt1(reRoll&~runVar.logScale) = runVar.pt(reRoll&~runVar.logScale)+prop(reRoll&~runVar.logScale);
	pt1(reRoll&runVar.logScale) = runVar.pt(reRoll&runVar.logScale).*exp(prop(reRoll&runVar.logScale));
	if ~isempty(runVar.bnd)
		reRoll = pt1<runVar.bnd(:,1)|pt1>runVar.bnd(:,2);
	else
		reRoll = false(size(pt1));
	end
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
pt0_M = runVar.pt;
pt0_M(runVar.logScale) = log(runVar.pt(runVar.logScale));
pt1_M = pt1;
pt1_M(runVar.logScale) = log(pt1(runVar.logScale));
bnd_M = runVar.bnd;
bnd_M(runVar.logScale,:) = log(runVar.bnd(runVar.logScale,:));

pdfBias = biasCalc(pt0_M,pt1_M,bnd_M,n,runVar.step);

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