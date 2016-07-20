function allResult = combineResults(varargin)

%combineResults merges all MCMC result structures (pts, logP and ptUn)
%   allResult = combineResults(result1,result2,....) All arguments of the
%   function are expected to be result structures. All result then combines
%   these together into a single structure. This function has no other
%   options or arguments.
%
%   A result structure has the following fields:
%       result.pts  - a list of parameter sets where each row is a set
%       result.logP - a list of corresponding -logP value for each
%       parameter set
%       result.ptUn - a list of logical flag that states whether the
%       corresponding parameter set is unique
%       result.T    - Tempering temperature of the result run
%       result.model- Parameter set that produced the best fit

allResult = varargin{1};

for ii = 2:length(varargin)
    tmpResult = varargin{2};
    allResult.pts = [allResult.pts;tmpResult.pts];
    allResult.logP = [allResult.logP;tmpResult.logP];
    allResult.ptUn = [allResult.ptUn;tmpResult.ptUn];
end

allResult.best = allResult.pts(allResult.logP==min(allResult.logP),:);

[allResult.logP,I] = sort(allResult.logP,'ascend');
allResult.pts = allResult.pts(I,:);
allResult.ptUn = allResult.ptUn(I);