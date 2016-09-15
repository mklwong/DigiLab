function hist3t(vals,t,nY)
%
% hist3t creates a heatmap of a set of temporal time courses.
%	hist3t(VALS,T,NY) with VALS being a NxM matrix where each M is the
%	number of sets of time profiles to be aggregated. T is a 1xM vector
%   and is the set of time points used for the time course. NY is an
%   integer and denotes the number of equal divisions the Y axis is divided
%   up into.

yVec = linspace(min(min(vals)),max(max(vals)),nY);
valVec = vals(:);
tVec = t(:,ones(1,size(vals,find(size(vals)~=length(t)))));
tVec = tVec(:);
hist3c([tVec valVec],{t,yVec})
shading interp
colormap([ones(1000,1) linspace(1,0,1000)' linspace(1,0,1000)'])