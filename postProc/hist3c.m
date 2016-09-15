function [n,x] = hist3c(varargin)
% coloured 3D histogram using a heapmap to represent the frequency of the
% distribution. Has the same inputs as hist3.
%
% See also: hist3

rmInd = [];

varargin(rmInd) = [];
hist3(varargin{:})
[a,x] = hist3(varargin{:});
n = max(max(a));
axesHandles = get(gca,'child');

if length(axesHandles)>1
	[~,hisHandleIndx] = intersect(get(axesHandles,'Type'),'surface');
else
	hisHandleIndx = 1;
end
set(axesHandles(hisHandleIndx),'FaceColor','interp','CDataMode','auto'); 
view(2)

end