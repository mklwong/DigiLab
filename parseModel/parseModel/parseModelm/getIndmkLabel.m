function [pFit,repInd,curInd] = getIndmkLabel(pFit,pVar,Bnd,repInd,text)

% Get Index Make Label
%	Support script for parseModelm and parseModelSBML

%% Get current index
curInd = find(isnan(pFit.lim(:,1)),1,'first'); %First empty element

%% Determine if this parameter is grouped
if pVar(1) < 0 %indicative that this parameter is grouped
	if isempty(repInd)
		indx = [];
	else
		indx = find(repInd(:,1)==pVar(1));
	end
	if ~isempty(indx) %group already exists, add to group
		curInd = repInd(indx,2);
	else			  %group does not exist, create group
		repInd(size(repInd,1)+1,:) = [pVar(1),curInd];
	end
end

%% Change description of parameter
if curInd == find(isnan(pFit.lim(:,1)),1,'first') 
%if curInd creates a new element in p (i.e. not reusing an already used p)
	pFit.desc{curInd,1} = text;
	if length(pVar) == 3
		pFit.lim(curInd,:) = pVar(2:3);
	else
		pFit.lim(curInd,:) = Bnd;
	end
else
	pFit.desc{curInd,1} = [pFit.desc{curInd,1} ' | ' text];
end