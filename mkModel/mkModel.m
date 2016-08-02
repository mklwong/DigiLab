function mkModel(name,clean)
%
% mkModel(name,clean)
%	Creates a new model for kinetic modelling from the model file template
%	with filename "name". Opens the file on run end.
%    
%   The clean option allows creation of a model file with no creation
%   guide.

existFile = exist([name '.m'],'file');
while existFile
	rePlc = input('Model already exists. Replace (y/n)? ','s');
	if strcmpi(rePlc,'y')
		rePlc = input('No backup will be made. Are you sure (y/n)? ','s');
	end
	existFile = ~(strcmpi(rePlc,'y') || strcmpi(rePlc,'n'));
	if existFile
		fprintf('Invalid entry. Please try again.\n')
	end
    if strcmpi(rePlc,'n')
		fprintf('Opening existing file\n')
        open([name '.m'])
        return
    end
end

curDir = which(mfilename);

if ~isempty(strfind(curDir,'\'))
	rmIndx = max(strfind(curDir,'\'));
elseif ~isempty(strfind(curDir,'/'))
	rmIndx = max(strfind(curDir,'/'));
else
	error('mkModel: How did I get here? Please contact developer!!!!!')
end
curDir(rmIndx+1:end) = [];

if nargin == 2
    curDir = [curDir 'modelTemplate_clean.m'];
else
    curDir = [curDir 'modelTemplate.m'];
end
copyfile(curDir,[name '.m'])

open([name '.m'])