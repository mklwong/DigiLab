function mkRules(name,clean)
%
% mkRule(name,clean)
%	Creates a new rule set for kinetic modelling from the model file 
%	template with filename "name". Opens the file on run end.
%    
%   The clean option allows creation of a ruleset file with no creation
%   guide.

existFile = exist([name '.m'],'file');
while existFile
	rePlc = input('Ruleset already exists. Replace (y/n)? ','s');
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
	error('mkRules: How did I get here? Please contact developer!!!!!')
end
curDir(rmIndx+1:end) = [];

if nargin == 2
    curDir = [curDir 'rulesTemplate_clean.m'];
else
    curDir = [curDir 'rulesTemplate.m'];
end
copyfile(curDir,[name '.m'])

open([name '.m'])