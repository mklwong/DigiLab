function mkModel(name)
%
% mkModel(name)
%	Creates a new model for kinetic modelling from the model file template
%	with filename "name". Opens the file on run end.

if exist([name '.m'],'file')
	fprintf('Model already exists. Opening file.\n')
else
	curDir = which(mfilename);
	strfind(curDir,'\');
	rmIndx = max(strfind(curDir,'\'));
	curDir(rmIndx+1:end) = [];
	curDir = [curDir 'template.m'];
	copyfile(curDir,[name '.m'])
end
open([name '.m'])