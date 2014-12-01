function model = parseModelSBML(curSBMod)

%Supported extensions
supExt = ['.xml ';
		  '.sbml'];

if nargin == 0
	[curSBMod,pathname] = uigetfile({'*.xml','SBML Files (.xml)'},'Selection SBML models to parse','MultiSelect','on');
end

% Remove file extensions from filename and get checksum
if iscell(curSBMod)
	curSBMod = curSBMod';
	for ii = 1:length(curSBMod)
		tmp = curSBMod{ii,1};
		for jj = 1:size(supExt,1);
			if strcmpi(deblank(supExt(jj,:)),tmp(end-length(deblank(supExt(jj,:)))+1:end))
				curSBMod{ii,1} = tmp(1:(end-length(deblank(supExt(jj,:)))));
				curSBMod{ii,2} = deblank(supExt(jj,:));
				rmAnnot([pathname curSBMod{ii,1} curSBMod{ii,2}]);
				curSBMod{ii,3} = md5([pathname curSBMod{ii,1} '_tmp' curSBMod{ii,2}]);
				break
			end
		end
	end
elseif ischar(curSBMod)
	tmp = curSBMod;
	clear curSBMod 
	for jj = 1:size(supExt,1);
		if strcmpi(deblank(supExt(jj,:)),tmp(end-length(deblank(supExt(jj,:)))+1:end))
			curSBMod{1,1} = tmp(1:(end-length(deblank(supExt(jj,:)))));
			curSBMod{1,2} = deblank(supExt(jj,:));
			rmAnnot([pathname curSBMod{1,1} curSBMod{1,2}]);
			curSBMod{1,3} = md5([pathname curSBMod{1,1} '_tmp' curSBMod{1,2}]);
			break
		end
	end
end

% If compare, then do a comparison, extracting models that no longer fit
% with a previous checksum or models whos SBML files are lost
modDel = [];
if nargin == 1
	[newSBML modDel] = checkSum;
else
	newSBML = curSBMod;
end

for ii = 1:size(modDel)
	
end

for ii = 1:size(newSBML,1)
	% Parse the SBML model and create structs and cells as required

	% Add the SBML file to the list of SBML models that makes up the matlab
	% model
	
end

% Remove the temporary SBML file made (with annotations removed)
for ii = 1:size(curSBMod,1)
	delete([curSBMod{ii,1} '_tmp' curSBMod{ii,2}])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Nested Functions
%% Check Sum
function [newSBML matModDel] = checkSum
	newSBML = {}; %To do list for parsing SBML models
	matModSum = matMod.sbmlFiles;
	matModDel = [];
	for jjj = 1:size(curSBMod,1)
		%Match file name of models in Matlab model with passed SBML models
		indx = find(cellfun(@(x)strcmpi(x,curSBMod{jjj,1}),matModSum(:,1)));
		if length(indx) > 1
			error('parseModelSBML:matModelCorrupted:Duplicate','Matlab model corrupted. Please recompile.')
		end
		%~~file name match
		if ~isempty(indx)     
			%~~model checksum not match
			if ~strcmpi(curSBMod{jjj,3},matModSum{indx,2}) 
				update = [];
				while isempty(update)
					update = questdlg(['The updated SBML model ' curSBMod{jjj,1} ' has been chosen. Update in the matlab model?'],'Update equivalent Model in Matlab','Yes','No','Yes');
				end
				%~~Matlab model prioritised. Do nothing
				
				%~~SBML model prioritised. Mark model indx in Matlab Model for removal. Add SBML file to "to do" list.
				if strcmp(update,'yes')  
					newSBML(end+1,:) = curSBMod(jjj,:);
					matModDel = [matModSum indx];
				end
			%~~model checksum matches. Do nothing.
			end
		%~~file name not match
		else 
			%Search for matching checkSum
			indx = find(cellfun(@(x)strcmpi(x,curSBMod{jjj,3}),matModSum{:,2}));
			%~~Matching checkSum
			if ~isempty(indx) 
				misnamedMatch = [];
				while isempty(misnamedMatch)
					misnamedMatch = questdlg(['The SBML model ''' curSBMod{jjj,1} ''' was found in the matlab model as ''' matModSum{indx,1} ''' What action would you like to take?'],['Rename Models as ' curSBMod{jjj,1}],['Rename Models as ' matModSum{indx,1}],'Leave as is.','Leave as is.');
				end
				if strcmpi(misnamedMatch,['Rename Models as' curSBMod{jjj,1}])
					matModSum{indx,1} = curSBMod{jjj,1};
				elseif strcmpi(misnamedMatch,['Rename Models as' matModSum{indx,1}])
					movefile([pathname curSBMod{jjj,1} curSBMod{jjj,2}],[matModSum{indx,1} curSBMod{jjj,2}])
					movefile([pathname curSBMod{jjj,1} '_tmp' curSBMod{jjj,2}],[matModSum{indx,1} '_tmp' curSBMod{jjj,2}])
					curSBMod{jjj,1} = matModSum{indx,1};
				end
			%~~No Matching checkSum
			else  
				%Add to "to do" list
				newSBML(end+1,:) = curSBMod(jjj,:);
			end
		end
	end

end

end

