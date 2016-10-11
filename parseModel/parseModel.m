function model = parseModel(modelList,varargin)

%   [model] = parseModel(modelname)
%   
%   Parse model into SigMat structure format to be used for simulation and
%   analysis of the model with:
%		- modelname as the name of the SigMat model as a text string or
%		function handle.
%		- model is the output of the SigMat text format as a structure
%		format.
%
%   The SigMat structure format looks like:
%	model
%		- .name     : original name of the model
%		- .modSpc   : Information for the species in the system
%			- .name   : Name of species
%			- .matVal : matrix value pre-matrix
%			- .pInd   : parameter index pre-matrix
%		- .modComp  : Information for compartments in the system
%			- .name   : Name of compartment
%			- .matVal : matrix value pre-matrix
%			- .pInd   : parameter index pre-matrix
%		- .pFit     : Information for free parameters in the system
%			- .desc : Name of parameter by row
%			- .lim  : limit of parameter by index (corresponding to each
%			          row)
%			- .npar : number of parameters
%			- .grp  : Grouping of parameters and the parameter index each
%			          group is assigned to.
%		- .param    : Parameter matrices used in simulation. 
%			- .name   : Name of matrices
%			- .matVal : matrix value pre-matrix
%			- .pInd   : parameter index pre-matrix
%		- .rxnRules : The rule file used to create the model structure
%
%   Model will contain a number of structures.
%

%% Function options
Names = ['expComp  '
	     'reparse  '];

% Default options
expComp = true;  % Explicit complex
reprse  = true;  % Reparse structures


%Parse optional parameters (if any)
for ii = 1:length(varargin)
	if ischar(varargin{ii}) %only enter loop if varargin{ii} is a parameter
		switch lower(deblank(varargin{ii}))
			case lower(deblank(Names(1,:))) %explicit modelling of complex or not
				expComp = varargin{ii+1};
			case lower(deblank(Names(2,:))) %explicit modelling of complex or not
				reprse = varargin{ii+1};
			case []
				error('Expecting Option String in input');
			otherwise
				error('Non-existent option selected. Check spelling.')
		end
	end
end

% Place modelname in cell as required
if ~iscell(modelList)
	modelList = {modelList};
end
modelNameList = modelList;

% Go through each cell and determine what type of data is contained and if
% they exist
for ii = 1:length(modelList)
	modelTest = modelList{ii};
	% Check what type of data each "modelList" is
	if ischar(modelTest)
		%String: either .m model or .sbml model
		modelNameList{ii} = 'str';
	elseif isa(modelTest,'function_handle')
		%Function handle: ode15s compatible function
		modelNameList{ii} = func2str(modelTest);
		if strcmpi(modelNameList{ii}(1),'@')
			modelNameList{ii} = 'func';
		else
			modelNameList{ii} = 'str';
		end
	elseif isstruct(modelTest)
		%Structure: Pre-parsed model
		if isfield(modelTest,'name') && isfield(modelTest,'rxnRules') && isfield(modelTest,'modSpc') && isfield(modelTest,'pFit') && isfield(modelTest,'param') && isfield(modelTest,'modComp')
			modelNameList{ii} = 'QSSA-m';
		else
			error('parseModel:unexpectedStruct','Unexpected structure detected. Make sure structure passed is a SigMat structure')
		end
	else
		%Other structure: Incompatible
		error('parseModel:unexpectedModelType',['Unexpected model type detected found in model ' num2str(ii) '. Only strings, function handles or model structures allow'])
	end
	
	% For the string inputs, check where the file locations are and work
	% out what type of file it is
	if strcmp(modelNameList{ii},'str')
		dotInd = strfind(modelNameList{ii},'.');
		if isempty(dotInd)
			% If no file extension is defined, we must check for both .m
			% and .sbml
			dotInd = length(modelList{ii})+1;
			% Try .m
			if exist([modelList{ii} '.m'],'file')
				is_m = true;
				modelList{ii} = [modelList{ii} '.m'];
			else
				is_m = false;
			end
			
			if exist([modelList{ii} '.sbml'],'file')
				is_sbml = true;
				modelList{ii} = [modelList{ii} '.sbml'];
			else
				is_sbml = false;
			end
			
			if is_m&&is_sbml
				error('parseModel:ambiguousModel',['Model ' modelList{ii}(1:(dotInd-1)) ' is both an ,.sbml and .m file. Please specify the extension.'])
			end
		end
		dotInd = dotInd(end);
		if strcmp(modelList{ii}(dotInd:end),'.m')
			modelNameList{ii} = 'QSSA-m';
		elseif strcmp(modelList{ii}(dotInd:end),'.sbml')
			modelNameList{ii} = 'QSSA-sbml';
		else
			error('parseModel:modelNotFound',['Model file ' modelList{ii} ' not found. Only .sbml or .m files accepted'])
		end
	end
end


%%%%%%%%%%%%%%%%%%
%% Enter Kernel %%
%%%%%%%%%%%%%%%%%%
% Check for case where there's only one structure and the switch is to not
% reparse it
if (length(modelList)==1 && isstruct(modelList{1})) && ~reprse
	model = modelList{1};
	return
end

% Compile each model into their reaction list form
modelRaw = cell(1,length(modelList));
for ii = 1:length(modelList)
	if strcmp(modelNameList{ii},'QSSA-m') 
		modelRaw{ii} = parseModelm(modelList{ii});
	elseif strcmp(modelNameList{ii},'QSSA-sbml')
		modelRaw{ii} = parseModelSBML(modelList{ii});
	else
		error('parseModel:badModelInput','Invalid model passed. Check inputs')
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Merge raw models and eliminate duplicate Entries %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:length(modelRaw)
	% Load the following pieces of data
	modName  = modelRaw{ii}.name;
	modComp  = modelRaw{ii}.modComp;
	modSpc   = modelRaw{ii}.modSpc;
	Bnd      = modelRaw{ii}.bnd;
	rxn      = modelRaw{ii}.rxn;
	rxnRules = modelRaw{ii}.rxnRules;
	spcMode  = modelRaw{ii}.spcMode;
	
	if ii == 1 % First model loaded is the "master model"
		modName_all  = modName;
		modComp_all  = modComp;
		modSpc_all   = modSpc;
		Bnd_all      = Bnd;
		rxn_all      = rxn;
		rxnRules_all = rxnRules;
		spcMode_all  = spcMode;
	else       % Subsequent models must be checked against "master model" for duplicates
		modName_all = [modName_all '|' modName];
		%% Check for duplicates in compartments
		% Check for clashes in parameter boundaries
		[jj,kk] = ismember(modComp(:,1),modComp_all(:,1));
		for mm = 1:length(kk)
			if jj(mm) % Verify if entries are the same
				misMatches = modComp{mm,2}~=modComp_all{kk(mm),2}; 
				if misMatches % Verify if not the same
					[newVal,strOut] = compareParam(modComp{mm,2},modComp_all{kk(mm),2});
					if isempty(newVal)
						error('parseModelm:UnresolvableParamClash',['Unresolvable clash found in ' modComp{mm,1} ' in model: ' modelname{ii}  ', ' stOut])
					else
						fprintf(['Mismatch in ' modComp{mm,1} ': ' strOut '\n'])
					end
					modComp_all{kk(mm),2} = newVal;
				end
			end
		end
		modComp(jj,:) = [];                   % Remove matching duplicates
		modComp_all = [modComp_all;modComp];  % Append to master model
		
		%% Check for duplicates in species
		[jj,kk] = ismember(modSpc(:,1),modSpc_all(:,1));
		for mm = 1:length(kk)
			if jj(mm) % Verify if entries are the same
				misMatches = modSpc{mm,3}~=modSpc_all{kk(mm),3}; 
				if misMatches % Verify if not the same
					[newVal,strOut] = compareParam(modSpc{mm,3},modSpc_all{kk(mm),3});
					if isempty(newVal)
						error('parseModelm:UnresolvableParamClash',['Unresolvable clash found in ' modSpc{mm,1} ' in model: ' modelname{ii}  ', ' stOut])
					else
						fprintf(['Mismatch in ' modSpc{mm,1} ': ' strOut '\n'])
					end
					modSpc_all{kk(mm),3} = newVal;
				end
			end
		end
		
		% Find mismatches in compartment name assigned to species
		jjIndx = find(jj);
		misMatches = ~ismember(modSpc(jjIndx,2),modSpc_all(kk(jj),2)); % Identify where mismatches are
		if any(misMatches)
			misMatchNames = modSpc(jjIndx(misMatches),1);
			misMatchNames(:,2) = {', '};
			misMatchNames(end,2) = {''};
			misMatchNames = misMatchNames';
			error(['The following species have conflicting compartments: ' horzcat(misMatchNames{:}) ', in model: ' modelname{ii} '. Please fix and ensure they are consistent.'])
		end
		
		modSpc(jj,:) = [];  % Remove matching duplicates
		modSpc_all = [modSpc_all;modSpc];
		
		%% Check for conflicts in default boundary
		bndNames = fieldnames(Bnd);
		bndAllNames = fieldnames(Bnd_all);
		for jj = 1:length(bndAllNames)
			[isMemb,membInd] = ismember(bndAllNames(jj),bndNames);
			if isMemb
				bndEqual = Bnd_all.(bndAllNames{jj}) == Bnd.(bndNames{membInd});
				if ~all(bndEqual)
					fprintf(['The following boundary: ' bndAllNames{jj} ', is inconsistent. Boundary encompassing both are set.\n'])
					Bnd_all.(bndAllNames{jj}) = [min([Bnd_all.(bndAllNames{jj})(1) Bnd.(bndNames{membInd})(1)]) max([Bnd_all.(bndAllNames{jj})(2) Bnd.(bndNames{membInd})(2)])];
				end
			end
		end
		% Check for boundaries that are in new model but not in aggregated
		% model
		isMemb = find(~ismember(bndNames,bndAllNames));
		for jj = 1:length(isMemb)
			Bnd_all.(bndNames(jj)) = Bnd.(bndNames(jj));
		end
		
		%% Merge reaction list
		rxn_all = [rxn_all rxn];

		%% Check for conflicting reaction rules abd
		if ~isequal(rxnRules_all,rxnRules)
			error(['Models to be combined are based on different rulesets. This is not allowed. Merged models must be based on the same ruleset. Please correct this before trying again.'])
		end
		if ~strcmp(spcMode_all,spcMode)
			error(['Models to be combined are based on different species mode. This is not allowed. Merged models must all either all be concentration or absolute amount based. Please correct this before trying again.'])
		end
	end
end

% Change label of global model descriptions
modComp  = modComp_all;
modSpc   = modSpc_all;
Bnd      = Bnd_all;
rxn      = rxn_all;

model.name    = modName_all;
%Save basic model attributes
if ~exist('rateLaw','var')
    model.rxnRules = rxnRules_all;
else
    if ischar(rateLaw)
        rateLaw = str2func(rateLaw);
    end
    model.rxnRules = rateLaw;
end

%Check that model rule file exists
if ~exist(func2str(model.rxnRules),'file')
    error('parseModelm:rxnRuleFileNotFound','Specified rate law file not found. Please check the containing folder is in the MATLAB Path')
end
param = model.rxnRules('ini');

%Convert all compartment and species arrays into structures
modSpc = cell2struct(modSpc,{'name','comp','matVal'},2);
modComp = cell2struct(modComp,{'name','matVal'},2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialise Output Structs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initalise parameter vector
pFit.desc = cell(100,1); pFit.lim  = nan(100,2); pFit.npar = 0; pFit.grp = zeros(0,2);
paramGrp = [-1,-1];

for ii = 1:length(param) % Create pInd for all params.
	param(ii).pInd = nan(size(param(ii).matVal));
end
param = expandTens(param);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cycle over compartments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:size(modComp)
	% Process parameter
    [modComp(ii),pInd,pFit] = testPar(modComp(ii),pFit);
	
	if ~isnan(pInd(2))
        parDesc = [num2str(pInd(2)) ' | Comp : ' modComp(ii).name];
		if ~isempty(pFit.desc{pInd(2)})
			parDesc(1:length(num2str(pInd(2)))) = ' ';
			oldLength = size(pFit.desc{pInd(2)},2);
			newLength = length(parDesc);
			pad(1:max(oldLength,newLength)) = ' ';
		else
			newLength = 0;
			oldLength = 0;
			pad = ' ';
		end
        pFit.desc{pInd(2)} = [pFit.desc{pInd(2)} pad(ones(size(pFit.desc{pInd(2)},1),1),1:(newLength-oldLength));
			                      parDesc        pad(1:(oldLength-newLength))];
		if isnan(pFit.lim(pInd(2),1)) %Insert default parameter range if custom value not defined
			pFit.lim(pInd(2),:) = Bnd.Comp;
		end
	end
	
    modComp(ii).pInd = pInd(2);
end		
% Constract modComp
tmpComp.name   = vertcat({modComp.name})';
tmpComp.matVal = vertcat(modComp.matVal);
tmpComp.pInd   = vertcat(modComp.pInd);

modComp = tmpComp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cycle over list of species
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:size(modSpc)
	% Process parameter
    [modSpc(ii),pInd,pFit] = testPar(modSpc(ii),pFit);
	
    if ~isnan(pInd(3))
        parDesc = [num2str(pInd(3)) ' | Conc : ' modSpc(ii).name];
		if ~isempty(pFit.desc{pInd(3)})
			parDesc(1:length(num2str(pInd(3)))) = ' ';
			oldLength = size(pFit.desc{pInd(3)},2);
			newLength = length(parDesc);
			pad(1:max(oldLength,newLength)) = ' ';
		else
			newLength = 0;
			oldLength = 0;
			pad = ' ';
		end
        pFit.desc{pInd(3)} = [pFit.desc{pInd(3)} pad(ones(size(pFit.desc{pInd(3)},1),1),1:(newLength-oldLength));
			                      parDesc        pad(1:(oldLength-newLength))];
		if isnan(pFit.lim(pInd(3),1)) %Insert default parameter range if custom value not defined
			pFit.lim(pInd(3),:) = Bnd.Conc;
		end
    end
    
    modSpc(ii).pInd = pInd(3);
end

% Contract modComp
tmpSpc.name   = vertcat({modSpc.name});
[~,tmpSpc.comp] = ismember(vertcat({modSpc.comp}),modComp.name);
tmpSpc.comp = tmpSpc.comp';
tmpSpc.comp = [tmpSpc.comp tmpSpc.comp*0]; %expand complex vector by 1 
                                           %column. For complex species 
										   %which have no defined reference
										   %species yet.
tmpSpc.matVal = vertcat(modSpc.matVal);
tmpSpc.pInd   = vertcat(modSpc.pInd);

modSpc = tmpSpc;

%Substitute in compartment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cycle over list of Reactions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

paramNames = {param.name}; %get list of parameter names
for ii=1:length(rxn)
	
	%% Determine if parameters are free or not
    [rxnPar,~,pFit,rxnTest] = testPar(rxn(ii),pFit);
	
	%% Turn reactions into maths using reaction rules
	[~,pMatBuilder,~,rxnDesc] = model.rxnRules('rxnRules',rxnTest,modSpc,expComp,ii); %Make all param values imaginary then use the matVal is produces to build parameter matrix
	[reqParam,matVal,modSpc] = model.rxnRules('rxnRules',rxnPar,modSpc,expComp,ii);
    
	% Find free parameters and their corresponding descriptions
	rxnFields = fieldnames(rxnTest);
	for jj = 1:length(rxnFields)
		if isnumeric(rxnTest.(rxnFields{jj})) && ~isempty(rxnTest.(rxnFields{jj}))
			if imag(rxnTest.(rxnFields{jj})) ~= 0
				pInd = imag(rxnTest.(rxnFields{jj}));
				if isempty(pFit.desc{pInd})
					pFit.desc{pInd} = [num2str(pInd) ' | ' rxnDesc.(rxnFields{jj}){2}];
				else
					% Appending label to already assigned parameter
					clear oldPad newPad
					dummyStr = num2str(pInd);
					dummyStr(1:end) = ' ';
					newStr = [dummyStr ' | ' rxnDesc.(rxnFields{jj}){2}];
					wNew = length(newStr);
					[hOld,wOld] = size(pFit.desc{pInd});
					oldPad(1:hOld,1:(max(0,wNew-wOld))) = ' '; 
					newPad(1,1:(max(0,wOld-wNew))) = ' '; 
					pFit.desc{pInd} = [pFit.desc{pInd} oldPad;newStr newPad];
				end
				if isnan(pFit.lim(pInd,:)) %is boundary still unfilled, then a default is requested
					pFit.lim(pInd,:) = Bnd.(rxnDesc.(rxnFields{jj}){1});
				end
			end
		end
	end
    
    % Loop through pre-matrices that require new additions and build the
    % pre-matrix
    for jj = 1:length(reqParam);
		%Make param builder
		pIndImag = imag(pMatBuilder{jj})==0;
		pMatBuilder{jj}(~pIndImag) = abs(imag(pMatBuilder{jj}(~pIndImag)));
		pMatBuilder{jj}(pIndImag) = NaN;
		
        %Append both into their required pre-matrices
        [~,paramInd,~] = intersect(paramNames,reqParam{jj});
        matInd = find(isnan(param(paramInd).matVal(:,1)),1,'first');
		param(paramInd).matVal(matInd:(matInd+size(matVal{jj},1)-1),:) = matVal{jj};
        param(paramInd).pInd(matInd:(matInd+size(matVal{jj},1)-1),:)   = pMatBuilder{jj};
	end
end

% remove NaN rows from generated tensors
for ii = 1:length(param)
	paramRmIndx = isnan(param(ii).matVal(:,1));
	param(ii) = contractTens(param(ii),paramRmIndx);
end
pFitRmIndx = isnan(pFit.lim(:,1));
pFit  = contractTens(pFit,pFitRmIndx);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Determine data and model relation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rmXData = [];  % List of dataSpcs to remove due to not having all species in the model.
% 
% % Loop through experimental species
% for ii = 1:size(dataSpc,1)
% 	% Loop through attached model state
% 	dataSpcTmp = zeros(size(dataSpc{ii,2}));
% 	for jj = 1:length(dataSpc{ii,2})
% 		
% 		% Determine if complex of model species needs to be included.
% 		% model species with start in them will have enzyme-substrate
% 		% complex included (flat -1).
% 		if strcmp(dataSpc{ii,2}{jj}(1),'*')
% 			incComp = -1;
% 			dataSpc{ii,2}{jj}(1) = [];
% 		else
% 			incComp = 1;
% 		end
% 		
% 		% Get species name index
% 		[~,~,indx]  = intersect(upper(dataSpc{ii,2}{jj}),upper(conc.name));
% 		
% 		% If required model species not in the simulation. Print warning and then exclude
%         if isempty(indx)
% 			storeError(model,[],[],[],['The state ' dataSpc{ii,2}{jj} ' required as an experimental equivalent state not found. This experimental state will be excluded from curve fitting'])
%             rmXData = ii;
% 			break
% 		end
% 		
% 		% Store index of model species including flag of whether to include
% 		% enzyme-substrate complex or not.
% 		dataSpcTmp(jj) = incComp*indx;
% 	end
% 	dataSpc{ii,2} = dataSpcTmp;
% end
% dataSpc(rmXData,:) = [];
% pFit.sim2dat = dataSpc;

%% Compile output
model.modComp = modComp;
model.modSpc = modSpc;
model.spcMode = spcMode_all;
model.param = param;
model.pFit = pFit;
model.raw.rxn = rxn;
model.raw.bnd = Bnd;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%End Main Function%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
function param = expandTens(param)
	paramFields = fieldnames(param);
	for ii = 1:length(param)
		for jj = 1:length(paramFields)
			if ~ismember(paramFields{jj},{'name'}) % Fields to not expand on contract
				arrayWidth = size(param(ii).(paramFields{jj}),2);
				param(ii).(paramFields{jj}) = [param(ii).(paramFields{jj}) ; nan(200,arrayWidth)];
			end
		end
	end
end

function param = contractTens(param,ind)
	paramFields = fieldnames(param);
	for jj = 1:length(paramFields)
		if ~ismember(paramFields{jj},{'name','npar','grp'}) % Fields to not expand or contract
			tmpTens = param.(paramFields{jj});
			tmpTens(ind,:) = [];
			param.(paramFields{jj}) = tmpTens;
		end
	end
end

function [valOut,strOut] = compareParam(val1,val2)
	% Make val1 the longer parameter description
	if length(val2) > length(val1)
		val = val2;
		val2 = val1;
		val1 = val;
	end
	
	% Compare the parameter descriptions
	if  length(val1)==4 && (length(val1)==length(val2)) %Both have custom grouped boundaries
		valOut = val1;
		valOut(3:4) = [min(valOut(3),val2(3)) max(valOut(4),val2(4))];
		strOut = 'Unequal boundary found. Largest custom boundary chosen.';
	elseif length(val1)==3 && (length(val1)==length(val2)) %Both have custom boundaries
		valOut = val1;
		valOut(2:3) = [min(valOut(2),val2(2)) max(valOut(3),val2(3))];
		strOut = 'Unequal boundary found. Largest custom boundary chosen.';
	elseif length(val1)==2 && (length(val1)==length(val2)) %Both have custom boundaries
		valOut = [];
		strOut = 'Clash found between dependence description. Please check.';
	elseif length(val1)==1 && (length(val1)==length(val2)) %Both have custom boundaries
		valOut = NaN;
		strOut = 'Parameter set as free parameter.';
	elseif length(val1)>length(val2)
		valOut = val1;
		if length(val2) == 3
			valOut(3:4) = [min(valOut(3),val2(2)) max(valOut(4),val2(3))];
			strOut = 'Reference parameter set. Largest custom boundary chosen.';
		else
			strOut = 'Unequal parameter description found. More complex option chosen.';
		end
	end
end