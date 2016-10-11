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
Names = ['expComp  '];

% Default options
flag = true;

%Parse optional parameters (if any)
for ii = 1:length(varargin)
	if ischar(varargin{ii}) %only enter loop if varargin{ii} is a parameter
		switch lower(deblank(varargin{ii}))
			case lower(deblank(Names(1,:))) %explicit modelling of complex or not
				flag(1) = varargin{ii+1};
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

% Go through each cell and determine what type of data is contained
for ii = 1:length(modelList)
	modelTest = modelList{ii};
	if ischar(modelTest)
		%String
	elseif isa(modelTest,'function_handle')
		%convert function handle to string
		modelNameList{ii} = func2str(modelTest);
	elseif isstruct(modelTest)
		%Previous model structure
		if isfield(modelTest,'name') && isfield(modelTest,'rxnRules') && isfield(modelTest,'modSpc') && isfield(modelTest,'pFit') && isfield(modelTest,'param') && isfield(modelTest,'modComp')
			modelNameList{ii} = 'struct';
		else
			error('parseModel:unexpectedStruct','Unexpected structure detected. Make sure structure passed is a SigMat structure')
		end
	else
		%Unusual structure
		error('parseModel:unexpectedModelType','Unexpected model type detected. Only strings, function handles or model structures allow')
	end
end


%% Kernel
modType = modelType(modelNameList);

if strcmp(modType,'QSSA-m') 
	for ii = 1:length(modelList)
		if isstruct(modelList{ii})
		elseif exist([modelList{ii}],'file')
			if strcmp(modelList{ii}(end-1:end),'.m')
				modelList{ii} = modelList{ii}(1:end-2);
			end
		else
			error('findTC:modelNotFound',['Model file ' modelList{ii} ' not found. Only .xml or .m files accepted'])
		end
	end
	model = parseModelm(modelList,flag);
elseif strcmp(modType,'QSSA-sbml')
	for ii = 1:length(modelname)
		if exist([func2str(modelname{ii})],'file')
			modelname{ii} = func2str(modelname{ii});
			if strcmp(modelname{ii}(end-4:end),'.sbml')
				modelname{ii} = modelname{ii}(1:end-5);
			end
		else
			error('findTC:modelNotFound',['Model file ' func2str(modelname{ii}) ' not found. Only .xml or .m files accepted'])
		end
	end
	model = parseModelSBML(modelname);
else
	error('modelObjective:badModelInput','Invalid model passed. Check inputs')
end
