function model = parseModel(modelname,varargin)

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
Names = ['expComp  ';
         'modelRule'];

% Default options
flag = true;
modelRules = @odeKinetic;

%Parse optional parameters (if any)
for ii = 1:length(varargin)
	if ischar(varargin{ii}) %only enter loop if varargin{ii} is a parameter
		switch lower(deblank(varargin{ii}))
			case lower(deblank(Names(1,:))) %explicit modelling of complex or not
				flag(1) = varargin{ii+1};
			case lower(deblank(Names(2,:))) %model Rules
				modelRules = varargin{ii+1};
			case []
				error('Expecting Option String in input');
			otherwise
				error('Non-existent option selected. Check spelling.')
		end
	end
end

if iscell(modelname)
	modelnameCell = modelname;
else
	modelnameCell{1} = modelname;
end

% Isolate model name and check that each are not pre-parsed structures
for ii = 1:length(modelnameCell)
	modelname = modelnameCell{ii};
	if ischar(modelname)
		%convert function handle to string
		modelname = str2func(modelname);
	elseif isstruct(modelname)
		if length(modelnameCell)>1
			error('parseModel:cannotCombineParsedModel','Currently proram not capable of combining parsed models. Feature will be incorporated in future versions')
		end
		model = modelname;
		if isfield(model,'name') && isfield(model,'rxnRules') && isfield(model,'modSpc') && isfield(model,'pFit') && isfield(model,'param') && isfield(model,'modComp')
			%Model structure detected
			modelname = str2func(model.name);
		else
			%Unusual structure
			error('parseModel:unexpectedModelType','Unexpected model type detected. Only strings, function handles or model structures allow')
		end
	end
	modelnameCell{ii} = modelname;
end
modelname = modelnameCell;

%% Kernel
modType = modelType(modelname);

if strcmp(modType,'ode15s')
	if exist('p','var')
		model = @(t,x) model(t,x,p);
	else
		model = @(t,x,p) model(t,x,p);
	end
elseif strcmp(modType,'QSSA-m') 
	for ii = 1:length(modelname)
		modelname{ii} = func2str(modelname{ii});
		if exist([modelname{ii}],'file')
			if strcmp(modelname{ii}(end-1:end),'.m')
				modelname{ii} = modelname{ii}(1:end-2);
			end
		else
			error('findTC:modelNotFound','Model file not found. Only .xml or .m files accepted')
		end
	end
	model = parseModelm(modelname,modelRules,flag);
elseif strcmp(modType,'QSSA-sbml')
	for ii = 1:length(modelname)
		if exist([modelname{ii}],'file')
			modelname{ii} = func2str(modelname{ii});
			if strcmp(modelname{ii}(end-4:end),'.sbml')
				modelname{ii} = modelname{ii}(1:end-5);
			end
		else
			error('findTC:modelNotFound','Model file not found. Only .xml or .m files accepted')
		end
	end
	model = parseModelSBML(modelname);
else
	error('modelObjective:badModelInput','Invalid model passed. Check inputs')
end
