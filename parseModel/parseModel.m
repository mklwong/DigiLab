function model = parseModel(modelname,varargin)

% First determine if ODE model or not. Currently done using error checking
% because not sure how to test if model is linked to a function or a
% script. The former implies an ode15 file while latter is either an SBML
% model or matlab QSSA model.

% Options

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

% Isolate model name
if ischar(modelname)
	%convert function handle to string
	modelname = str2func(modelname);
elseif isstruct(modelname)
	model = modelname;
	if isfield(model,'name') && isfield(model,'rxnRules') && isfield(model,'modSpc') && isfield(model,'pFit') && isfield(model,'param') && isfield(model,'modComp')
		%Model structure detected
		modelname = str2func(model.name);
	else
		%Unusual structure
		error('parseModel:unexpectedModelType','Unexpected model type detected. Only strings, function handles or model structures allow')
	end
end
	 
%% Kernel

modType = modelType(modelname);

if strcmp(modType,'ode15s')
	if exist('p','var')
		model = @(t,x) model(t,x,p);
	else
		model = @(t,x,p) model(t,x,p);
	end
elseif strcmp(modType,'QSSA-m') || strcmp(modType,'QSSA-sbml')
	
	modelname = func2str(modelname);
	if ~exist('model','var')
		% If model not created, parsing model
		if strcmp(modelname((end-1):end),'.m')			% check for file extension
			model = parseModelm(modelname,modelRules,flag);
		elseif strcmp(modelname((end-3):end),'.xml')   % check for file extension
			model = parseModelSBML(modelname);
		elseif exist([modelname '.m'],'file')    % Test if .m file exists
			model = parseModelm(modelname,modelRules,flag);
		elseif exist([modelname '.xml'],'file') % Test if .xml file exists
			model = parseModelSBML(modelname);
		else
			error('findTC:modelNotFound','Model file not found. Only .xml or .m files accepted')
		end
	end
else
	error('modelObjective:badModelInput','Invalid model passed. Check inputs')
end
