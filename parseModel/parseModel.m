function model = parseModel(modelname,p,varargin)

% First determine if ODE model or not. Currently done using error checking
% because not sure how to test if model is linked to a function or a
% script. The former implies an ode15 file while latter is either an SBML
% model or matlab QSSA model.

% Options

%% Function options
Names = ['expComp ';
         'paramDef'];

% Default options
expComp = true;
if ischar(modelname)
	model.name = modelname;
	model.rxnRules = @odeKinetic;
elseif isa(modelname,'function_handle')
	model.name = func2str(modelname);
	model.rxnRules = @odeKinetic;
else
	model = modelname;
end
	 
if ~exist('p')
	p = [];
end

%Parse optional parameters (if any)
for ii = 1:length(varargin)
	if ischar(varargin{ii}) %only enter loop if varargin{ii} is a parameter
		switch lower(deblank(varargin{ii}))
			case lower(deblank(Names(1,:)))
				expComp = varargin{ii+1};
			case lower(deblank(Names(2,:)))
				model.rxnRules = varargin{ii+1};
			case []
				error('Expecting Option String in input');
			otherwise
				error('Non-existent option selected. Check spelling.')
		end
	end
end

%% Kernel

modType = modelType(model.name);

if strcmp(modType,'ode15s')
	odeFile = model;
	if nargin == 2
		model = @(t,x) odeFile(t,x,p);
	else
		model = @(t,x,p) odeFile(t,x,p);
	end
elseif strcmp(modType,'QSSA')
	if isstruct(model) %if model is already a structure, then do some cursory checks to see if model is already built. If so then quit
		if isfield(model,'name') && isfield(model,'rxnRules') && isfield(model,'conc') && isfield(model,'pFit') && isfield(model,'param')
			return
		end
	end

	% Parsing models
	if ischar(model.name)
		if strcmp(model.name((end-1):end),'.m')
			model = parseModelm(model,expComp,p);
		elseif exist([model.name '.m'],'file')
			model = parseModelm(model,expComp,p);
		elseif strcmp(model.name((end-3):end),'.xml')
			model = parseModelSBML(model.name);
		elseif exist([model.name '.xml'],'file')
			model = parseModelSBML(model.name);
		else
			error('findTC:modelNotFound','Model file not found. Only SBML or .m files accepted')
		end
	else
		error('findTC:modelClassUnknown','Unable to process model. Check model type')
	end
	
	% Impose passed parameter on reaction parameters, else use default in
	% tensor
	if ~isempty(p)
		if isrow(p)
			p = p';
		end
		model = model.rxnRules('insParam',model,p);
	end
else
	error('modelObjective:badModelInput','Invalid model passed. Check inputs')
end
