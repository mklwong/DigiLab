function model = parseModel(model,p,varargin)

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
rxnRules = @odeKinetic;
	 
%Parse optional parameters (if any)
for ii = 1:length(varargin)
	if ischar(varargin{ii}) %only enter loop if varargin{ii} is a parameter
		switch lower(deblank(varargin{ii}))
			case lower(deblank(Names(1,:)))
				expComp = varargin{ii+1};
			case lower(deblank(Names(2,:)))
				rxnRules = varargin{ii+1};
			case []
				error('Expecting Option String in input');
			otherwise
				error('Non-existent option selected. Check spelling.')
		end
	end
end

%% Kernel

modType = modelType(model);

if strcmp(modType,'ode15s')
	odeFile = model;
	if nargin == 2
		mod = @(t,x) odeFile(t,x,p);
	else
		mod = @(t,x,p) odeFile(t,x,p);
	end
elseif strcmp(modType,'QSSA')
	if isa(model,'function_handle')
		model = func2str(model);
	end

	% Parsing models
	if ischar(model)
		if strcmp(model((end-1):end),'.m')
			mod = parseModelm(model,rxnRules,expComp,p);
		elseif exist([model '.m'],'file')
			mod = parseModelm(model,rxnRules,expComp,p);
		elseif strcmp(model((end-3):end),'.xml')
			mod = parseModelSBML(model);
		elseif exist([model '.xml'],'file')
			mod = parseModelSBML(model);
		else
			error('findTC:modelNotFound','Model file not found. Only SBML or .m files accepted')
		end
		mod.name = model;
	elseif isstruct(model)
		mod = model;
	else
		error('findTC:modelClassUnknown','Unable to process model. Check model type')
	end
	
	model = mod;
	
	% Impose passed parameter on reaction parameters, else use default in
	% tensor
	if nargin == 2
		if isrow(p)
			p = p';
		end
		model = rxnRules('insParam',model,p);
	end
else
	error('modelObjective:badModelInput','Invalid model passed. Check inputs')
end
