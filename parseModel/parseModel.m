function mod = parseModel(model,p)

% First determine if ODE model or not. Currently done using error checking
% because not sure how to test if model is linked to a function or a
% script. The former implies an ode15 file while latter is either an SBML
% model or matlab QSSA model.

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
			mod = parseModelm(model);
		elseif exist([model '.m'],'file')
			mod = parseModelm(model);
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
	
	% Impose passed parameter on reaction parameters, else use default in
	% tensor
	if nargin == 2
		if isrow(p)
			p = p';
		end
		mod.k0.tens(mod.k0.pInd(:,1),2) = (p(abs(mod.k0.pInd(:,2))).*mod.k0.sign.*mod.k0.factor).^sign(mod.k0.pInd(:,2));
		mod.k1.tens(mod.k1.pInd(:,1),3) = (p(abs(mod.k1.pInd(:,2))).*mod.k1.sign.*mod.k1.factor).^sign(mod.k1.pInd(:,2));
		mod.k2.tens(mod.k2.pInd(:,1),4) = (p(abs(mod.k2.pInd(:,2))).*mod.k2.sign.*mod.k2.factor).^sign(mod.k2.pInd(:,2));
		mod.G.tens(mod.G.pInd(:,1),4)  = (p(abs(mod.G.pInd(:,2))).*mod.G.sign.*mod.G.factor).^sign(mod.G.pInd(:,2));
		mod.x.tens(~isnan(mod.x.pInd)) = p(mod.x.pInd(~isnan(mod.x.pInd))).*mod.x.factor(~isnan(mod.x.pInd));
		%EXCEPTION: Checking for NaN tensor values (should all be done)
		if [find(isnan(mod.k0.tens(:,end)))' find(isnan(mod.k1.tens(:,end)))' find(isnan(mod.k2.tens(:,end)))' find(isnan(mod.G.tens(:,end)))']
			error('parseModel:UndefinedFreeTensorValues','Some tensor values are still NaN. Free parameters were undefined.')
		end
	end
else
	error('modelObjective:badModelInput','Invalid model passed. Check inputs')
end
