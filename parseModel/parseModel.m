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
		mod.k0.tens = addp(p,mod.k0,2);
        mod.k1.tens = addp(p,mod.k1,3);
        mod.k2.tens = addp(p,mod.k2,4);
        mod.G.tens = addp(p,mod.G,4);
		mod.x.tens = addp(p,mod.x,1);
		%EXCEPTION: Checking for NaN tensor values (should all be done)
		if [find(isnan(mod.k0.tens(:,end)))' find(isnan(mod.k1.tens(:,end)))' find(isnan(mod.k2.tens(:,end)))' find(isnan(mod.G.tens(:,end)))']
			error('parseModel:UndefinedFreeTensorValues','Some tensor values are still NaN. Free parameters were undefined.')
		end
	end
else
	error('modelObjective:badModelInput','Invalid model passed. Check inputs')
end

end
function tens = addp(p,tensCat,ind)
    toUpdate = ~isnan(tensCat.pInd);
    maxp = max(tensCat.pInd);
    if maxp > length(p)
        error('parseModel:insufficientParams','Insufficient parameters parsed. Please check parameter input')
    end
    tensCat.tens(toUpdate,ind) = p(tensCat.pInd(toUpdate)).*tensCat.tens(toUpdate,ind);
    tens = tensCat.tens;
end