function modType = modelType(model)
modType = [];
for ii = 1:length(model)
	if ischar(model(ii))
		model = str2func(model(ii));
	elseif isstruct(model(ii))
		modType = 'QSSA-m';
		continue
	end

	% Test if function handle by trying it as a function handle of 3 inputs.
	try 
		modelTest = model{ii};
		modelTest([],[],[]); %an ODE model should accept t, x and p
		if ~strcmp(modType,'ode15s') || isempty(modType)
			error('modelType:ModelsOfDifferentTypesPassed','Models of different types passed. Algorithm cannot currently combine different types. Please recheck')
		end
		modType = 'ode15s';
		continue
	catch msg
	end
end

if strcmp(msg.identifier,'MATLAB:scriptNotAFunction')
	% A script exists with the name parsed but it is not a function. Likely
	% to be a QSSA with matlab file
	if ~strcmp(modType,'QSSA-m') && ~isempty(modType)
		error('modelType:ModelsOfDifferentTypesPassed','Models of different types passed. Algorithm cannot currently combine different types. Please recheck')
	end
	modType = 'QSSA-m';
elseif strcmp(msg.identifier,'MATLAB:UndefinedFunction')
	% A script exists with the name parsed but it is not a function. Likely
	% to be a QSSA
	if ~strcmp(modType,'QSSA-sbml') && ~isempty(modType)
		error('modelType:ModelsOfDifferentTypesPassed','Models of different types passed. Algorithm cannot currently combine different types. Please recheck')
	end
	modType = 'QSSA-sbml';
elseif strcmp(msg.identifier,'MATLAB:badsubscript') || strcmp(msg.identifier,'MATLAB:maxrhs') || strcmp(msg.identifier,'MATLAB:TooManyInputs')
	fprintf('In correct number of parameters for ode15s function. It must have input arguments (t,x,p) to valid for us in parseModel')
	msg.throw
else
	msg.throw
end