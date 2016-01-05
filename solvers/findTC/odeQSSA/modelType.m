function modType = modelType(model)

try 
	model([],[],[]); %an ODE model should accept t, x and p
	msg.identifier = '';
catch msg
end

if (strcmp(msg.identifier,'MATLAB:scriptNotAFunction') || strcmp(msg.identifier,'MATLAB:UndefinedFunction')) || ischar(model) || isstruct(model)
	if ischar(model)
		if ~exist(model)
			error('modeType:fileNotFound','Model/file does not exist.')
		end
	end
	modType = 'QSSA';
elseif strcmp(msg.identifier,'MATLAB:badsubscript') || strcmp(msg.identifier,'MATLAB:maxrhs') || strcmp(msg.identifier,'MATLAB:TooManyInputs')
	modType = 'ode15s';
elseif strcmp(msg.identifier,'MATLAB:dispatcher:InexactCaseMatch')
	error(msg.message)
end