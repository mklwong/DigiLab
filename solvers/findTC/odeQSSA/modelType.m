function modType = modelType(model)

if ischar(model)
	model = str2func(model);
elseif isstruct(model)
	modType = 'QSSA-m';
	return
end
% Test if function handle by trying it as a function handle of 3 inputs.
try 
	model([],[],[]); %an ODE model should accept t, x and p
	modType = 'ode15s';
	return
catch msg
end

if strcmp(msg.identifier,'MATLAB:scriptNotAFunction')
	% A script exists with the name parsed but it is not a function. Likely
	% to be a QSSA with matlab file
	modType = 'QSSA-m';
elseif strcmp(msg.identifier,'MATLAB:UndefinedFunction')
	% A script exists with the name parsed but it is not a function. Likely
	% to be a QSSA
	modType = 'QSSA-sbml';
elseif strcmp(msg.identifier,'MATLAB:badsubscript') || strcmp(msg.identifier,'MATLAB:maxrhs') || strcmp(msg.identifier,'MATLAB:TooManyInputs')
	fprintf('In correct number of parameters for ode15s function. It must have input arguments (t,x,p) to valid for us in parseModel')
	msg.throw
else
	msg.throw
end