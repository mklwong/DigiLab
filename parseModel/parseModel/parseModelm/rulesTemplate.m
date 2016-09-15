function varargout = [rxnRules](method,varargin)

switch lower(deblank(method))
	
case 'ini'
%<Insert parameters>
varargout = {param};

case 'rxnrules'
[rxn,modSpc,flag,ii] = varargin{:};
%<Insert reaction rules>
varargout = {reqTens,matVal,modSpc,rxn};

case 'compile'
if length(varargin) == 2
	[model,tspan] = varargin{:};
elseif length(varargin) == 1
	model = varargin{:};
	tspan = [0 1]; %no non-dimensionalise
end
%<Insert pre-compilation steps>
varargout = {model};

case 'dyneqn'
[t,x,model] = varargin{:};
%<Insert compilation steps and dynamic equation calculation>
varargout = {dx_dt};
end