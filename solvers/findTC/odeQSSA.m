function [t,Y,YComp,model] = odeQSSA(modelRaw,tspan,varargin)
%
% odeQSSA
%
%	Take parameter inputs from findTC.
%	Default 

%% To Implement
% Solver choice in options

%% Code to handle options
% Option names
Names = ['p       '
         'inp     '
         '-r      '
		 'y0      '
		 'odeopts '
	     'errDir  '
		 '-b      '];
     
%Initialise potental options
p    = [];
x0   = [];
tmpInp  = [];
noRamp = false;
errDir = false;
noBasal = false;

detectOpts = true;
%Parse optional parameters (if any)
for ii = 1:length(varargin)
	if detectOpts %only enter loop if varargin{ii} is a parameter
		detectOpts = false;
		switch lower(deblank(varargin{ii}))
			case lower(deblank(Names(1,:)))   %Parameters
				p = varargin{ii+1};
			case lower(deblank(Names(2,:)))   %Input
				tmpInp = varargin{ii+1};
			case lower(deblank(Names(3,:)))   %No initial ramping
				noRamp = true;
				detectOpts = true;
			case lower(deblank(Names(4,:)))   %Initial conditions
				x0 = varargin{ii+1};
			case lower(deblank(Names(5,:)))   %Ode options
				options = varargin{ii+1};
			case lower(deblank(Names(6,:)))   %Set directory for error output
				errDir = varargin{ii+1};
			case lower(deblank(Names(7,:)))   %No basal
				noBasal = true;
			case []
				error('Expecting Option String in input');
			otherwise
				error('Non-existent option selected. Check spelling.')
		end
	else
	        detectOpts = true;
	end
end

%% Compile model if new one inserted
if isrow(p)
	p = p';
end

model = parseModel(modelRaw,p);

% %Correct dimension of x0 and tspan
if isrow(x0)
    x0 = x0';
end
if isrow(tspan)
    tspan = tspan';
end

% Create initial concentration vector
if isempty(x0)
	x0 = model.conc.tens;
elseif size(x0,1) ~= size(model.conc.tens)
	x0(length(model.conc.tens)) = 0;
end

inpConst = zeros(length(x0),1);
inpFun = @(t)zeros(length(x0),1);
%Input values: make all into either function handles or vectors
% This component looks at the experiment-simulation name pair, then
% compares the experiment name with the name given in the 
if iscell(tmpInp) %state name-val pair
	protList = model.conc.name;
	inpFunInd = [];
	inpConstInd = [];
	for ii = 1:size(tmpInp,1)
		if ischar(tmpInp{ii,1})
			[~,stateInd] = intersect(upper(protList),upper(tmpInp{ii,1})); % Match input state
			tmpInp{ii,1} = stateInd; %Replace the name with index
			%Separate spikes and gradual inputs
			if ~isempty(stateInd)
				if isa(tmpInp{ii,2},'function_handle')
					inpFunInd = [inpFunInd ii];
				else
					inpConstInd = [inpConstInd ii];
				end
			end
		end
	end
	
	%Insert constant values
	inpConst(vertcat(tmpInp{inpConstInd,1})) = vertcat(tmpInp{inpConstInd,2});
	inpFunCell = cell(length(x0),1);
	inpFunCell(:) = {@(t)0};
	inpFunCell(tmpInp{inpFunInd,1}) = tmpInp(inpFunInd,2);
	inpFun = @(t) cellfun(@(f) f(t),inpFunCell);
	
elseif isa(tmpInp,'function_handle') %vector of function handles
	[a,b] = size(tmpInp(1));
	if b>1
		if a~=1
			error('findTC:inpfunDimWrong','Dimensions of input function handle is wrong.')
		end
		if b <= length(x0)
			inpFun = @(t) [tmpInp(t)';zeros(length(x0)-b,1)];
		else
			error('findTC:tooManyInpState','Too many input states in vector-val method.')
		end
	elseif b==1
		if a <= length(x0)
			inpFun = @(t) [tmpInp(t);zeros(length(x0)-a,1)];
		else
			error('findTC:tooManyInpState','Too many input states in vector-val method.')
		end
	else
		error('findTC:inpfunDimWrong','Dimensions of input function handle is wrong.')
	end
elseif size(tmpInp,2)==2             %Ind val pair
	inpConst(tmpInp(:,1)) = tmpInp(:,2);
elseif min(size(tmpInp))==1          %vector of spiked final concentration
	[a,b] = size(tmpInp);
	if a == 1
		tmpInp = tmpInp';
		a = b;
	end
	if a > length(x0)
		error('findTC:tooManyInpState','Too many input states in vector-val method.')
	end
	inpConst = [tmpInp;zeros(length(x0)-a,1)];
elseif size(tmpInp,2)>2
	error('odeQSSA:inpArrayDimWrong','Dimension of system input incorrect. Check your inputs')
end
	
%% From here on in, all time is non-dimensionalised, as ode will be run from t = 0 to 1.
%Normalisation factor which makes sure the ramp in will ramp in the full
%amount after the run time.
normInp = @(t) inpFun(t*(tspan(end)-tspan(1))+tspan(1))*(tspan(end)-tspan(1)); %time shift and non-dimensionalise inp;

model = model.rxnRules('nondim',model,tspan,normInp);

%ODE Solver options and warning
if ~exist('options','var')
	options = odeset('relTol',9e-6,'NonNegative',ones(size(x0)));
end
warnstate('error')

try
%% Solving
% Ramping
if ~noRamp
	modelRamp = model;
	if noBasal
		%Do not basal the time course. Set all rate parameters to zero.
		modelRamp.k0 = @(t) (inpConst+x0)*2*normpdf(t,0,0.2);
		modelRamp = model.rxnRules('ramp',modelRamp);
	else
		%Initial basal
		%Enter the single input at the end of the time course
		modelRamp.k0 = @(t) x0*2*normpdf(t,0,0.2)+inpConst*2*normpdf(t,1,1e-6);
	end
	dx_dt = @(t,x) model.rxnRules('dynEqn',t,x,modelRamp);
	[t,Y] = ode15s(dx_dt,[0 (1-1e-6) 1],x0*0,options);
	y0 = Y(end,:)';
	y0(y0<0) = 0;
end

%Run
model.k0 = model.fullSigma;
dx_dt = @(t,x) model.rxnRules('dynEqn',t,x,model);
[t,Y] = ode15s(dx_dt,[0 1],y0,options);

t = t*(tspan(end)-tspan(1))+tspan(1); %Restore to original units
if length(tspan)>2
	Y = interp1(t,Y,tspan);
	t = tspan;
end
YComp  = Y;
catch errMsg
%% Error catching
	YComp = inf(length(tspan),length(x0));
	t = tspan;
	if errDir 
		storeError(modelRaw,x0,p,errMsg,errMsg.message,errDir)
	else
		storeError(modelRaw,x0,p,errMsg,errMsg.message)
	end
end
Y = compDis(model,YComp);      %dissociate complex
warnstate('on') %Switch warnings back to warnings
end

% ODE Solve warning messages. Turn them into errors so the error catcher
% can catch them.
function warnstate(state)
    warning(state,'MATLAB:illConditionedMatrix');
    warning(state,'MATLAB:ode15s:IntegrationTolNotMet');
    warning(state,'MATLAB:nearlySingularMatrix');
    warning(state,'MATLAB:singularMatrix');
end