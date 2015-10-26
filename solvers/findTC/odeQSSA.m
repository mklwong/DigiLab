function [t,Y,YComp] = odeQSSA(model,tspan,varargin)
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
         'r       '
		 'y0      '
		 'odeopts '
		 'rd      '
		 '-b      '];
     
%Initialise potental options
p    = [];
x0   = [];
tmpInp  = [];
rampOnly = false;
rampDebug = false;
ramp = true;
basal = true;

%Parse optional parameters (if any)
for ii = 1:length(varargin)
	if ischar(varargin{ii}) %only enter loop if varargin{ii} is a parameter
		switch lower(deblank(varargin{ii}))
			case lower(deblank(Names(1,:)))   %Parameters
				p = varargin{ii+1};
			case lower(deblank(Names(2,:)))   %Input
				tmpInp = varargin{ii+1};
			case lower(deblank(Names(3,:)))   %No initial ramping
				ramp = false;
			case lower(deblank(Names(4,:)))   %Ramping in initial condition
				rampOnly = true;
				rampDebug = false;
			case lower(deblank(Names(5,:)))
				x0 = varargin{ii+1};
			case lower(deblank(Names(6,:)))
				options = varargin{ii+1};
			case lower(deblank(Names(7,:)))   %Ramping in initial condition and debug (see ramping time course)
				rampOnly = true;
				rampDebug = true;
			case lower(deblank(Names(8,:)))   %No basal (for running faster)
				basal = false;
			case []
				error('Expecting Option String in input');
			otherwise
				error('Non-existent option selected. Check spelling.')
		end
	end
end

%% Compile model if new one inserted
if isempty(p)
	[k0,k1,k2,G,model] = prepareTensor(model);
else
	[k0,k1,k2,G,model] = prepareTensor(model,p);
end

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
% Non-dimensionalise.
k0 = k0*(tspan(end)-tspan(1));
k1 = k1*(tspan(end)-tspan(1));
k2(:,end) = k2(:,end)*(tspan(end)-tspan(1));

normFac = 1/trapz(linspace(0,1,10000),normpdf(linspace(0,1,10000),0,0.2));
%Normalisation factor which makes sure the ramp in will ramp in the full
%amount after the run time.
normInp = @(t) [inpFun(t*(tspan(end)-tspan(1))+tspan(1))*(tspan(end)-tspan(1))]; %time shift and non-dimensionalise inp;

basalSigma = @(t) k0*ones(1,length(t)); 
fullSigma = @(t) normInp(t) + k0*ones(1,length(t)); 

%ODE Solver options and warning
if ~exist('options','var')
	options = odeset('relTol',9e-6,'NonNegative',ones(size(x0)));
end
warnstate('error')

%% Solving
% Ramping
if ramp && basal
	rampSigma = @(t) x0*normFac*normpdf(t,0,0.2);
	dx_dt = @(t,x) dynEqn(t,x,G,rampSigma,k1*0,zeros(0,4));
	[t,Y] = ode45(dx_dt,[0 1],x0*0,options);
	x0 = Y(end,:)';
end

% Solve Basal Condition
if basal && ~rampOnly
	%Run
	dx_dt = @(t,x) dynEqn(t,x,G,basalSigma,k1,k2);
	converge = false;
	ii = 1;
	while ~converge
		[t,Y] = ode15s(dx_dt,[0 2^ii],x0,options);
		ii = ii + 1;
		converge = (sqrt(sum(abs((Y(end,:)-Y(end-1,:))./(Y(end,:)+eps)).^2))/length(Y(end,:)))<1e-4;
		x0 = Y(end,:)';
	end
end

% Ramping
if ramp
	rampSigma = @(t) (inpConst + x0)*normFac*normpdf(t,0,0.2);
	dx_dt = @(t,x) dynEqn(t,x,G,rampSigma,k1*0,zeros(0,4));
	[t,Y] = ode45(dx_dt,[0 1],x0*0,options);
	x0 = Y(end,:);
end

% Solve ODE
if ~rampOnly
	%Run
	dx_dt = @(t,x) dynEqn(t,x,G,fullSigma,k1,k2);
	[t,Y] = ode15s(dx_dt,[0 1],x0,options);
end
t = t*(tspan(end)-tspan(1))+tspan(1); %Restore to original units
if length(tspan)>2
	Y = interp1(t,Y,tspan);
	t = tspan;
end

YComp  = Y;
Y = compDis(model,Y);      %dissociate complex

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