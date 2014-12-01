function [y0,y0c,status] = findSS(model,dose,state,varargin)
% findSS(model,dose,state,varargin)
%	required inputs are:
%		- model : either 
%				* ode15s model (passed as function handle), 
%				* QSSA .m script file (passed as function handle or string)
%				* SBML .m file (passed as function handle or string)
%		- dose  : 1xn vector of dose to spike in. Note doses are ALWAYS
%				  injected instantaneously at the start (i.e. step profile).
%		- states: 1xn cell array of either protein names or integers per cell,
%				  corresponding to the state associated with the dose
%				  marked 
%			e.g. dose = [1 ; states = {'insulin';    OR states =  {12;
%						 10;           'PDGF'   ;                  15;
%						 3];           'VEGF'   };                 16};
%				 this means [insulin] = 1; [PDGF] = 10 and [VEFG] = 3;
%			     with indicies ind(insulin) = 12, ind(PDGF) = 15 and ind(VEGF) = 16
%
%				
%	Optional inputs are (pass using options:value pairs)
%		- tEnd   : the time span for each steady state run
%		- x0     : initial condition (excluding the dose input). If not
%			  	   passed will search from ode file
%		- p      : vector of parameters. If not passed will search from ode
%				   file.
%		- odeopts: Only relevant for ODE15s models. The odeopts to be used
%		           when using the ode solver. Defaults are as ode function
%		           defaults
%		- tol    : relative difference between initial and final state of
%				   each run before steady state is said ot be achieved.
%				   Defaults at 1e-3;
%
% Status codes:
% -1 = steady state produces inf
% -2 = steady state run is taking too long

%Defaults
tspan = [0 1800];
varargpass = cell(0);
tol = 1e-3;

%% Parse optional parameters (if any)
% Options
Names = ['tEnd   '
         'p      '
		 'x0     '
		 'tol    '
		 'odeopts'
		 '-b     '];

odeopts = [];
p = [];
basal = true;

% Processing of options
for ii = 1:length(varargin)
	if ischar(varargin{ii}) %only enter loop if varargin{ii} is a parameter
		switch lower(deblank(varargin{ii}))
			case lower(deblank(Names(1,:)))   %tspan
				tspan(2) = varargin{ii+1};
			case lower(deblank(Names(2,:)))   %Parameters
				varargpass{end+1} = varargin{ii};
				varargpass{end+1} = varargin{ii+1};
				varargOde15s{3} = varargin{ii+1};
				p = varargin{ii+1};
			case lower(deblank(Names(3,:)))   %Initial condition
				varargpass{end+1} = varargin{ii};
				varargpass{end+1} = varargin{ii+1};
				varargOde15s{1} = varargin{ii+1};
				x0Indx = ii+1;
			case lower(deblank(Names(4,:)))   %Criterion for reaching steady state
				tol = varargin{ii+1};
			case lower(deblank(Names(5,:)))   %odeopts
				odeopts = varargin{ii+1};
				varargOde15s{2} = odeopts;
			case lower(deblank(Names(6,:)))   %No basal
				basal = false;
			case []
				error('findSS:expectingVal',['Expecting Option String after ''' varargin{ii} ''' input']);
			otherwise
				error('findSS:unknownOption',['Non-existent option: ''' varargin{ii} ''' selected. Check spelling.'])
		end
	end
end
inp = cell(length(dose),2);
% Create inp cell array
if iscell(state)
	for ii = 1:length(dose)
		inp{ii,1} = state{ii};
		inp{ii,2} = dose(ii);
	end
else
	inp{1,1} = state;
	inp{1,2} = dose;
end
modType = modelType(model);
if strcmp(modType,'ode15s')
	if ischar(inp{1,1})
		error('findSS:inpStateStr','Input state must be converted to an x index is using an ODE15s model')
	end
	inp = cell2mat(inp);
	x0 = varargOde15s{1};
	x0(inp(:,1)) = x0(inp(:,1))+inp(:,2);
	varargOde15s{1} = x0;
	varargpass = varargOde15s;
else
	varargpass{end+1} = 'inp';
	varargpass{end+1} = inp;
end

t2 = tic;
if basal
	[t,Y,YComp,status] = findTC(model,tspan,varargpass{:});
else
	[t,Y,YComp,status] = findTC(model,tspan,varargpass{:},'-b');
end
y0_fst = Y(1,:);
y0     = Y(end,:);
y0c    = YComp(end,:);
varargOde15s{1} = y0c;
if exist('x0Indx','var')
	varargpass{x0Indx} = y0c;
else
	varargpass{end+1} = 'x0';
	varargpass{end+1} = y0c;
	x0Indx = length(varargpass);
end

while tolCal(Y(end-1,:),Y(end,:),tol)~=0
	tspan = tspan*5+tspan(2);
	[t,Y,YComp,status] = findTC(model,tspan,varargpass{:},'-r','-b');
	y0_fst = Y(1,:);
	y0     = Y(end,:);
	y0c    = YComp(end,:);
	varargpass{x0Indx} = y0c;
	varargOde15s{1} = y0c;
	% Exit triggers
	if isinf(y0)   
		status = -1;
		return
	end
	%if (toc(t2)) > 5
%		status = -2;
%		return
%	end
end

% Error documentation (only if parameter given. If not then assume passed
% at a higher level)
if exist('p','var')
	if status == -1  
		storeError(model,y0,p,msg,'Initial steady state produced Inf')
	elseif status == -2
		storeError(model,y0,p,msg,'Initial steady state taking too long')
	end
else
	if status == -1  
		error('findSS:InfResult','Initial steady state produced Inf\n')
	elseif status == -2
		error('findSS:ExcessRunTime','Initial steady state taking too long')
	end
end

function resid = tolCal(a,b,tol)
	resid = floor(abs(a-b)./max([a;b])/tol)*tol;
	resid(max([a;b])==0) = 0;
	resid = max(resid);