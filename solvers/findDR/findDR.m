function [dspan,Y,YComp] = findDR(model,dspan,state,varargin)

% doseResponse(model,dose,state,varargin)
%	required inputs are:
%		- model : either 
%				* ode15s model (passed as function handle), 
%				* QSSA .m script file (passed as function handle or string)
%				* SBML .xml file (passed as function handle or string)
%		- dose  : 1xn vector of dose to spike in. Note doses are ALWAYS
%				  injected instantaneously at the start (i.e. step profile).
%		- states: 1xn cell of either protein names or integers per cell,
%				  corresponding to the state associated with the dose
%				  marked 
%			e.g. dose = {[ 1  2  3]; states = {'insulin';  OR states =  {12;
%						 [10 20 30];           'PDGF'   ;                15;
%						 [ 3  3  4]};          'VEGF'   };               16};
%				 this means [insulin] = [1 2 3]; 
%							[PDGF] = [10 20 30]; and 
%							[VEFG] = [3 3 4];
%			     with indicies ind(insulin) = 12, ind(PDGF) = 15 and ind(VEGF) = 16
%
%				
%	Optional inputs are (pass using options:value pairs)
%		- tEnd : the time span for each steady state run
%		- p    : vector of parameters. If not passed will search from ode
%				 file.
%		- x0   : initial condition (excluding the dose input). If not
%				 passed will search from ode file
%		- tol  : relative difference between initial and final state of
%				 each run before steady state is said ot be achieved.
%				 Defaults at 1e-3;
%
% Status codes:
% -1 = steady state produces inf
% -2 = steady state run is taking too long

model = parseModel(model);

% Determining number of doses to do
if iscell(dspan)
	lenD = length(dspan{1}(:));
	doses = zeros(length(dspan),lenD);
	for ii = 1:length(dspan)
		doses(ii,:) = dspan{ii}(:)';
	end
else
	lenD = length(dspan(:));
	doses = dspan(:)';
end

for ii = 1:lenD
	[Ytmp,YCtmp] = findSS(model,doses(:,ii),state,varargin{:});
	if ii == 1
		Y = nan(lenD,length(Ytmp));
		YComp = nan(lenD,length(YCtmp));
	end
	Y(ii,:) = Ytmp;
	YComp(ii,:) = YCtmp;
end

if isrow(dspan)
	dspan = dspan';
end