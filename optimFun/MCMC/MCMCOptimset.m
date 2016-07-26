function varargout = MCMCOptimset(varargin)

% MCMCOptimset sets the options required for the MCMC function
%   OPTIONS = MCMCoptimset('PARAM1',VALUE1,'PARAM2',VALUE2,...) creates a
%   structure in OPTIONS that can be inputted into the MCMC algorithm.
%
%   PARAM is not case sensitive.
%
%   The options available are:
%       - T   : Tempering temperature.                           (Default = 1)
%       - PRIR: A prior set of points the MCMC will select from. (Empty by default)
%       - PMIN: Minimum probability the MCMC accepted point must be before
%               it's included in the posterior. 0 <= Pmin <= 1   (Default = 0)
%       - MODE: MCMC run mode. Defines whether workers will communicate
%               their progress with each other. Value is boolean (Default = 0)
%       - VIS : Whether the MCMC will show extra outputs and plot a
%               representative figure of it's progress in parameter space.
%               Value must be boolean.                           (Default = false) 
%
%   OPTIONS = MCMCOptimset outputs the default option.

%========================%
%== Option Definitions ==%
%========================%

Names = ['PtNo       ' %number of points to store    | defaults at 10,000
		 'T          ' %value greater than zero      | defaults to 1
         'prior      ' %struct                       | defaults to empty
         'Pmin       ' %value less than one          | defaults to 0
		 'AcceptRatio' %Acceptance ratio before      | defaults to 0.5
		               %changing step size
         'display    ' %debug, text, terminal or off | defaults to text
		 'DispInt    ' %display interval: number of  | defaults to 10
		               %times the program prints     |
				 	   %a completion percentage      | 
		 'Walltime   ' %max run time in mins         | defaults to 60
		 'Stepi      ' %starting average step size   | defaults to 0.25
         'adaptFun   ' %function handle of how       | defaults to "adaptStep.m"
                       %MCMC step is adapted         |
                       %function of (scenario,pt0,pt1,step)
		 'PropDis    ' %function handle as function  | defaults to "propDis.m"
		               %of (p,bounds,stepsize)       |
		 'ParMode    ' %true, number or false        | defaults to true
		 'PassNo     ' %Number of points stored      | defaults at 500
		               %before passing.              |
		 'Resample   ' %Number of points stepped     | defaults at 50
		               %before attempting a new      |
					   %point (new point chosen      |
				  	   %subject to metropolis        |
					   %algorithm).                  |
		 'maxStep    ' %Maximum step size            | defaults at 1
		 'dir        ' %String for output location   | defaults at blank
		 'seed       ' %Initialisation seed          | defaults empty
		 ];         
     
%===============================%
%== create default opts struc ==%
%===============================%

opts.ptNo      = 10000;
opts.T         = 1;
	 prior.pts  = [];
	 prior.logP = [];
	 prior.T    = 1;
opts.prior     = prior;
opts.Pmin      = 0;
opts.rjtRto    = 0.5;
opts.disp      = 'Terminal';
opts.dispInt   = 10;
opts.walltime  = 60;
opts.stepi      = 0.25;
opts.propDis   = @(runVar) propDis(runVar);
opts.parMode   = true;
opts.passNo    = 500;
opts.resample  = 50;
opts.adaptFun  = @adaptStep;
opts.maxStep   = 1;
opts.dir       = '.';
opts.seed      = [];

%========================%
%== Construct new opts ==%
%========================%
detectOpts = true;
for ii = 1:length(varargin)
	if ii == 1
		if isstruct(varargin{1})
			opts = varargin{1};
		end
	end
	if ischar(varargin{ii}) && detectOpts
		detectOpts = false;
		switch lower(varargin{ii})
			case lower(deblank(Names(1,:)))    %ptNo
				opts.ptNo = varargin{ii+1};
			case lower(deblank(Names(2,:)))    %T
				opts.T = varargin{ii+1};
			case lower(deblank(Names(3,:)))    %prir
				opts.prior = varargin{ii+1};
			case lower(deblank(Names(4,:)))    %Pmin
				opts.Pmin = varargin{ii+1};
			case lower(deblank(Names(5,:)))    %actpRto
				opts.rjtRto = varargin{ii+1};
			case lower(deblank(Names(6,:)))    %display
				opts.disp = varargin{ii+1};
			case lower(deblank(Names(7,:)))    %display int
				opts.dispInt = varargin{ii+1};
			case lower(deblank(Names(8,:)))    %wall time
				opts.walltime = varargin{ii+1};
			case lower(deblank(Names(9,:)))    %step size
				opts.stepi = varargin{ii+1};
			case lower(deblank(Names(10,:)))   %adaptStep
				opts.adaptStep = varargin{ii+1};
			case lower(deblank(Names(11,:)))   %Proposal distribution
				opts.propDis = varargin{ii+1};
			case lower(deblank(Names(12,:)))   %parallel mode
				opts.parMode = varargin{ii+1};
			case lower(deblank(Names(13,:)))   %passNo mode
				opts.passNo = varargin{ii+1};
			case lower(deblank(Names(14,:)))   %resample
				opts.resample = varargin{ii+1};
			case lower(deblank(Names(15,:)))   %max step
				opts.maxStep = varargin{ii+1};
			case lower(deblank(Names(16,:)))   %storage directory
				opts.dir = varargin{ii+1};
			case lower(deblank(Names(17,:)))   %seed
				opts.seed = varargin{ii+1};
            otherwise
				warning('MCMCoptimset:unknowninput',['Input options ' num2str(ii) ' is non-existent. Check spelling.'])
					
		end
	else
		detectOpts = true;
	end
end

varargout{1} = opts;