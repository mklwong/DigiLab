function varargout = MCMCOptimset(varargin)

% MCMCOptimset sets the options required for the MCMC function
%   OPTIONS = MCMCoptimset('PARAM1',VALUE1,'PARAM2',VALUE2,...) creates a
%   structure in OPTIONS that can be inputted into the MCMC algorithm.
%
%   PARAM is not case sensitive.
%
%   The options available are:
%       BASIC OPTIONS
%		- PTNO : Number of points to store. (Default = 10,000)
%
%       - T    : Tempering temperature.     (Default = 1)
%
%       - PRIOR: A prior set of points the MCMC will select from. (Empty = default)
%			     This must contain the fields:
%                       - prior.pts : List of points generated from previous
%                                     MCMC run
%                       - prior.logP: List of logP, corresponding to the
%                                     points in prior.pts
%                       - prior.ptUn: Logical vector of whether the
%                                     corresponding points in prior.pts is
%                                     unique or not.
%                These will have been outputted from a previous MCMC run.
%
%       - PMIN : Minimum probability the MCMC accepted point must be before
%                it's included in the posterior. 0 <= Pmin <= 1   (Default = 0)
%
%       - PARMODE : MCMC parallel mode. Defines how many cores will be used.
%                   If false, will not use parallel mode, if true will use 
%                   all cores possible. (Default = true)
%
%		- PASSNO: Number of points secondary labs will store before passing
%		          them to the primary lab. (Default = 500)
%
%       DISPLAY OPTIONS
%       - DISPLAY : Defines the verbosity of MCMC outputs. Options are:
%                        - Debug: All outputs from all labs are printed 
%                          into one text file per lab. Data transmission
%                          between labs are recorded.
%                        - Text: Outputs from lab 1 is printed into a text
%                          file. Only progress percentage is recorded.
%                        - Terminal: Same as text, but output is printed on
%                          the terminal.
%                        - Off: No display is shown.                 
%                          (Default = Terminal)
%
%		- DISPINT : Number of progress points to output through the run. If
%		            10, then program will print at 10%, 20%,...,90% progress. 
%                   (Default = 10)
%
%		- DIR: Directory to print run outputs to (Default = empty)
%
%       RUNTIME OPTIONS
%       - STALLTIME : Max time (in mins) allowed by the algorithm to find
%                     no new points before termination occurs (default = 5)
%
%       - WALLTIME  : Max time (in mins) the algorithm is allowed to run.
%                     (default = 60)
%
%		ADVANCED MCMC BEHAVIOUR OPTIONS
%		- ACCEPTRATIO: Target acceptance ratio (default = 0.5)
%
%		- STEPI : Starting step size normalised to parameter dimension size
%                 (default = 0.25)
%
%		- MAXSTEP: Max step size the program can evolve the system to.
%				   (default = 1)
%
%		- RESAMPLE: Number of points to be found before program will
%		            randomly find a new seed point. (Default = 50)
%
%		- PROPDIS: Function handle point to the function to be used for the
%		           proposal distribution (Default = @propDis)
%
%		- SEED: Random generator stream which the MCMC run will start from.
%		        (Default = empty)
%
%   OPTIONS = MCMCOptimset outputs the default options.

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
         'Stalltime  ' %max time between stored      | defaults to 5
                       %points in mins         
		 'Stepi      ' %starting average step size   | defaults to 0.25
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
opts.stalltime = 5;
opts.stepi     = 0.25;
opts.propDis   = @propDis;
opts.parMode   = true;
opts.passNo    = 500;
opts.resample  = 50;
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
			case lower(deblank(Names(9,:)))    %stall time
				opts.stalltime = varargin{ii+1};
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