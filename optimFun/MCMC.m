function [pts,logP,status] = MCMC(objfun,bnd,opts)

% [pstrir,P] = MCMC(objfun,bnd,opts)
%
% Note the outputted probability is -log(P) of the real probability. Also
% it is the probability of the tempered landscape.
%
% MCMC chain with tempering.
%
% model = scenario file for model
% ptNo = number of points that will form the posterior
% U    = experimental training data
%
% Opts is a struct that contains optional parameters
% T    = annealing temperature
% prir = the results prior from a previous saved MCMCSeriel run
% Pmin = Minimum tempered probability required for acceptance into
%        posterior
%
% If prior is empty, then the MCMC generates from scratch. Else, the 
% program will randomly sample from the prior every so often.

%% Kernel
%=== Create or access program options ===%
if nargin == 2
    opts = MCMCOptimset;
elseif nargin == 3
	if ~isstruct(opts)
		error('MCMC:InvalidInput','Invalid options given for MCMC')
	end
end

%== Integrity Check ==%
if sum(double(bnd(:,1)>bnd(:,2)))>0
	error('MCMC:LowerBoundLargerThanUpperBound','Lower bound is set to be larger than the upper bound')
end

%== Process displays ==%
if ~strcmpi(opts.disp,'off')
    if strcmpi(opts.disp,'full') && ~opts.parMode
        clf
    end
	startTime = clock;
	fprintf('Start time is %2.0f:%2.0f:%2.2f (%2.0f-%2.0f-%4.0f)\n',startTime(4),startTime(5),startTime(6),startTime(3),startTime(2),startTime(1))
end

if opts.parMode
    parComp('open');
end

%% Process prior based on boundary
if ~isempty(opts.prir)
	prir = opts.prir;
	[prir_n,~] = size(prir.pts);
	% Remove prior points that are outside the boundary
	if prir_n>0
		rmIndx = ~logical(prod(double(prir.pts>(ones(prir_n,1)*bnd(:,1)') & prir.pts<(ones(prir_n,1)*bnd(:,2)')),2));
		prir.logP(rmIndx) = [];
		prir.pts(rmIndx,:) = [];
	end
	opts.prir = prir;
end

%% MCMC Kernel split into multi core parallel and single core mode
if opts.parMode
	spmd
		[ptRaw,logPRaw,status] = MCMCRun(objfun,bnd,opts);
	end
else
    [ptRaw,logPRaw,status] = MCMCRun(objfun,bnd,opts);
end

%% Run Completion
% Store lab one's compiled work
status = status{1};
pts = ptRaw{1};
logP   = logPRaw{1};

% Remove unused data slots
pts(isnan(logP),:) = [];
logP(isnan(logP)) = [];  

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% MCMC Function %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ptLocal,logPLocal,status] = MCMCRun(objfun,bnd,opts)
% Parallel mode protocols
% Packet tags are:
%    1: Sending and receiving points and data
%    2: Communication of run status
%
% Program status:
%    1: Enter and continue MCMC run
%    0: MCMC exited with no errors
%   -1: MCMC exited as stall time exceeded.

t1 = tic; %t1 is the MCMC begin time (in global terms)
t2 = tic; %t2 is the time to last completion checkpoint 
%this is also the counter used for run termination.

%Existence test for parallel computing parameters
if ~exist('numlabs','builtin')
    numlab = 1;
    labindx = 1;
else
    numlab = numlabs();
    labindx = labindex;
end

%Check that boundary is correctly supplied
if size(bnd,1) ~= 2 && size(bnd,2) ~= 2
	error('MCMC:BoundaryInvalid','Boundary not supplied correctly. Must be [lb,ub] where lb and ub are 1xN vectors where N is the number of free parameters')
end

% Reinitialise the random number stream after spmd generated stream (which
% is always the same
mystream = RandStream.create('mrg32k3a','seed',sum(clock*100),'NumStreams',numlab,'StreamIndices',labindx); 
RandStream.setGlobalStream(mystream);

if labindx == 1
    ptNoMax = opts.ptNo;
else
    ptNoMax = opts.passNo;
end

% Initialise variables for storing points
ptLocal = zeros(ptNoMax,size(bnd,1));    
logPLocal  = nan(ptNoMax,1);    


% Initialise counters and markers
stallWarn = 0;  %number of stall cycles triggered (program stops when this hits 10)
pt_n   = 0;
stp_n  = 0;
logP      = Inf;  %Set this to enter point selection and acceptance loop
prir_n = size(opts.prir.pts,1);
pltHndl    = []; % for initialising the program into the first plot of the run
adaptFun = opts.adaptFun;
opts.step = ones(size(bnd(:,1)))*opts.stepi;
opts    = adaptFun('initial',bnd(:,1),[],opts);
dumPt = (bnd(:,1)+min([bnd(:,2),bnd(:,1)*1e2],[],2))/2;
[~,~,logScl] = opts.propDis(dumPt',bnd,opts);
nprogress = 1;
acptCnt = int8(rand(1,20)>opts.rjtRto);
pt_uniQ_n = 0;

% Markov Chain started below.
status = 1;
%% Start loop
while status == 1
	
    %% New active point selection
    if mod(stp_n,opts.resample)==0 || ~exist('pt','var')
		if (strcmpi(opts.disp,'full') || strcmpi(opts.disp,'text')) && labindx == 1
            fprintf('New starting point...\n')
            if exist('pt','var')
                fprintf('Existing point:  ')
                fprintf('%4.2e, ',pt(1:end-1))
                fprintf('%4.2e\n',pt(end))
			end
		end
        if prir_n > 1
            rngPt = ceil(rand(1)*prir_n);
            ptTest = opts.prir.pts(rngPt,:);
            logPNew  = opts.prir.logP(rngPt);
		elseif prir_n == 1
			% If only a single prior point is given, do not jump around the
			% parameter space by reseeding. Also do not put boundaries on
			% the fitting.
			ptTest = opts.prir.pts;
			logPNew  = opts.prir.logP;
			opts.resample = Inf;
			% Remove boundaries of the run, because with only one seed
			% point, the run is assumed to be exploratory so should be
			% unbounded.
			bnd(:,1) = 0;
			bnd(:,2) = Inf;
		else
			if sum(bnd(:,1)==0 | isinf(bnd(:,2)))
				error('mcmc:unboundNoPrior','Cannot be run with no boundary when no prior is given')
			end
            ptTest = 10.^(rand(size(bnd,1),1).*(log10(bnd(:,2))-log10(bnd(:,1)))+log10(bnd(:,1)));
            pt2 = rand(size(bnd,1),1).*(bnd(:,2)-bnd(:,1))+bnd(:,1);
            ptTest(~logScl) = pt2(~logScl);
            logPNew  = objfun(ptTest);
            opts.resample = Inf;
        end % New candidate point
		
        if (strcmpi(opts.disp,'full') || strcmpi(opts.disp,'text')) && labindx == 1
            fprintf('Candidate point: ')
            fprintf('%4.2e, ',ptTest(1:end-1))
            fprintf('%4.2e\n',ptTest(end))
        end % Active visualisation: initial new point   
		
        % Metropolis acceptance criteria for new point
        thres = rand(1);
		a     = min([1 exp(-(logPNew-logP)/opts.T)]);
        if a > thres
            pt = ptTest;
            logP = logPNew;
            % Reinitialise MCMC parameters
            opts.step = ones(size(bnd(:,1)))*opts.stepi;
            if (strcmpi(opts.disp,'full') || strcmpi(opts.disp,'text')) && labindx == 1
                fprintf('chosen... (%7.1fs)\n',toc(t1))
            end % Active visualisation: new point chosen
        else
            if (strcmpi(opts.disp,'full') || strcmpi(opts.disp,'text'))  && labindx == 1
                fprintf('rejected... (%7.1fs)\n',toc(t1))
            end % Active visualisation: old point retained
        end
    end

    %% MCMC evolution of active point
    [pt,logP,testRes,delPt,opts] = MCMCKernel(objfun,pt,logP,bnd,opts);
    acptCnt = [acptCnt(2:end) testRes];
	opts = adaptFun(acptCnt,pt,delPt,opts);

    %% Intermediate plotting of points (full display, only at single core mode)
    if ~opts.parMode && strcmpi(opts.disp,'full')
        if  isempty(pltHndl)
            subplot(3,1,1)
			% Draw the boundaries
            pltHndl = plot(pt(1),pt(2),'x',...
            pt(1),pt(2),'o',...
            [bnd(1,1) bnd(1,2)],[bnd(2,1) bnd(2,1)],...
            [bnd(1,2) bnd(1,2)],[bnd(2,1) bnd(2,2)],...
            [bnd(1,2) bnd(1,1)],[bnd(2,2) bnd(2,2)],...
            [bnd(1,1) bnd(1,1)],[bnd(2,2) bnd(2,1)]);
            if logScl(1)
                set(gca,'XScale','log')
            end
            if logScl(2)
                set(gca,'YScale','log')
            end
            subplot(3,1,2)
            pltHndl2 = semilogy(1,logP/opts.T,[0 1],-log10([opts.Pmin opts.Pmin]),':');
			subplot(3,1,3)
			pltHndl3 = plot(1,sqrt(sum(opts.step.^2)));
        else
			subplot(3,1,1)
			xlabel('Param 1')
			ylabel('Param 2')
			title(['Pts saved = ' num2str(pt_uniQ_n,'%d')])
            set(pltHndl(1),'XData',[get(pltHndl(1),'XData') pt(1)],'YData',[get(pltHndl(1),'YData') pt(2)])
            n   = get(pltHndl2(1),'XData');
            set(pltHndl(2),'XData',pt(1),'YData',pt(2))
			
			set(pltHndl2(1),'XData',[n n(end)+1],'YData',[get(pltHndl2(1),'YData') logP/opts.T])
			set(pltHndl2(2),'XData',[0 n(end)+1])
			set(pltHndl3,'XData',[n n(end)+1],'YData',[get(pltHndl3,'YData') sqrt(sum(opts.step.^2))])
			drawnow
		end
	end
	
	%% Point storage
    %Acceptance criteria for storing of active point
    if logP/opts.T <= -log(opts.Pmin)
        if testRes
            pt_uniQ_n = pt_uniQ_n + 1;
            if (strcmpi(opts.disp,'full') || strcmpi(opts.disp,'text')) && labindx==1
                disp(pt_uniQ_n+1)
            end
        end
        pt_n = pt_n + 1;
        logPLocal(pt_n) = logP;
        ptLocal(pt_n,:) = pt;
        t2 = tic;
        stallWarn = 0;
	end
	
	%% Program escape
    % Exit clause for other labs
    if opts.parMode && labindx ~= 1
        if labProbe(1,2)
            if (strcmpi(opts.disp,'full') || strcmpi(opts.disp,'text')) && labindex==1
                fprintf('Exit signal received. Quitting.\n\r');
            end
            status = labReceive(1,2);
        end
	end
	
    %% Parallel mode packet send and receive
    if opts.parMode
        if labindx == 1 
            while labProbe('any',1)
                dat = labReceive();
                ptNew    = dat{1};
                logPNew  = dat{2};
                slave_pt_uniQ_n = dat{3};
                pt_uniQ_n = pt_uniQ_n + slave_pt_uniQ_n;
                pt_n_New = length(logPNew);
                pt_n_Get = min([ptNoMax-pt_n pt_n_New]);
                ptLocal((pt_n+1):(pt_n+pt_n_Get),:) = ptNew(1:pt_n_Get,:);
                logPLocal((pt_n+1):(pt_n+pt_n_Get)) = logPNew(1:pt_n_Get,:);
                pt_n = pt_n+pt_n_Get;
                if (strcmpi(opts.disp,'full') || strcmpi(opts.disp,'text')) && labindex==1
                    fprintf('Data packet received (pt_n = %d).\n\r',slave_pt_uniQ_n);
                end
            end
        elseif labindx > 1 && pt_uniQ_n == ptNoMax
            labSend({ptLocal,logPLocal,pt_uniQ_n},1,1);
            logPLocal = nan(size(logPLocal));
            ptLocal = nan(size(ptLocal));
            pt_n = 0;
            if (strcmpi(opts.disp,'full') || strcmpi(opts.disp,'text'))
                fprintf('Data packet sent.\n\r');
            end
        end %when packet full, send
    end

    %% Print progress report for running

    if pt_uniQ_n >= nprogress*opts.ptNo/opts.dispInt && labindx == 1
        nprogress = nprogress + 1;
        fprintf('%3.0f%% done after %7.1f seconds.\n',(pt_uniQ_n/ptNoMax*100),toc(t1))
        if pt_uniQ_n/ptNoMax >= 1
            fprintf('Run complete. Quitting. \n')
            if opts.parMode
                labSend(0,2:numlab,2) %Send stop signal
                if (strcmpi(opts.disp,'full') || strcmpi(opts.disp,'text')) && labindex==1
                    fprintf('Exit signal sent.\n\r');
                end
            end
            status = 0;
        end
    end

    %% Stall handling
    if labindx == 1 && toc(t2)>((stallWarn+1)*opts.walltime*60/10)
        stallWarn = stallWarn + 1;
        if labindx == 1
            fprintf('Program still running, but stuck in low probability area (%2.2f). (Last:%7.1fs|Tot:%7.1fs)\n',stallWarn,toc(t1),toc(t2))
            if strcmpi(opts.disp,'full') || strcmpi(opts.disp,'text')
                fprintf('Program still running, but stuck in low probability area (%2.2f). (Last:%7.1fs|Tot:%7.1fs)\r\n',stallWarn,toc(t1),toc(t2));
            end
        end
        if toc(t2)>opts.walltime*60
            fprintf('Lab stop triggered due to taking too long...')
            if strcmpi(opts.disp,'full') || strcmpi(opts.disp,'text')
                fprintf('Lab stop triggered due to taking too long...');
            end
            if opts.parMode
                if strcmpi(opts.disp,'full') || strcmpi(opts.disp,'text')
                fprintf('Exit signal sent. Quitting.\n\r');
                end
                labSend(-1,2:numlab,2) %Send stop signal
            end
            status = -1;
        end 
    end
end

%% Run completion
if opts.parMode
	if strcmpi(opts.disp,'full')
		fclose(h);
	end
labBarrier
end

if labindx == 1
    fprintf('Done!\n')
end

if ~opts.parMode
    status = {status};
    ptLocal = {ptLocal};
    logPLocal = {logPLocal};
end %Make results into cells for numlabs = 1 so finished results can be handled the same way regardless of parallel run or not
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Sub function starts below %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pt,logP,ptTest,pt1,opts] = MCMCKernel(objfun,pt0,logP0,bnd,opts)

% MCMCKernel generates and tests the new point point using the hasting
% ratio for MCMC. Point proposal is conducted

% Generate next point
[pt1,pdfBias] = opts.propDis(pt0,bnd,opts);

%Check point still in boundary
if sum(pt1 < bnd(:,1)) || sum(pt1 > bnd(:,2))
	error('MCMC:boundaryBreached','The program ran outside the boundary. Likely to be a bug in the proposal distribution')
end

% Try next point
logP1 = objfun(pt1);      % Obtain likelihood of proposed step

% Metropolis Algorithm
test = min([1 exp((logP0-logP1)/opts.T)/pdfBias]); % Calculate hastings term.
										  % To make metropolis, set pdfBias
										  % to one in your proposal
										  % distribution function.

thres = rand(1);   % Randomise threshold

ptTest = thres < test;

if ptTest
	pt = pt1;
	logP = logP1;
else
	pt = pt0;
	logP = logP0;
end

end