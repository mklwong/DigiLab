function [pts,logP,status] = MCMC(objfun,pt0,bnd,opts)

% [pstrir,P] = MCMC(objfun,pt0,bnd,opts)
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
% prior = the results prior from a previous saved MCMCSeriel run
% Pmin = Minimum tempered probability required for acceptance into
%        posterior
%
% If prior is empty, then the MCMC generates from scratch. Else, the 
% program will randomly sample from the prior every so often.

%% Kernel
%=== Create (if no passed) or access program options ===%
if nargin == 3
    opts = MCMCOptimset;
elseif nargin == 4
	if ~isstruct(opts)
		error('MCMC:InvalidInput','Invalid options given for MCMC')
	end
end

%== Integrity Check ==%
bnd = sort(bnd,2);

%== Initialise displays ==%
if ~strcmpi(opts.disp,'off')
    if strcmpi(opts.disp,'full') && ~opts.parMode
        clf
    end
	startTime = clock;
	fprintf('Start time is %2.0f:%2.0f:%2.2f (%2.0f-%2.0f-%4.0f)\n',startTime(4),startTime(5),startTime(6),startTime(3),startTime(2),startTime(1))
end

%== Process parallel computing ==%
if opts.parMode
    if ~islogical(opts.parMode)
        parRes = parComp('open',opts.parMode);
    else
        parRes = parComp('open');
    end
	if parRes == -1
		opts.parMode = false;
	end
end

%== Add initial points and functions into optimisation settings ==%
if ~isempty(pt0)
	if isrow(pt0)
		opts.pt0 = pt0;
	else
		opts.pt0 = pt0';
	end
	if sum(double(opts.pt0<bnd(:,1)' | opts.pt0>bnd(:,2)'))
		find(double(opts.pt0<bnd(:,1)' | opts.pt0>bnd(:,2)'))
		error('MCMC:InitialPointOutOfBound','Initial point is outside the required boundary.')
	end
end
runVar.obj = objfun;
runVar.bnd = bnd;

%% MCMC Kernel split into multi core parallel and single core mode
if opts.parMode
	spmd
		[ptRaw,logPRaw,status] = MCMCRun(runVar,opts);
	end
else
    [ptRaw,logPRaw,status] = MCMCRun(runVar,opts);
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

function [ptLocal,logPLocal,status] = MCMCRun(runVar,opts)
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

%% Existence test for parallel computing parameters
if ~exist('numlabs','builtin')
    numlab = 1;
    labindx = 1;
else
    numlab = numlabs();
    labindx = labindex;
end

%% Process prior
if ~isempty(opts.prior.pts)
	prior = opts.prior;
	[prir_n,~] = size(prior.pts);
	% Remove prior points that are outside the boundary
	if prir_n>0
		rmIndx = ~logical(prod(double(prior.pts>(ones(prir_n,1)*runVar.bnd(:,1)') & prior.pts<(ones(prir_n,1)*runVar.bnd(:,2)')),2));
		prior.logP(rmIndx) = [];
		prior.pts(rmIndx,:) = [];
	end
	opts.prior = prior;
	clear prior rmIndx prir_n
end

if ~isempty(opts.prior)
	% Extract priors
	priorPts  = opts.prior.pts;
	priorLogP = opts.prior.logP;
	% Remove duplicate probabilities
	dupIndx = find(opts.prior.logP(1:end-1)-opts.prior.logP(2:end)==0);
	priorLogP(dupIndx)  = [];
	priorPts(dupIndx,:) = [];
	% Turn objective score into PDF and then CDF for prior
	[priorP,I] = sort(exp(-priorLogP));
	priorPts = priorPts(I,:);
	priorLogP = priorLogP(I);
	priorP = cumsum(priorP);
	% Remove low probability points
	rmPts = priorP<(1e-4*max(priorP));	
	priorP(rmPts) = [];
	priorLogP(rmPts) = [];
	priorPts(rmPts,:) = [];
	% Normalise CDF to 1
    runVar.prior.logP = priorLogP;
    runVar.prior.pts = priorPts;
	runVar.priorP = priorP/max(priorP);
	clear dupIndx rmPts priorLogP priorPts priorP
end

%% Initialisation
% Reinitialise the random number stream after spmd generated stream (which
% is always the same
mystream = RandStream.create('mrg32k3a','seed',sum(clock*100),'NumStreams',numlab,'StreamIndices',labindx); 
RandStream.setGlobalStream(mystream);

if labindx == 1
    ptNoMax = opts.ptNo;
else
    ptNoMax = opts.passNo;
end

if ~isfield(opts,'pt0')
	varNo = size(runVar.bnd,1);
else
	varNo = length(runVar.pt0);
end
% Initialise variables for storing points
ptLocal = zeros(ptNoMax,varNo);    
logPLocal  = nan(ptNoMax,1);    

% Initialise counters and markers
stallWarn = 0;  %number of stall cycles triggered (program stops when this hits 10)
pt_n   = 0;
runVar.logP   = Inf;  %Set this to enter point selection and acceptance loop

pltHndl    = []; % for initialising the program into the first plot of the run
adaptFun = opts.adaptFun;
dumVar.p0 = varNo; %Dummy runVar for adaptFun
dumVar.dp = [];
opts    = adaptFun('initial',dumVar,opts);
nprogress = 1; %Blocks of progress passed (each block is printed as a debug message)
acptCnt = int8(rand(1,20)>opts.rjtRto);
pt_uniQ_n = 0; %Counter for unique points.

% Markov Chain started below.
status = 1;
%% Start loop
while status == 1
	
    %% New active point selection
    if mod(pt_n,opts.resample)==0 || ~isfield(runVar,'pt')
		if (strcmpi(opts.disp,'full') || strcmpi(opts.disp,'text')) && labindx == 1
            fprintf('New starting point...\n')
            if isfield(runVar,'pt')
                fprintf('Existing point:  ')
                fprintf('%4.2e\n',runVar.logP)
			end
		end
        if ~isempty(opts.prior.pts)
			% Select new point from prior based on goodness of fit of the
			% prior
            if (strcmpi(opts.disp,'full') || strcmpi(opts.disp,'text')) && labindx == 1
                fprintf('Selected randomly from prior\n')
            end
            rngPt = rand(1);
			newPtInd = ceil(interp1([0;runVar.priorP],0:length(runVar.priorP),rngPt));
            ptTest = opts.prior.pts(newPtInd,:);
			logPNew  = runVar.obj(ptTest);
		elseif isfield(opts,'pt0')
			% If only a single prior point is given, do not jump around the
			% parameter space by reseeding. Also do not put boundaries on
			% the fitting.
            if (strcmpi(opts.disp,'full') || strcmpi(opts.disp,'text')) && labindx == 1
                fprintf('Selected from initial start point\n')
            end
			ptTest = opts.pt0;
			logPNew  = runVar.obj(ptTest);
			opts.resample = Inf;
			% Remove boundaries of the run, because with only one seed
			% point, the run is assumed to be exploratory so should be
			% unbounded.
			if isempty(runVar.bnd)
				runVar.bnd = ones(varNo,2);
				runVar.bnd(:,1) = -Inf;
				runVar.bnd(:,2) = Inf;
			end
		else
			if sum(runVar.bnd(:,1)==0 | isinf(runVar.bnd(:,2)))
				error('mcmc:unboundNoPrior','Cannot be run with no boundary when no prior is given')
            end
            if (strcmpi(opts.disp,'full') || strcmpi(opts.disp,'text')) && labindx == 1
                fprintf('Selected randomly from within boundary\n')
            end
            ptTest = rand(size(runVar.bnd,1),1).*(runVar.bnd(:,2)-runVar.bnd(:,1))+runVar.bnd(:,1);
            logPNew  = runVar.obj(ptTest);
            opts.resample = Inf;
        end % New candidate point
		
        if (strcmpi(opts.disp,'full') || strcmpi(opts.disp,'text')) && labindx == 1
            fprintf('Candidate point: ')
            fprintf('%4.2e \n',logPNew)
        end % Active visualisation: initial new point   

        % Metropolis acceptance criteria for new point
        thres = rand(1);
		a     = min([1 exp(-(logPNew-runVar.logP)/opts.T)]);
        if a > thres
            runVar.pt = ptTest;
            runVar.logP = logPNew;
            % Reinitialise MCMC parameters
            opts.step = ones(size(runVar.bnd(:,1)))*opts.stepi;
            if (strcmpi(opts.disp,'full') || strcmpi(opts.disp,'text')) && labindx == 1
                fprintf('chosen... (%7.1fs)\n',toc(t1))
				pt_n = pt_n + 1;
            end % Active visualisation: new point chosen
        else
            if (strcmpi(opts.disp,'full') || strcmpi(opts.disp,'text'))  && labindx == 1
                fprintf('rejected... (%7.1fs)\n',toc(t1))
            end % Active visualisation: old point retained
		end
		
    end

    %% MCMC evolution of active point
    [runVar,opts] = MCMCKernel(runVar,opts);
    acptCnt = [acptCnt(2:end) runVar.ptTest];
	opts = adaptFun(acptCnt,runVar,opts);

    %% Intermediate plotting of points (full display, only at single core mode)
    if ~opts.parMode && strcmpi(opts.disp,'full')
		% Test Scale
		logTest = log(runVar.bnd(:,2))-log(runVar.bnd(:,1));
		logScale = logTest>1&imag(logTest)==0;
		if ~logScale(1)
			fordDirX = runVar.pt(1)+opts.basis(1)*opts.step(1);
			RevDirX = runVar.pt(1)-opts.basis(1)*opts.step(1);
		else
			fordDirX = runVar.pt(1)*opts.basis(1)*opts.step(1);
			RevDirX = runVar.pt(1)/(opts.basis(1)*opts.step(1));
		end
		
		if ~logScale(2)
			fordDirY = runVar.pt(2)+opts.basis(2)*opts.step(2);
			RevDirY = runVar.pt(2)-opts.basis(2)*opts.step(2);
		else
			fordDirY = runVar.pt(2)*opts.basis(2)*opts.step(2);
			RevDirY = runVar.pt(2)/(opts.basis(2)*opts.step(2));
		end
		
        if isempty(pltHndl)
            subplot(2,2,[1 3])
			

			% Draw the boundaries
			curBasis = [opts.basis null(opts.basis')];%%%
			if isfield(opts,'pt0')
				pltHndl = plot(runVar.pt(1),runVar.pt(2),'x',...
				runVar.pt(1),runVar.pt(2),'o',...
				[fordDirX RevDirX],[fordDirY RevDirY]);%%%
			else
				pltHndl = plot(runVar.pt(1),runVar.pt(2),'x',...
				runVar.pt(1),runVar.pt(2),'o',...
				[runVar.pt(1) runVar.pt(1)+opts.basis(1)*opts.step(1)],[runVar.pt(2) runVar.pt(2)+opts.basis(2)*opts.step(1)],...%%%%
				[runVar.bnd(1,1) runVar.bnd(1,2)],[runVar.bnd(2,1) runVar.bnd(2,1)],...
				[runVar.bnd(1,2) runVar.bnd(1,2)],[runVar.bnd(2,1) runVar.bnd(2,2)],...
				[runVar.bnd(1,2) runVar.bnd(1,1)],[runVar.bnd(2,2) runVar.bnd(2,2)],...
				[runVar.bnd(1,1) runVar.bnd(1,1)],[runVar.bnd(2,2) runVar.bnd(2,1)]);
			end
			if logScale(1)
				set(gca,'XScale','log')
			end
			if logScale(2)
				set(gca,'YScale','log')
			end
			% Other diagnostic plots
            subplot(2,2,2)
            pltHndl2 = semilogy(1,runVar.logP/opts.T,[0 1],-log10([opts.Pmin opts.Pmin]),':');
 			subplot(2,2,4)
			pltHndl3 = plot(1,sqrt(sum(opts.step.^2)));
        else
			subplot(2,2,[1 3])
			xlabel('Param 1')
			ylabel('Param 2')
			title(['Pts saved = ' num2str(pt_uniQ_n,'%d')])
            set(pltHndl(1),'XData',[get(pltHndl(1),'XData') runVar.pt(1)],'YData',[get(pltHndl(1),'YData') runVar.pt(2)])
            n   = get(pltHndl3(1),'XData');
            set(pltHndl(2),'XData',runVar.pt(1),'YData',runVar.pt(2))
			curBasis = [opts.basis null(opts.basis')];
			set(pltHndl(3),'XData',[fordDirX RevDirX],'YData',[fordDirY RevDirY])
			set(pltHndl2(1),'XData',[n n(end)+1],'YData',[get(pltHndl2(1),'YData') runVar.logP/opts.T])
			set(pltHndl2(2),'XData',[0 n(end)+1])
 			set(pltHndl3,'XData',[n n(end)+1],'YData',[get(pltHndl3,'YData') sqrt(sum(opts.step.^2))])		
			drawnow
		end
	end
	
	%% Point storage
    %Acceptance criteria for storing of active point
    if runVar.logP/opts.T <= -log(opts.Pmin)
        if runVar.ptTest
            pt_uniQ_n = pt_uniQ_n + 1;
%             if (strcmpi(opts.disp,'full') || strcmpi(opts.disp,'text')) && labindx==1
%                 disp(pt_uniQ_n+1)
%             end
		end
        pt_n = pt_n + 1;
        logPLocal(pt_n) = runVar.logP;
        ptLocal(pt_n,:) = runVar.pt;
        t2 = tic;
        stallWarn = 0;
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
			pt_uniQ_n = 0;
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

function [runVar,opts] = MCMCKernel(runVar,opts)

% MCMCKernel generates and tests the new point point using the hasting
% ratio for MCMC. Point proposal is conducted

pt0 = runVar.pt;
logP0 = runVar.logP;

% Generate next point
[pt1,pdfBias] = opts.propDis(pt0,runVar.bnd,opts);

%Check point still in boundary
if sum(pt1 < runVar.bnd(:,1)) || sum(pt1 > runVar.bnd(:,2))
	error('MCMC:boundaryBreached','The program ran outside the boundary. Likely to be a bug in the proposal distribution')
end

% Try next point
if ~isrow(pt1)
	pt1 = pt1';
end
if ~isrow(pt0)
	pt0 = pt0';
end
logP1 = runVar.obj(pt1);      % Obtain likelihood of proposed step

% Metropolis Algorithm
test = min([1 exp((logP0-logP1)/opts.T)/pdfBias]); % Calculate hastings term.
														  % To make metropolis, set pdfBias
														  % to one in your proposal
														  % distribution function.

thres = rand(1);   % Randomise threshold

runVar.ptTest = thres < test;

%ABC Rejection when less than threshold

absthres = 1;
if logP1 < absthres
	runVar.ptTest = true;
end

if runVar.ptTest
	runVar.pt = pt1;
	runVar.logP = logP1;
else
	runVar.pt = pt0;
	runVar.logP = logP0;
end
try
	runVar.delPt = pt1-pt0;
catch
	keyboard
end
end