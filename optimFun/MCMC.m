function [pts,logP,ptUnique,status] = MCMC(objfun,pt0,bndry,opts)

% [pstrir,P] = MCMC(objfun,pt0,bndry,opts)
%
% Note the outputed probability is value of the objective function. It is
% assumed that the relationship between the objective function and the
% likelihood is P = exp[-objfun(p)]
%
% MCMC chain with tempering.
%
% objfun = Objective function to be explored.
% pt0    = Starting point for the MCMC chain. Can be left blank for a
%          random start point to be chosen automatically.
% bndry  = Boundry within which the MCMC algorithm will be restricted
%          within. This is an N x 2 matrix where N is the number of
%          parameters to be fitted. The first column contains the lower
%          bounds while the second column contains the upper bounds.
% Opts   = A struct that contains MCMC settings. It is advised that this be
%          generated using MCMCOptimset
%
% pt0 and bndry are both optional inputs, but AT LEAST ONE of these input
% arguments must be entered. This is because the randomised start point
% requires a finite space to sample from. Conversely an unbounded search
% must have a defined starting point.

%% Input Parameter Integrity Check
%== Parameter 3: Boundary ==%
bndry = sort(bndry,2); %Make sure boundary is sorted as lower bound in first column and upper bound in second column

%== Parameter 2: Initial point ==%
if ~isempty(pt0)
	if isrow(pt0)
		pt0 = pt0';
	end
	if ~isempty(bndry)
		if any(pt0<bndry(:,1) | pt0>bndry(:,2))
			find(pt0<bndry(:,1) | pt0>bndry(:,2))
			error('MCMC:InitialPointOutOfBound','Initial point is outside the required boundary.')
		end
	end
end

%=== Parameter 4: Options ===%
if nargin == 3         % When options not passed, used defaults.
    opts = MCMCOptimset;
elseif nargin == 4     % When options passed, check that the option class is correct
	if ~isstruct(opts)
		error('MCMC:InvalidInput','Invalid options given for MCMC')
	end
end

%% Condensing input parameters for passing into MCMC Kernel
runVar.pt = pt0;
runVar.obj = objfun;
runVar.bnd = bndry;

%% Preprocessing Prior such that only highest probabilities points are considered
if ~isempty(opts.prior.pts)
	prior = opts.prior;
	[prir_n,~] = size(prior.pts);
	% Remove prior points that are outside the boundary
	if prir_n>0
		rmIndx = ~all( (prior.pts>(ones(prir_n,1)*runVar.bnd(:,1)')) & (prior.pts<(ones(prir_n,1)*runVar.bnd(:,2)')) ,2);
		prior.logP(rmIndx) = [];
		prior.pts(rmIndx,:) = [];
	end
	% Extract priors
	priorPts  = prior.pts(prior.ptUn,:);
	priorLogP = prior.logP(prior.ptUn)/opts.T;
	% Turn objective score into PDF and then CDF for prior
	[priorP,I] = sort(exp(-priorLogP));
	priorPts = priorPts(I,:);
	priorLogP = priorLogP(I);
	priorP = cumsum(priorP);
	% Remove points that are 1000x less likely to be true compared to the
	% best point
	rmPts = priorP<(1e-3*max(priorP));	
	priorMin = max(priorP(rmPts));
	priorP(rmPts) = [];
	priorLogP(rmPts) = [];
	priorPts(rmPts,:) = [];
	% Normalise CDF to 1
	runVar.priorP = (priorP-priorMin)/(max(priorP)-priorMin);
    opts.prior.pts = priorPts;
    opts.prior.logP = priorLogP;
	clear dupIndx rmPts priorLogP priorPts priorP prior rmIndx prir_n
else
	runVar.priorP = [];
end

%% Check for and open parallel computing
if opts.parMode
	% Utilises parComp to check computer's parallel computing capacity, and
	% open the cluster if available.
	if islogical(opts.parMode)
		parRes = parComp('open');
	else
		parRes = parComp('open',opts.parMode);
	end
	% If cluster not available, revert to non-parallel mode
    if parRes == -1
		opts.parMode = false;
    end
end

%% Check for whether the working folder exists
if ~exist(opts.dir,'dir')
	mkdir(opts.dir);
end

%% MCMC Kernel split into multi core parallel and single core mode
if opts.parMode
	spmd
        warningSwitch('off')
		[ptRaw,logPRaw,ptUniqueRaw,status] = MCMCKernel(runVar,opts);
		warningSwitch('on')
	end
else
    warningSwitch('off')
    [ptRaw,logPRaw,ptUniqueRaw,status] = MCMCKernel(runVar,opts);
	warningSwitch('on')
end

%% Run Completion
% Store lab one's compiled work
status = status{1};
pts    = ptRaw{1};
logP   = logPRaw{1};
ptUnique = ptUniqueRaw{1};

% Remove unused data slots
pts(isnan(logP),:) = [];
ptUnique(isnan(logP)) = [];
logP(isnan(logP))  = [];  

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% MCMC Function %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ptLocal,logPLocal,ptUniqLocal,status] = MCMCKernel(runVar,opts)

%
% Parallel mode protocols
% Packet tags are:
%    1: Sending and receiving points and data
%    2: Communication of run status
%
% Program status:
%    1: Enter and continue MCMC run
%    0: MCMC exited with no errors
%   -1: MCMC exited as stall time exceeded.

%%
% ==========================================
% =================Preamble=================
% ==========================================

%       -----------------------------
% ----- Initialise parallel computing ------
%       -----------------------------
% Compatibility with versions without parallel computing toolbox
if ~exist('numlabs','builtin')
    numlab = 1;
    labindx = 1;
else
    numlab = numlabs();
    labindx = labindex;
end

% Distinguish between master and slave in terms of total points that need to be saved
if labindx == 1
    ptNoMax = opts.ptNo;
else
    ptNoMax = opts.passNo;
end

% Set up random number stream
if isempty(opts.seed) % If no random seed given. Generate and save
	mystream = RandStream.create('mrg32k3a','seed',sum(clock*100),'NumStreams',numlab,'StreamIndices',labindx); % Initialise
	save([opts.dir '/Seed-T_' num2str(opts.T) '-Slave ' num2str(labindx)],'runVar','opts','mystream');
else                  % Else use the existing seed
	mystream = opts.seed{labindx};
end
RandStream.setGlobalStream(mystream); %set the seed

%      ---------------------------
% ----- Initialise Display Outputs ------
%      ---------------------------

% Start run timers
t1   = tic;   % t1 is the MCMC begin time (in global terms)
t2   = tic;   % t2 is the time to last completion checkpoint 
tNow = clock; % Current time

% Create output files if necessary
runVar.outputName = [];
if strcmpi(opts.disp,'debug')
    runVar.outputName = [opts.dir '/Output-T_' num2str(opts.T) '-Slave ' num2str(labindx) '.txt'];
    runVar.outFileHandle = fopen(runVar.outputName,'wt');
elseif strcmpi(opts.disp,'text') && labindx == 1
	runVar.outputName = [opts.dir '/Output-T_' num2str(opts.T) '.txt' ];
    runVar.outFileHandle = fopen(runVar.outputName,'wt'); 
elseif labindx~=1 || strcmpi(opts.disp,'off')
	runVar.outFileHandle = []; %If outhandle is [], does not print
else 
    runVar.outFileHandle = 1;  %If outhandle is 1, prints to terminal
end
fprintf_cust(runVar.outFileHandle,'-------------------------------\n');
fprintf_cust(runVar.outFileHandle,'--------Run Information--------\n');
fprintf_cust(runVar.outFileHandle,'-------------------------------\n');
fprintf_cust(runVar.outFileHandle,'T = %1.2f | ',opts.T);
fprintf_cust(runVar.outFileHandle,'Run Begins at %2.0f:%2.0f:%2.0f (%2.0f-%2.0f-%4.0f) \n',tNow([4:6 3:-1:1]));

%      -------------------------
% ----- Initialise Run Variables ------
%      -------------------------
% Choose start point
if ~isempty(runVar.pt)
	varNo = length(runVar.pt);
	opts.resample = Inf; %no resampling when picking one start point
elseif ~isempty(runVar.priorP)
	varNo = size(opts.prior.pts(1,:),2);
	rngPt = rand(1);
	newPtInd = ceil(interp1([0;runVar.priorP],0:length(runVar.priorP),rngPt));
    runVar.pt    = opts.prior.pts(newPtInd,:)';
elseif isempty(runVar.bnd)
	error('mcmc:unboundNoPrior','MCMC cannot be run with no boundary when no prior is given')
else
	varNo = size(runVar.bnd,1);
	runVar.pt   = seedPt(runVar);
end
runVar.logP = runVar.obj(runVar.pt);

% Initialise run functions and parameters
runVar = opts.propDis('adapt',runVar,opts);
runVar.ptTest = mod(1:20,2); % Vector of result of past acceptance tests. As an initial start point, we assume suggests and failures were alternating evenly.
if isempty(runVar.bnd)
	runVar.logScale = false(size(runVar.pt));
else
	[~,runVar.logScale] = seedPt(runVar);
end

% Initialise storage variables (store 10 times more than necessary in case
% of duplicates)
ptLocal = zeros(ptNoMax*10,varNo);
ptUniqLocal = false(ptNoMax*10,1);   
logPLocal   = nan(ptNoMax*10,1);    

% Initialise run monitors
stepCount = 1;   % Total number of steps (independent of Pmin threshold)
rjtCount  = 0;   % Continuous rejection count
stallWarn = 0;   % Number of stall cycles triggered (program stops when this hits 10)
pltHndl   = [];  % For initialising the program into the first plot of the run
nprogress = 1;   % Blocks of progress passed (each block is printed as a debug message)

% ==========================================
% ============== MCMC Start ================
% ==========================================

% Initialise lab status
if opts.parMode  
 	otherLabStat = zeros(1,numlabs);
end

status = 1; % Enter MCMC loop

%%
while status == 1
	
printCheckpoint('1',runVar.outputName,opts.disp);	

%      --------------------
% ----- Find new seed point ------
%      --------------------
if mod(stepCount,opts.resample)==0
	if ~isempty(runVar.priorP)
		% Select new point from prior based on goodness of fit of the
		% prior
		rngPt = rand(1);
		newPtInd = ceil(interp1([0;runVar.priorP],0:length(runVar.priorP),rngPt));
		ptTest = opts.prior.pts(newPtInd,:)';
		logPNew  = runVar.obj(ptTest);
	else
		ptTest = seedPt(runVar);
		logPNew  = runVar.obj(ptTest);
	end % New candidate point

	% Metropolis acceptance criteria for new point
	thres    = rand(1);
	testProb = min([1 exp(-(logPNew-runVar.logP)/opts.T)]);
	if (testProb > thres) || rjtCount > opts.resample
		runVar.pt = ptTest;
		runVar.logP = logPNew;
		% Reinitialise MCMC parameters
		opts.step = ones(size(runVar.bnd(:,1)))*opts.stepi;
		rjtCount = 0;
	end
end

   
%      --------------------
% ----- Evolve Active Point ------
%      --------------------
[runVar,opts] = MCMCEvolve(runVar,opts);
stepCount = stepCount + 1;

if mod(floor(toc(t1))/60,10) == 0 && strcmpi(opts.disp,'text')
	tNow = clock;
	fprintf_cust(runVar.outFileHandle,'Time elapsed - %1.0f | Real time - %2.0f:%2.0f:%2.0f \n',floor(toc(t1)/60),tNow(4:6));
end
printCheckpoint('3',runVar.outputName,opts.disp);

%      ----------------------------
% ----- Plot Distribution of Points ------
%      ----------------------------
if ~opts.parMode && strcmpi(opts.disp,'full')
	% Test Scale
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
		title(['Pts saved = ' num2str(sum(ptUniqLocal),'%d')])
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
printCheckpoint('4',runVar.outputName,opts.disp);

%Debug workspace saving
if strcmpi(opts.disp,'debug')
	save([opts.dir '/DebugWorkspace-T_' num2str(opts.T) '-Slave' num2str(labindx) '.mat']) %Save entire workspace
end

%       ---------------------
% ----- Store point in worker ------
%       ----------------------
if runVar.logP/opts.T <= -log(opts.Pmin)
	ptUniq = false;
	if runVar.ptTest(end)
		ptUniq = true;
		rjtCount = 0;
    else
        rjtCount = rjtCount + 1;
	end
	pt_n = find(isnan(logPLocal),1,'first');
	if isempty(pt_n) %Increase size of storage matrices is necessary
		ptLocal     = [ptLocal;zeros(ptNoMax*10,varNo)];
		ptUniqLocal = [ptUniqLocal;false(ptNoMax*10,1)]; 
		logPLocal   = [logPLocal;nan(ptNoMax*10,1)];
		pt_n = find(isnan(logPLocal),1,'first');
	end
	ptUniqLocal(pt_n) = ptUniq;
	logPLocal(pt_n)   = runVar.logP;
	ptLocal(pt_n,:)   = runVar.pt;
	t2 = tic;
	stallWarn = 0;
else
	rjtCount = rjtCount + 1;
end
printCheckpoint('5',runVar.outputName,opts.disp);
    
%       -----------------------------
% ----- Primary worker progress check ------
%       -----------------------------
if (sum(ptUniqLocal) >= nprogress*opts.ptNo/opts.dispInt) && labindx == 1
	nprogress = nprogress + 1;
	tNow = clock;
	fprintf_cust(runVar.outFileHandle,'%3.0f%% done after %7.1f seconds. | (%2.0f:%2.0f:%2.0f) \n',(sum(ptUniqLocal)/ptNoMax*100),toc(t1),tNow(4:6));
    if sum(ptUniqLocal)/ptNoMax >= 1
		fprintf_cust(runVar.outFileHandle,'Current block complete. | (%2.0f:%2.0f:%2.0f) \n',tNow(4:6));
		status = 0;
		if opts.parMode
			labSend(status,2:numlab,2) %Send stop signal
			fprintf_cust(runVar.outFileHandle,'Exit signal sent. | (%2.0f:%2.0f:%2.0f) \n',tNow(4:6));
		end
	end
end
printCheckpoint('6',runVar.outputName,opts.disp);

%       ------------------------------
% ----- Stall check (lack of progress) ------
%       ------------------------------
if labindx == 1 && toc(t2)>((stallWarn+1)*opts.stalltime*60/10) && status~=0
	stallWarn = stallWarn + 1;
	tNow = clock;
	fprintf_cust(runVar.outFileHandle,'Program still running, but stuck in low probability area (%2.2f). (Last:%7.1fs|Tot:%7.1fs)  | (%2.0f:%2.0f:%2.0f) \n',stallWarn,toc(t1),toc(t2),tNow(4:6));
	if toc(t2)>opts.stalltime*60
		fprintf_cust(runVar.outFileHandle,'Lab stop triggered due to taking too long (Stalltime)... | (%2.0f:%2.0f:%2.0f) \n',tNow(4:6));
		status = -1;
		if opts.parMode
			labSend(status,2:numlab,2) %Send stop signal
			fprintf_cust(runVar.outFileHandle,'Exit signal sent. Quitting.\n | (%2.0f:%2.0f:%2.0f) \n',tNow(4:6));
		end
	end 
end

%       ------------------------------
% ----- Walltime check (lack of progress) ------
%       ------------------------------
if toc(t2)>opts.walltime*60
    fprintf_cust(runVar.outFileHandle,'Lab stop triggered due to taking too long (Walltime)... | (%2.0f:%2.0f:%2.0f) \n',tNow(4:6));
    status = -2;
    if opts.parMode
        labSend(status,2:numlab,2) %Send stop signal
        fprintf_cust(runVar.outFileHandle,'Exit signal sent. Quitting.\n | (%2.0f:%2.0f:%2.0f) \n',tNow(4:6));
    end
end 
printCheckpoint('7',runVar.outputName,opts.disp);

%       ----------------------------------
% ----- Secondary worker exit signal check ------
%       ----------------------------------
if opts.parMode && labindx ~= 1
	if labProbe(1,2)
		tNow = clock;
		fprintf_cust(runVar.outFileHandle,'Exit signal received. Quitting. | (%2.0f:%2.0f:%2.0f) \n',tNow(4:6));
		status = labReceive(1,2);
	end
end
printCheckpoint('8',runVar.outputName,opts.disp);

%       -----------------------------------------------
% ----- Send data from secondary work to primary worker ------
%       -----------------------------------------------
if labindx == 1 && opts.parMode && status == 1
	while labProbe('any',1)
		[dat,srcIndx] = labReceive('any');
		if isnumeric(dat)
			otherLabStat(srcIndx) = dat;
		else
			ptNew     = dat{1};
			logPNew   = dat{2};
			ptUniqNew = dat{3};
			pt_n = find(isnan(logPLocal),1,'first');
			pt_n_New = length(logPNew);
			ptLocal((pt_n+1):(pt_n+pt_n_New),:) = ptNew(1:pt_n_New,:);
			ptUniqLocal((pt_n+1):(pt_n+pt_n_New)) = ptUniqNew(1:pt_n_New,:);
			logPLocal((pt_n+1):(pt_n+pt_n_New)) = logPNew(1:pt_n_New,:);
			if ~strcmpi(opts.disp,'off')
				tNow = clock;
				fprintf_cust(runVar.outFileHandle,'Data packet received (pt_n = %d) from slave %d. | (%2.0f:%2.0f:%2.0f) \n',nansum(ptUniqNew),srcIndx,tNow(4:6));
            end
		end
	end
elseif (labindx > 1 && nansum(ptUniqLocal) == ptNoMax) && status == 1
    ptLast = find(~isnan(logPLocal),1,'last');
	labSend({ptLocal(1:ptLast,:),logPLocal(1:ptLast),ptUniqLocal(1:ptLast)},1,1);
	ptUniqLocal = false(size(logPLocal));
	logPLocal = nan(size(logPLocal));
	ptLocal = zeros(size(ptLocal));
	tNow = clock;
	fprintf_cust(runVar.outFileHandle,'Data packet sent. | (%2.0f:%2.0f:%2.0f) \n',tNow(4:6));
end
printCheckpoint('9',runVar.outputName,opts.disp)

end
%% Run completion
printCheckpoint('10',runVar.outputName,opts.disp)

%       -----------------------------------------------------
% ----- Cycle through each secondary worker and look for data ------
%       -----------------------------------------------------

fprintf_cust(runVar.outFileHandle,'Cycling through slave-labs to complete outstanding data transfers.\n')
t3 = tic; %end timer
if labindx == 1 && opts.parMode
	otherLabStat(1) = 11; %Set successful exit for lab 1
	while any(otherLabStat~=11) %If any lab still hasn't successfully exited run, keep trying to receive data. If any "otherLabStat" is not 11, then there is still at least 1 labSend left to go.
        if labProbe('any')
			[dat,srcIndx] = labReceive('any');
			if isnumeric(dat)
				otherLabStat(srcIndx) = dat;
			else
				ptNew     = dat{1};
				logPNew   = dat{2};
				ptUniqNew = dat{3};
				pt_n = find(isnan(logPLocal),1,'first');
				pt_n_New = length(logPNew);
				ptLocal((pt_n+1):(pt_n+pt_n_New),:) = ptNew(1:pt_n_New,:);
				ptUniqLocal((pt_n+1):(pt_n+pt_n_New)) = ptUniqNew(1:pt_n_New,:);
				logPLocal((pt_n+1):(pt_n+pt_n_New)) = logPNew(1:pt_n_New,:);
				if ~strcmpi(opts.disp,'off')
					tNow = clock;
					fprintf_cust(runVar.outFileHandle,'Data packet received (pt_n = %d) from slave %d. | (%2.0f:%2.0f:%2.0f) \n',nansum(ptUniqNew),srcIndx,tNow(4:6));
				end
			end
		end
		pause(1)
		if toc(t3)> 600
			fprintf_cust(runVar.outFileHandle,'Waiting too long, quiting...\n');
			break
		elseif mod(toc(t3),30)>0 && mod(toc(t3),30) <1
			fprintf_cust(runVar.outFileHandle,'Still waiting on labs...');
			fprintf_cust(runVar.outFileHandle,':%d',find(otherLabStat~=11));
			fprintf_cust(runVar.outFileHandle,'\n');
			labSend(status,2:numlab,2)
		end
	end
elseif opts.parMode
	labSend(11,1); %Send successful exit message back to lab 1.
end

printCheckpoint('11',runVar.outputName,opts.disp)
    
fprintf_cust(runVar.outFileHandle,'--------------------------------\n');
fprintf_cust(runVar.outFileHandle,'------------Run Ended-----------\n');
fprintf_cust(runVar.outFileHandle,'--------------------------------\n\n\n');

if ~isempty(runVar.outFileHandle) && runVar.outFileHandle~=1
    fclose(runVar.outFileHandle);
end

if ~opts.parMode
    status = {status};
    ptLocal = {ptLocal};
	ptUniqLocal = {ptUniqLocal};
    logPLocal = {logPLocal};
end %Make results into cells for numlabs = 1 so finished results can be handled the same way regardless of parallel run or not
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Sub function starts below %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [runVar,opts] = MCMCEvolve(runVar,opts)

% MCMCKernel generates and tests the new point point using the hasting
% ratio for MCMC. Point proposal is conducted

% Take out before previous 
pt0   = runVar.pt;
logP0 = runVar.logP;

printCheckpoint('2.1.1',runVar.outputName,opts.disp);
% Generate next point
[pt1,pdfBias] = opts.propDis('newpt',runVar);
printCheckpoint('2.1.2',runVar.outputName,opts.disp);
%Check point still in boundary
if ~isempty(runVar.bnd)
	if sum(pt1 < runVar.bnd(:,1)) || sum(pt1 > runVar.bnd(:,2))
		error('MCMC:boundaryBreached','The program ran outside the boundary. Likely to be a bug in the proposal distribution')
	end
end
printCheckpoint('2.1.3',runVar.outputName,opts.disp);
% Obtain likelihood of proposed step
logP1 = runVar.obj(pt1);      
if isnan(logP1)
	logP1 = Inf;
end
printCheckpoint('2.1.4',runVar.outputName,opts.disp);
% Metropolis Algorithm andABC Rejection when less than threshold 
test = min([1 exp((logP0-logP1)/opts.T)/... % Calculate hastings term to make metropolic, set
	                    pdfBias]);          % pdfbias to one in your proposal distribution function
thres = rand(1);   % Randomise threshold
runVar.ptTest(end+1) = thres < test;
runVar.ptTest(1) = [];
% absthres = 1;
% if logP1 < absthres
% 	runVar.ptTest = true;
% end
printCheckpoint('2.1.5',runVar.outputName,opts.disp);
% Replace current point
if runVar.ptTest(end)
	runVar.pt = pt1;
	runVar.logP = logP1;
else
	runVar.pt = pt0;
	runVar.logP = logP0;
end
runVar.delPt = pt1-pt0;

runVar = opts.propDis('adapt',runVar,opts);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pt,logScale] = seedPt(runVar)

% Determine whether to use logarithmic scale to sample boundary or to use
% linear
bndRng = (runVar.bnd(:,2)-runVar.bnd(:,1));
bndRng(isinf(bndRng)) = 1; %Unbounded ones are set to 1

%Determine scale
logRng = log10(runVar.bnd(:,2))-log10(runVar.bnd(:,1));
%          - magnitude of range > 1
%          - magnitude of range non-imaginary (i.e. boundary crosses zero)
%          - magnitude of range not infinity (i.e. one boundary IS zero)
logScale = ((logRng>1&imag(logRng)==0)&(~isinf(logRng)));

ptMidLin = ((runVar.bnd(:,2)-runVar.bnd(:,1))/2+runVar.bnd(:,1));
ptMidLog = 10.^((log10(runVar.bnd(:,2))-log10(runVar.bnd(:,1)))/2+log10(runVar.bnd(:,1)));

% Generate random variable. Undirected
prop = (rand(size(runVar.bnd(:,1)))-0.5);
pt(~logScale,1) = ptMidLin(~logScale)+prop(~logScale).*bndRng(~logScale);
pt( logScale,1) = ptMidLog(logScale).*(10.^(prop(logScale).*logRng(logScale)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function warningSwitch(trip)

% Place list of warnings to switch below
warning(trip,'parseModel:PreparsedModel');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function printCheckpoint(num,outputName,disOpts)

if strcmpi(disOpts,'debug') && ~isempty(outputName)
    chkPtFileHndl = fopen([outputName(1:(end-4)) 'checkPoint.txt'],'wt');
    fprintf(chkPtFileHndl,'%s',num);
    fclose(chkPtFileHndl); 
end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = fprintf_cust(outFileHandle,text,varargin)
	if ~isempty(outFileHandle)
        pause(0.1)
		out = fprintf(outFileHandle,text,varargin{:});
	end
end