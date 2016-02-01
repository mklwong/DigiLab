% MCMCSeriel Documentation.
%
% This script utilises the MCMC algorithm in seriel mode to perform
% parameter fitting of the experimental data contained within DAT to the
% simulated data with topology given by MODEL.
%
% The program will run an MCMC from a user defined starting tempering
% temperature, incrementing down exponentially down to a temperature of 1. 
% The number of increments is also user defined.
%
% It is possible to provide a
% prior for the program to run on, but only as a series of points
% distributed based on the prior distribution, with the subsequent MCMC 
% accounting only for the likelihood, not the posterior. As a result, care
% must be taken when using this program with a prior. A way to ensure this
% problem is not 

%% Kernel
% if exist(data,'file')
% 	load(data);  % Load objective
% elseif exist(['objectiveData\' data],'file')
% 	load(['objectiveData\' data])
% else
% 	error('Objective file not found')
% end

% Get model name
if ischar(model)
    modName = model;
elseif isa(model,'function_handle')
    modName = func2str(model);
end

%Check for optional entries. Otherwise, impose as required.
fprintf('-------------------------------\n')
fprintf('--------Run Information--------\n')
fprintf('-------------------------------\n')
% Get model locatoin
modelLoc = which(modName);
if isempty(modelLoc)
    fprintf('Model not found. Quitting...\n')
    return
else
    fprintf('Testing Model: %s\n',modName)
end
modelLoc(end-length(modName)-1:end) = [];

% Load prior file (if exists)
prirDat = strrep(prirDat,'/','\'); %Change forward slash to backslash.
if isempty(prirDat)
    fprintf('No prior given.\n')
elseif exist([modelLoc prirDat],'file') == 2
    load([modelLoc prirDat]); % Load prior file
    fprintf('prior file Loaded.\n')
else
	error('Specified prior not found')
end

% Generate random string for run ID and create unique folder for run
% results
[~,sysName] = system('hostname');   % Get hostname for number generator
sysName(end) = [];
symChar  = ['a':'z' 'A':'Z' '0':'9'];        %possible symbols
c = clock;
rng('shuffle'); %shuffle the random generator

runID = [num2str(c(3)) '-' num2str(c(2)) '_' num2str(c(4)) '-' num2str(c(5)) '_' sysName '_' symChar(randi(length(symChar),1,5))]; %random string of 5
if ~exist([modelLoc '/' runID],'file')
	mkdir(modelLoc,runID);                       %Make subdirectory
end

% If no prior, initialise.
if ~exist('result','var')
    prir = struct();
    prir.logP = [];
    prir.pts  = [];
    prir.T    = 1;
else
    prir = result;
	fprintf('prior entered.\n')
end

%Determine Scheduler by power factor x^n. The stps-1 is because the last value
%needs to be one, which corresponds to x^0. T vector goes from largest to
%smallest
if stps > 1
	dstp = stps/(stps-1);
	fac = exp(log(TSpan)/(stps));
	T = fac.^(linspace(stps,0,stps));
	fileLocs = cell(1,length(T));
else
	T = TSpan;
end

% Print tempering schedule
fprintf('Tempering Schedule:\n')
fprintf('%6.2f  ',T)
fprintf('\n')

% Parse Model
model = parseModel(model);
if ~exist('objFun')
	objFun = @(p) modelObjective(model,p,U);
end
if ~exist('opts','var')
	opts = MCMCOptimset();
end
opts = MCMCOptimset(opts,'Pmin',Pmin,'parmode',ParMode,'PtNo',ptNo,'dir',[modelLoc,runID],'display',dispOpt);
%%=======================================================================%%
%%=======================START RUN=======================================%%
%%=======================================================================%%

ii = 1;  %Temperature schedule counter

while ii ~= length(T)+1
    
    %=======    Timer setup     =======%
    t1 = tic;

    %=======     Run MCMC       =======%
    fprintf('\n**** MCMC Start ****\n')
    fprintf('Tempering at T = %6.2f\n',T(ii))
    opts = MCMCOptimset(opts,'T',T(ii)','Prior',prir);
    [pts,logP,status] = MCMC(objFun,[],model.pFit.lim,opts);
    
    fprintf('Temperature %6.2f done after %7.1f seconds.\n',T(ii),toc(t1))
    
    %======= Save the posterior =======%
    fileName = ['T=' num2str(T(ii),'%6.2f') '.mat'];
	if exist([modelLoc '/' runID '/' fileName],'file')
        fprintf('Existing results file found for this temperature. Load...\n')
        load([modelLoc '/' runID '/' fileName]);
        fprintf('Appending new results...\n')
        logPnew = [result.logP; logP];
        ptsNew  = [result.pts ; pts ];
    else
        logPnew = logP;
        ptsNew = pts;
	end
	result.logP  = logPnew;
	result.pts   = ptsNew;
	result.T     = T(ii);
	result.model = model;
	result.best  = pts(logP==min(logP),:);
    if exist('U','var')
        result.obj   = U;
    end
    save([modelLoc runID '/' fileName],'result');
    fileLocs{ii} = fileName;

    %======= Decide next Temperature =======%
    % Decide the next step based on whether the MCMC completed or not.
    % Order in the below loop is:
    %   - Failed at highest temperature (go higher)
    %   - Failed at other temperature   (try between 1 higher and current)
    %   - Succeed                       (progress to next temperature)
    if status == -1
        if ii == 1
            fprintf('MCMC failed at highest temperature. Going higher!\n')
            T = sort([T fac^(log(max(T))/log(fac)+dstp)],'descend');
            fileLocs = [fileLocs(1:ii-1) {[]} fileLocs(ii:end)];  %make space for new temperature in list of save files
            fprintf('Tempering Schedule:\n')
            fprintf('%6.2f  ',T)
            fprintf('\n')
            status = 0;
        else
            T = sort([T fac.^((log(T(ii-1))/log(fac)+log(T(ii))/log(fac))/2)],'descend');
            if (T(ii-1)-T(ii))/T(ii) < 0.05
                fprintf('MCMC failed. New temperature not significantly different to last sucessful temperature.\n')
                fprintf('QUITTING\n')
                break
            end
            fprintf('MCMC failed. Trying a a higher temerature T = %f.\n',T(ii))
            fileLocs = [fileLocs(1:ii-1) {[]} fileLocs(ii:end)];  %make space for new temperature in list of save files
            fprintf('Tempering Schedule:\n')
            fprintf('%6.2f  ',T)
            fprintf('\n')
            status = 0;
            % Load a prior if it exists
            if ~isempty(fileLocs{ii-1})
                load([modelLoc runID '/' fileLocs{ii-1}]);
                prir = result;
            end
        end
    else
        ii = ii + 1;
        prir = result;
    end
    fprintf('****Run Complete****\n')
end

fprintf('Seriel Mode Complete.\n')
% close parallel computing pool
parComp('close')