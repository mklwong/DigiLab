function [ptRej,PRej,status] = nestSample(model,nIter,T,U,varargin)

% nestSample explores the parameter space of objFun using the nested 
% sampling algorithm.
%
%	[pstrir,P,status] = nestSample(objFun,nIter);;

% Run settings
nActive = nIter;	%number of active points
nMCMC	= 20;		%iterations for Monte-Carlo Markov-Chain
nDisp   = 100;		%Number of outputs to print over the run
acptRto	= 0.5;      %The program will adjust its stepsize to retain this 
                    %acceptance ratio

% This program is coded using -log(P) to avoid floating point errors.
                    
%% Kernel %%
% Start run timer
t1 = tic;

% Setup %
if ischar(model)
    model = str2func(model);
end

% Get parameter space info
[~,pTot,bnd,~] = model('info');
lb     = bnd{1};
ub     = bnd{2};
logScl = bnd{3};

% Processing of MCMC run parameters
sprdi   = exp(-log(nActive)/pTot);      %Starting width of proposal distribution function
rjtMax = round(1/acptRto);     %Calculate maximum number of rejections from 
                               %acceptance ratio. Rounded for integer

% Randomly create active points
ptAct = 10.^(rand(nActive,pTot).*(ones(nActive,1)*(log10(ub)-log10(lb)))+(ones(nActive,1)*log10(lb)));
pt2 = rand(nActive,pTot).*(ones(nActive,1)*(ub-lb))+ones(nActive,1)*lb;
ptAct(:,~logScl) = pt2(:,~logScl);

% Determine likelihood of initial active points
spmd
    labPtNo   = floor(nActive/numlabs);
    PActStore = nan(labPtNo,1);
    dispCnt = 0;
    for ii = 1:length(PActStore)
        t2 = tic;
        PActStore(ii) = modelObjective(model,[],ptAct((labindex-1)*labPtNo+ii,:),U);
        if toc(t2)>5
            storeError(model,[],ptAct((labindex-1)*labPtNo+ii,:),'nestedSampling taking too long at this point')
            error('nestSample:lagAtPt','nestedSampling taking too long at this point')
        end
        if mod(ii,ceil(labPtNo/nDisp))==0
            dispCnt = dispCnt + 1;
            fprintf('Active points %3.0f%% created after %7.1f seconds.\n',dispCnt/nDisp*100,toc(t1))
        end
    end
    if numlabs == 1
         PActStore = {PActStore};
    end
end

PAct = nan(nActive,1);
i1 = 1;
i2 = 0;
for ii = 1:length(PActStore)
    i2 = i2 + length(PActStore{ii});
    PAct(i1:i2) = PActStore{ii};
    i1 = i1 + length(PActStore{ii});
end

PAct = PAct/T;

fprintf('Active point likelihoods found.\n')

t1 = tic;
% Initialise rejected pool
PRej	= nan(nIter,1);
ptRej	= nan(nIter,pTot);
dispCnt = 0;
        
% Begin random sampling %
for ii = 1:nIter
    % Find and store smallest likelihood from active points
	PRej(ii) = max(PAct);
    minIndx = find(PRej(ii)==PAct,1,'first');
    ptRej(ii,:) = ptAct(minIndx,:);
    sprd = sprdi;
	
    % Remove found smallest likelihood from active points and replace with
    % NaN for now.
    ptAct(minIndx:end,:) = [ptAct(minIndx+1:end,:);NaN*ptRej(ii,:)];
    PAct(minIndx:end) = [PAct(minIndx+1:end);NaN];
    ptSeed = NaN*ptRej(ii,:);
    % Pick random active point to seed next active point (reroll if NaN
    % accidentally found)
    while prod(double(isnan(ptSeed)))       
        apSeed = ceil(rand(1)*(nActive-1));
        ptSeed = ptAct(apSeed,:);
    end
    
    PSeed  = modelObjective(model,[],ptSeed,U)/T;
    
    % MCMC Step %
    for jj = 1:nMCMC
        while PSeed > min(PRej)
            [ptSeed PSeed sprd] = MCMCKernel(model,ptSeed,PSeed,sprd,lb,ub,logScl,U,T,rjtMax);
            if toc(t2) > 60
                keyboard
            end
        end
    end
    t2 = tic;

    % Add new point to active point
    ptAct(end,:) = ptSeed;
    PAct(end) = PSeed;
	
	if mod(ii,ceil(nIter/nDisp))==0
		dispCnt = dispCnt + 1;
		fprintf('%3.0f%% done after %7.1f seconds.\n',dispCnt/nDisp*100,toc(t1))
	end
end

status = 0;