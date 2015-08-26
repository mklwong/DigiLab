function [resid status]= modelObjective(model,p,U,debug)

% resid = modelObjective(model,x0,p,U)
%
% This function finds the residual for topology file model, given parameter
% set p and inputs and objectives given by U.
%
% The output is log10 of (sum of the ((error divided by std)squared)).
%
% Obj = 0.5 is a pretty good objective for this scheme

% Dose response has no error checking mechanism yet.
% Also work to unify time courses and dose response cases not done yet.

%% Compile model
modType = modelType(model);
model = parseModel(model,p);

%% Unify data if not already done
% To be filled in


%% Get steady state
[~,y0_default,status] = findSS(model,0,{1});

if status < 0
	resid = Inf;
	return
end

%% Begin Fitting
ChiSq = zeros(1,length(U));  %Initialise time course Objective vector
normExp = [];                %the vector that stores normalisation factors

for ii = 1:length(U)
    
    %% Find and match simulated and experimental states   
    objProt = U(ii).states;      % Obtain protein names from objective
    modOut = model.pFit.sim2dat; % Get labels to match to experimental states
	modProt = model.conc.name;      % Get protein names in simulation
    [matches,simIndx,expIndx] = intersect(upper(modOut(:,1)),upper(objProt));
	
	%expIndx is the row index of the experimental data that is fed in U
	%simIndx is the row index of the label cell array associated with the
	%    expIndx state.
	%
	%In other words:
	%     modOut = {'State1',[simStates];      objProt = {'Name1'; data = [row1; 
	%               'State2',[simStates];                 'Name2';         row2; 
	%               'State3',[simStates];      expIndx -> 'Name3';         row3; 
	%     simIndx-> 'State4',[simStates];                 'Name4';         row4; 
	%               'State5',[simStates]};                'Name5'};        row5]; 
	% If State4 is the same as Name3
	
    %% Run simulation  
    % Simulate
	x = U(ii).x;
	
	% Get steady state if an initial condition has been passed
	if ~isempty(U(ii).x0)
		y0 = model('rxn',p,U(ii).x0);
		[y0,status] = findSS(model,[0 max(vertcat(U.x))],y0,p);
		if status < 0
			resid = Inf;
			return
		end
	else
		y0 = y0_default;
	end

	if strcmpi(U(ii).type,'time')
        inp = U(ii).input;

		try
			[xSim,YSim,YSim2] = findTC(model,[0 max(x)-x(1)],'x0',y0,'inp',inp,'-b');
		catch msg
			storeError(model,y0,p,msg,[msg.identifier '; ' msg.message]);
			xSim  = [x(1) max(x)];
			YSim  = Inf*ones(2,length(y0));
			YSim2 = YSim;
			YSim  = YSim(:,model.x.comp>=0);  %remove complex states from Y
		end

    elseif strcmpi(U(ii).type,'dose')
		if size(U(ii).x,1)~= size(inpIndx)
			dspan = U(ii).x';
			if size(inpIndx) == 1
				inpIndx = {inpIndx};
			else
				inpIndx = mat2cell(inpIndx);
			end
		else
			dspan = U(ii).x;
		end
        [xSim,YSim,YSim2] = findDR(model,dspan,inpIndx,'x0',y0);
	end
	
    %% Extract states
    YSimTest = extractSim(modOut(simIndx,2),YSim,YSim2);
    [YExpMean YExpStd] = extractExp(expIndx,U(ii).mean,U(ii).std);
    
    %% Create + extract normalising factors
    if strcmpi(U(ii).scale{1},'rel')
        normExp = getNormExp(normExp,U(ii).scale{2},YSimTest,YExpMean,matches);
        [~,normIndx] = intersect(normExp{1,U(ii).scale{2}},matches);
        normNow = normExp{U(ii).scale{2},2};
        normNow = normNow(normIndx);
    elseif strcmp(U(ii).scale{1},'abs')
        normNow = 1;
    end

    %% Determine objective value using normalising factor
    ChiSq(ii) = objValueExp(interp1(xSim,YSimTest,x),YExpMean,YExpStd,normNow);
	%Penalise normalisation factor if it's too small.
	if ChiSq(ii) == 0 && (1+sum(1./(normNow)))==Inf
		ChiSq(ii) = Inf;
	else
		ChiSq(ii) = ChiSq(ii)*(1+sum(1e-4./(normNow)));
	end
    
    %% Visualisation
	if exist('debug','var')
		if debug
			YSim = extractSim(modOut(simIndx,2),YSim,YSim2);
			visObj(xSim,YSim,U(ii).x,YExpMean,YExpStd,normNow,U(ii).type,matches,U(ii).input);
		end
	end
end

%% Sum Final Objective
pFinal = sum(ChiSq);
resid = pFinal;

end