function model = parseModelm(model,rxnRules,expComp)

%   out = parseModelm(model,debug)
%

% parseRxn mostly complete. Now completing inserting the parseRxn generated
% tensor into the full tensor. Also need to consider how pInd works in the
% case of rate parameter tensors. Try and unify it to how it works in the x
% case.
%
% Biggest concerns seems to be the dynamic variable name on the left hand
% side of assignments.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Program Internal Description %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Each parameter has six potential formats:
%	1) [val]        : Known parameter
%	2) [NaN]        : Unknown parameter, default range used
%	3) [NaN lb ub]  : Unknown parameter, custom range

%   4) [NaN grp]       : Unknown parameter that is grouped. Default range used.
%                     Multiplicative factor is assumed to be one.
%   5) [NaN grp lb ub] : Unknown parameter with other parameters with the same
%                     value, custom range.
%   6) [NaN grp+factor*i]: Unknown parameter that is a multiplicative factor of
%                     another parameter. Parameter of this format is the
%                     dependent.
% grp is a positive integer. All parameters with grp that is the same
% integer are in the same group. To distinguish it from "val", the first
% element of the vector is negative is it indicates an unknown parameter
% that is part of a group.
%
% This program processes this format by going through all parameters and
% assignment as necessary.
%
% While assigning parameters to their correct tensor locations, a set of
% support vectors are generated linking a hypothetical set of parameters to
% their locations in tensors, their sign (based on reaction rules, since 
% some tensor values need to be negative, e.g. consumption) and their index
% from the inputted parameter vector. A negative index is used to denote
% division by the parameter. E.g. Indx = 1 is el_(ij) = p(1) while Indx =
% -3 indicates el_(ij) = 1/p(3);
%
% As an example, p(1) might be a Michaelis Constant, so this generates 6
% tensor values. So the generated vectors will be:
%	 G_ID = [...;   G = [       ...          ;   G_Sign = [...;
%            -1;         S_indx,S_indx,E_indx,0;              1 ;
%            -1;         S_indx,E_indx,S_indx,0;              1 ;
%            -1;         E_indx,S_indx,E_indx,0;              1 ;
%            -1;         E_indx,E_indx,S_indx,0;              1 ;
%            -1;         C_indx,S_indx,E_indx,0;             -1 ;
%            -1];        C_indx,E_indx,S_indx,0];            -1];
%
% Finally the tensors are generated from the parameters by going:
%	G(:,4) = (p(abs(G_ID)).*(G_Sign)).^(sign(G_ID)); Which is then passed into the ode solver.
%
% When a parameter with reserved index is encountered, the pID number will
% still go up in order, but it will store the reservation index and pair it
% with the pID generated. For instance.
%	resID = [ ...  ;
%            -3,10];
%
% Now everytime reservation index -3 is encountered, the pID will not be
% incremented, instead it will set the ID of that as 10.
%
% Compartmentalisation is not yet implemented. 
%
% %%%%%%%%%%%%%Output Groupings%%%%%%%%%%%%
% The output of this function are in the form of structs.
%
% The params struct will contain a struct with 6 fields, all of which are
% themselves structs:
%	k3  : (bimolecular type tensor)
%	k2  : (unimolecular type tensors)
%	k1  : (zeroth order reaction type tensors)
%	G   : (QSS type tensors)
%	x   : (concentration related tensors)
%	pFit: (parameter fitting related tensors. This includes parameter
%	       description, limits, and relationship between model states and
%	       experimental states)
%
% The k3,k2,k1 and G structs have the following fields
%	- tens  : The tensor itself
%	- sign  : The sign of the parameter when placed in the tensor
%	- pInd  : The index of the p vector which the free tensor elements are
%	          associated with
%	- factor: What the parameter is multiplied by. Useful in the case of
%	          grouped parameters that proportional but not the same as each
%	          other (e.g. p(2) = 4*p(1), then p(2) is given a factor of
%	          four and made to equal p(1))
%
% The x struct has the following fields:
%	- names  : name of species
%	- tens   : the vector of initial concentrations
%	- comp   : compartment size of each species
%	- selfLoc: indexes within the x.tens vector that are free variables
%	- pInd   : The index of the p vector which the free tensor elements are
%	        associated with
%	- factor : For a dependent parameter, the multiplicative factor the
%	           initial condition is wrt the parent parameter.
%
% The pFit struct has the following fields:
%	- desc   : description of each parameter
%	- lim    : limits of each parameter
%	- sim2dat: Linking simulation species with experimental species
	 
%%%%%%%%%%%%%%%%%%%%%
%% Import model file
%%%%%%%%%%%%%%%%%%%%%
if ischar(model)
	modelname = model;
	clear model
	model.name = modelname;
end

% Initialise parameters
model.rxnRules = rxnRules;
param = model.rxnRules('ini');
rxn = sigRxnList();
v = rxn; %Legacy code. For backward compatibility.

run(modelname); 
%loads the following 
%	- spcComp: compartment info for model.
%	- modSpc:  species infor for model
%	- dataSpc: how data readout and model species are related
%	- Bnd:	   default boundaries for parameter types
%	- rxn:	   list of reactions of the model

rxn(1) = []; % remove 
% Legacy code
if size(v) > 1
    error('Please change all v in your model file into rxn. v no longer recognised as reaction network variable')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialise Output Structs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initalise parameter vector
pFit.desc = cell(100,1); pFit.lim  = nan(100,2);
curParInd = 0;
paramGrp = [-1,-1];

% Initialise chemical species
[a,~]     = size(spcComp);
comp.tens    = nan(a,1);    % Compartment size
comp.name    = cell(a,1);   % Compartment name
comp.pInd    = nan(a,1);    % Vector showing the parameter index a free comprtment will use

% Initialise chemical species
[a,~]     = size(modSpc);
conc.tens    = nan(a,1);    % Initial concentration
conc.name    = cell(a,1);   % Species name
conc.comp    = ones(a,1);   % Compartment Index
conc.pInd    = nan(a,1);    % Vector showing the parameter index a free state will use

for ii = 1:length(param) % Create pInd for all params.
	param(ii).pInd = nan;
end
param = expandTens(param);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cycle over compartments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:size(spcComp,1)
	% Save Name
	comp.name(ii) = spcComp(ii,1);
	
	% Process parameter
    [val,freeParam,grp] = testPar(spcComp{ii,2});
	comp.tens(ii) = val;
	
	parDesc = ['Comp : ' spcComp{ii,1}];
	if freeParam
		if length(spcComp{ii,2})==3
			custBnd = spcComp{ii,2}(2:3);
		else
			custBnd = [];
		end
		[pFit,conc,paramGrp,curParInd,putParInd] = procFreeParam(pFit,curParInd,parDesc,conc,custBnd,Bnd.Conc,paramGrp,grp);
		comp.pInd(ii) = putParInd; % store p index of the compartment is assigned to
	end
end		

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cycle over list of species
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:a
    % Save Name
	conc.name{ii} = modSpc{ii,1};
    
    % Get compartment size for each species
	[~,comptIndx] = intersect(upper(spcComp(:,1)),upper(modSpc{ii,2}));
	conc.comp(ii) = comptIndx;

    % Process parameter
    [val,freeParam,grp] = testPar(modSpc{ii,3});
    
	conc.tens(ii) = val; % Save either parameter value of parameter multiplicative factor
    parDesc = ['Conc : ' modSpc{ii,1}];
	if freeParam
		if length(modSpc{ii,3})==3
			custBnd = modSpc{ii,3}(2:3);
		else
			custBnd = [];
		end
		[pFit,conc,paramGrp,curParInd,putParInd] = procFreeParam(pFit,curParInd,parDesc,conc,custBnd,Bnd.Conc,paramGrp,grp);
		conc.pInd(ii) = putParInd; % store p index of the x0 is assigned to
	end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cycle over list of Reactions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii=1:length(rxn)
	
	%% Determine if parameters are free or not
    % Test parameter for 'k'
    [rxn(ii).k,freeParam(1),grp(1),bnd{1}] = testPar(rxn(ii).k);
	
    % Test parameter for 'Km'
    [rxn(ii).Km,freeParam(2),grp(2),bnd{2}] = testPar(rxn(ii).Km);
	
	% Test parameter for 'n'
    [rxn(ii).n,freeParam(3),grp(3),bnd{3}] = testPar(rxn(ii).n);
	
	% Test parameter for 'A'
    [rxn(ii).r,freeParam(4),grp(4),bnd{4}] = testPar(rxn(ii).r);
	
	%% Turn reactions into maths using reaction rules
    [reqTens,tensVal,parDesc,conc,comptMod] = model.rxnRules('rxnRules',rxn(ii),conc,comp,expComp,ii);
	
	tensNames = {param.name};
	
	% Loop through new additions 
    for jj = 1:length(reqTens);
		
		[~,reqIndx] = intersect(upper(tensNames),upper(reqTens{jj}));
		
		% Expand tensor as required
		tensLength  = find(isnan(param(reqIndx).tens(:,1)),1,'first'); %Length of current tensor
		appndLength = size(tensVal{jj},1)-1;                      %Length to be appended
		appndIndx   = tensLength:(tensLength+appndLength);        %Indicies in current tensor new values to be appended to
		if appndLength + tensLength + 10 > length(param(reqIndx).pInd)
			param(reqIndx) = expand(param(reqIndx));
		end
		
		% Insert parameter indicies for free parameters
		if freeParam(jj)
			[pFit,param,paramGrp,curParInd,putParInd] = procFreeParam(pFit,curParInd,parDesc{jj},param,bnd{jj},Bnd.(reqTens{jj}),paramGrp,grp(jj));
			param(reqIndx).pInd(appndIndx) = putParInd;               %Append parameter index
		else
			param(reqIndx).pInd(appndIndx) = NaN;                     %Mark no parameter needed
		end
		param(reqIndx).tens(appndIndx,:) = tensVal{jj};               %Append tensor
    end
end

% remove NaN rows from generated tensors
for ii = 1:length(param)
	paramRmIndx = isnan(param(ii).tens(:,1));
	param(ii) = contractTens(param(ii),paramRmIndx);
end
pFitRmIndx = isnan(pFit.lim(:,1));
pFit  = contractTens(pFit,pFitRmIndx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine data and model relation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmXData = [];  % List of dataSpcs to remove due to not having all species in the model.

% Loop through experimental species
for ii = 1:size(dataSpc,1)
	% Loop through attached model state
	dataSpcTmp = zeros(size(dataSpc{ii,2}));
	for jj = 1:length(dataSpc{ii,2})
		
		% Determine if complex of model species needs to be included.
		% model species with start in them will have enzyme-substrate
		% complex included (flat -1).
		if strcmp(dataSpc{ii,2}{jj}(1),'*')
			incComp = -1;
			dataSpc{ii,2}{jj}(1) = [];
		else
			incComp = 1;
		end
		
		% Get species name index
		[~,~,indx]  = intersect(upper(dataSpc{ii,2}{jj}),upper(conc.name));
		
		% If required model species not in the simulation. Print warning and then exclude
        if isempty(indx)
			storeError(model,[],[],[],['The state ' dataSpc{ii,2}{jj} ' required as an experimental equivalent state not found. This experimental state will be excluded from curve fitting'])
            rmXData = ii;
			break
		end
		
		% Store index of model species including flag of whether to include
		% enzyme-substrate complex or not.
		dataSpcTmp(jj) = incComp*indx;
	end
	dataSpc{ii,2} = dataSpcTmp;
end
dataSpc(rmXData,:) = [];
pFit.sim2dat = dataSpc;

%% Compile output
model.conc = conc;
model.pFit = pFit;
model.param = param;
model.comp = comp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%End Main Function%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
function param = expandTens(param)
	paramFields = fieldnames(param);
	for ii = 1:length(param)
		for jj = 1:length(paramFields)
			if ~ismember(paramFields{jj},{'name'}) % Fields to not expand on contract
				arrayWidth = size(param(ii).(paramFields{jj}),2);
				param(ii).(paramFields{jj}) = [param(ii).(paramFields{jj}) ; nan(200,arrayWidth)];
			end
		end
	end
end

function param = contractTens(param,ind)
	paramFields = fieldnames(param);
	for jj = 1:length(paramFields)
		if ~ismember(paramFields{jj},{'name'}) % Fields to not expand on contract
			tmpTens = param.(paramFields{jj});
			tmpTens(ind,:) = [];
			param.(paramFields{jj}) = tmpTens;
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%
function [pFit,param,paramGrp,curParInd,putParInd] = procFreeParam(pFit,curParInd,parDesc,param,custBnd,defBnd,paramGrp,grp)

if ~all(grp~=paramGrp(:,1))            %%existing group of free parameters
	putParInd = paramGrp(grp==paramGrp(:,1),2);
else                             %% ungrouped or new group of free parameter
	curParInd = curParInd + 1; % New parameter required

	%Make new group
	if grp ~= 0                  
		paramGrp(end+1,:) = [grp curParInd];
	end

	% Set parameter boundary
	if ~isempty(custBnd)       %% If there is custom set boundary
		pFit.lim(curParInd,:) = custBnd;
	else                             %% Else use default boundary
		pFit.lim(curParInd,:) = defBnd;
	end
	pFit.desc{curParInd}	= parDesc;
	putParInd = curParInd;
end
	
end