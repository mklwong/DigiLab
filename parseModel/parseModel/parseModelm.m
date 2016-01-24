function model = parseModelm(model,rxnRules,flags)

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

rxn(1) = []; %When rxn is initialised it always creates an empty reaction. Remove this.

run(modelname); 
%loads the following 
%	- modComp: compartment info for model.
%	- modSpc:  species infor for model
%	- dataSpc: how data readout and model species are related
%	- Bnd:	   default boundaries for parameter types
%	- rxn:	   list of reactions of the model

% Legacy code
if size(v) > 1
    error('Please change all v in your model file into rxn. v no longer recognised as reaction network variable')
end

%Convert all compartment and species arrays into structures
modSpc = cell2struct(modSpc,{'name','comp','matVal'},2);
modComp = cell2struct(modComp,{'name','matVal'},2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialise Output Structs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initalise parameter vector
pFit.desc = cell(100,1); pFit.lim  = nan(100,2); pFit.npar = 0; pFit.grp = zeros(0,2);
paramGrp = [-1,-1];

for ii = 1:length(param) % Create pInd for all params.
	param(ii).pInd = nan(size(param(ii).matVal));
end
param = expandTens(param);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cycle over compartments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:size(modComp)
	% Process parameter
    [modComp(ii),pInd,pFit] = testPar2(modComp(ii),pFit);
	
    if ~isnan(pInd(2))
        parDesc = ['Comp : ' modComp(ii).name];
        pFit.desc{pInd(2)} = parDesc;
        pFit.lim(pInd(2),:) = Bnd.Comp;
    end
    
    modComp(ii).pInd = pInd(2);
end		
% Constract modComp
tmpComp.name   = vertcat({modComp.name})';
tmpComp.matVal = vertcat(modComp.matVal);
tmpComp.pInd   = vertcat(modComp.pInd);

modComp = tmpComp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cycle over list of species
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:size(modSpc)
	% Process parameter
    [modSpc(ii),pInd,pFit] = testPar2(modSpc(ii),pFit);
	
    if ~isnan(pInd(3))
        parDesc = ['Conc : ' modSpc(ii).name];
        pFit.desc{pInd(3)} = parDesc;
        pFit.lim(pInd(3),:) = Bnd.Conc;
    end
    
    modSpc(ii).pInd = pInd(3);
end
% Constract modComp
tmpSpc.name   = vertcat({modSpc.name});
[~,tmpSpc.comp] = ismember(vertcat({modSpc.comp}),modComp.name);
tmpSpc.comp = tmpSpc.comp';
tmpSpc.matVal = vertcat(modSpc.matVal);
tmpSpc.pInd   = vertcat(modSpc.pInd);

modSpc = tmpSpc;

%Substitute in compartment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cycle over list of Reactions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

paramNames = {param.name}; %get list of parameter names
BndCell = [fieldnames(Bnd) struct2cell(Bnd)];
for ii=1:length(rxn)
	
	%% Determine if parameters are free or not
    [rxn(ii),pInd,pFit] = testPar2(rxn(ii),pFit);
	
	%% Turn reactions into maths using reaction rules
    [reqParam,matVal,pBuild,parDesc,modSpc] = model.rxnRules('rxnRules',rxn(ii),modSpc,flags,ii);
	
    % Get definition of parameters
    [~,pList2Desc] = ismember(parDesc(:,1),fieldnames(rxn(ii))); % parDesc to pInd mapping
     
    % Remove non-free parameters from parDesc and insert definition
    parDesc(isnan(pInd(pList2Desc)),:) = [];    %Remove parameters that are not free from parDesc
    pList2Desc(isnan(pInd(pList2Desc))) = [];   %Remove equivalent indices from pList2Desc
    % Account for grouped parameters that have popped up again. They need
    % to be skipped so "leading description" is not replaced.
    for jj = 1:length(pList2Desc)
        pFit.desc{pInd(pList2Desc(jj))} = [pFit.desc{pInd(pList2Desc(jj))} ' | ' parDesc{jj,3}]; %Insert definition
    end
    
    % Loop through pre-matrices that require new additions and build the
    % pre-matrix
    for jj = 1:length(reqParam);
		%Make param builder
        pIndBuild = nan(size(pBuild{jj}));
		[~,pBuild2Params] = ismember(pBuild{jj},fieldnames(rxn(ii))); % parDesc to pInd mapping
		pIndBuild(pBuild2Params~=0) = 0;
        [~,pBuild2pInd] = ismember(pBuild{jj},parDesc(:,1)); % parDesc to pInd mapping
        pIndBuild(pBuild2pInd~=0) = pInd(pList2Desc((pBuild2pInd(pBuild2pInd~=0))));
        pIndBuild = pIndBuild(ones(1,size(matVal{jj},1)),:);

        %Append both into their required pre-matrices
        [~,paramInd,~] = intersect(paramNames,reqParam{jj});
        matInd = find(isnan(param(paramInd).matVal(:,1)),1,'first');
        param(paramInd).matVal(matInd:(matInd+size(matVal{jj},1)-1),:) = matVal{jj};
        param(paramInd).pInd(matInd:(matInd+size(matVal{jj},1)-1),:)   = pIndBuild;
    end
    
    % Get boundary definition for parameter. Note this comes after making
    % pre-matrices because non-free parameters have been removed from
    % parDesc. But in the next step we will be removing more than that. So
    % compiling pre-matrix is most convenient before the boundary
    % incorporation step.
    
    %1) Check back to pFit.lim to see if any boundaries of mapped free
    %parameters are already defined (i.e. they are custom).  Remove these
    %indices.
    parDesc(~isnan(pFit.lim(pInd(pList2Desc),1)),:) = [];
    pList2Desc(~isnan(pFit.lim(pInd(pList2Desc),1))) = [];
    
    % Get the field from Bnd which corresponds to the free parameter in
    % question
    if ~isempty(pList2Desc)
        [~,bnd2Desc] = ismember(parDesc(:,2),fieldnames(Bnd)); % parDesc to pInd mapping
        pFit.lim(pInd(pList2Desc),:) = vertcat(BndCell{bnd2Desc,2});
    end
end

% remove NaN rows from generated tensors
for ii = 1:length(param)
	paramRmIndx = isnan(param(ii).matVal(:,1));
	param(ii) = contractTens(param(ii),paramRmIndx);
end
pFitRmIndx = isnan(pFit.lim(:,1));
pFit  = contractTens(pFit,pFitRmIndx);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Determine data and model relation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rmXData = [];  % List of dataSpcs to remove due to not having all species in the model.
% 
% % Loop through experimental species
% for ii = 1:size(dataSpc,1)
% 	% Loop through attached model state
% 	dataSpcTmp = zeros(size(dataSpc{ii,2}));
% 	for jj = 1:length(dataSpc{ii,2})
% 		
% 		% Determine if complex of model species needs to be included.
% 		% model species with start in them will have enzyme-substrate
% 		% complex included (flat -1).
% 		if strcmp(dataSpc{ii,2}{jj}(1),'*')
% 			incComp = -1;
% 			dataSpc{ii,2}{jj}(1) = [];
% 		else
% 			incComp = 1;
% 		end
% 		
% 		% Get species name index
% 		[~,~,indx]  = intersect(upper(dataSpc{ii,2}{jj}),upper(conc.name));
% 		
% 		% If required model species not in the simulation. Print warning and then exclude
%         if isempty(indx)
% 			storeError(model,[],[],[],['The state ' dataSpc{ii,2}{jj} ' required as an experimental equivalent state not found. This experimental state will be excluded from curve fitting'])
%             rmXData = ii;
% 			break
% 		end
% 		
% 		% Store index of model species including flag of whether to include
% 		% enzyme-substrate complex or not.
% 		dataSpcTmp(jj) = incComp*indx;
% 	end
% 	dataSpc{ii,2} = dataSpcTmp;
% end
% dataSpc(rmXData,:) = [];
% pFit.sim2dat = dataSpc;

%% Compile output
model.modSpc = modSpc;
model.pFit = pFit;
model.param = param;
model.modComp = modComp;
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
		if ~ismember(paramFields{jj},{'name','npar','grp'}) % Fields to not expand or contract
			tmpTens = param.(paramFields{jj});
			tmpTens(ind,:) = [];
			param.(paramFields{jj}) = tmpTens;
		end
	end
end