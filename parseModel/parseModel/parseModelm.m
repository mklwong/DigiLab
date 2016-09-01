function model = parseModelm(model,flags)

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
%   6) [factor grp]: Unknown parameter that is a multiplicative factor of
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
modelname = model;
clear model
model.name = '';

for ii = 1:length(modelname)
	model.name = [model.name modelname{ii} '|'];
	% Initialise parameters
	rxn = sigRxnList();
	rxn(1) = []; %When rxn is initialised it always creates an empty reaction. Remove this.
	run(modelname{ii}); 

	%loads the following 
	%	- modComp: compartment info for model.
	%	- modSpc:  species infor for model
	%	- dataSpc: how data readout and model species are related
	%	- Bnd:	   default boundaries for parameter types
	%	- rxn:	   list of reactions of the model
	
	if ii == 1
		modComp_all  = modComp;
		modSpc_all   = modSpc;
		Bnd_all      = Bnd;
		rxn_all      = rxn;
		rxnRules_all = rxnRules;
	else
		% Check for same compartment entry between models
		[jj,kk] = ismember(modComp(:,1),modComp_all(:,1));
		for mm = 1:length(kk)
			if jj(mm) % Verify if entries are the same
				misMatches = modComp{mm,2}~=modComp_all{kk(mm),2}; 
				if misMatches % Verify if not the same
					[newVal,strOut] = compareParam(modComp{mm,2},modComp_all{kk(mm),2});
					if isempty(newVal)
						error('parseModelm:UnresolvableParamClash',['Unresolvable clash found in ' modComp{mm,1} ' in model: ' modelname{ii}  ', ' stOut])
					else
						warning(['Mismatch in ' modComp{mm,1} ': ' strOut])
					end
					modComp_all{kk(mm),2} = newVal;
				end
			end
		end
		modComp(jj,:) = [];  % Remove matching duplicates
		modComp_all = [modComp_all;modComp];
		
		% Check for conflicts between models in model species initial
		% condition
		[jj,kk] = ismember(modSpc(:,1),modSpc_all(:,1));
		for mm = 1:length(kk)
			if jj(mm) % Verify if entries are the same
				misMatches = modSpc{mm,3}~=modSpc{kk(mm),3}; 
				if misMatches % Verify if not the same
					[newVal,strOut] = compareParam(modSpc{mm,3},modSpc{kk(mm),2});
					if isempty(newVal)
						error('parseModelm:UnresolvableParamClash',['Unresolvable clash found in ' modSpc{mm,1} ' in model: ' modelname{ii}  ', ' stOut])
					else
						warning(['Mismatch in ' modSpc{mm,1} ': ' strOut])
					end
					modSpc_all{kk(mm),3} = newVal;
				end
			end
		end
		
		% Find mismatches in compartment name
		jjIndx = find(jj);
		misMatches = ~ismember(modSpc(jjIndx,2),modSpc_all(kk(jj),2)); % Identify where mismatches are
		if any(misMatches)
			misMatchNames = modSpc(jjIndx(misMatches),1);
			misMatchNames(:,2) = {', '};
			misMatchNames(end,2) = {''};
			misMatchNames = misMatchNames';
			error(['The following species have conflicting compartments: ' horzcat(misMatchNames{:}) ', in model: ' modelname{ii} '. Please fix and ensure they are consistent.'])
		end
		
		modSpc(jj,:) = [];  % Remove matching duplicates
		modSpc_all = [modSpc_all;modSpc];
		
		% Check for conflicts between models in boundary conditions
		bndNames = fieldnames(Bnd);
		bndAllNames = fieldnames(Bnd_all);
		for jj = 1:length(bndAllNames)
			[isMemb,membInd] = ismember(bndAllNames(jj),bndNames);
			if isMemb
				bndEqual = Bnd_all.(bndAllNames{jj}) == Bnd.(bndNames{membInd});
				if ~all(bndEqual)
					warning(['The following boundary: ' bndAllNames{jj} ', is inconsistent. Boundary encompassing both are set.'])
					Bnd_all.(bndAllNames{jj}) = [min([Bnd_all.(bndAllNames{jj})(1) Bnd.(bndNames{membInd})(1)]) max([Bnd_all.(bndAllNames{jj})(2) Bnd.(bndNames{membInd})(2)])];
				end
			end
		end
		% Check for boundaries that are in new model but not in aggregated
		% model
		isMemb = find(~ismember(bndNames,bndAllNames));
		for jj = 1:length(isMemb)
			Bnd_all.(bndNames(jj)) = Bnd.(bndNames(jj));
		end
		
		% Merge reaction list
		rxn_all = [rxn_all rxn];
		
		% Check for conflicting reaction rules
		if ~isequal(rxnRules_all,rxnRules)
			error(['Models to be combined are based on different rulesets. This is not allowed. Merged models must be based on the same ruleset. Please correct this before trying again.'])
		end
	end
end
modComp = modComp_all;
modSpc  = modSpc_all;
Bnd     = Bnd_all;
rxn     = rxn_all;
model.name = model.name(1:end-1);
rxnRules   = rxnRules_all;

%Load Model Rules
if ~exist('rateLaw','var')
    model.rxnRules = rxnRules;
else
    if ischar(rateLaw)
        rateLaw = str2func(rateLaw);
    end
    model.rxnRules = rateLaw;
end

%Check that model rule file exists
if ~exist(func2str(model.rxnRules),'file')
    error('parseModelm:rxnRuleFileNotFound','Specified rate law file not found. Please check the containing folder is in the MATLAB Path')
end
param = model.rxnRules('ini');

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
    [modComp(ii),pInd,pFit] = testPar(modComp(ii),pFit);
	
	if ~isnan(pInd(2))
        parDesc = [num2str(pInd(2)) ' | Comp : ' modComp(ii).name];
		if ~isempty(pFit.desc{pInd(2)})
			parDesc(1:length(num2str(pInd(2)))) = ' ';
			oldLength = size(pFit.desc{pInd(2)},2);
			newLength = length(parDesc);
			pad(1:max(oldLength,newLength)) = ' ';
		else
			newLength = 0;
			oldLength = 0;
			pad = ' ';
		end
        pFit.desc{pInd(2)} = [pFit.desc{pInd(2)} pad(ones(size(pFit.desc{pInd(2)},1),1),1:(newLength-oldLength));
			                      parDesc        pad(1:(oldLength-newLength))];
		if isnan(pFit.lim(pInd(2),1)) %Insert default parameter range if custom value not defined
			pFit.lim(pInd(2),:) = Bnd.Comp;
		end
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
model.spcMode = spcMode;
for ii = 1:size(modSpc)
	% Process parameter
    [modSpc(ii),pInd,pFit] = testPar(modSpc(ii),pFit);
	
    if ~isnan(pInd(3))
        parDesc = [num2str(pInd(3)) ' | Conc : ' modSpc(ii).name];
		if ~isempty(pFit.desc{pInd(3)})
			parDesc(1:length(num2str(pInd(3)))) = ' ';
			oldLength = size(pFit.desc{pInd(3)},2);
			newLength = length(parDesc);
			pad(1:max(oldLength,newLength)) = ' ';
		else
			newLength = 0;
			oldLength = 0;
			pad = ' ';
		end
        pFit.desc{pInd(3)} = [pFit.desc{pInd(3)} pad(ones(size(pFit.desc{pInd(3)},1),1),1:(newLength-oldLength));
			                      parDesc        pad(1:(oldLength-newLength))];
		if isnan(pFit.lim(pInd(3),1)) %Insert default parameter range if custom value not defined
			pFit.lim(pInd(3),:) = Bnd.Conc;
		end
    end
    
    modSpc(ii).pInd = pInd(3);
end

% Contract modComp
tmpSpc.name   = vertcat({modSpc.name});
[~,tmpSpc.comp] = ismember(vertcat({modSpc.comp}),modComp.name);
tmpSpc.comp = tmpSpc.comp';
tmpSpc.comp = [tmpSpc.comp tmpSpc.comp*0]; %expand complex vector by 1 
                                           %column. For complex species 
										   %which have no defined reference
										   %species yet.
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
    [rxn(ii),pInd,pFit] = testPar(rxn(ii),pFit);
	
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
		if isempty(pFit.desc{pInd(pList2Desc(jj))})
			pFit.desc{pInd(pList2Desc(jj))} = [num2str(pInd(pList2Desc(jj))) ' | ' parDesc{jj,3}]; %Insert definition
		else
			clear oldPad newPad
			dummyStr = num2str(pInd(pList2Desc(jj)));
			dummyStr(1:end) = ' ';
			newStr = [dummyStr ' | ' parDesc{jj,3}];
			wNew = length(newStr);
			[hOld,wOld] = size(pFit.desc{pInd(pList2Desc(jj))});
			oldPad(1:hOld,1:(max(0,wNew-wOld))) = ' '; 
			newPad(1,1:(max(0,wOld-wNew))) = ' '; 
			pFit.desc{pInd(pList2Desc(jj))} = [pFit.desc{pInd(pList2Desc(jj))} oldPad;newStr newPad];
		end
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
		try
			param(paramInd).matVal(matInd:(matInd+size(matVal{jj},1)-1),:) = matVal{jj};
		catch msg
			keyboard
		end
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

function [valOut,strOut] = compareParam(val1,val2)
	% Make val1 the longer parameter description
	if length(val2) > length(val1)
		val = val2;
		val2 = val1;
		val1 = val;
	end
	
	% Compare the parameter descriptions
	if  length(val1)==4 && (length(val1)==length(val2)) %Both have custom grouped boundaries
		valOut = val1;
		valOut(3:4) = [min(valOut(3),val2(3)) max(valOut(4),val2(4))];
		strOut = 'Unequal boundary found. Largest custom boundary chosen.';
	elseif length(val1)==3 && (length(val1)==length(val2)) %Both have custom boundaries
		valOut = val1;
		valOut(2:3) = [min(valOut(2),val2(2)) max(valOut(3),val2(3))];
		strOut = 'Unequal boundary found. Largest custom boundary chosen.';
	elseif length(val1)==2 && (length(val1)==length(val2)) %Both have custom boundaries
		valOut = [];
		strOut = 'Clash found between dependence description. Please check.';
	elseif length(val1)==1 && (length(val1)==length(val2)) %Both have custom boundaries
		valOut = NaN;
		strOut = 'Parameter set as free parameter.';
	elseif length(val1)>length(val2)
		valOut = val1;
		if length(val2) == 3
			valOut(3:4) = [min(valOut(3),val2(2)) max(valOut(4),val2(3))];
			strOut = 'Reference parameter set. Largest custom boundary chosen.';
		else
			strOut = 'Unequal parameter description found. More complex option chosen.';
		end
	end
end