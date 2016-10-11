function model = parseModelm(model)

%   out = parseModelm(model)
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
	
if ischar(model)
	rxn = sigRxnList();
	rxn(1) = []; %When rxn is initialised it always creates an empty reaction. Remove this.
	run(model(1:(end-2))); 
	modelName = model;
	model = struct();
	model.name = modelName;
elseif isstruct(model)
	grpsChecked = []; % Parameter checking
	rxn = model.raw.rxn;       % Putting in reactions
	Bnd = model.raw.bnd;       % Putting in default boundaries
	rxnRules = model.rxnRules; % Putting used reaction rules

	% De-parsing model compartments
	modComp = model.modComp.name;
	modComp = [modComp num2cell(model.modComp.matVal)];
	% Apply the parameter definition and boundaries for free species
	for jj = find(~isnan(model.modComp.pInd))'
		if any(model.pFit.grp(:,2)==model.modComp.pInd(jj))
			% Check for grouped parameters
			curGrp = model.pFit.grp(:,model.pFit.grp(:,2)==model.modComp.pInd(jj),1);
			if any(curGrp==grpsChecked)
				%If group encountered has already been found
				modComp(jj,2) = {[curGrp model.modComp.matVal(jj)]};
			else
				%If group encountered has not yet been found
				modComp(jj,2) = {[NaN curGrp model.pFit.lim(model.modComp.pInd(jj),:)]};
				grpsChecked = [grpsChecked curGrp];
			end
		else
			% For ungrouped parameters (we will ignore case where the
			% default boundary is used).
			modComp(jj,2) = {[NaN model.pFit.lim(model.modComp.pInd(jj),:)]};
		end
	end

	% De-parsing model species
	modSpc = model.modSpc.name';
	modSpc = [modSpc model.modComp.name(model.modSpc.comp(:,1))];
	modSpc = [modSpc num2cell(model.modSpc.matVal)];
	% Apply the parameter definition and boundaries for free species
	for jj = find(~isnan(model.modSpc.pInd))'
		if any(model.pFit.grp(:,2)==model.modSpc.pInd(jj))
			% Check for grouped parameters
			curGrp = model.pFit.grp(:,model.pFit.grp(:,2)==model.modSpc.pInd(jj),1);
			if any(curGrp==grpsChecked)
				%If group encountered has already been found
				modSpc(jj,3) = {[curGrp model.modSpc.matVal(jj)]};
			else
				%If group encountered has not yet been found
				modSpc(jj,3) = {[NaN curGrp model.pFit.lim(model.modSpc.pInd(jj),:)]};
				grpsChecked = [grpsChecked curGrp];
			end
		else
			% For ungrouped parameters (we will ignore case where the
			% default boundary is used).
			modSpc(jj,3) = {[NaN model.pFit.lim(model.modSpc.pInd(jj),:)]};
		end
	end
	modSpc(model.modSpc.comp(:,2)~=0,:)=[];   % Remove custom created species
	spcMode = model.spcMode;
	modelname = model.name;	
	model = struct();
	model.name = modelname;
end

model.rxnRules = rxnRules;
model.modComp = modComp;
model.modSpc = modSpc;
model.spcMode = spcMode;
model.rxn = rxn;
model.bnd = Bnd;
