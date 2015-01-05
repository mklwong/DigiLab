function out = parseModelm(model,varargin)

%   out = parseModelm(model,debug)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Program Internal Description %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Each parameter has six potential formats:
%	1) [val]        : Known parameter
%	2) [NaN]        : Unknown parameter, default range used
%   3) [NaN grp]       : Unknown parameter that is grouped. Default range used.
%                     Multiplicative factor is assumed to be one.
%	4) [NaN lb ub]  : Unknown parameter, custom range
%   5) [NaN grp lb ub] : Unknown parameter with other parameters with the same
%                     value, custom range.
%   6) [NaN grp factor]: Unknown parameter that is a multiplicative factor of
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
% The base struct will contain a struct with 6 fields, all of which are
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

%%%%%%%%%%%%%%%%%%%
%% Function options
%%%%%%%%%%%%%%%%%%%
% Default options
expComp = 1;

Names = 'expComp ';

%Parse optional parameters (if any)
for ii = 1:length(varargin)
	if ischar(varargin{ii}) %only enter loop if varargin{ii} is a parameter
		switch lower(deblank(varargin{ii}))
			case lower(deblank(Names(1,:)))
				expComp = varargin{ii+1};
			case []
				error('Expecting Option String in input');
			otherwise
				error('Non-existent option selected. Check spelling.')
		end
	end
end
	 
%%%%%%%%%%%%%%%%%%%%%
%% Import model file
%%%%%%%%%%%%%%%%%%%%%
rxn = struct();
rxn(end).label = [];
    rxn(end).sub = [];  
    rxn(end).prod= []; 
    rxn(end).enz = [];
	rxn(end).Km = [];
    rxn(end).k  = [];  
    
v = rxn; %Legacy code. For backward compatibility.

if isa(model,'function_handle')
	model = func2str(model);
end

out.name = model;
xMod = [];
run(model); %loads xMod and v, xData and boundaries from model

rxn(1) = [];
% Legacy code
if size(v) > 1
    error('Please change all v in your model file into rxn. v no longer recognised as reaction network variable')
end

% Convenience functions
notFixed = @(val) (isnan(val) || val < 0);
	%determine if a parameter is defined as fixed or not. This is
	%essentially determining if it's passed as NaN or negative, both of
	%which indicate free variables.
tall = @(struct) min([find(isnan(struct.tens(:,1)),1,'first')-1 length(struct.tens(:,1))]); 
	%NaN is used to indicate empty slots in a tensor. This is so numbers 
	%filled can be easily determined by counting non-NaN elements.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialise Output Structs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pFit.desc = cell(100,1); pFit.lim  = nan(100,2);
curParamIndx = 0;
storeGrp = [-1,-1];

repInd  = []; %repeatedly used parameters (in multiple contexts)

G.tens  = nan(100,4); G.sign  = nan(100,1); G.pInd  = nan(100,2); G.factor   = nan(100,1);
k2.tens = nan(100,4); k2.sign = nan(100,1); k2.pInd = nan(100,2); k2.factor  = nan(100,1);
k1.tens = nan(100,3); k1.sign = nan(100,1); k1.pInd = nan(100,2); k1.factor  = nan(100,1);
k0.tens = nan(100,2); k0.sign = nan(100,1); k0.pInd = nan(100,2); k0.factor  = nan(100,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cycle over list of species
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[a,~]     = size(xMod);
x.tens    = nan(a,1);    % Initial concentration
x.name    = cell(a,1);   % Species name
x.comp    = ones(a,1);   % Compartment Index
x.pInd    = nan(a,1);    % Vector showing the parameter index a free state will use

for ii = 1:a
    % Save Name
	x.name{ii} = xMod{ii,1};
    
    % Get compartment size for each species
	[~,comptIndx] = intersect(upper(xComp(:,1)),upper(xMod{ii,2}));
	x.comp(ii) = xComp{comptIndx,2};

    % Process parameter
    [val,freeParam,grp] = testPar(xMod{ii,3});
    
	x.tens(ii) = val; % Save either parameter value of parameter multiplicative factor
    
    if freeParam
        if ~all(grp~=storeGrp)            %%existing group of free parameters
            x.pInd(ii) = storeGrp(grp==storeGrp,2);
        else                             %% ungrouped or new group of free parameter
            curParamIndx = curParamIndx + 1; % New parameter required

            %Make new group
            if grp ~= 0                  
                storeGrp(end+1,:) = [grp curParamIndx];
            end

            % Set parameter boundary
            if length(xMod{ii,3}) == 3       %% If there is custom set boundary
                pFit.lim(curParamIndx,:) = xMod{ii,3}(2:3);
            else                             %% Else use default boundary
                pFit.lim(curParamIndx,:) = BndConc;
            end
            pFit.desc{curParamIndx}	= ['Conc - ' xMod{ii,1}];
            x.pInd(ii) = curParamIndx; % store p index the variable x0 represents
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cycle over list of Reactions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii=1:length(rxn)
	
	%% Reaction classifier (based on number of substrates)
    [subIndx,prodIndx,rxnType] = classifyReaction(rxn(ii),x);
    
	%% Determine if parameters are free or not
    % Test 'k'
    testVal = rxn.k;
    [val(1),freeParam(1),factor(1)] = testPar(testVal);

    % Test 'Km'
    testVal = rxn.Km;
    [val(2),freeParam(2),factor(2)] = testPar(testVal);
	
	%% Set tensor targets and values based on reaction type
    parseRxn(rxnType,prodIndx,subIndx,val,x)
	
	% Insert values into tensors. Expand tensors are necessary
    for jj = 1:length(reqTens);
        if freeParam(jj)
            [pFit,repInd,pInd] = getIndmkLabel(pFit,rxn(ii).(rateParam{jj}),bnd{jj},repInd,paramDesc{jj});
            if (find(isnan(pFit.lim(:,1)),1,'first')+10)>size(pFit.lim,1)
                pFit.desc = [pFit.desc;cell(100,1)]; 
                pFit.lim  = [pFit.lim;nan(100,2)];
            end
        else
            pInd = [];
        end
		eval([reqTens{jj} '= appendTens(' reqTens{jj} ',tensVal{jj},sign{jj},tensInd{jj},pInd,factor(jj));'])
        if eval(['(tall(' reqTens{jj} ')+10)>size(' reqTens{jj} '.tens,1)'])
            eval([reqTens{jj} '.tens = [' reqTens{jj} '.tens;nan(100,size(' reqTens{jj} ',2)];']) 
            eval([reqTens{jj} '.sign = [' reqTens{jj} '.sign;nan(100,1)];']) 
            eval([reqTens{jj} '.pInd = [' reqTens{jj} '.pInd;nan(100,1)];']) 
        end
    end
end
G.pInd(:,2) = -G.pInd(:,2); % Make the G.pInd negative since they are all inserted as denominators.

%% Determining how Model states interacts with experimental states
rmXData = [];
for ii = 1:size(xData,1)
	for jj = 1:length(xData{ii,2})
		if strcmp(xData{ii,2}{jj}(1),'*')
			tmpSign = -1;
			xData{ii,2}{jj}(1) = [];
		else
			tmpSign = 1;
		end
		[~,~,indx]  = intersect(upper(xData{ii,2}{jj}),upper(x.name));
        if isempty(indx)
			storeError(model,[],[],[],['The state ' xData{ii,2}{jj} ' required as an experimental equivalent state not found. This experimental state will be excluded from curve fitting'])
            rmXData = ii;
			break
		end
		indx_rxn(jj) = tmpSign*indx;
	end
	xData{ii,2} = indx_rxn ;
end
xData(rmXData,:) = [];
pFit.sim2dat = xData;

%% Final processing
% Find clamped species and remove them from tensors so they won't change
% themselves due to internal events, but will induce change
%
% Also remove clamp rows by turning their tens values into 0 to make them
% have no influence. These rows can't be removed because the pInd to  
% tensInd reference pairs to link the parameter vector to tensor matrix
% have already been fixed at this point. It will be too difficult to
% recalculate. So it's easier to just leave them there. This will also
% retain the model connection in the model file for clarity.
clampSpc = find(isinf(abs(x.comp)));
if ~isempty(clampSpc)
	for ii = clampSpc'
		k0Indx = (k0.tens(:,1)==ii);
		k0.tens(k0Indx,end) = 0;
		k0.pInd(k0Indx,:) = NaN;
		k0.sign(k0Indx,:) = NaN;
		k0.factor(k0Indx,:) = NaN;
		
		k1Indx = (k1.tens(:,1)==ii);
		k1.tens(k1Indx,end) = 0;
		k1.pInd(k1Indx,:) = NaN;
		k1.sign(k1Indx,:) = NaN;
		k1.factor(k1Indx,:) = NaN;
		
		k2Indx = (k2.tens(:,1)==ii);
		k2.tens(k2Indx,end) = 0;
		k2.pInd(k2Indx,:) = NaN;
		k2.sign(k2Indx,:) = NaN;
		k2.factor(k2Indx,:) = NaN;

		GIndx = (G.tens(:,1)==ii);
		G.tens(GIndx,end) = 0;
		G.pInd(GIndx,:) = NaN;
		G.sign(GIndx,:) = NaN;
		G.factor(GIndx,:) = NaN;
	end
end
% remove NaN rows
k0.tens(isnan(k0.tens(:,1)),:) = [];
k0.pInd(isnan(k0.pInd(:,1)),:) = [];
k0.sign(isnan(k0.sign(:,1)),:) = [];
k0.factor(isnan(k0.factor(:,1)),:) = [];

k1.tens(isnan(k1.tens(:,1)),:) = [];
k1.pInd(isnan(k1.pInd(:,1)),:) = [];
k1.sign(isnan(k1.sign(:,1)),:) = [];
k1.factor(isnan(k1.factor(:,1)),:) = [];

k2.tens(isnan(k2.tens(:,1)),:) = [];
k2.pInd(isnan(k2.pInd(:,1)),:) = [];
k2.sign(isnan(k2.sign(:,1)),:) = [];
k2.factor(isnan(k2.factor(:,1)),:) = [];

G.tens(isnan(G.tens(:,1)),:) = [];
G.pInd(isnan(G.pInd(:,1)),:) = [];
G.sign(isnan(G.sign(:,1)),:) = [];
G.factor(isnan(G.factor(:,1)),:) = [];
	
limIndx = isnan(pFit.lim(:,1));
pFit.desc(limIndx,:) = [];
pFit.lim(limIndx,:) = [];

%% Compile output
out.x = x;
out.pFit = pFit;
out.k0 = k0;
out.k1 = k1;
out.k2 = k2;
out.G  = G;
end

%%%%%%%%%%%%%%%%%%%%%
function tens = appendTens(tens,vals,signs,tensInd,pInd,fact)
n = size(vals,2);
try
	tens.tens(tensInd,1:n) = vals;
catch
	keyboard
end
if ~isempty(pInd)
    tens.pInd(tensInd,:) = [tensInd' pInd*ones(length(tensInd),1)];
    tens.sign(tensInd) = signs;
    tens.factor(tensInd) = fact;
end
end