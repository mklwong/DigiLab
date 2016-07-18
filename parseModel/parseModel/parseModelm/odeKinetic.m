function varargout = odeKinetic(method,varargin)

% Compartmentalisation has been implemented. For each reaction, the larger
% volume appears in the numerator. Thus, the reaction rate is with respect
% to the larger volume of the substrates (this is because flux is kept 
% constant, but with different volumes, this results in different rate of 
% change of concentration. So one of them has to be fixed. It turns out for 
% universal agreement of relationship between flux and reaction rate, we 
% just have to fix whether to use the larger volume or smaller volume as 
% standard when adjusting reaction rates for volume differences.

switch lower(deblank(method))
	
case 'ini'
	% Put your tensors into the parameter "param" in any for you want. Note
	% tensors will be populated down the row. The "name" field is required
	% and must contain a string.
	%
	% These are what will appear in the kinetic equation
	
	param(1).name = 'Km';
	param(1).matVal = nan(1,5);
	%[rate_indx,rate_R1_indx,R2_indx,geo_pindx,Km_pindx]
	
	param(2).name = 'k2';
	param(2).matVal = nan(1,5);
	%[prod_indx,R1_indx,R2_indx,geo_pindx,k2_pindx]
    
	param(3).name = 'k1';
	param(3).matVal = nan(1,4);
	%[prod_indx,R1_indx,geo_pindx,k1_pindx]
    
	param(4).name = 'k0';
	param(4).matVal = nan(1,2);
	%[prod_indx,k0_pindx]
    
	param(5).name = 'Hill';
	param(5).matVal = nan(1,7);
	%[prod_indx,sub_indx,enz_indx,geo_pindx,k1_pindx,Km_pindx,n_pindx]
	varargout = {param};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'rxnrules'
% Require inputs are: [rxn,x,flags]
% Required outputs are: [reqTens,tensVal,paramDesc,x*]	
%	* x is used when new species need to be created. It is passed in the
%	output to inform the rest of the program of the species change it
%	implements
%
% Ensure that when collapsing the tensors, that the tensor collapses into
% the first column. I.e. the index of the actual derivative in the tensor
% equation must match the index of the first row of the index based tensor.
%
% E.g. index    i        j         k         value
%     Param = [ 1        3         4          104  ];
%
% Rate is dx_dt = param*x*x. The indices implies dx_dt(1) = 104*x(3)*x(4)

[rxn,modSpc,flag,ii] = varargin{:};

% Convert all non-cell species references to one that has cell
if ~iscell(rxn.sub)
	if isempty(rxn.sub)
		rxn.sub = '';
	end
	rxn.sub = {rxn.sub};
end
if ~iscell(rxn.prod)
	rxn.prod = {rxn.prod};
end
if ~iscell(rxn.enz)
	rxn.enz = {rxn.enz};
end

% Match species in reaction with species master list
[~,testSub,subIndx]  = intersect(upper(rxn.sub) ,upper(modSpc.name));
[~,testProd,prodIndx] = intersect(upper(rxn.prod),upper(modSpc.name));
[~,testEnz,enzIndx]  = intersect(upper(rxn.enz) ,upper(modSpc.name));

% Check that there are no unmatched species in the reaction. If any
% unmatched then throw an error
if length(testSub) ~= length(rxn.sub)  && ~isempty(rxn.sub)
	rxn.sub(testSub) = [];
	rxn.sub(2,:) = {', '};
	rxn.sub(end,end) = {'"'};
	error('odeKinetic:SubstrateNotInSpeciesList',['Substrates "' horzcat(rxn.sub{:,:}) ' missing in reaction with description: "' rxn.desc '"'])
elseif length(testProd) ~= length(rxn.prod) 
	rxn.prod(testProd) = [];
	rxn.prod(2,:) = {', '};
	rxn.prod(end,end) = {'"'};
	error('odeKinetic:ProductNotInSpeciesList',['Product "' horzcat(rxn.prod{:,:}) ' missing in reaction with description: "' rxn.desc '"'])
elseif length(testEnz) ~= length(rxn.enz)
	rxn.enz(testEnz) = [];
	rxn.enz(2,:) = {', '};
	rxn.enz(end,end) = {'"'};
	error('odeKinetic:EnzymeNotInSpeciesList',['Enzyme "' horzcat(rxn.enz{:,:}) ' missing in reaction with description: "' rxn.desc '"'])
end

% Make 1x0 vector if matrix is empty
if isempty(subIndx)
	subIndx = zeros(0,1);
end
if isempty(prodIndx)
	prodIndx = zeros(0,1);
end
if isempty(enzIndx)
	enzIndx = zeros(0,1);
end

% Change vector direction for legacy matlab versions
if isrow(subIndx)
	subIndx = subIndx';
end
if isrow(prodIndx)
	prodIndx = prodIndx';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Classifier%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nSub = length(subIndx);
if nSub ==0
    if ~isempty(enzIndx) 
		if isempty(rxn.Km)
			subIndx  = [enzIndx;subIndx];
			prodIndx = [enzIndx;prodIndx];
			rxnType = 'uni';
		end
    else
        rxnType = 'syn';
    end
elseif nSub == 1
    if ~isempty(enzIndx) 
        if isempty(rxn.Km)
            subIndx   = [enzIndx;subIndx];
            prodIndx  = [enzIndx;prodIndx];
            rxnType = 'bi';
		elseif ~isempty(rxn.n)
			rxnType = 'hillFun';
        elseif ~isempty(rxn.Km)
            rxnType = 'enzQSSA';
        else
            error('parseModelm:exceptionClassifier1','Unknown error');
        end
    else
        rxnType = 'uni';
    end
elseif nSub == 2
	rxnType = 'bi';
else
    error('parseModelm:TooManySubstrates',['There are too many substrates in reaction ' num2str(ii)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rules for classified reaction types
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k  = rxn.k;
Km = rxn.Km;
n  = rxn.n;
r = rxn.r;
if isempty(r)
	r = 1;
end

subVec  = ones(size(subIndx));
prodVec = ones(size(prodIndx));

	switch rxnType
		case 'syn'
			% Header
			reqTens   = {'k0'};
			% Synthesis
			paramDesc = {'k','k0',['k    : phi -> ' modSpc.name{prodIndx}]};  

			% Maths
			matVal   = {[prodIndx k*prodVec]};
            pBuild   = {{'','k'}};
		case 'uni'
			% Header
			reqTens   = {'k1'};
			if isempty(prodIndx)
				% Degradation
                rxnFormula = [modSpc.name{subIndx} ' -> phi'];
                paramDesc = {'r','r',['VolOverlap: ' rxnFormula];
                             'k','k1',['k: ' rxnFormula]};
			elseif subIndx(1)==prodIndx(1)
				% Enzyme mediated Synthesis
				rxnFormula = [modSpc.name{prodIndx(2:end)} ' | ' modSpc.name{subIndx}];
                paramDesc = {'r','r',['VolOverlap: ' rxnFormula];
                             'k','k1',['k/Km: ' rxnFormula]};
			else
				% Unimolecular conversion/dissociation
				rxnFormula = [modSpc.name{subIndx} ' -> ' modSpc.name{prodIndx}];
                paramDesc = {'r','r',['VolOverlap: ' rxnFormula];
                             'k','k1',['k: ' rxnFormula]};
            end

                         
			% Maths
				% Everything else
				matVal   = {[subIndx         subIndx               subVec*[r -k];
							 prodIndx        subIndx'*prodVec     prodVec*[r  k]]};
                pBuild   = {{'','','r','k'}};
                            
		case 'bi'
			% Header
			reqTens   = {'k2'};
			if isempty(rxn.enz)
				% Association
                rxnFormula = [modSpc.name{subIndx} ' -> ' modSpc.name{prodIndx}];
				paramDesc = {'r','r',['volOverlap: ' rxnFormula];
                             'k','k2',['k: ' rxnFormula];};
			elseif (prodIndx(1)==subIndx(1)) && length(prodIndx)==1
				% Simplified enzyme mediated degradation
                rxnFormula = [modSpc.name{subIndx(2)} ' ->  phi | ' modSpc.name{subIndx(1)}];
				paramDesc = {'r','r',['volOverlap: ' rxnFormula];
                             'k','k2',['kc/Km: ' rxnFormula]};
            elseif (prodIndx(1)==subIndx(1))
				% Simplified enzyme kinetics
                rxnFormula = [modSpc.name{subIndx} ' -> ' modSpc.name{prodIndx}];
				paramDesc = {'r','r',['volOverlap: ' rxnFormula];
                             'k','k2',['kc/Km: ' rxnFormula]};
			end

			% Maths
			matVal   = {[[subIndx  subVec*subIndx'   subVec*[r -k]];
						  [prodIndx prodVec*subIndx' prodVec*[r  k]]]};
            pBuild   = {{'','','','r','k'}};
                      
		case 'enzQSSA'
			%Determine if overlap volume given. If not, use substrate
			%volume

			if flag(1)
				%Make new complex species
				compIndx = length(modSpc.name)+1;
				modSpc.name{compIndx} = [modSpc.name{subIndx} '-' modSpc.name{enzIndx}];
				modSpc.comp(compIndx,:) = modSpc.comp([subIndx enzIndx]); %Assume complex formed is in same compartment as substrate
				modSpc.matVal(compIndx) = 0;
				modSpc.pInd(compIndx) = NaN;

				%Insert tensor values
				reqTens   = {'k1','Km'};
                rxnFormula = [modSpc.name{subIndx} ' -> ' modSpc.name{prodIndx} ' [' modSpc.name{enzIndx} ']'];
				paramDesc = {'r','r',['VolOverlap: ' rxnFormula];
                             'k','k1',['kc: ' rxnFormula];
							 'Km','Km',['Km: ' rxnFormula]};
				matVal   = {[subIndx     compIndx          [r -k];
							  prodIndx prodVec*compIndx  prodVec*[r k]];
							 [compIndx  subIndx enzIndx r -Km;
							  compIndx  enzIndx subIndx r -Km;
							  subIndx   subIndx enzIndx r  Km;
							  subIndx   enzIndx subIndx r  Km;
							  enzIndx   subIndx enzIndx r  Km;
							  enzIndx   enzIndx subIndx r  Km;
							  ]};
                pBuild   = {{'','','r','k'};
                            {'','','','r','Km'}};
			else
				reqTens   = {'k2','Km'};
                rxnFormula = [modSpc.name{subIndx} ' -> ' modSpc.name{prodIndx} ' [' modSpc.name{enzIndx} ']'];
				paramDesc = {'r','r',['VolOverlap : ' rxnFormula];
                             'k','k2',['kc/Km : ' rxnFormula];
							 'Km','Km',['Km    : ' rxnFormula]};
						 
				matVal   = {[subIndx      subIndx enzIndx        subVec*[r -k];
							  prodIndx prodVec*[subIndx enzIndx] prodVec*[r  k]];
							  [subIndx subIndx enzIndx r Km;
							   subIndx enzIndx subIndx r Km;
							   enzIndx subIndx enzIndx r Km;
							   enzIndx enzIndx subIndx r Km]};
                pBuild   = {{'','','','r','k'};
                            {'','','','r','Km'}};
			end
		case 'hillFun'
			reqTens   = {'Hill'};
                rxnFormula = [modSpc.name{subIndx} ' -> ' modSpc.name{prodIndx} ' [' modSpc.name{enzIndx} ']'];
				paramDesc = {'r','r',['VolOverlap: ' rxnFormula]
                             'k','k1',['kc_Hill: ' rxnFormula];
							 'Km','Km',['Km_Hill: ' rxnFormula];
							 'n','n',['n_Hill: ' rxnFormula]};
				matVal   = {[subIndx  subIndx enzIndx  subVec*[r -k Km n];
							 prodIndx  subIndx enzIndx prodVec*[r  k Km n]]
                             };
                pBuild   = {{'','','','r','k','Km','n'};
                            {'','','','r','k','Km','n'}};
	end
	varargout = {reqTens,matVal,pBuild,paramDesc,modSpc};
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'compile'
% Require inputs are:   [model,tspan]
% Required outputs are: [model]	
if length(varargin) == 2
	[model,tspan] = varargin{:};
elseif length(varargin) == 1
	model = varargin{:};
	tspan = [0 1]; %no non-dimensionalise
end

model.matInd = []; %Initialise matrix index order (to be used in dynEqn component).
for ii = 1:length(model.param)
	switch lower(model.param(ii).name)
		case 'k0'
			model.param(ii).matVal = nonDim(model.param(ii).matVal,tspan,2);
            tmpVec = zeros(length(model.modSpc),1);
            tmpVec(model.param(ii).matVal(:,1)) = model.param(ii).matVal(:,2);
            model.param(ii).matVal = tmpVec;
		case 'k1'
			model.param(ii).matVal = nonDim(model.param(ii).matVal,tspan,4);	
            model.param(ii).matVal(:,3) = model.comp(model.param(ii).matVal(:,2)).*model.param(ii).matVal(:,4).*model.param(ii).matVal(:,3);
            model.param(ii).matVal(:,4) = [];
            model.param(ii).matVal = sparse(model.param(ii).matVal(:,1),model.param(ii).matVal(:,2),model.param(ii).matVal(:,3));
            if ~all(size(model.param(ii).matVal)==length(model.modSpc))
                model.param(ii).matVal(length(model.modSpc),length(model.modSpc))=0;
            end
            model.matInd(1) = ii;
		case 'k2'
			model.param(ii).matVal = nonDim(model.param(ii).matVal,tspan,5);
            compUsed = min(model.comp(model.param(ii).matVal(:,[2 3])),[],2);
            model.param(ii).matVal(:,4) = compUsed.*model.param(ii).matVal(:,5).*model.param(ii).matVal(:,4);
            model.param(ii).matVal(:,5) = [];
            model.matInd(2) = ii;
		case 'km'
			compUsed = min(model.comp(model.param(ii).matVal(:,[2 3])),[],2);
            model.param(ii).matVal(:,4) = model.comp(model.param(ii).matVal(:,1)).*model.param(ii).matVal(:,5)./(model.param(ii).matVal(:,4).*compUsed);
            model.param(ii).matVal(:,5) = [];
            model.matInd(3) = ii;
		case 'hill'
            % go from [ind1 ind2 ind3 r k Km n] to [ind1 ind2 ind3 r*k*V Km n]
			model.param(ii).matVal = nonDim(model.param(ii).matVal,tspan,5);
            compUsed = min(model.comp(model.param(ii).matVal(:,[2 3])),[],2);
            model.param(ii).matVal(:,4) = compUsed.*model.param(ii).matVal(:,4).*model.param(ii).matVal(:,5);
            model.param(ii).matVal(:,5) = [];
            model.matInd(4) = ii;
		end
end

varargout = {model};

case 'dyneqn'
%% Dynamic Equation construction
%Required inputs are [t,x,model]
%Required outputs are [dx_dt]

[t,x,model] = varargin{:};

% ~~Initialise matrices~~
G = zeros(length(model.modSpc));
V = G;
HillTerm = G;

% ~~Construction of k1 matrix~~
WInd = model.matInd(1);
W = model.param(WInd).matVal;

% ~~Construction of k2 matrix~~
VInd = model.matInd(2);
VTmp = sparse(model.param(VInd).matVal(:,1),model.param(VInd).matVal(:,2),x(model.param(VInd).matVal(:,3)).*model.param(VInd).matVal(:,4));
[a,b] = size(VTmp);
V(1:a,1:b) = VTmp;

% ~~Construction of G matrix~~
GInd = model.matInd(3);
GTmp = sparse(model.param(GInd).matVal(:,1),model.param(GInd).matVal(:,2),x(model.param(GInd).matVal(:,3))./model.param(GInd).matVal(:,4));
[a,b] = size(GTmp);
G(1:a,1:b) = GTmp;

% ~~Construction of HillFun matrix~~
HillInd = model.matInd(4);
HillRatio = (x(model.param(HillInd).matVal(:,3))./model.param(HillInd).matVal(:,5)).^model.param(HillInd).matVal(:,6);
HillVal = model.param(HillInd).matVal(:,4).*HillRatio./(HillRatio+1);
HillTmp = sparse(model.param(HillInd).matVal(:,1),model.param(HillInd).matVal(:,2),HillVal);
[a,b] = size(HillTmp);
HillTerm(1:a,1:b) = HillTmp;

% Solve
varargout{1} = (eye(length(x))+G)\((V*x+(W+HillTerm)*x+model.sigma(t).*model.comp)./model.comp);
end