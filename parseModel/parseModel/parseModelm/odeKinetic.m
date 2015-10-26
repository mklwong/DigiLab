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
	
	param(1).name = 'Km';
	param(1).tens = nan(1,4);
	
	param(2).name = 'k2';
	param(2).tens = nan(1,4);
	
	param(3).name = 'k1';
	param(3).tens = nan(1,3);
	
	param(4).name = 'k0';
	param(4).tens = nan(1,2);
	
	varargout = {param};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'rxnrules'
% Require inputs are: [rxn,x,expComp]
% Required outputs are: [reqTens,tensVal,paramDesc,x]	
%
% Ensure that when collapsing the tensors, that the tensor collapses into
% the first column. I.e. the index of the actual derivative in the tensor
% equation must match the index of the first row of the index based tensor.
%
% E.g. index    i        j         k         value
%     Param = [ 1        3         4          104  ];
%
% Derivative is dx_dt = param*x*x;
%
% We require dx_dt(1) = 104*x(3)*x(4)

[rxn,x,expComp,ii] = varargin{:};

[~,~,subIndx]  = intersect(upper(rxn.sub) ,upper(x.name));
[~,~,prodIndx] = intersect(upper(rxn.prod),upper(x.name));
[~,~,enzIndx]  = intersect(upper(rxn.enz) ,upper(x.name));

%Change for legacy matlab versions
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
        subIndx  = [enzIndx;subIndx];
        prodIndx = [enzIndx;prodIndx];
        rxnType = 'uni';
    else
        rxnType = 'syn';
    end
elseif nSub == 1
    if ~isempty(enzIndx) 
        if isempty(rxn.Km)
            subIndx   = [enzIndx;subIndx];
            prodIndx  = [enzIndx;prodIndx];
            rxnType = 'bi';
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

overlap = rxn.A;

if isempty(overlap)
	overlap = 1;
end

prodComp = x.comp(prodIndx);
subComp  = x.comp(subIndx);
subVec  = ones(size(subIndx));
prodVec = ones(size(prodIndx));

switch rxnType
	case 'syn'
		% Header
		reqTens   = {'k0'};
		% Synthesis
		paramDesc = {['k    : phi -> ' x.name{prodIndx}]};  
		
		% Maths
		tensVal   = {[prodIndx k*prodVec]};
	case 'uni'
		% Header
		reqTens   = {'k1'};
		if isempty(prodIndx)
			% Degradation
			paramDesc = {['k    : ' x.name{subIndx} ' -> phi'];};
		elseif subIndx(1)==prodIndx(1)
			% Enzyme mediated Synthesis
			paramDesc = {['k    : phi -> ' x.name{prodIndx(2:end)} ' | ' x.name{subIndx}];};
		else
			% Unimolecular conversion/dissociation
			paramDesc = {['k    : ' x.name{subIndx} ' -> ' x.name{prodIndx}];};
		end
		
		% Maths
			% Everything else
			tensVal   = {[subIndx         subIndx              -k*overlap*subComp/subComp;
						  prodIndx        subIndx'*prodVec      k*overlap*subComp*(1./prodComp)]};
	case 'bi'
		% Header
		reqTens   = {'k2'};
		if isempty(rxn.enz)
			% Association
			paramDesc = {['k    : ' x.name{subIndx} ' -> ' x.name{prodIndx}];};
		elseif (prodIndx(1)==subIndx(1)) && length(prodIndx)==1
			% Simplified enzyme mediated degradation
			paramDesc = {['kc/Km: ' x.name{subIndx(2)} ' ->  phi | ' x.name{subIndx(1)}]};
		elseif (prodIndx(1)==subIndx(1))
			% Simplified enzyme kinetics
			paramDesc = {['kc/Km: ' x.name{subIndx(2)} ' -> ' x.name{prodIndx(2)} '| ' x.name{subIndx(1)}]};
		end
		
		if isempty(rxn.A)
			overlap = min(x.comp(subIndx));
		end
		
		% Maths
		try
		tensVal   = {[[subIndx  subVec*subIndx'  -k*overlap./x.comp(subIndx)];
					  [prodIndx prodVec*subIndx'  k*overlap*(1./x.comp(prodIndx))]]};
		catch
			keyboard
		end
	case 'enzQSSA'
		if expComp
			%Make new complex species
			compIndx = length(x.name)+1;
			x.name{compIndx} = [x.name{subIndx} '-' x.name{enzIndx}];
			x.comp(compIndx) = x.comp(subIndx); %Assume complex formed is in same compartment as substrate
			x.tens(compIndx) = 0;
			x.pInd(compIndx) = NaN;
			try
			tensVal   = {[subIndx     compIndx      -k*overlap*(1./x.comp(subIndx));
						  prodIndx prodVec*compIndx  k*overlap*(1./x.comp(prodIndx))];
						 [subIndx   subIndx enzIndx  Km;
						  subIndx   enzIndx subIndx  Km;
						  enzIndx   subIndx enzIndx  Km;
						  enzIndx   enzIndx subIndx  Km;
						  compIndx  subIndx enzIndx -Km;
						  compIndx  enzIndx subIndx -Km]};
			catch msg
				printErr(msg)
				keyboard
			end
			reqTens   = {'k1','Km'};
			paramDesc = {['kc   : ' x.name{subIndx} ' -> ' x.name{prodIndx} ' [' x.name{enzIndx} ']'];
						 ['Km   : ' x.name{subIndx} ' -> ' x.name{prodIndx} ' [' x.name{enzIndx} ']']};
		else
			tensVal   = {[subIndx      subIndx enzIndx       -k*overlap*(1./x.comp(subIndx));
						  prodIndx prodVec*[subIndx enzIndx]  k*overlap*(1./x.comp(prodIndx))];
						  [subIndx subIndx enzIndx  Km;
						   subIndx enzIndx subIndx  Km;
						   enzIndx subIndx enzIndx  Km;
						   enzIndx enzIndx subIndx  Km]};
			reqTens   = {'k2','Km'};
			paramDesc = {['kc/Km : ' x.name{subIndx} ' -> ' x.name{prodIndx} ' [' x.name{enzIndx} ']'];
						 ['Km    : ' x.name{subIndx} ' -> ' x.name{prodIndx} ' [' x.name{enzIndx} ']']};
		end
end
varargout = {reqTens,tensVal,paramDesc,x};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'insparam'
% Require inputs are:   [model,p]
% Required outputs are: [model]	
[model,p] = varargin{:};

% Putting param into tensors
for ii = 1:length(model.param)
	freeInd = find(~isnan(model.param(ii).pInd));
	if ~isempty(freeInd)
		model.param(ii).tens(freeInd,end) = model.param(ii).tens(freeInd,end).*p(model.param(ii).pInd(freeInd));
	end
end

% Putting concentration into tensors
freeInd = find(~isnan(model.conc.pInd));
if ~isempty(freeInd)
	model.conc.tens(freeInd,end) = model.conc.tens(freeInd,end).*p(model.conc.pInd(freeInd));
end

varargout = {model};
end