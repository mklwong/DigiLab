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
	param(1).tens = nan(1,4);
	
	param(2).name = 'k2';
	param(2).tens = nan(1,4);
	
	param(3).name = 'k1';
	param(3).tens = nan(1,3);
	
	param(4).name = 'k0';
	param(4).tens = nan(1,2);
	
	param(5).name = 'k1MM';
	param(5).tens = nan(1,4);
	
	param(6).name = 'KmMM';
	param(6).tens = nan(1,4);
	
	param(7).name = 'Hill_n';
	param(7).tens = nan(1,4);
	
	% All possible fields in a reaction. These are what will appear in
	% 'rxnrules'
	rxn = struct();
	rxn(end).label = [];
    rxn(end).sub = [];  
    rxn(end).prod= []; 
    rxn(end).enz = [];
	rxn(end).Km = [];
    rxn(end).k  = [];  
	rxn(end).A  = [];
	rxn(end).n  = [];
	
	varargout = {param,rxn};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'rxnrules'
% Require inputs are: [rxn,x,expComp]
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

overlap = rxn.A;

if isempty(overlap)
	overlap = 1;
end

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
			tensVal   = {[subIndx         subIndx              -k*overlap*subVec;
						  prodIndx        subIndx'*prodVec      k*overlap*prodVec]};
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
		
		% Maths
		try
		tensVal   = {[[subIndx  subVec*subIndx'  -k*overlap*subVec];
					  [prodIndx prodVec*subIndx'  k*overlap*prodVec]]};
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

			%Insert tensor values
			reqTens   = {'k1','Km'};
			paramDesc = {['kc   : ' x.name{subIndx} ' -> ' x.name{prodIndx} ' [' x.name{enzIndx} ']'];
						 ['Km   : ' x.name{subIndx} ' -> ' x.name{prodIndx} ' [' x.name{enzIndx} ']']};
			tensVal   = {[subIndx     compIndx      -k*overlap;
						  prodIndx prodVec*compIndx  k*overlap];
						 [subIndx   subIndx enzIndx  Km;
						  subIndx   enzIndx subIndx  Km;
						  enzIndx   subIndx enzIndx  Km;
						  enzIndx   enzIndx subIndx  Km;
						  compIndx  subIndx enzIndx -Km;
						  compIndx  enzIndx subIndx -Km]};
		else
			reqTens   = {'k2','Km'};
			paramDesc = {['kc/Km : ' x.name{subIndx} ' -> ' x.name{prodIndx} ' [' x.name{enzIndx} ']'];
						 ['Km    : ' x.name{subIndx} ' -> ' x.name{prodIndx} ' [' x.name{enzIndx} ']']};
			tensVal   = {[subIndx      subIndx enzIndx       -k*overlap*subVec;
						  prodIndx prodVec*[subIndx enzIndx]  k*overlap*prodVec];
						  [subIndx subIndx enzIndx  Km;
						   subIndx enzIndx subIndx  Km;
						   enzIndx subIndx enzIndx  Km;
						   enzIndx enzIndx subIndx  Km]};
		end
	case 'hillFun'
		reqTens   = {'k1MM','KmMM','Hill_n'};
			paramDesc = {['kc_Hill : ' x.name{subIndx} ' -> ' x.name{prodIndx} ' [' x.name{enzIndx} ']'];
						 ['Km_Hill : ' x.name{subIndx} ' -> ' x.name{prodIndx} ' [' x.name{enzIndx} ']'];
						 ['n_Hill : ' x.name{subIndx} ' -> ' x.name{prodIndx} ' [' x.name{enzIndx} ']']};
			tensVal   = {[subIndx  subIndx enzIndx -k*overlap*subVec;
						 prodIndx  subIndx enzIndx  k*overlap*prodVec];
						 [subIndx  subIndx enzIndx Km;
						 prodIndx  subIndx enzIndx Km];
						  [subIndx subIndx enzIndx  n;
						 prodIndx  subIndx enzIndx  n];};
end
varargout = {reqTens,tensVal,paramDesc,x};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'insparam'
% Require inputs are:   [model,p]
% Required outputs are: [model]	
[model,p] = varargin{:};

% Putting param into tensors
x0 = model.conc.tens;
for ii = 1:length(model.param)
	freeInd = find(~isnan(model.param(ii).pInd));
	if ~isempty(freeInd)
		try
		model.param(ii).tens(freeInd,end) = model.param(ii).tens(freeInd,end).*p(model.param(ii).pInd(freeInd));
		catch
			keyboard
		end
	end
	% Pre-generate full matrix if possible
	if size(model.param(ii).tens,2)==2
		model.tensor.k0 = full(sparse(model.param(ii).tens(:,1),ones(size(model.param(ii).tens(:,1))),model.param(ii).tens(:,2)));
		if length(model.tensor.k0)~= length(x0)
			model.tensor.k0(size(x0,1),1) = 0;
		end
	elseif size(model.param(ii).tens,2)==3
		model.tensor.k1 = full(sparse(model.param(ii).tens(:,1),model.param(ii).tens(:,2),model.param(ii).tens(:,3)));
		if (size(model.tensor.k1,1)~= length(x0) || size(model.tensor.k1,2)~= length(x0))
			model.tensor.k1(size(x0,1),size(x0,1)) = 0;
		end
	end
end

% Putting concentration into tensors
freeInd = find(~isnan(model.conc.pInd));
if ~isempty(freeInd)
	model.conc.tens(freeInd,end) = model.conc.tens(freeInd,end).*p(model.conc.pInd(freeInd));
end

varargout = {model};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 'nondim'
%Required inputs are [model,tspan,normInp]
%Required outputs are [model]

[model,tspan,normInp] = varargin{:};

model.tensor.k0 = model.tensor.k0*(tspan(end)-tspan(1));
model.tensor.k1 = model.tensor.k1*(tspan(end)-tspan(1));
model.param(2).tens(:,end) = model.param(2).tens(:,end)*(tspan(end)-tspan(1));
model.param(5).tens(:,end) = model.param(5).tens(:,end)*(tspan(end)-tspan(1));

model.basalSigma = @(t) model.tensor.k0*ones(1,length(t)); 
model.fullSigma = @(t) normInp(t) + model.tensor.k0*ones(1,length(t)); 

varargout{1} = model;

case 'ramp'
% This part of the script modifies the reaction rates such that no reaction
% occurs during the ramping phase.
model = varargin{1};
	model.tensor.k0 = model.tensor.k0*0;
	model.tensor.k1 = model.tensor.k1*0;
	model.param(2).tens(:,end) = 0;
	model.param(5).tens(:,end) = 0;
varargout{1} = model;

case 'dyneqn'
%Required inputs are [t,x,model]
%Required outputs are [dx_dt]

[t,x,model] = varargin{:};
x(x<0) = 0; %sometimes the system goes to less than zero. When it wants to do this, set it to zero

M = zeros(size(model.tensor.k1));
L = M;
hillTerm = M;
compVal = model.conc.comp;

% Make compartment size never infinite. If it is, make it 1 million times
% larger.
indx = isinf(compVal);
compVal(indx) = 0;
compVal(indx) = max(compVal)*1e6;

% Compartment size correction for unimolecular type reactions.
sourceCompk1 = compVal';
sourceCompk1 = sourceCompk1(ones(1,length(compVal)),:);

% Compartment size correction for bimolecular reaction. Take smaller of two
% compartments
sourceCompk2 = compVal';
sourceCompk2 = sourceCompk2(ones(1,length(compVal)),:,[1 1]);
sourceCompk2(:,:,2) = sourceCompk2(:,:,1)';
sourceCompk2 = min(sourceCompk2,[],3);

%% Tensor construction
model.param(5).tens(:,4) = min(compVal(model.param(5).tens(:,2:3)),[],2).*model.param(5).tens(:,4).*(x(model.param(7).tens(:,3)).^model.param(7).tens(:,4))./(model.param(6).tens(:,4)+x(model.param(7).tens(:,3)).^model.param(7).tens(:,4));

MTmp = sparse(model.param(1).tens(:,1),model.param(1).tens(:,3),x(model.param(1).tens(:,2))./model.param(1).tens(:,4));
[a,b] = size(MTmp);
M(1:a,1:b) = MTmp;

LTmp = sparse(model.param(2).tens(:,1),model.param(2).tens(:,2),model.param(2).tens(:,4).*x(model.param(2).tens(:,3)).*min(compVal(model.param(2).tens(:,2:3)),[],2));
[a,b] = size(LTmp);
L(1:a,1:b) = LTmp;

hillTmp = sparse(model.param(5).tens(:,1),model.param(5).tens(:,2),model.param(5).tens(:,4));
[a,b] = size(hillTmp);
hillTerm(1:a,1:b) = hillTmp;

%% Solve
varargout{1} = ((eye(length(x))+M)\(L*x+(model.tensor.k1.*sourceCompk1)*x+model.k0(t).*compVal+hillTerm*x))./compVal;

end