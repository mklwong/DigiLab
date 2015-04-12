function [reqTens,tensInd,tensVal,pow,bnd,paramDesc,x] = parseRxn(rxnType,subIndx,prodIndx,enzIndx,val,x,Bnd)

% parseRxnImp
%   This reaction parser uses the dQSSA but calculates the complex
%   concentration explicitly.

prodIndx = prodIndx';
prodVec = ones(size(prodIndx));

%% Create comma separated Substrate and Product List
subList = cell(1,2*length(subIndx));
subList(1:2*length(subIndx)) = {', '};
subList(1:2:2*length(subIndx)) = x.name(subIndx);
if ~isempty(subList)
    subList(end) = [];
end

prodList = cell(1,2*length(prodIndx));
prodList(1:2*length(prodIndx)) = {', '};
prodList(1:2:2*length(prodIndx)) = x.name(prodIndx);
if ~isempty(prodList)
    prodList(end) = [];
end

%% Types
	switch rxnType
		case 'syn'
			tensInd   = {1:length(prodIndx)};
			tensVal   = {[prodIndx val(1)*prodVec]};
			pow      = {prodVec};
			bnd       = {Bnd.k0};
			reqTens   = {'k0'};
			paramDesc = {['k    : phi -> ' x.name{prodIndx}]};
		case 'uni' %or degradation
			tensInd   = {(1:length(prodIndx)+1)};
			tensVal   = {[subIndx         subIndx        -val(1)*ones(size(subIndx));
				          prodIndx     subIndx'*prodVec    val(1)*prodVec]};
			pow      = {[1;prodVec]};
			bnd       = {Bnd.k1};
			reqTens   = {'k1'};
			if isempty(prodIndx)
				paramDesc = {['k    : ' x.name{subIndx} ' -> phi'];};
			elseif subIndx(1)==prodIndx(1)
				paramDesc = {['k    : phi -> ' x.name{prodIndx(2:end)} ' | ' x.name{subIndx}];};
            else
                paramDesc = {['k    : ' x.name{subIndx} ' -> ' prodList{:}];};
			end
		case 'bi'  %or mass action enzyme kinetics
			tensInd   = {(1:length(prodIndx)+2)};
			tensVal   = {[[subIndx(1) subIndx(1) subIndx(2)     -val(1);
				           subIndx(2) subIndx(1) subIndx(2)     -val(1)];
						  [prodIndx     prodVec*subIndx     prodVec*val(1)]]};
			pow      = {[1;1;prodVec]};
			bnd       = {Bnd.k2};
			reqTens   = {'k2'};
			if length(prodIndx)==1
                if subIndx(2)==prodIndx(1)
                    paramDesc = {['kc/Km: ' x.name{subIndx(2)} ' ->  phi | ' x.name{subIndx(1)}]};
                else
                    paramDesc = {['k    : ' subList{:} ' -> ' prodList{:}]};
                end
			elseif length(prodIndx)==2
                if subIndx(2)==prodIndx(2)
                    paramDesc = {['kc/Km: ' x.name{subIndx(2)} ' -> ' x.name{prodIndx(2)} '| ' x.name{subIndx(1)}]};
                else
                    paramDesc = {['k    : ' subList{:} ' -> ' prodList{:}]};
                end
            else
                paramDesc = {['k    : ' subList{:} ' -> ' prodList{:}]};
			end
		case 'enzQSSA'
				%Make new complex species
				comIndx = length(x.name)+1;
                x.name{comIndx} = [x.name{subIndx} '-' x.name{enzIndx}];
                x.comp(comIndx) = -x.comp(subIndx);
                x.tens(comIndx) = 0;
				
				tensInd   = {(1:length(prodIndx)+1);
							 (1:6)};
				tensVal   = {[subIndx  comIndx -val(1);
					          prodIndx' ones(size(prodIndx',1),1)*[comIndx  val(1)]];
							 [subIndx subIndx enzIndx  val(2);
							  subIndx enzIndx subIndx  val(2);
							  enzIndx subIndx enzIndx  val(2);
							  enzIndx enzIndx subIndx  val(2);
							  comIndx subIndx enzIndx -val(2);
							  comIndx enzIndx subIndx -val(2)]};
				pow      = {[1 1],[-1 -1 -1 -1 -1 -1]};
				bnd       = {Bnd.k1,Bnd.Km};
				reqTens   = {'k1','G'};
				paramDesc = {['kc   : ' subList{:} ' -> ' prodList{:} ' [' x.name{enzIndx} ']'];
					         ['Km   : ' subList{:} ' -> ' prodList{:} ' [' x.name{enzIndx} ']']};
	end