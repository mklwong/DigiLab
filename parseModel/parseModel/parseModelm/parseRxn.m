function parseRxn(rxnType,prodIndx,subIndx,val,x)

prodIndx = prodIndx';
    prodVec = ones(size(prodIndx));
	switch rxnType
		case 'syn'
			tensInd   = {1:length(prodIndx)};
			tensVal   = {[prodIndx val(1)*proVec]};
			sign      = {proVec};
			bnd       = {Bndk0};
			reqTens   = {'k0'};
			paramDesc = {['k    : phi -> ' x.name{prodIndx}]};
		case 'uni' %or degradation
			tensInd   = {tall(k1)+(1:length(prodIndx)+1)};
			tensVal   = {[subIndx         subIndx              -val(1);
				          prodIndx subIndx'*proVec val(1)*size(prodIndx')]};
			sign      = {[-1;0*prodIndx'+1]};
			bnd       = {Bndk1};
			reqTens   = {'k1'};
			if (length(subIndx)==length(prodIndx))==1
				paramDesc = {['k    : ' x.name{subIndx} ' -> ' x.name{prodIndx}];};
			elseif isempty(prodIndx)
				paramDesc = {['k    : ' x.name{subIndx} ' -> phi'];};
			elseif subIndx(1)==prodIndx(1)
				paramDesc = {['k    : phi -> ' x.name{prodIndx(2:end)} ' | ' x.name{subIndx}];};
			end
		case 'bi'  %or mass action enzyme kinetics
			tensInd   = {(1:length(prodIndx)+2)};
			tensVal   = {[[subIndx(1) subIndx(1) subIndx(2) -val(1);
				           subIndx(2) subIndx(1) subIndx(2) -val(1)];
						  [prodIndx prodVec*subIndx(1)  prodVec*val(1)];
                          [prodIndx prodVec*subIndx(2)  prodVec*val(1)]]};
			sign      = {[-1;-1;prodVec]};
			bnd       = {Bndk2};
			reqTens   = {'k2'};
			if isempty(rxn(ii).enz)
				paramDesc = {['k    : ' x.name{subIndx} ' -> ' x.name{prodIndx}];};
			elseif length(prodIndx)==2
				paramDesc = {['kc/Km: ' x.name{subIndx(2)} ' -> ' x.name{prodIndx(2)} '| ' x.name{subIndx(1)}]};
			elseif length(prodIndx)==1
				paramDesc = {['kc/Km: ' x.name{subIndx(2)} ' ->  phi | ' x.name{subIndx(1)}]};
			end
		case 'enzQSSA'
			if expComp
				%Make new complex species
				comIndx = tall(x)+1;
				x.name{comIndx} = [x.name{subIndx} '-' x.name{enzIndx}];
				x.comp(comIndx) = -x.comp(subIndx); %Assume complex formed is in same compartment as substrate
				x.tens(comIndx) = 0;
				x.pInd(comIndx)    = NaN;
				x.selfLoc(comIndx) = NaN;
				
				tensInd   = {tall(k1)+(1:length(prodIndx)+1);
							 tall(G)+(1:6)};
				tensVal   = {[subIndx  comIndx -val(1);
					          prodIndx' ones(size(prodIndx',1),1)*[comIndx  val(1)]];
							 [subIndx subIndx enzIndx  1/val(2);
							  subIndx enzIndx subIndx  1/val(2);
							  enzIndx subIndx enzIndx  1/val(2);
							  enzIndx enzIndx subIndx  1/val(2);
							  comIndx subIndx enzIndx -1/val(2);
							  comIndx enzIndx subIndx -1/val(2)]};
				sign      = {[-1 1],[1 1 1 1 -1 -1]};
				bnd       = {Bndk1,BndKm};
				reqTens   = {'k1','G'};
				paramDesc = {['kc   : ' x.name{subIndx} ' -> ' x.name{prodIndx} ' [' x.name{enzIndx} ']'];
					         ['Km   : ' x.name{subIndx} ' -> ' x.name{prodIndx} ' [' x.name{enzIndx} ']']};
			else
				tensInd   = {tall(k1)+(1:length(prodIndx)+1);
							 tall(G)+1:4};
				tensVal   = {[subIndx  subIndz enzIndx -val(1);
					          prodIndx' ones(size(prodIndx',1),1)*[subIndz enzIndx  val(1)]];
							  [subIndx subIndx enzIndx  1/val(2);
							   subIndx enzIndx subIndx  1/val(2);
							   enzIndx subIndx enzIndx  1/val(2);
							   enzIndx enzIndx subIndx  1/val(2)]};
				sign      = {[-1 1];
					         [1 1 1 1]};
				bnd       = {Bndk2;pBndG};
				reqTens   = {'k2','G'};
				paramDesc = {['kc/Km : ' x.name{subIndx} ' -> ' x.name{prodIndx} ' [' x.name{enzIndx} ']'];
					         ['Km    : ' x.name{subIndx} ' -> ' x.name{prodIndx} ' [' x.name{enzIndx} ']']};
			end
	end