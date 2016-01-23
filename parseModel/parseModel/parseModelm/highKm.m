%% Define postProc states
modComp = {'Cytosol',1;
			'PM',NaN;
			'Outside',1000};

modSpc ={	'IR'		,'PM'		,	NaN	; %1
			'pIR'		,'PM'		,	0	; %2
			'piIR'		,'PM'		,	0	; %3
			'AKT'		,'PM'		,	NaN	; %4
			'pAKT'		,'PM'		,	0	; %5
			'AKTSub'	,'Cytosol'	,	NaN	; %6
			'pAKTSub'	,'Cytosol'	,	0	; %7
			'TSC2'		,'Cytosol'	,	NaN	; %8
			'piTSC2'	,'Cytosol'	,	0	; %9
			'mTORC1'	,'Cytosol'	,	NaN	; %10
			'pimTORC1'	,'Cytosol'	,	0	; %11
			'Insulin'	,'Outside'	,	0	; %12
      };

 dataSpc = {'p473AKT',{'*pAKT'};...
          'Insulin',{'*Insulin'}};
	   


% Features of default parameters
Bnd.k0 = [1e01 1e04];
Bnd.k1  = [5e-5 0.5];
Bnd.k2   = [5e-5 5e-1];
Bnd.Km   = [1e-2 1e02];
Bnd.Conc = [1e-1 1e1];
Bnd.Comp = [1e-1 1e1];
Bnd.r = [0 1];

% IR
rxn(end+1).desc = 'IR -> pIR | Insulin';
    rxn(end).sub = {'IR'};  
    rxn(end).prod= 'pIR'; 
	rxn(end).enz = 'Insulin';
    rxn(end).k  = NaN;  
    rxn(end).Km  = NaN; 
rxn(end+1).desc = 'pIR -> IR';
    rxn(end).sub='pIR';  
    rxn(end).prod='IR';
    rxn(end).k = [NaN 1 3]; 
rxn(end+1).desc = 'IR -> piIR | mTORC1 (2)';
    rxn(end).sub = 'IR';  
    rxn(end).prod= 'piIR';  
    rxn(end).enz = 'mTORC1';
    rxn(end).k = [NaN 1 2 4]; 
	rxn(end).Km = [NaN 1];
rxn(end+1).desc = 'piIR -> IR';
    rxn(end).sub = 'piIR';  
    rxn(end).prod= 'IR';  
    rxn(end).k = [3 1]; 

%AKT
rxn(end+1).desc = 'AKT -> pAKT(1)';
    rxn(end).sub = 'AKT';
    rxn(end).prod= 'pAKT';
    rxn(end).enz = 'pIR';
    rxn(end).k = NaN; 
    rxn(end).Km = NaN;
rxn(end+1).desc = 'pAKT -> AKT';
    rxn(end).sub ='pAKT';  
    rxn(end).prod='AKT'; 
    rxn(end).k = NaN; 

%Akt Substrates
rxn(end+1).desc = 'AKTSub->pAKTSub | pAKT (4)';
    rxn(end).sub = 'AKTSub';  
    rxn(end).prod= 'pAKTSub';
    rxn(end).enz = 'pAKT';
    rxn(end).k = NaN; 
    rxn(end).Km = NaN; 
rxn(end+1).desc = 'pAS160 -> AS160';
    rxn(end).sub ='pAKTSub';  
    rxn(end).prod='AKTSub';
    rxn(end).k = NaN; 

%TSC2
rxn(end+1).desc = 'TSC2->piTSC2  | pAKT (7)';
    rxn(end).sub = 'TSC2'; 
    rxn(end).prod= 'piTSC2';
    rxn(end).enz = 'pAKT';
    rxn(end).k = NaN;  
	rxn(end).Km = NaN;
rxn(end+1).desc = 'piTSC2 -> TSC2';
    rxn(end).sub = 'piTSC2';
    rxn(end).prod= 'TSC2';
    rxn(end).k   = NaN; 

%mTORC1
rxn(end+1).desc = 'mTORC1->pimTORC1 | TSC2 (8)';
    rxn(end).sub  = 'mTORC1';
    rxn(end).prod = 'pimTORC1';
    rxn(end).enz  = 'TSC2';
    rxn(end).k = NaN;  
	rxn(end).Km = NaN;  
rxn(end+1).desc = 'pimTORC1 -> mTORC1';
    rxn(end).sub = 'pimTORC1'; 
    rxn(end).prod= 'mTORC1';
    rxn(end).k = NaN; 
