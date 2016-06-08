%% Synthesis

spcMode = 'a';

%% Compartment definition
%
% spcComp = {'Compartment name', relative size};
%
modComp = {'Cyto1', NaN;
           'Cyto2', NaN};

%% Model species definition
%
% modSpc ={'State name', 'Compatment'  , conc/param};

modSpc = {'I'         ,'Cyto2'  , 0
		  'A'         ,'Cyto1'  , 1
	      'pA'        ,'Cyto1'  , 0
	      'B'         ,'Cyto1'  , 1
		  'p1B'       ,'Cyto1'  , 0
		  'p2B'       ,'Cyto1'  , 0
		  'C'         ,'Cyto2'  , 1
		  'D'         ,'Cyto2'  , 1
		  'CD'        ,'Cyto2'  , 0
		  'pCD'       ,'Cyto2'  , 0
		  };

%% Relationship between simulation state and model state association

% dataSpc = {'Exp State Name',{'Sim State Name 1','Sim State Name 2'}};

 dataSpc = {};

%% Features of default parameters
% Bnd* = [lb ub]
% * is the parameter type
%
Bnd.k0  = [1e01 1e04];
Bnd.k1  = [5e-5 5e-1];
Bnd.k2   = [5e-5 5e-1];
Bnd.Km   = [1e-2 1e02];
Bnd.Conc = [1e-1 1e1];
Bnd.n   = [1 4];
Bnd.Comp = [0 1];
Bnd.r    = [0 1];

%% Reactions
% Reactions are stored in the variable rxn. Each reaction has
% between 3 to 6 fields depending on the reaction to be created. (c) marks
% fields that are compulsory. The fields are:
%   - desc (c): Identifier of the reaction. Used for descling outputs.
%   - k     (c): Reaction rate. Is a parameter.
%   - Sub      : List of substrates written as a cell array.
%   - Prod     : List of products written as a cell array.
%   - Enz      : Mediating enzyme
%   - Km       : Michaelis Constant for enzymatic reaction. Is a parameter.
%   

% All enzyme kinetic reactions are reversible (as per the model in paper)

rxn(end+1).desc = 'B -> p1B | A';
    rxn(end).sub = 'B';  
    rxn(end).prod= 'p1B';
	rxn(end).enz = 'A';
    rxn(end).k   = NaN; 
	rxn(end).Km   = NaN;
rxn(end+1).desc = 'p1B -> B | A';
    rxn(end).sub = 'p1B';  
    rxn(end).prod= 'B';
	rxn(end).enz = 'A';
    rxn(end).k   = NaN; 
	rxn(end).Km   = NaN;
	
rxn(end+1).desc = 'p1B -> B';
    rxn(end).sub = 'p1B';  
    rxn(end).prod= 'B';
    rxn(end).k   = NaN; 
	
rxn(end+1).desc = 'C + D -> CD';
    rxn(end).sub = {'C','D'};  
    rxn(end).prod= 'CD';
    rxn(end).k   = NaN;
	
rxn(end+1).desc = 'CD -> C + D';
    rxn(end).sub = 'CD';  
    rxn(end).prod= {'C','D'};
    rxn(end).k   = NaN; 
	
rxn(end+1).desc = 'CD -> pCD | I';
    rxn(end).sub = 'CD';  
    rxn(end).prod= 'pCD';
	rxn(end).enz = 'I';
    rxn(end).k   = NaN; 
	rxn(end).Km  = NaN;
rxn(end+1).desc = 'pCD -> CD | I';
    rxn(end).sub = 'pCD';  
    rxn(end).prod= 'CD';
	rxn(end).enz = 'I';
    rxn(end).k   = NaN; 
	rxn(end).Km  = NaN;
	
rxn(end+1).desc = 'pCD -> CD';
    rxn(end).sub = 'pCD';  
    rxn(end).prod= 'CD';
    rxn(end).k   = NaN; 
	
rxn(end+1).desc = 'A -> pA | pCD';
    rxn(end).sub = 'A';  
    rxn(end).prod= 'pA';
	rxn(end).enz = 'pCD';
    rxn(end).k   = NaN;
	rxn(end).Km  = NaN;
rxn(end+1).desc = 'pA -> A | pCD';
    rxn(end).sub = 'pA';  
    rxn(end).prod= 'A';
	rxn(end).enz = 'pCD';
    rxn(end).k   = NaN;
	rxn(end).Km  = NaN;
	
rxn(end+1).desc = 'pA -> A';
    rxn(end).sub = 'pA';  
    rxn(end).prod= 'A';
    rxn(end).k   = NaN; 
	
rxn(end+1).desc = 'B -> p2B | pCD';
    rxn(end).sub = 'B';  
    rxn(end).prod= 'p2B';
	rxn(end).enz = 'pCD';
    rxn(end).k   = NaN; 
	rxn(end).Km  = NaN;
rxn(end+1).desc = 'p2B -> B | pCD';
    rxn(end).sub = 'p2B';  
    rxn(end).prod= 'B';
	rxn(end).enz = 'pCD';
    rxn(end).k   = NaN; 
	rxn(end).Km  = NaN;
	
rxn(end+1).desc = 'p2B -> B | p1B';
    rxn(end).sub = 'p2B';  
    rxn(end).prod= 'B';
	rxn(end).enz = 'p1B';
    rxn(end).k   = NaN; 
	rxn(end).Km  = NaN;
rxn(end+1).desc = 'B -> p2B | p1B';
    rxn(end).sub = 'B';  
    rxn(end).prod= 'p2B';
	rxn(end).enz = 'p1B';
    rxn(end).k   = NaN; 
	rxn(end).Km  = NaN;
	
rxn(end+1).desc = 'p2B -> B';
    rxn(end).sub = 'p2B';  
    rxn(end).prod= 'B';
    rxn(end).k   = NaN; 
