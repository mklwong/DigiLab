%% Enzyme Kinetic (dQSSA)

rxnRules = @odeKinetic;
spcMode = 'c';

%% Compartment definition
%
% spcComp = {'Compartment name', relative size};
%
modComp = {'Cyto1', NaN;
           'Cyto2', NaN;
		   'Cyto3', NaN;
		   'Cyto4', NaN};

%% Model species definition
%
% modSpc ={'State name', 'Compatment'  , conc/param};

modSpc = {'A'         ,'Cyto1' , 1;
		  'C'         ,'Cyto2' , 1;
          'B'         ,'Cyto3' , 0;
          'D'         ,'Cyto4' , 1};

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
Bnd.Comp = [0 1];
Bnd.r    = [0 1];
Bnd.n   = [1 4];

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

%% AKT Translocation Mechanics
rxn(end+1).desc = 'A -> B | C';
    rxn(end).sub = 'A';
    rxn(end).prod = 'B'; 
    rxn(end).enz = 'C';
    rxn(end).k   = NaN; 
    rxn(end).Km  = NaN;
rxn(end+1).desc = 'B -> A | D';
    rxn(end).sub = 'B';  
    rxn(end).prod = 'A'; 
    rxn(end).enz = 'D';
    rxn(end).k   = NaN; 
    rxn(end).Km  = NaN; 
rxn(end+1).desc = 'B -> A';
    rxn(end).sub = 'B';  
    rxn(end).prod = 'A'; 
    rxn(end).k   = NaN;