%% Synthesis

%% Compartment definition
%
% spcComp = {'Compartment name', relative size};
%
spcComp = {'Cyto', NaN;
           'Cyto2',NaN};

%% Model species definition
%
% modSpc ={'State name', 'Compatment'  , conc/param};

modSpc = {'A'         ,'Cyto'  , 1;
          'B'         ,'Cyto'  , 0;
          'C'         ,'Cyto'  , 1};

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
Bnd.k1MM  = [5e-5 5e-1];
Bnd.KmMM   = [5e-5 5e-1];
Bnd.Hill_n   = [1e-2 1e02];
Bnd.Conc = [1e-1 1e1];

%% Reactions
% Reactions are stored in the variable rxn. Each reaction has
% between 3 to 6 fields depending on the reaction to be created. (c) marks
% fields that are compulsory. The fields are:
%   - Label (c): Identifier of the reaction. Used for labelling outputs.
%   - k     (c): Reaction rate. Is a parameter.
%   - Sub      : List of substrates written as a cell array.
%   - Prod     : List of products written as a cell array.
%   - Enz      : Mediating enzyme
%   - Km       : Michaelis Constant for enzymatic reaction. Is a parameter.
%   

%% AKT Translocation Mechanics
rxn(end+1).label = 'A -> B | C (Hillfun)';
    rxn(end).sub = 'A';  
    rxn(end).prod= 'B';
	rxn(end).enz = 'C';
    rxn(end).k   = NaN; 
	rxn(end).Km  = NaN;
	rxn(end).n   = NaN;
rxn(end+1).label = 'B -> A';
    rxn(end).sub = 'B';  
    rxn(end).prod= 'A';
    rxn(end).k   = NaN; 