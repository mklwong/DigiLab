%% Forward: Parameters
% A common feature in these model files is the idea of a "parameter". Any
% variable that is a parameter is defined using six potential formats:
%	1) [val]        : Known and fixed parameter
%	2) [NaN]        : Unknown parameter, default range used
%   3) [-grp]       : Unknown parameter that is grouped. Default range used.
%                     Multiplicative factor is assumed to be one.
%	4) [NaN lb ub]  : Unknown parameter, custom range
%   5) [-grp lb ub] : Unknown parameter with other parameters with the same
%                     value, custom range.
%   6) [-grp factor]: Unknown parameter that is a multiplicative factor of
%                     another parameter. Parameter of this format is the
%                     dependent.
% grp is a positive integer. All parameters with grp that is the same
% integer are in the same group. To distinguish it from "val", the first
% element of the vector is negative is it indicates an unknown parameter
% that is part of a group.

%% Compartment definition
%
% xComp = {'Compartment name', relative size};
%
xComp = {'Cyto', 1;
         'PM'  , 0.05;
         'ES'  , Inf;
        };

%% Model species definition
%
% xMod ={'State name', 'Compatment'  , conc/param};

xMod = {'AKT'     ,'Cyto', NaN;
        'mAKT'    ,'PM'  , 0;
        'p473mAKT','PM'  , 0;
        'p473AKT' ,'Cyto', 0;
        'mTORC2'  ,'PM'  , 1};

%% Relationship between simulation state and model state association

% xData = {'Exp State Name',{'Sim State Name 1','Sim State Name 2'}};

 xData = {'p473AKT',{'*p473AKT','*p473mAKT'}};

%% Features of default parameters
% Bnd* = [lb ub]
% * is the parameter type
%
Bndk0  = [1e01 1e04];
Bndk1  = [5e-5 0.5];
Bndk2   = [5e-5 5e-1];
BndKm   = [1e-2 1e02];
BndConc = [1e-1 1e1];

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

rxn(end+1).label = 'mAKT -> p473mAKT | mTORC2';
    rxn(end).sub = 'mAKT';  
    rxn(end).prod= 'p473mAKT'; 
    rxn(end).enz = 'mTORC2';
	rxn(end).Km  = NaN; 
    rxn(end).k   = 0.1; 
