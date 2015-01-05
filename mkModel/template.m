%% Define postProc states
% This template contains comments that guide the creation of a model. These
% are marked by lines with a simple '%' at the front. These can be deleted
% once you are familiar with how to create a model. Leave the lines with
% '%%'.

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
% ~~~Guide~~~
% "Compartment name" is the name of the compartment in the system
%
% Relative size is the relative size of the compartment. Can be infinite
% (in which case the concentration of species in that compartment is
% fixed).
%
% ~~~Example~~~
% We are modelling the a cell with only a cytoplasm (cyto) and plasma
% membrane (PM). The extracellular space (ES) is considered far larger than 
% either compartments which we consider infinite in size. This would give 
% the following.
%
xComp = {'Cyto', 1;
         'PM'  , 0.05;
         'ES'  , Inf;
        };

%% Model species definition
%
% xMod ={'State name', 'Compatment'  , conc/param};
%
% ~~~Guide~~~
% "State" is the name of the state to be included in the simulation
%
% Compartment name of the compartment as defined in the compartment part.
%
% Conc/param is either the concentration of the state. This is a parameter.
%
% ~~~Example~~~
% We are modelling the protein AKT. It can either be cytosolic or plasma
% membrane bound. It has two phosphorylation sites and can be dual
% phosphorylated. AKT initially exists unphosphorylated and in the cytosol.
% mTORC2 is also included to phosphorylation AKT.
% We would define this by:
%
<<<<<<< HEAD
xMod = {'AKT'     ,'Cyto',NaN;
        'mAKT'    ,'PM',  0;
        'p473mAKT','PM',  0;
        'p473AKT' ,'Cyto',0;};
=======
xMod = {'AKT'     ,'Cyto', NaN;
        'mAKT'    ,'PM'  , 0;
        'p473mAKT','PM'  , 0;
        'p473AKT' ,'Cyto', 0;
        'mTORC2'  ,'PM'  , 1};
>>>>>>> master

%% Relationship between simulation state and model state association

% xData = {'Exp State Name',{'Sim State Name 1','Sim State Name 2'}};

% ~~~Guide~~~
%
% Overview: The way this is designed is each experiment state is
% "constructed" from the simulated states. Thus one experimental states
% maps to multiple simulated states.
%
% Specifics:
% "Exp Name" is the name of the state as labelled in the experimental
% data.
%
% "Sim State Name 1" and "Sim State Name 2" are the state names as labelled
% in the simulation. The states included are summed together.
%
% An asterix (*) at the front of simulation state names implies the
% implicit complexes for that state (calculated by not created explicitly
% in xMod, such as ES complexes) are included as well.
%
%
% ~~~Example~~~
% Phosphorylation of p473AKT is experimental state to be compared. The
% experiment is a Western Blot without cell fractionation, hence plasma
% membrane and cytosolic AKT are captured. In our model the two pools of
% AKT are discrete, hence the group of simulation states corresponding to
% the experimental state are p473AKT and p473mAKT. Additionally, enzyme
% substrate complexes for both the simulation states are included, so they
% need to be asterixed.
%
% The final row for xData is hence:
 xData = {'p473AKT',{'*p473AKT','*p473mAKT'}};

%% Features of default parameters
%
% Bnd* = [lb ub]
% * is the parameter type
%
% ~~~Guide~~~
% The above gives the default boundaries for each type of parameter. There
% are currently five parameter types:
%   k0   - Zeroth order synthesis type
%   k1   - First order degradataion type
%   k2   - Second order associaition type
%   Km   - Association/Dissociation/Michaelis constant type
%   Conc - Concentration type
%
Bndk0  = [1e01 1e04];
Bndk1  = [5e-5 0.5];
Bndk2   = [5e-5 5e-1];
BndKm   = [1e-2 1e02];
BndConc = [1e-1 1e1];

%% Reactions
%
% ~~~Guide~~~
%
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
% Begin a new reaction by placing (end+1) after rxn. Continue a reaction by
% placing (end) after rxn [The reason for this is due to programming. Just
% follow this rule and you will not run into trouble].

% The combination used for each reaction determines the rules used in the
% backend of the program. This will be explained more later.
%
% For example, with an unknown Km and a known reaction rate of 0.1mol/s,
% the phosphorylation of AKT at the membrane by mTORC2 can be written as:

rxn(end+1).label = 'mAKT -> p473mAKT | mTORC2';
    rxn(end).sub = 'mAKT';  
    rxn(end).prod= 'p473mAKT'; 
    rxn(end).enz = 'mTORC2';
	rxn(end).Km  = NaN; 
    rxn(end).k   = 0.1; 

 % More information on types of reactions:
 % If the fields (excluding label and k) are:
 %      - 1 Sub: Degradataion
 %      - 1 Sub, 1 Enz: Enzyme mediated degradataion
 %      - 1 Sub, 1 Prod : Transformation
 %      - 1 Sub, 2+ Prod: dissociation
 %      - 2 Sub, 1 Prod : Association 
 %      - 1 Prod: Constant synthesis
 %      - 1 Prod, 1 Enz: Enzyme mediated synthesis
 %      - 1 Sub, 1 Enz, Prod: Enzymatic reaction that can produce any 
 %        number of products
