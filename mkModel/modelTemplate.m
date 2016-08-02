%% Define postProc states
% This template contains comments that guide the creation of a model. These
% are marked by lines with a simple '%' at the front. These can be deleted
% once you are familiar with how to create a model. Leave the lines with
% '%%'.

%% Forward: Parameters
% A common feature in these model files is the idea of a "parameter". Any
% variable that is a parameter is defined using six potential formats:
%	1) [val]           : Known and fixed parameter
%	2) [NaN]           : Unknown parameter, default range used
%	3) [NaN lb ub]     : Unknown parameter, custom range
%   5) [grp*i]         : Unknown parameter with other parameters with the same
%                        value, custom range. Note the "i" is iota, i.e. imaginary.
%   6) [factor + grp*i]: Unknown parameter that is a multiplicative factor of
%                        another parameter. Parameter of this format is
%                        dependent. Note the "i" is iota, i.e. imaginary.
% grp is a positive integer. All parameters with grp that is the same
% integer are in the same group.

%% Concentration Setting
%
% spcMode = 'a';
%
% ~~~Guide~~~
% "species mode" is a switch that determines how quantification in the
% model are handled. It can either be 'a' or 'c' for 'absolute' or
% 'concentration'. For 'absolute', all quantifications are considered in
% amounts, that is molarity (for instance). For 'concentration', all
% quantifications are considered in concentration, that is molar value (for
% instance).
% 
% The impact on this is in whether the quantification is divided by the
% compartment size or not (is divided if absolute, is not otherwise).

spcMode = 'a'; 


%% Compartment definition
%
% modComp = {'Compartment name', relative size};
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
modComp = {'Cyto', 1;
           'PM'  , 0.05;
           'ES'  , Inf;
           };

%% Model species definition
%
% modSpc ={'State name', 'Compatment'  , conc/param};
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
modSpc = {'AKT'     ,'Cyto', 0;
          'mAKT'    ,'PM'  , NaN;
          'p473mAKT','PM'  , 0;
          'p473AKT' ,'Cyto', 0;
          'mTORC2'  ,'PM'  , 1};

%% Features of default parameters
%
% Bnd.(param) = [lb ub]
% (param) is the parameter type
%
% ~~~Guide~~~
% The above gives the default boundaries for each type of parameter. There
% are currently five parameter types:
%   k0   - Zeroth order synthesis type
%   k1   - First order degradataion type
%   k2   - Second order associaition type
%   Km   - Association/Dissociation/Michaelis constant type
%   Conc - Concentration type
%   n    - Hill coefficient
%   r    - Geometry factor (defines effective fraction of overlap area that
%                            is reactive in a cross compartment reaction)
%   comp - Compartment size
%
Bnd.k0   = [1e01 1e04];
Bnd.k1   = [5e-5 0.5];
Bnd.k2   = [5e-5 5e-1];
Bnd.Km   = [1e-2 1e02];
Bnd.Conc = [1e-1 1e1];
Bnd.n    = [1 4];
Bnd.r    = [0 1];
Bnd.Comp = [0 1];

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
%   - r        : Geometry factor
%   - n        : Hill coefficient
%   
% Begin a new reaction by placing (end+1) after rxn. Continue a reaction by
% placing (end) after rxn [The reason for this is due to programming. Just
% follow this rule and you will not run into trouble].

% The combination used for each reaction determines the rules used in the
% backend of the program. This will be explained more later.
%
% For example the following reaction:
% - Describes the enzymatic conversation of mAKT to p473mAKT by mTORC2
% - The reaction is modelled with a hill function
% - Reaction has an unknown Km and a known reaction rate of 0.1mol/s,
% - "r" is an unknown parameter but is between 0.1 and 1.
% - The hill coefficient is 3x the value of the geometry factor "r". So "n"
%   and "r" are grouped together as parameter group 1.

rxn(end+1).desc = 'mAKT -> p473mAKT | mTORC2';
    rxn(end).sub = 'mAKT';  
    rxn(end).prod= 'p473mAKT'; 
    rxn(end).enz = 'mTORC2';
	rxn(end).Km  = NaN; 
    rxn(end).k   = 0.1; 
	rxn(end).r   = [NaN 1 0.1 1]; 
	rxn(end).n   = [3 1]; 

 % More information on types of reactions:
 % If the fields (excluding label and k) are:
 %      - 1 Sub               : Degradataion
 %      - 1 Sub, 1 Enz        : Enzyme mediated degradataion
 %      - 1 Sub, 1 Prod       : Transformation
 %      - 1 Sub, 2+ Prod      : dissociation
 %      - 2 Sub, 1 Prod       : Association 
 %      - 1 Prod              : Constant synthesis
 %      - 1 Prod, 1 Enz       : Enzyme mediated synthesis
 %      - 1 Sub, 1 Enz, Prod  : Enzymatic reaction that can produce any 
 %                              number of products
 %		- n, Km and k included: Hill function