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
% We would define this by:
%
xMod = {'AKT'     ,'Cyto',NaN;
        'mAKT'    ,'PM',  0;
        'p473mAKT','PM',  0;
        'p473AKT' ,'Cyto',0;};



%% Relationship between simulation state and model state association

xData = {'Exp Name',{'*Sim 1','*Sim 2'}};
% ~~~Guide~~~
% "Exp Name" is the name of the state as labelled in the experimental
% data.
%
% "Sim 1" and "Sim 2" are the state names as labelled in the simulation.
% Each row relates the experimental state to a set of equivalent states
% in the simulation. The results of the states are added together and
% then compared to the experimental state.
%
% An asterix (*) at the front of simulation state names implies the
% enzyme-substate complexes for that state is included as well.
%
%
% ~~~Example~~~
% Phosphorylation of p473AKT is experimental state to be compared. The
% experiment is a Western Blot without cell fractionation, hence plasma
% membrane and cytosolic AKT are captured. In our model the two pools of
% AKT are discrete, hence the group of simulation states corresponding to
% the experimental state is p473AKT and p473mAKT. Additionally, enzyme
% substrate complexes for both the simulation states are included, so they
% need to be asterixed.
%
% The final row for xData is hence:
%	- {'p473AKT',{'*p473AKT',p473mAKT'}}
%       ^string    ^asterix
%     of name of    marks the
%     name of exp   E-S complex
%     state         is included

%% Features of default parameters
Bndk0  = [1e01 1e04];
Bndk1  = [5e-5 0.5];
Bndk2   = [5e-5 5e-1];
BndKm   = [1e-2 1e02];
BndConc = [1e-1 1e1];
% ~~~Guide~~~
% The above gives the default boundaries for each type of parameter. The
% variable names above are "bound - param type". So Bndk0 is the boundary
% of k0 type parameters. k0 type are zeroth order types.
%
% All that needs to be done is to place the lower and upper bound enclosed
% in square brackets.

%% Reactions
v(end+1).label = 'S->P | E';
    v(end).sub = 'S';  
    v(end).prod= 'P'; 
    v(end).enz = 'E';
	v(end).Km  = NaN; 
    v(end).k  = [-1]; 
% ~~~Guide~~~
% The above defines a reaction. The first line is alway:
%	>  v(end+1).label = 
%	A string (i.e. letters enclosed in ' ). This is simply used to label
%	the parameters on output. So pick a label that describes the reaction
% The following lines define the participating species and reaction
% parameters. For reactants, begin the line with:
%	> v(end).sub =
%	Note: Maximum of two reactants are allowed.
% For products begin the line with:
%	> v(end).prod =
% For enzymes, begin the line with:
%	> v(end).enz =
%   Note: maximum of 1 enzyme is allowed.
% For each of the above, list the species within braces, and enclose each
% species in inverted commas separated by commas. E.g: {'AKT','PIP3'}. The
% braces can be ignored if only one species is used.
%
% The number of reactants, product and enzymes used will determine the
% reaction type. This is classified by the program.
%
% For parameters, two parameter types are allowed:
%	> v(end).k  = 
%	This is a normal rate parameters
%	> v(end).Km  = 
%	This is a quasi-steady state parameter (includes Michaelis constant and
%	dissociation constants).
%
% Again the definition of the parameters defines the way the reaction is
% implemented. E.g. only k given models by mass action, only Km given
% models as quasi steady state, both k and Km given models as enzyme
% reaction using dQSSA.