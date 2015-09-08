function [val,freeParam,grp,bnd] = testPar(testVal)

% This function tests three value parameter set for model parameter.
%
% It tests for:
%	- Whether parameter is free or not.
%	- If free, whether it is grouped or not. 
%		- If grouped, which group it belongs to.
%		- If grouped, places the multiplicative factor in first value of
%		triad for further processing.
%		- If ungrouped, sets flag "ungrouped parameter" and sets 
%		  multiplicative factor as 1.
%	- If not free, sets flag "not a free parameter".

% Default values
freeParam = false;
grp = 0;
bnd = [];

% Parameter not used in this reaction. Skip.
if isempty(testVal)     
	val = [];
    return;
end

cond = testVal(1);

if isnan(cond)   % Free parameter.
    val = 1; % Multiplicative factor set to 1.
    freeParam = true;
	if length(testVal) == 2 || length(testVal) == 4     % Grouped
		grp = real(testVal(2)); % Save the group
		val = imag(testVal(2)); % Save the multiplicative factor
		freeParam = true;
		if val == 0 %if by laziness no real component give, set it as 1.
			val = 1;
		end
		if length(testVal) == 4
			bnd = testVal(3:4);
		end
	elseif length(testVal) == 3 % With boundaries
		bnd = testVal(2:3);
	end
elseif isreal(cond)  % Is a fixed parameter.
    val = testVal(1);
end