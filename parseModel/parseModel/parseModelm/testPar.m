function [val,freeParam,grp] = testPar(testVal)

% Default values
val = 0;   % Val is either the parameter value (if fixed) or the multiplicative factor (if unfixed).
freeParam = false;
grp = 0;

if isempty(testVal)     % Parameter not used in this reaction. Skip.
    return;
end

cond = testVal(1);


if isnan(cond)   % Is an ungroup free parameter.
    val = 1; % Multiplicative factor set to 1.
    freeParam = true;
elseif isreal(cond)  % Is a fixed parameter.
    val = testVal(1);
elseif ~isreal(cond)
% Is a grouped free parameter. 
% Imaginary part defines group.
% Real part defines multiplicative factor.      
%
% In a grouped parameter, sometimes the real part is not given is the
% multiplicative factor is one. We account for this below
    if testVal(1) == 0
        testVal(1) = 1;
    end
    val = real(testVal(1)); % Save the multiplicative factor 
    freeParam = true;
    grp = imag(testVal(1));
end