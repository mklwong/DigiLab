function [k,kSSA,kSSA_r,xi,inp] = makeParam(tInjStrt,tInjEnd)

k = 10.^((rand(1,36)-0.5)*2);
xiv = 10.^((rand(1,36)-0.5)*2);
%k([1:7 17:end]) = 0;

% Concentration rules
xiv(4) = xiv(3)*(rand(1)+0.5); % Enforce concentration relationship between C and D
xiv(5) = 100;				   % Enforce concentration relationship between C and D

% Set dominate direction in enzymatic reaction and enforce rapid
% equilibrium assumption.
k = adjustK(k);

% Determine quasi-steady state parameters (i.e. Michaelis Constants)
[kSSA kSSA_r] = getQSSA(k);

% Create initial condition vector
xi = [0 xiv(1) 0 xiv(2) 0 0 xiv(3) xiv(4) 0 0];

if nargin > 0
	% Time profile for I injection
	inpAmp = [xiv(5);zeros(19,1)];
	inp{1} = @(t)inpAmp                    *normpdf(t,(tInjEnd+tInjStrt)/2,(tInjEnd-tInjStrt)/10);
	inp{2} = @(t)[inpAmp(1:10);zeros(10,1)]*normpdf(t,(tInjEnd+tInjStrt)/2,(tInjEnd-tInjStrt)/10);
else
	inp = @(t) 0;
end
