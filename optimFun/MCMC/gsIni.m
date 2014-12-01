function dAdt = gsIni(t,mu,sig)

% dAdt = gsIni(t,mu,sig)
%
% Generate a gaussian function. t is the time, mu is the time when dAdt is 
% greatest. Sig is the standard deviation.

dAdt = 1/sqrt(2*pi())/sig*exp(-(t-mu).^2/2/(sig)^2);