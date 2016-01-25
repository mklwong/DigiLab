function [kSSA kSSA_r] = getQSSA(k) 

%Create Michaelis constant for QSSA models
getQSS = @(x) [x(3) (x(2)+x(3)+eps)/x(1)]; 
kSSA = [getQSS(k(1:3)),...   %Rnx 1f
	    getQSS(k(4:6)),...   %Rnx 1r
		k(7:9),...           %Rnx 2-4
		getQSS(k(10:12)),... %Rnx 5f
		getQSS(k(13:15)),... %Rnx 5r
	    k(16),...            %Rnx 6f
		getQSS(k(17:19)),... %Rnx 7f
		getQSS(k(20:22)),... %Rnx 7r
		k(23),...            %Rnx 8
		getQSS(k(24:26)),... %Rnx 9f
		getQSS(k(27:29)),... %Rnx 9r
	    getQSS(k(30:32)),... %Rnx 10f
	    getQSS(k(33:35)),... %Rnx 10r
		k(36)];              %Rnx 11

% Paramters for generating initial conditions (setting all non-Michaelis
% constant to zero).
kSSA_r = kSSA;
kSSA_r([1 3 5 6 7 8 10 12 13 15 17 18 20 22 24 26]) = 0; 