function YOut = extractSim(outCombo,YSim,YSim2)

% This function creates the customs states from the definitions as provided
% by the model file.
%
% outCombo might be something like this:
%   {[1 3];[4 -6]}
%
% In this example, YOut will have 2 states, the first state be the sum of
% YSim(:,1) and YSim(:,3). The second state will be YSim(:,4)+YSim2(:,6).

% Initialise YTmp, based on number of time points and number of states 
% needed to be generated.
YOut = zeros(length(YSim(:,1)),length(outCombo));

for ii = 1:length(outCombo)
    indx = outCombo{ii};
    YOut(:,ii) = sum(YSim(:,indx(indx>0)),2)+sum(YSim2(:,-indx(indx<0)),2);
end