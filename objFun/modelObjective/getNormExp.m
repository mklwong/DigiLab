function normExp = getNormExp(normExp,expGrp,YSim,YExp,matches)

% Determine if normalising factors for this group exists
if length(normExp)<=expGrp   %Case where it doesn't
    normExp{expGrp,1} = matches;
    tmp               = exp(-max(YSim)./max(YExp)*1e6)+max(YSim)./max(YExp);
    normExp{expGrp,2} = tmp;
else                        %Case where it does
    [~,svdIndx] = intersect(upper(matches),upper(normExp{1,expGrp})); %Determine which normfacs already saved
    
    %Delete duplicates
    YSim(:,svdIndx)  = [];
    YExp(:,svdIndx)  = [];
    matches(svdIndx) = [];
    
    %Append match names
    normExp{expGrp,1} = [normExp{expGrp,1};matches];
    
    % Get normalising factors
    tmp               = exp(-max(YSim)./max(YExp)*1e6)+max(YSim)./max(YExp);
    normExp{expGrp,2} = [normExp{expGrp,2} tmp];
end