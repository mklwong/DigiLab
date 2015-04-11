function [subIndx,prodIndx,enzIndx,rxnType] = classifyReaction(rxn,x)

[~,~,subIndx]  = intersect(upper(rxn.sub) ,upper(x.name));
[~,~,prodIndx] = intersect(upper(rxn.prod),upper(x.name));
[~,~,enzIndx]  = intersect(upper(rxn.enz) ,upper(x.name));

nSub = length(subIndx);
if nSub ==0
    if ~isempty(rxn.enz) 
        subIndx  = [subIndx enzIndx];
        prodIndx = [prodIndx enzIndx];
        rxnType = 'uni';
    else
        rxnType = 'syn';
    end
elseif nSub == 1
    if ~isempty(rxn.enz) 
        if isempty(rxn.Km)
            subIndx  = [enzIndx subIndx];
            prodIndx = [enzIndx prodIndx];
            rxnType = 'bi';
        elseif ~isempty(rxn.Km)
            rxnType = 'enzQSSA';
        else
            error('parseModelm:exceptionClassifier1','Unknown error');
        end
    else
        rxnType = 'uni';
    end
elseif nSub == 2
    rxnType = 'bi';
else
    error('parseModelm:TooManySubstrates',['There are too many substrates in reaction ' num2str(ii)])
end