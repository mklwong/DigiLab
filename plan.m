Define parameter
%rateVal - [prod_indx,R1_indx,geo_pindx,k1_pindx] (odeKinetic)
%bndType – {boundary type of each row above} (odeKinetic) – for indexes these are obviously not relevant, so it nan.

What the definition does is the boundary type of all required parameters are created. This list of parameter boundary types are then cross checked against the list of boundaries given by the input model to ensure that all default boundaries are specified in the model. 

Create p_index field for all parameter that mirrors the rateVal field
% XX.pInd = XX.rateVal (parseModelm)

Next we move onto the reaction classification

Parse Reaction
Test parameters in reaction and assign parameter index
% [rxn,pInd,pList,grp] = testPar(rxn,)
-	Rxn are modified to give the multiplicative factor if is free parameter. Gives the nominal value is not free parameter.
-	pInd = NaN if not free param, is the pIndex it is assigned if free parameter.
-	pList is the list of parameters, with new parameters created
Classify reaction
% [model.rxnRules('rxnRules',rxn(ii),conc,comp,expComp,ii);

Next we move onto the reaction classification

Fill in gaps in knowledge about parameter
We need to fill in
-	Boundary: based on the classification and hence what each free parameter now means
-	Parameter description. Again we know what each parameter means.
-	
