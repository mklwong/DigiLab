Define parameter
%rateVal - [prod_indx,R1_indx,geo_pindx,k1_pindx] (odeKinetic)
%bndType – {boundary type of each row above} (odeKinetic) – for indexes these are obviously not relevant, so it nan.

What the definition does is the boundary type of all required parameters are created. 
This list of parameter boundary types are then cross checked against the list of boundaries 
given by the input model to ensure that all default boundaries are specified in the model. 

Create p_index field for all parameter that mirrors the rateVal field
% XX.pInd = XX.rateVal (parseModelm)

Next we move onto the reaction classification



Parse Reaction
Test parameters in reaction and assign parameter index
% [rxn,pInd,pList] = testPar(rxn,pList,defBnd)
Input
- obj(ii): The object to be parsed
Output
- obj(ii): Spits the object back out with original numeric elements modified to give the 
multiplicative factor if is free parameter. Gives the nominal value is not free parameter.
- pInd: NaN if not free param, is the pIndex it is assigned if free parameter.
Is a list of all the parameter fields
- pList: is the list of free parameters in the model, with new parameters created. 
One field is the number of free parameters. Another lists the parameter groups and
the parameter index the group corresponds to, and finally the 

Classify reaction
[reqParam,rateVal,geoVal,parDesc,conc] = model.rxnRules('rxnRules',rxn(ii),conc,comp,expComp,ii);
Inputs
- 'rxnRules': flag
- rxn(ii): input the entire ii-th reaction for classification
- conc: vector of species
- flags: miscellaneous other important flags 
Outputs
- reqParam: list of parameter fields that need to be modified
- rateVal: list of arrays that need to be appended to each paramater field previously specified
- parDesc: descriptions of the parameters
- conc: spat back out in case new species have been added

Put Parameter in pre-matrix

Fill in gaps in knowledge about parameter
We need to fill in
-	Boundary: based on the classification and hence what each free parameter now means
-	Parameter description. Again we know what each parameter means.
-	
