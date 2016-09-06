function [strct,pInd,pList,testStrct] = testPar(strct,pList)

% Description of outputs
% - strct is the final reaction structure after the parameter values are
% changed
% - testStrct is the reaction structure for building the parameter matrix
% (which has all parameter values set to NaN).
% - pInd gives the parameter index which will be applied for each parameter
% - pList updated pList with grouped parameter sets updated

testStrct = strct;

% Loop through each field and initialise pInd from that
flds = fieldnames(strct);
pInd = NaN(1,length(flds));

for ii = 1:length(flds)
    %Parameters are only numeric
    if isnumeric(strct.(flds{ii}))
        [strct.(flds{ii}),custBnd,grp] = classPar(strct.(flds{ii}));
        
        %Determine if free parameter necessary
        if ~all(isnan(custBnd)) || (all(custBnd==0)&&~isnan(grp))
            % Assume group found. Find the parameter index for that group
            if ~isnan(grp)
                % Match group number to existing groups. If doesn't exist 
                pInd(ii) = max([pList.grp(pList.grp(:,1)==grp,2) NaN]);
            end
           
            % If parameter index not inserted, then that means a new one
            % needs to be made
            if isnan(pInd(ii))
                pList.npar = pList.npar+1;
                pInd(ii) = pList.npar;
                if ~isnan(grp)
                    pList.grp = [pList.grp;[grp pInd(ii)]];
                end
            end
            
            % Enter custom boundaries
            if ~all(custBnd==0)
                pList.lim(pInd(ii),:) = custBnd;
			end
			testStrct.(flds{ii}) = pInd(ii)*1i;
		else
			testStrct.(flds{ii}) = -abs(testStrct.(flds{ii}));
		end
    end
end

end

%%%%%%%%%%%%%
%%%%%%%%%%%%%

function [matVal,bnd,grp] = classPar(vec)
% matVal is the values that will be entered into the matrix
% bnd    is the boundary. If
%           - 1x2 array, then that is a custom boundary
%           - []       , default/group boundary used
%           - NaN      , fixed and no boundary required
% grp    is the group of this parameter. NaN means not in group
if isempty(vec)
    matVal = [];
    bnd = nan;
    grp = nan;
elseif length(vec)==1
    grp = NaN;
    if ~isnan(vec(1))
        bnd     = NaN;
        matVal  = vec(1);
    else
        bnd     = 0;
        matVal  = 1;
    end
elseif length(vec) == 2
    grp     = vec(2);
    bnd     = [];
    if ~isnan(vec(1))
        matVal = vec(1);
    else
        matVal = 1;
    end
elseif length(vec) == 3
    if ~isnan(vec(1))
        %Leading numeric ambigious
        error('SIGMAT.testPar:AmbigVal','First value in 1x3 parameter array is not NaN. Interpretation is ambigious.')
    else
        bnd = vec(2:3);
        grp = NaN;
        matVal = 1;
    end
elseif length(vec) == 4
    matVal = 1;
    bnd = vec(3:4);
    grp = vec(2);
else
    %Incorrect number of values passed in parameter
    error('SIGMAT.testPar:WrongNumVal','Incorrect number of values passed in parameter')
end
end