classdef sigRxnList
	properties
		desc = '';
		sub = cell(0,0);  
    	prod= cell(0,0); 
    	enz = cell(0,0);
		Km = [];
    	k  = [];  
		r  = []; %Is 1 by default
		n  = [];
	end
	methods
		%% Pre checks of class attribute types
		% String, species types
		function obj = set.desc(obj,descIn)
			if ~ischar(descIn)
				error('sigMatMod:DescAssignmentWrong','Error entering description. Can only be a string.')
			end
			obj.desc = descIn;
		end
		function obj = set.sub(obj,subIn)
			if ischar(subIn) && ~isempty(subIn)
				subIn = {subIn};
			elseif isempty(subIn)
				subIn = cell(0,0);
			elseif ~iscell(subIn)
				error('sigMatMod:SubAssignmentWrong','Error entering substrate list. Can only be a string or cell.')
			end
			obj.sub = subIn;
				end
		function obj = set.prod(obj,prodIn) 
			if ischar(prodIn) && ~isempty(prodIn)
				prodIn = {prodIn};
			elseif isempty(prodIn)
				prodIn = cell(0,0);
			elseif ~iscell(prodIn)
				error('sigMatMod:ProdAssignmentWrong','Error entering product list. Can only be a string or cell.')
			end
			obj.prod = prodIn;
		end
		function obj = set.enz(obj,enzIn)
			if ischar(enzIn)&& ~isempty(enzIn)
				enzIn = {enzIn};
			elseif isempty(enzIn)
				enzIn = cell(0,0);
			elseif ~iscell(enzIn)
				error('sigMatMod:EnzAssignmentWrong','Error entering enzyme list. Can only be a string or cell.')
			end
			obj.enz = enzIn;
		end
		
		% Numerical, parameter types
		function obj = set.Km(obj,KmIn)
			if ~isnumeric(KmIn)
				error('sigMatMod:KmAssignmentWrong','Error entering Equilibrium Constant. Can only be a numerical array.')
			end
			obj.Km = KmIn;
		end
		function obj = set.k(obj,kIn)
			if ~isnumeric(kIn)
				error('sigMatMod:kAssignmentWrong','Error entering rate constant. Can only be a numerical array.')
			end
			obj.k = kIn;
		end
		function obj = set.r(obj,rIn)
			if ~isnumeric(rIn)
				error('sigMatMod:rAssignmentWrong','Error entering geometry factor. Can only be a numerical array.')
			end
			obj.r = rIn;
				end
		function obj = set.n(obj,nIn)
			if ~isnumeric(nIn)
				error('sigMatMod:nAssignmentWrong','Error entering hill coefficient. Can only be a numerical array.')
			end
			obj.n = nIn;
		end
		
		%% Other methods
	end
end