function Y = compDis(model,Y)

% Y2 = compDis(Y,G)
x = model.conc.tens;
G = model.param(1).tens;

if ~isempty(G)
	% Remove the terms which correlate with change in enzyme/substrate
	% concentration
	G(G(:,1)==G(:,2) | G(:,1)==G(:,3),:)=[];

	% Column 1 now contains all the complex indexes. We use only column 2 to
	% identify enzymes and substrates as column 3 is a duplicate. Column 4 is
	% the Km's so we remove them.
	G(:,[3 4]) = [];

	% Sort G by sub/enz index
	[~,I] = sort(G(:,2)); 
	G = G(I,:);

	M = eye(size(x,1));
	scale = model.comp.tens(model.conc.comp);
	scale = scale(:,ones(1,length(scale)));
	M((G(:,2)-1)*max(G(:,1))+G(:,1))=1;
	Y = Y*(M.*scale);
	Y = Y./(scale(:,ones(1,size(Y,1))))';
	
	% Remove complexes
	rmIndx = unique(G(:,1));
	Y(:,rmIndx) = [];
end