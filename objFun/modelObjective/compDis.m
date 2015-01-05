function Y = compDis(model,Y)

% Y2 = compDis(Y,G)
x = model.x.tens;
G = model.G.tens;

if ~isempty(G)
	% Remove the terms which correlate with change in enzyme/substrate
	% concentration
	G(G(:,1)==G(:,2) | G(:,1)==G(:,3),:)=[];

	% Column 1 now contains all the complex indexes. We use onl column 2 to
	% identify enzymes and substrates as column 3 is a duplicate. Column 4 is
	% the Km's so we remove them.
	G(:,[3 4]) = [];

	% Sort G by sub/enz index
	[~,I] = sort(G(:,2)); 
	G = G(I,:);

	M = eye(size(x,1));
	M((G(:,2)-1)*max(G(:,1))+G(:,1))=1;
	Y = Y*M;
end