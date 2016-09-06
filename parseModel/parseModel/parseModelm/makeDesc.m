function outDesc = makeDesc(sub,prod,enz,param)

% Add +'s and spaces to sub list as necessary
if isempty(sub)
	sub = 'Phi';
end
if ischar(sub)
	sub = {sub};
end

if size(sub,1)~=1
	sub = sub';
end
appendSub = cell(size(sub));
appendSub(:) = {' + '};
appendSub(end) = {''};
sub = [sub;appendSub];

% Add +'s and spaces to prod list as necessary
if isempty(prod)
	prod = 'Phi';
end
if ischar(prod)
	prod = {prod};
end
if size(prod,1)~=1
	prod = prod';
end
appendprod = cell(size(prod));
appendprod(:) = {' + '};
appendprod(end) = {''};
prod = [prod;appendprod];

if ischar(enz)
	enz = {enz};
end

% Compile the description
if isempty(enz)
	outDesc = [param ' : ' horzcat(sub{:}) ' -> ' horzcat(prod{:})];
else
	outDesc = [param ' : ' horzcat(sub{:}) ' -> ' horzcat(prod{:}) ' | ' enz{:}];
end