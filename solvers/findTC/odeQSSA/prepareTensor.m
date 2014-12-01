function [k0,k1,k2,G,model] = prepareTensor(model,p)

if nargin == 1
	model = parseModel(model);
elseif nargin == 2
	model = parseModel(model,p);
end

% Get x0 to know the required tensor size
x0 = model.x.tens;

% Generate tensors using sparse function
k0 = full(sparse(model.k0.tens(:,1),model.k0.tens(:,1)*0+1,model.k0.tens(:,2)));
k1 = full(sparse(model.k1.tens(:,1),model.k1.tens(:,2),model.k1.tens(:,3)));
k2 = model.k2.tens;
G  = model.G.tens;

% Make up matrices to correct dimension
if length(k0)~= length(x0)
	k0(size(x0,1),1) = 0;
end
if (size(k1,1)~= length(x0) || size(k1,2)~= length(x0))
	k1(size(x0,1),size(x0,1)) = 0;
end

