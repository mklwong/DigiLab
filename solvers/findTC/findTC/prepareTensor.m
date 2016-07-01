function model = prepareTensor(model,p)

if nargin == 1
	model = parseModel(model);
elseif nargin == 2
	model = parseModel(model,p);
end

% Get x0 to know the required tensor size
x0 = model.conc.tens;

% Generate tensors using sparse function
model.tensor.k0 = full(sparse(model.param(4).tens(:,1),ones(size(model.param(4).tens(:,1))),model.param(4).tens(:,2)));
model.tensor.k1 = full(sparse(model.param(3).tens(:,1),model.param(3).tens(:,2),model.param(3).tens(:,3)));
model.tensor.k2 = model.param(2).tens;
model.tensor.G  = model.param(1).tens;

% Make up matrices to correct dimension
if length(model.tensor.k0)~= length(x0)
	model.tensor.k0(size(x0,1),1) = 0;
end
if (size(model.tensor.k1,1)~= length(x0) || size(model.tensor.k1,2)~= length(x0))
	model.tensor.k1(size(x0,1),size(x0,1)) = 0;
end

