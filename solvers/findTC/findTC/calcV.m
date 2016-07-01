function dx_dt = calcV(model,t,x,p)

if isrow(x)
	x = x';
end

[k0,k1,k2,G] = prepareTensor(model,p);
k0 = @(t) k0;

dx_dt = dynEqn(t,x,G,k0,k1,k2);