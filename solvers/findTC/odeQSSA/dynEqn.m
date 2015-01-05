function dx_dt = dynEqn(t,x,G,k0,k1,k2)
%
% Evaluates the tensor equation dx_dt = Inv(I - G*x)\(k2*x*x + k1*x + k0(t))

x(x<0) = 0; %sometimes the system goes to less than zero. When it wants to do this, set it to zero

M = zeros(size(k1));
L = M;

MTmp = sparse(G(:,1),G(:,3),G(:,4).*x(G(:,2)));
[a,b] = size(MTmp);
M(1:a,1:b) = MTmp;

LTmp = sparse(k2(:,1),k2(:,2),k2(:,4).*x(k2(:,3)));
[a,b] = size(LTmp);
L(1:a,1:b) = LTmp;

dx_dt = (eye(length(x))+M)\(L*x+k1*x+k0(t));

end