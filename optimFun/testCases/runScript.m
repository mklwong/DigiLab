opts = MCMCOptimset('display','terminal','T',10,'parMode',2,'ptNo',10000,'Pmin',0);

% Rosenbrock function
opts = MCMCOptimset(opts,'dir','./Rosenbrock');
obj = @(X,Y) (10*(Y-X.^2)).^2+(X-1).^2;
bnd = [-10 10; -10 10];
pts = optimFuncTest(@(p)obj(p(1),p(2)),bnd,opts);
[X,Y,Z] = surfFuncTest(obj,bnd);
figure(1)
subplot(2,1,1);
hist3c(pts,[100 100])
shading interp
xlim(bnd(1,:))
ylim(bnd(2,:))
subplot(2,1,2);
surf(X,Y,exp(-Z))
shading interp
view(2)

% Ackley's function
opts = MCMCOptimset(opts,'dir','./Ackleys');
obj = @(X,Y) -5*log10(1-((-20*exp(-0.2*sqrt(0.5*(X.^2+Y.^2)))-exp(0.5*(cos(2*pi()*X)+cos(2*pi()*Y)))+exp(1)+20)/14.4));
bnd = [-5 5; -5 5];
pts = optimFuncTest(@(p)obj(p(1),p(2)),bnd,opts);
[X,Y,Z] = surfFuncTest(obj,bnd);
figure(2)
subplot(2,1,1);
hist3c(pts,[100 100])
xlim(bnd(1,:))
ylim(bnd(2,:))
shading interp
subplot(2,1,2);
surf(X,Y,(-Z))
shading interp
view(2)

% Eggholder function
opts = MCMCOptimset(opts,'dir','./Eggholder');
obj = @(X,Y) -10*log10(1-(-(Y+47).*sin((abs(X/2+Y+47)).^0.5)-X.*sin((abs(X-Y-47)).^0.5)+959.6408)/2100);
bnd = [-512 512; -512 512];
pts = optimFuncTest(@(p)obj(p(1),p(2)),bnd,opts);
[X,Y,Z] = surfFuncTest(obj,bnd);
figure(3)
subplot(2,1,1);
hist3c(pts,[100 100])
xlim(bnd(1,:))
ylim(bnd(2,:))
shading interp
subplot(2,1,2);
surf(X,Y,exp(-Z))
shading interp
view(2)