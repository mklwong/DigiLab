opts = MCMCOptimset('display','terminal','T',10,'parMode',2,'ptNo',70000,'Pmin',0);

% Rosenbrock function
fprintf('-------------------------------------\n')
fprintf('-----Testing Rosenbrock Function-----\n')
fprintf('-------------------------------------\n')
opts = MCMCOptimset(opts,'dir','./Rosenbrock4');
obj = @(X,Y) (10*(Y-X.^2)).^2+(X-1).^2;
bnd = [-10 10; -10 10];
result = optimFuncTest(@(p)obj(p(1),p(2)),bnd,opts);
[X,Y,Z] = surfFuncTest(obj,bnd);
figure(1)
subplot(2,1,1);
hist3c(result.pts,[100 100])
shading interp
xlim(bnd(1,:))
ylim(bnd(2,:))
subplot(2,1,2);
surf(X,Y,10.^(-Z))
shading interp
view(2)

% Ackley's function
fprintf('-------------------------------------\n')
fprintf('------Testing Ackley''s Function------\n')
fprintf('-------------------------------------\n')
opts = MCMCOptimset(opts,'dir','./Ackleys');
obj = @(X,Y) -log10((1-(((-20*exp(-0.2*sqrt(0.5*(X.^2+Y.^2)))-exp(0.5*(cos(2*pi()*X)+cos(2*pi()*Y)))+exp(1)+20)-0.10168)/14.5)));
bnd = [-5 5; -5 5];
result = optimFuncTest(@(p)obj(p(1),p(2)),bnd,opts);
[X,Y,Z] = surfFuncTest(obj,bnd);
figure(2)
subplot(2,1,1);
hist3c(result.pts,[100 100])
xlim(bnd(1,:))
ylim(bnd(2,:))
shading interp
subplot(2,1,2);
surf(X,Y,10.^(-Z))
shading interp
view(2)

% Eggholder function
fprintf('-------------------------------------\n')
fprintf('-----Testing Eggholder Function------\n')
fprintf('-------------------------------------\n')
opts = MCMCOptimset(opts,'dir','./Eggholder');
obj = @(X,Y) -10*log10(1-(-(Y+47).*sin((abs(X/2+Y+47)).^0.5)-X.*sin((abs(X-Y-47)).^0.5)+959.6408)/2100);
bnd = [-512 512; -512 512];
result = optimFuncTest(@(p)obj(p(1),p(2)),bnd,opts);
[X,Y,Z] = surfFuncTest(obj,bnd);
figure(3)
subplot(2,1,1);
hist3c(result.pts,[100 100])
xlim(bnd(1,:))
ylim(bnd(2,:))
shading interp
subplot(2,1,2);
surf(X,Y,10.^(-Z))
shading interp
view(2)