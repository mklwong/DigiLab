%%%testCase%%%

% This program tests whether the odeQSSA program is working correctly.
clear resid

v = rand(1,4);
pars = [10.^((rand(1,12)-0.5)*4) rand(1,1)*3+1];

%% Synthesis
method = 'syn';

x0 = [0 0 0 0 0 0 0];
k = [pars(1)];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k),[0 10],x0);
[tMod,~,YMod] = findTC(@SynProg,tReal,'p',[v(1:2) pars(1)],'-b');

resid(1) = sum((YReal(:,1)-YMod(:,1))./(YReal(:,1)+eps));

%% Uni
method = 'uni';
%Dissociation

x0 = [pars(1) 0 0 0 0 0 0];
k = [pars(2) 0 0 0];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k,v),[0 10],x0);
[tMod,~,YMod] = findTC(@UniDis,tReal,'p',[v(1:2) pars(2)],'y0',x0(1:3),'-b');
resid(2) = sum(sum((YReal(:,1:3)-YMod(:,1:3))./ (ones(length(tReal),1)*max(YReal(:,1:3))) ));

%Interconversion
x0 = [pars(1) 0 0 0 0 0 0];
k = [0 pars(2) 0 0];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k,v),[0 10],x0);
[tMod,~,YMod] = findTC(@UniCon,tReal,'p',[v(1:2) pars(2)],'y0',x0(1:2),'-b');
resid(3) = sum(sum((YReal(:,[1 2])-YMod(:,[1 2]))./(ones(length(tReal),1)*max(YReal(:,1:2))) ));

%Degradation
x0 = [pars(1) 0 0 0 0 0 0];
k = [0 0 pars(2) 0];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k,v),[0 10],x0);
[tMod,~,YMod] = findTC(@UniDeg,tReal,'p',[v(1:2) pars(2)],'y0',x0(1:3),'-b');
resid(4) = sum((YReal(:,1)-YMod(:,1))./(ones(length(tReal),1)*max(YReal(:,1)) ));

%SynMM
x0 = [0 pars(1) 0 0 0 0 0];
k = [0 0 0 pars(2)];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k,v),[0 10],x0);
[tMod,~,YMod] = findTC(@SynMM,tReal,'p',[v(1:2) pars(2)],'y0',x0(1:3),'-b');
resid(5) = sum((YReal(:,1)-YMod(:,1))./(ones(length(tReal),1)*max(YReal(:,1)) ));

%% Bi
method = 'bi';
% Association
x0 = [pars(1) pars(2) 0 0 0 0 0];
k = [pars(3) 0 0];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k,v),[0 10],x0);
[tMod,~,YMod] = findTC(@BiAss,tReal,'p',[v(1:2) k(1)],'y0',x0(1:3),'-b');
resid(6) = sum(sum((YReal(:,1:3)-YMod(:,1:3))./(ones(length(tReal),1)*max(YReal(:,1:3))) ));

% MM Enzyme Kinetics
x0 = [pars(1) pars(2) 0 0 0 0 0];
k = [0 pars(3) 0];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k,v),[0 10],x0);
[tMod,~,YMod] = findTC(@BiMM,tReal,'p',[v(1:2) k(2)],'y0',x0(1:3),'-r','-b');
resid(7) = sum(sum((YReal(:,1:3)-YMod(:,1:3))./(ones(length(tReal),1)*max(YReal(:,1:3))) ));

% MM Enzyme Kinetics Deg
x0 = [pars(1) pars(2) 0 0 0 0 0];
k = [0 0 pars(3)];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k,v),[0 10],x0);
[tMod,~,YMod] = findTC(@BiMMDeg,tReal,'p',[v(1:2) k(3)],'y0',x0(1:3),'-b');
resid(8) = sum(sum((YReal(:,1)-YMod(:,1))./(ones(length(tReal),1)*max(YReal(:,1))) ));

%% Enzyme Kinetic
%% Single
method = 'enzQSSA';
x0 = [pars(1) pars(2) 0 0 0 0 0];
k = [pars(3:5) 0 0 0 pars(7)];
kQSSA = [k(3) (k(2)+k(3)+eps)/k(1) k(6) (k(5)+k(6)+eps)/k(4)];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k,v),[0 200],x0);
[tMod,~,YMod] = findTC(@enzKinetic,tReal,'p',[v kQSSA k(end)],'y0',x0(1:4),'-b');
resid(9) = sum((YReal(end,[1:3 5])-YMod(end,[1:3 5]))./YReal(end,[1:3 5]) );

%% Reversible
method = 'enzQSSA';
x0 = [pars(1) pars(2) pars(10) pars(11) 0 0 0];
k = [pars(3:8) 0];
kQSSA = [k(3) (k(2)+k(3))/k(1) k(6) (k(5)+k(6))/k(4)];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k,v),[0 200],x0);
[tMod,~,YMod] = findTC(@enzKinetic,tReal,'p',[v kQSSA k(end)],'y0',x0(1:4),'-b');
resid(10) = sum((YReal(end,1:6)-YMod(end,1:6))./YReal(end,1:6) );

%% Hill Function
method = 'hillFun';
x0 = [pars(1) 0 pars(2) 0 0 0 0];
k = [pars(3:7) pars(12) pars(8)];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k,v),[0 200],x0);
[tMod,~,YMod] = findTC(@hillFun,tReal,'p',[v(1:2) k],'y0',x0(1:3),'-b');
resid(11) = sum((YReal(:,1)-YMod(:,1))./(ones(length(tReal),1)*max(YReal(:,1))+eps));

%% Complex Model
[k,kSSA,kSSA_r,xi,inp] = makeParam(0,1);
xi(1) = 1;
[tReal,YReal] = ode15s(@(t,x) ComplexMA(t,x,k,v(1:2)),[0 500],[xi zeros(1,10)]);
[tMod,Y,YMod] = findTC(@ComplexdQSSA,tReal,'p',[v(1:2) kSSA],'y0',xi,'-b');
resid(12) = 10.^(-mean(-log10(abs(YReal(end,:)-YMod(end,:))./max([YReal(end,:);YMod(end,:)]))));
