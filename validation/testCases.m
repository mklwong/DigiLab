%%%testCase%%%

% This program tests whether the odeQSSA program is working correctly.

v = [1 2];

%% Synthesis
method = 'syn';

x0 = [0 0 0 0 0 0 0];
k = [1];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k),[0 10],x0);
[tMod,~,YMod] = findTC(@SynProg,tReal,'p',[v 1],'-b');

resid(1) = sum(YReal(:,1)-YMod(:,1));

%% Uni
method = 'uni';
%Dissociation
x0 = [1 0 0 0 0 0 0];
k = [1 0 0 0];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k,v),[0 10],x0);
[tMod,~,YMod] = findTC(@UniDis,tReal,'p',[v 1],'-b');
resid(2) = sum(sum(YReal(:,1:3)-YMod(:,1:3)));

%Interconversion
x0 = [1 0 0 0 0 0 0];
k = [0 1 0 0];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k,v),[0 10],x0);
[tMod,~,YMod] = findTC(@UniCon,tReal,'p',[v 1],'-b');
resid(3) = sum(sum(YReal(:,[1 2])-YMod(:,[1 2])));

%Degradation
x0 = [1 0 0 0 0 0 0];
k = [0 0 1 0];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k,v),[0 10],x0);
[tMod,~,YMod] = findTC(@UniDeg,tReal,'p',[v 1],'-b');
resid(4) = sum(YReal(:,1)-YMod(:,1));

%SynMM
x0 = [0 1 0 0 0 0 0];
k = [0 0 0 1];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k,v),[0 10],x0);
[tMod,~,YMod] = findTC(@SynMM,tReal,'p',[v 1],'-b');
resid(5) = sum(YReal(:,1)-YMod(:,1));

%% Bi
method = 'bi';
% Association
x0 = [1 1 0 0 0 0 0];
k = [1 0 0];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k,v),[0 10],x0);
[tMod,~,YMod] = findTC(@BiAss,tReal,'p',[v 1],'-b');
resid(6) = sum(sum(YReal(:,1:3)-YMod(:,1:3)));

% MM Enzyme Kinetics
x0 = [1 1 0 0 0 0 0];
k = [0 1 0];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k,v),[0 10],x0);
[tMod,~,YMod] = findTC(@BiMM,tReal,'p',[v 1],'-r','-b');
resid(7) = sum(sum(YReal(:,1:3)-YMod(:,1:3)));

% MM Enzyme Kinetics Deg
x0 = [1 1 0 0 0 0 0];
k = [0 0 1];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k,v),[0 10],x0);
[tMod,~,YMod] = findTC(@BiMMDeg,tReal,'p',[v 1],'-b');
resid(8) = sum(sum(YReal(:,1:2)-YMod(:,1:2)));

%% Enzyme Kinetic
method = 'enzQSSA';
% Single
x0 = [1 1 1 0 0 0 0];
k = [1 1 1 1 1 1 1];
kQSSA = [k(3) (k(2)+k(3))/k(1) k(6) (k(5)+k(6))/k(4)];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k,v),[0 30],x0);
[tMod,~,YMod] = findTC(@enzKinetic,tReal,'p',[v kQSSA k(end)],'y0',x0(1:4),'-b');
resid(9) = sum(YReal(end,1:6)-YMod(end,1:6));

% Reversible
x0 = [1 1 0 1 0 0 0];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k,v),[0 30],x0);
[tMod,~,YMod] = findTC(@enzKinetic,tReal,'p',[v kQSSA k(end)],'y0',x0(1:4),'-b');
resid(10) = sum(YReal(end,1:6)-YMod(end,1:6));

% Hill Function
method = 'hillFun';
x0 = [1 0 10 0 0 0 0];
k = [1 1 1 1 1 1 ];
v = [1 1];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k,v),[0 10],x0);
[tMod,~,YMod] = findTC(@hillFun,tReal,'p',[v k],'y0',x0(1:3),'-b');
resid(11) = sum(YReal(end,1)-YMod(end,1));