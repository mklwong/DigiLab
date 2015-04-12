%%%testCase%%%

% This program tests whether the odeQSSA program is working correctly.

%% Synthesis
method = 'syn';

x0 = [0 0 0 0 0 0 0];
k = [1];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k),[0 10],x0);
[~,~,YMod] = findTC(@SynProg,tReal,'p',1,'-b');

resid(1) = sum(YReal(:,1)-YMod(:,1));

%% Uni
method = 'uni';
%Dissociation
x0 = [1 0 0 0 0 0 0];
k = [1 1 1 0];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k),[0 10],x0);
[tMod,~,YMod] = findTC(@UniDis,tReal,'p',1,'-b');
resid(2) = sum(YReal(:,1)-YMod(:,1));

%Interconversion
x0 = [1 0 0 0 0 0 0];
k = [1 0 1 0];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k),[0 10],x0);
[tMod,~,YMod] = findTC(@UniCon,tReal,'p',1,'-b');
resid(3) = sum(YReal(:,1)-YMod(:,1));

%Degradation
x0 = [1 0 0 0 0 0 0];
k = [1 1 0 0];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k),[0 10],x0);
[tMod,~,YMod] = findTC(@UniDeg,tReal,'p',1,'-b');
resid(4) = sum(YReal(:,1)-YMod(:,1));

%SynMM
x0 = [1 0 0 0 0 0 0];
k = [1 1 1 1];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k),[0 10],x0);
[tMod,~,YMod] = findTC(@SynMM,tReal,'p',1,'-b');
resid(5) = sum(YReal(:,2)-YMod(:,1));

%% Bi
method = 'bi';
% Association
x0 = [1 1 0 0 0 0 0];
k = [1 1 0];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k),[0 10],x0);
[tMod,~,YMod] = findTC(@BiAss,tReal,'p',1,'-b');
resid(6) = sum(YReal(:,1)-YMod(:,1));

% MM Enzyme Kinetics
x0 = [1 1 0 0 0 0 0];
k = [1 1 1];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k),[0 10],x0);
[tMod,~,YMod] = findTC(@BiMM,tReal,'p',1,'-b');
resid(7) = sum(YReal(:,1)-YMod(:,1));

% MM Enzyme Kinetics Deg
x0 = [1 1 0 0 0 0 0];
k = [1 0 1];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k),[0 10],x0);
[tMod,~,YMod] = findTC(@BiMMDeg,tReal,'p',1,'-b');
resid(8) = sum(YReal(:,1)-YMod(:,1));

%% Enzyme Kinetic
method = 'enzQSSA';
% Single
x0 = [1 1 0 0 1 0 0];
k = [1 1 1 0 1 1];
kQSSA = [1 2 1 Inf];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k),[0 30],x0);
[tMod,~,YMod] = findTC(@enzKinetic,tReal,'p',kQSSA,'-b');
resid(9) = sum(YReal(:,1)-YMod(:,1));

% Reversible
x0 = [1 1 0 0 1 0 0];
k = [1 1 1 1 1 1 ];
kQSSA = [1 2 1 2];
[tReal,YReal] = ode15s(@(t,x) testReactions(t,x,method,k),[0 10],x0);
[tMod,~,YMod] = findTC(@enzKinetic,tReal,'p',kQSSA,'-b');
resid(10) = sum(YReal(:,1)-YMod(:,1));