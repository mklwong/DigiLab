function [X,Y,Z] = surfFuncTest(obj,bnd)

[X,Y] = meshgrid(linspace(bnd(1,1),bnd(1,2),250),linspace(bnd(2,1),bnd(2,2),250));

Z = obj(X,Y);

