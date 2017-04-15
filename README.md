DigiLab - The Digital Cellular Laboratory in MATLAB

=========

The Digital Cellular Laboratory in MATLAB is a MATLAB package that is designed to 
easily model large biochemical networks with complex rate equations. It achieves this
by allowing modularisation of biochemical networks and abstracting away any modelling 
component that is more complex that designing the network connections and state of the
network (i.e. parameter and initial condition selection). However, the deeper mathematics
can be edited by changing the underlying default rate equations used by the package.

To use this package, either download the package, unzip the contents into
any folder. Then add all contents of the folder into the MATLAB search
path. Alternately you can create a fork of the repository (following [this
guide](https://help.github.com/articles/fork-a-repo/) ), and then adding the
path of your local fork to the MATLAB search path. The functions within the
package should then work.

The basic usage steps are:

1. modelName = mkModel('filename') - this creates a barebones sigMat model to build from.
2. [t,Y] = findTC(modelName,t) - this simulates the model name.
3. posterior = MCMC(objFun,pt0,bnd) - this performs parameter fitting on 
                                         the objective function objFun.

Type '>> help [function]' for more information about usage of these functions.

You can find a more comprehensive walkthrough of this package in the 
documentation for this package in documentation.pdf

Licensing
=========

This work is copyright. Apart from any use permitted under the 
Copyright Act 1968, no part may be reproduced by any process, 
nor may any other exclusive right be exercised, without the permission
of Martin Wong, Sydney University, 2016.



This program is free software: you can redistribute it and/or modify 
it under the terms of the GNU General Public License as published by 
the Free Software Foundation, either version 3 of the License, or 
(at our option) any later version.





This program is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.



You should have received a copy of the GNU General Public License 
along with this program.  If not, see <http://www.gnu.org/licenses/>.
