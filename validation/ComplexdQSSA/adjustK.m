function k = adjustK(k)

% Make the catalytic rate larger in the favoured direction
reorder = @(x) [max(x) min(x)]; 
k([3 6]) = reorder(k([3 6]));     %Rxn 1
k([12 15]) = reorder(k([12 15])); %Rxn 5
k([19 22]) = reorder(k([19 22])); %Rxn 7
k([26 29]) = reorder(k([26 29])); %Rxn 9
k([32 35]) = reorder(k([32 35])); %Rxn 10

% Enforce rapid equilibrium assumption
quikDis = @(x) [x(1)*100 max(x([3 6]))*100 x(3) x(4)*100 max(x([3 6]))*100 x(6)];
k(1:6) = quikDis(k(1:6));
k(10:15) = quikDis(k(10:15));
k(17:22) = quikDis(k(17:22));
k(24:29) = quikDis(k(24:29));
k(30:35) = quikDis(k(30:35));