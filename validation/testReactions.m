function dx_dt = testReactions(t,x,type,k)

switch type
    case 'syn'
        % Synthesis k = [1]
        dx_dt(1) = k(1);
    case 'uni'
        % Dissociation: k = [1 1 1 0]
        % Interconversion k = [1 0 1 0]
        % Degradation:  k = [1 1 0 0]
        % Enzyme kinetic Syn (MM): k = [1 1 1 1]
        dx_dt(1) = -k(1)*x(1) + k(1)*k(4)*x(1);
        dx_dt(2) =  k(1)*k(3)*x(1);
        dx_dt(3) =  k(1)*k(2)*k(3)*x(1);
    case 'bi'
        % Association:             k = [1 1 0]
        % MM Enzyme Kinetic:       k = [1 1 1]
        % MM enzyme Kinetic Deg:   k = [1 0 1]
        dx_dt(1) = -k(1)*x(1)*x(2);
        dx_dt(2) = -k(1)*x(1)*x(2) + k(1)*k(3)*x(1)*x(2);
        dx_dt(3) =  k(1)*k(2)*x(1)*x(2);
    case 'enzQSSA'
        %k(1-3) is the forward. k(4-6) is the reverse.
        dx_dt(1) = -k(1)*x(1)*x(2)+k(2)*x(3)+k(4)*x(6);
        dx_dt(2) = -k(1)*x(1)*x(2)+k(2)*x(3)+k(3)*x(3);
        dx_dt(3) =  k(1)*x(1)*x(2)-k(2)*x(3)-k(3)*x(3);
        dx_dt(4) = -k(4)*x(4)*x(5)+k(5)*x(6)+k(3)*x(3);
        dx_dt(5) = -k(4)*x(4)*x(5)+k(5)*x(6)+k(6)*x(6);
        dx_dt(6) =  k(4)*x(4)*x(5)-k(5)*x(6)-k(6)*x(6);
        
end
dx_dt(7) = 0;
dx_dt = dx_dt';