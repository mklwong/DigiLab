function dx_dt = testReactions(t,x,type,k,v)

switch type
    case 'syn'
        % Synthesis k = [1]
        dx_dt(1) = k(1);
    case 'uni'
        % Dissociation: k = [1 0 0 0]
        % Interconversion k = [0 1 0 0]
        % Degradation:  k = [0 0 1 0]
		% SynMM: k = [0 0 0 1]
        dx_dt(1) = (-(k(1)+k(2)+k(3))*x(1)+k(4)*x(2))*v(1)/v(1);
        dx_dt(2) =  (k(1)+k(2))*x(1)*v(1)/v(2);
        dx_dt(3) =  k(1)*x(1)*v(1)/v(1);
    case 'bi'
        % Association:             k = [1 0 0]
        % MM Enzyme Kinetic:       k = [0 1 0]
        % MM enzyme Kinetic Deg:   k = [0 0 1]
		vOverlap = min(v);
        dx_dt(1) = -k(1)*x(1)*x(2)*vOverlap/v(1) - k(2)*x(1)*x(2)*vOverlap/v(1) - k(3)*x(1)*x(2)*vOverlap/v(1);
        dx_dt(2) = -k(1)*x(1)*x(2)*vOverlap/v(2);
        dx_dt(3) =  k(1)*x(1)*x(2)*vOverlap/v(2) + k(2)*x(1)*x(2)*vOverlap/v(2);
    case 'enzQSSA'
        %k(1-3) is the forward. k(4-6) is the reverse.
		v12 = min(v(1),v(2));
		v34 = min(v(3),v(4));
        dx_dt(1) =  (-k(1)*x(1)*x(2)*v12 +k(2)*x(5)*v12 +k(6)*x(6)*v34      +k(7)*x(3)*v(3))/v(1);
        dx_dt(2) =  (-k(1)*x(1)*x(2)*v12 +k(2)*x(5)*v12 +k(3)*x(5)*v12                      )/v(2);
        
        dx_dt(3) =  (-k(4)*x(3)*x(4)*v34 +k(5)*x(6)*v34 +k(3)*x(5)*v12      -k(7)*x(3)*v(3))/v(3);
        dx_dt(4) =  (-k(4)*x(3)*x(4)*v34 +k(5)*x(6)*v34 +k(6)*x(6)*v34                      )/v(4);
		
		dx_dt(5) =  k(1)*x(1)*x(2)-k(2)*x(5)-k(3)*x(5);
        dx_dt(6) =  k(4)*x(3)*x(4)-k(5)*x(6)-k(6)*x(6);
		
	case 'hillFun'
		dx_dt(1) = -k(1)*x(1)*x(3).^k(3)/(k(2).^k(3)+x(3).^k(3))+k(4)*x(2)*v(2)/v(1);
        dx_dt(2) =  k(1)*x(1)*x(3).^k(3)/(k(2).^k(3)+x(3).^k(3))*v(1)/v(2)-k(4)*x(2);
		dx_dt(3) = 0;
end

dx_dt(7) = 0;
dx_dt = dx_dt';