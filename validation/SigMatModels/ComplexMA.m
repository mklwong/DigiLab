function dx_dt = ComplexMA(t,x,k,v)

%% Definition of states
%x(1) =  I
%x(2) =  A
%x(3) =  pA
%x(4) =  B
%x(5) =  p1B
%x(6) =  p2B
%x(7) =  C
%x(8) =  D
%x(9) =  CD
%x(10) = pCD
%x(11) = A-B
%x(12) = A-p1B
%x(13) = I-CD
%x(14) = I-pCD
%x(15) = pCD-A
%x(16) = pCD-pA
%x(17) = pCD-B
%x(18) = pCD-pB
%x(19) = p1B-B
%x(20) = p1B-p2B

v12 = min(v(1),v(2));
%% Reaction velocities
rxn1_1f = k(1)*x(2)*x(4)*v(1);	%A + B -> A-B
rxn1_2f = k(2)*x(11)*v(1);		%A-B -> A + B
rxn1_3f = k(3)*x(11)*v(1);		%A-B -> A-p1B
rxn1_1r = k(4)*x(2)*x(5)*v(1);	%A + p1B -> A-p1B
rxn1_2r = k(5)*x(12)*v(1);		%A-p1B -> p1B + A
rxn1_3r = k(6)*x(12)*v(1);	    %A-p1B -> B + A

rxn2 = k(7)*x(5)*v(1);          %p1B -> B
rxn3 = k(8)*x(7)*x(8)*v(2);     %C+D -> CD
rxn4 = k(9)*x(9)*v(2);			%CD -> C + D

rxn5_1f = k(10)*x(1)*x(9)*v(2);  %CD + I -> I-CD
rxn5_2f = k(11)*x(13)*v(2);      %I-CD -> CD + I
rxn5_3f = k(12)*x(13)*v(2);      %I-CD -> I-pCD
rxn5_1r = k(13)*x(1)*x(10)*v(2); %pCD + I -> I-pCD
rxn5_2r = k(14)*x(14)*v(2);      %I-pCD -> pCD + I
rxn5_3r = k(15)*x(14)*v(2);      %I-pCD -> I-CD

rxn6 = k(16)*x(10)*v(2);		 %pCD -> CD

rxn7_1f = k(17)*x(2)*x(10)*v12;	%pCD + A -> pCD-A
rxn7_2f = k(18)*x(15)*v12;      %pCD-A -> A + pCD
rxn7_3f = k(19)*x(15)*v12;      %pCD-A -> pCD-pA
rxn7_1r = k(20)*x(3)*x(10)*v12; %pCD + pA -> pCD-pA
rxn7_2r = k(21)*x(16)*v12;		%pCD-pA -> pCD + pA
rxn7_3r = k(22)*x(16)*v12;		%pCD-pA -> pCD-A

rxn8 = k(23)*x(3)*v(1);			%pA -> A

rxn9_1f = k(24)*x(4)*x(10)*v12; %pCD + B -> pCD-B
rxn9_2f = k(25)*x(17)*v12;		%pCD-B -> B + pCD
rxn9_3f = k(26)*x(17)*v12;      %pCD-B -> pCD-pB
rxn9_1r = k(27)*x(6)*x(10)*v12; %pCD + pB -> pCD-pB
rxn9_2r = k(28)*x(18)*v12;		%pCD-pB -> pCD + pB
rxn9_3r = k(29)*x(18)*v12;		%pCD-pB -> pCD-B

rxn10_1f = k(30)*x(5)*x(6)*v(1); %p1B + p2B -> p1B-p2B
rxn10_2f = k(31)*x(19)*v(1);	 %p1B-p2B -> p1B + p2B
rxn10_3f = k(32)*x(19)*v(1);	 %p1B-p2B -> p1B-B
rxn10_1r = k(33)*x(5)*x(4)*v(1); %p1B + B -> p1B-B
rxn10_2r = k(34)*x(20)*v(1);     %p1B-B -> p1B + B
rxn10_3r = k(35)*x(20)*v(1);     %p1B-B -> p1B-p2B

rxn11 = k(36)*x(6)*v(1);		 %p2B -> pB


% Constructing ODEs
dx_dt(1) = (-rxn5_1f +rxn5_2f -rxn5_1r +rxn5_2r)/v(2);
dx_dt(2) = (-rxn1_1f +rxn1_2f -rxn1_1r +rxn1_2r -rxn7_1f +rxn7_2f +rxn8)/v(1);
dx_dt(3) = (-rxn7_1r +rxn7_2r -rxn8)/v(1);

dx_dt(4) = (-rxn1_1f +rxn1_2f  +rxn2 -rxn9_1f +rxn9_2f  -rxn10_1r +rxn10_2r +rxn11)/v(1);
dx_dt(5) = (-rxn1_1r +rxn1_2r  -rxn2 -rxn10_1f +rxn10_2f -rxn10_1r +rxn10_2r)/v(1);
dx_dt(6) = (-rxn9_1r +rxn9_2r -rxn11 -rxn10_1f +rxn10_2f)/v(1);
dx_dt(7) = (-rxn3    +rxn4)/v(2);
dx_dt(8) = (-rxn3    +rxn4)/v(2);
dx_dt(9) = ( rxn3    -rxn4    -rxn5_1f +rxn5_2f +rxn6)/v(2);
dx_dt(10)= (                  -rxn5_1r +rxn5_2r -rxn6  -rxn9_1f +rxn9_2f -rxn9_1r +rxn9_2r  -rxn7_1f + rxn7_2f -rxn7_1r +rxn7_2r)/v(2);

dx_dt(11)=  (rxn1_1f  -rxn1_2f  -rxn1_3f +rxn1_3r)/v(1);
dx_dt(12)=  (rxn1_1r  -rxn1_2r  -rxn1_3r +rxn1_3f)/v(1);
dx_dt(13)=  (rxn5_1f  -rxn5_2f  -rxn5_3f +rxn5_3r)/v(2);
dx_dt(14)=  (rxn5_1r  -rxn5_2r  -rxn5_3r +rxn5_3f)/v(2);
dx_dt(15)=  (rxn7_1f  -rxn7_2f  -rxn7_3f +rxn7_3r)/v12;
dx_dt(16)=  (rxn7_1r  -rxn7_2r  -rxn7_3r +rxn7_3f)/v12;
dx_dt(17)=  (rxn9_1f  -rxn9_2f  -rxn9_3f +rxn9_3r)/v12;
dx_dt(18)=  (rxn9_1r  -rxn9_2r  -rxn9_3r +rxn9_3f)/v12;
dx_dt(19)=  (rxn10_1f  -rxn10_2f  -rxn10_3f +rxn10_3r)/v(1);
dx_dt(20)=  (rxn10_1r  -rxn10_2r  -rxn10_3r +rxn10_3f)/v(1);

dx_dt = dx_dt';