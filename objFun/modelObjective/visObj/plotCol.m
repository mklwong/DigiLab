function rgbOut = plotCol(n,ntot,scheme)

Names(1,:) = 'rainbow';

switch scheme
    case lower(deblank(Names(1,:)))
          %initial, 
        col_i = [1 0 0];
        r     = [0    -1 0    0  ];
        g     = [0.75  0 0 -0.75];
        b     = [0     0 1    0  ];
end

Col = [0 r;0 g;0 b];

if ntot ~= 1
	%Space in colour space n goes
	nCol = (n-1)/(ntot-1)*length(r);

	rgbOut = (col_i'+sum(Col(:,1:(floor(nCol)+1)),2)+mod(nCol,1)*Col(:,(ceil(nCol)+1)))';
else
	rgbOut = Col(:,2);
end