residAggregate = nan(1000,12);
for ii = 1:1000;
	testCases;
	residAggregate(ii,:) = resid;
	fprintf('.')
	if mod(ii,20) == 0
		fprintf('|%d\n',ii)
	end
end