residAggregate = nan(100,12);
for ii = 1:100;
	testCases;
	residAggregate(ii,:) = resid;
	fprintf('.')
	if mod(ii,20) == 0
		fprintf('|%d\n',ii)
	end
end