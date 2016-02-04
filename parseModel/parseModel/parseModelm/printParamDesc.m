function printParamDesc(model,p)

for ii = 1:length(model.pFit.desc); 
	for jj = 1:size(model.pFit.desc{ii},1)
		if jj == 1 && nargin == 2
			fprintf('%s - %s \n',model.pFit.desc{ii}(jj,:),p(ii)); 
		else
			fprintf('%s \n',model.pFit.desc{ii}(jj,:)); 
		end
	end
end