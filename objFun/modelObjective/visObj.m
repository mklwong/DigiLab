function visObj(xSim,YSim,xExp,YExpMean,YExpStd,normNow,type,matches,input)

% Get number of x steps and states
[nx,ny] = size(YExpMean);

% Rescale YExp to Simulation units
YExpMean = YExpMean.*(ones(nx,1)*normNow);
YExpStd = YExpStd.*(ones(nx,1)*normNow);

if isa(input{2},'function_handle')
	input{2} = 'Custom Function ';
end

figure
clf
hold on
% Cycle through each state assigning a different colour
for ii = 1:ny
    if strcmpi(type,'time')
        h(ii) = plot(xSim,YSim(:,ii),'color',plotCol(ii,ny,'rainbow'));
        plot(xExp,YExpMean(:,ii),'LineStyle','none','MarkerEdgeColor',plotCol(ii,ny,'rainbow'),'Marker','o')
        plot(xExp,YExpMean(:,ii)-YExpStd(:,ii),'LineStyle','none','MarkerEdgeColor',plotCol(ii,ny,'rainbow'),'Marker','^')
        plot(xExp,YExpMean(:,ii)+YExpStd(:,ii),'LineStyle','none','MarkerEdgeColor',plotCol(ii,ny,'rainbow'),'Marker','v')
        title(['Dose = ' num2str(input{2}) '{\mu}M'])
        %set(gca,'YScale','log')
    elseif strcmpi(type,'dose')
        h(ii) = semilogx(xSim,YSim(:,ii),'color',plotCol(ii,ny,'rainbow'));
		semilogx(xExp,YExpMean(:,ii),'o','LineStyle','none','MarkerEdgeColor',plotCol(ii,ny,'rainbow'))
		semilogx(xExp,YExpMean(:,ii)-YExpStd(:,ii),'^','LineStyle','none','MarkerEdgeColor',plotCol(ii,ny,'rainbow'))
		semilogx(xExp,YExpMean(:,ii)+YExpStd(:,ii),'v','LineStyle','none','MarkerEdgeColor',plotCol(ii,ny,'rainbow'));
		set(gca,'XScale','log')
    end
end
legend(h,matches)
hold off