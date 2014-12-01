function [lb ub] = reBound(res,varargin)

if ~isempty(intersect({'bounds'},fieldnames(res)))
	lbCur = res.bounds{1};
	ubCur = res.bounds{2};
	scl = res.bounds{3};
else
	model = varargin{1};
	[~,~,bnds] = model('info');
	lbCur = bnds{1};
	ubCur = bnds{2};
	scl   = bnds{3};
end

midPt = res.pstrir(res.P==min(res.P),:);

logRng = log10(ubCur)-log10(lbCur);
lbLog = 10.^(round(log10(midPt)*2-logRng)/2);
ubLog = 10.^(round(log10(midPt)*2+logRng)/2);

linRng = ubCur-lbCur;
lbLin = midPt+linRng/2;
ubLin = midPt-linRng/2;
ubLin = max([ubLin;zeros(size(ubLin))]); %make sure lower bound non-zero

lb = lbLin;
lb(scl) = lbLog(scl);
ub = ubLin;
ub(scl) = ubLog(scl);

fprintf('lbOther = [');  fprintf('%3.2e ',lb); fprintf('];\r');
fprintf('ubOther = [');  fprintf('%3.2e ',ub); fprintf('];\r');
fprintf('sclOther = ['); fprintf('%3.2e ',scl); fprintf('];\r');