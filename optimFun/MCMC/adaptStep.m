function opts = adaptStep(acptCnt,p0,p1,opts)

if isfield(opts,'step')
	if isrow(opts.step)
		opts.step = opts.step';
	end
end

% Initial step vector initiation
if ischar(acptCnt)
	if strcmpi(acptCnt,'initial')
		opts.step = opts.stepi+p0*0;
	end
	return
end

% Check step ratio
rjtRto = sum(acptCnt)/length(acptCnt);


if rjtRto <= opts.rjtRto
    opts.step = max([1e-4+0*opts.step opts.step/1.1],[],2);
elseif rjtRto > opts.rjtRto
    opts.step = min([opts.maxStep+0*opts.step opts.step*1.1],[],2);
end