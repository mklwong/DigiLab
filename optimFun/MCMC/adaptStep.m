function opts = adaptStep(acptCnt,runVar,opts)

if isfield(opts,'step')
	if isrow(opts.step)
		opts.step = opts.step';
	end
end

% Initial step vector initiation
if ischar(acptCnt)
	if strcmpi(acptCnt,'initial')
		opts.step  = opts.stepi;
		opts.basis = [1;zeros(runVar.p0-1,1)];
		opts.basisSkew = 0;
	end
	return
end

% Check step ratio
rjtRto = sum(acptCnt)/length(acptCnt);

% Adapt step
%% Normal undirected adaptation
if rjtRto <= opts.rjtRto  %Success
    opts.step = max([1e-2+0*opts.step opts.step/1.1],[],2);
elseif rjtRto > opts.rjtRto % Fail
    opts.step = min([opts.maxStep+0*opts.step opts.step*1.1],[],2);
end

%% Test directed adaptation
% curBasis = opts.basis(:,1);
% testBasis = (runVar.dp')/norm(runVar.dp',2);
% amount = dot(testBasis,curBasis);
% testBasis = sign(amount)*testBasis;
% if rjtRto <= opts.rjtRto    %Fail
% 	opts.basisSkew = max([0 opts.basisSkew-abs(amount)]);
%     opts.step      = max([1e-2+0*opts.step opts.step/1.01],[],2);
% elseif rjtRto > opts.rjtRto %Success
% 	newBasis = curBasis+(opts.step(1)/norm(runVar.dp',2))/9*testBasis;
% % 	if abs(amount) < 0.2
% % 		keyboard
% % 	end
% 	opts.basis     = newBasis/norm(newBasis,2);
% 	if isnan(norm(newBasis,2))
% 		keyboard
% 	end
% 	opts.basisSkew = opts.basisSkew+abs(amount);
%     opts.step = min([opts.maxStep+0*opts.step opts.step*(1+0.1*abs(amount))],[],2);
% end