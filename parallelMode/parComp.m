function varargout = parComp(mode,numWorkers)

modeOpts = ['open '
            'close'];
        
if nargin == 1
    numWorkers = Inf;
end

switch lower(deblank(mode))
    case lower(deblank(modeOpts(1,:)))
        if exist('parpool','file') || exist('matlabpool','file')
            clustInfo = parcluster('local');
            maxWorkers = clustInfo.NumWorkers;
            numWorkers = min([maxWorkers numWorkers]);
            if exist('parpool','file')
                parObj = gcp('nocreate');
                if isempty(parObj)
                    parpool(numWorkers);
                end
            elseif exist('matlabpool','file')
                matlabpool('open',numWorkers);
			end
			varargout{1} = 0;
        else
            warning('parComp:noParallelToolbox','This machine does not support parallel computing. Reverting to single core mode')
            varargout{1} = -1;
        end

    case lower(deblank(modeOpts(2,:)))
        if exist('parpool','file')
            parObj = gcp('nocreate');
            delete(parObj)
        elseif exist('matlabpool','file')
            matlabpool('close')
        end
        varargout{1} = 0;
	otherwise
		error('parComp:incorrectMode','Unknown mode selected. Must either be open or close')
end