function varargout = parComp(mode,numWorkers)

%parComp Cross version toggle of parallel computer pools
% this function can open and close a pool of workers for parallel computing
% purposes. It checks matlab version and availability of the parallel
% computing toolbox to determine if parallel computing is possible and
% what function (parpool vs matlabpool) to use.
%
% status = parComp('open',numWorkers) attempts to open a pool. If the 
% parallel computing toolbar isn't installed, parComp returns -1. Else
% parComp attempts to open the number of workers specified by numWorkers. 
% parComp will only open up to the maximum number of workers available in 
% the pool.
%
% status = parComp('close') attempts to close an open pool.

modeOpts = ['open '
            'close'];
        
if nargin == 1
    numWorkers = Inf;
end

switch lower(deblank(mode))
    case lower(deblank(modeOpts(1,:))) %Opening cluster
        a = ver; %get matlab toolbox availability
        parBoxExist =  ~isempty(intersect({a.Name},'Parallel Computing Toolbox'));
        if parBoxExist
            matVer = version('-release');
	    matSubVer = matVer(5);
            matVer = str2double(matVer(1:4));
            if matVer>=2014 || (matVer==2013 && strcmpi(matSubVer,'b'))
                clustInfo = parcluster('local');
                maxWorkers = clustInfo.NumWorkers;
                numWorkers = min([maxWorkers numWorkers]); %Make sure required workers is not more than max workers
                parObj = gcp('nocreate');
                if isempty(parObj)
                    parpool(numWorkers);
                else
                    fprintf('Cluster already open\n')
                end
                varargout{1} = 0;
			else
                clustInfo = findResource;
                maxWorkers = clustInfo.ClusterSize;
                numWorkers = min([maxWorkers numWorkers]); %Make sure required workers is not more than max workers
                if matlabpool('size')==0
                    matlabpool('open');
                else
                    fprintf('Cluster already open\n')
                end
                varargout{1} = 0;
            end
        else %No cluster
            warning('parComp:noParallelToolbox','This machine does not support parallel computing. Reverting to single core mode')
            varargout{1} = -1;
        end
    case lower(deblank(modeOpts(2,:))) %Closing cluster
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
