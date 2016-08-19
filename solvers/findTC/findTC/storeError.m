function storeError(model,x0,p,matmsg,msg,Errdir)

% storeError(model,p,msg)
%
% This function stores some information for debugging purposes when an
% error is generated. The stores:
% - The model (if custom, will fall under custom).
% - The parameter set attempted
% - The error generated
% - The time of generation

%% Extract model name
if nargin == 6
	dirLoc = Errdir;
else

if isfield(model,'name')
	modName = model.name;
	dirLoc = which(modName);
	dirLoc(end-1:end) = [];
elseif ischar(model)
	dirLoc = ['./' model];
else
	dirLoc = './custom';
end
end
dirLoc = [dirLoc '-Errors'];

%% Check if directory for storage exists
if ~exist(dirLoc,'dir')
    mkdir(dirLoc);
end

if exist('labindex');
	labVal = labindex;
else
	labVal = 1;
end

%% Generate error identity code
errorID = floor(rand(1)*1000000);
errorID = ['e' num2str(errorID)];

t = clock();

% %% Print the output into the error report text file
% h = fopen([dirLoc '/runError.txt'],'a');
% fprintf(h,'\r\n%s: %s\r\n',modName,msg);
% fprintf(h,'Error found at %2.0f-%2.0f-%4.0f, %2.0f:%2.0f\r\n',t(3),t(2),t(1),t(4:5));
% 
% %Print parameter
% if isempty(p)
% 	fprintf(h,'No p vector given. \n\r');
% elseif isnumeric(p)
%     fprintf(h,'p = [%4.4e',p(1));
%     fprintf(h,',%4.4e',p(2:end));
%     fprintf(h,']\r\n');
% else
%     fprintf(h,'p is the wrong format. \n\r');
% end
% %Print initial condition
% if isempty(x0)
% 	fprintf(h,'no x0 vector given. \n\r');
% elseif isnumeric(x0)
% fprintf(h,'x0 = [%4.4e',x0(1));
% fprintf(h,',%4.4e',x0(2:end));
% fprintf(h,']\r\n');
% else
% 	fprintf(h,'x0 is the wrong format. \n\r');
% end
% 
% fprintf(h,'Files saved with filename %s.mat',errorID);
% fprintf(h,'\r\n-----------------------------------------');
% fclose(h);

%% Save state into error file
if exist([dirLoc '/errors-lab-' labVal '.mat'],'file')
    load([dirLoc '/errors-lab-' labVal '.mat']);
else
	errs = struct([]);
end
if isempty(errs)
	nerrs = 1;
else
	nerrs = length(errs)+1;
end
errs(nerrs).p   = p;
errs(nerrs).x0    = x0;
errs(nerrs).model = model;
errs(nerrs).errID = errorID;
errs(nerrs).msg   = matmsg;
errs(nerrs).time  = [num2str(t(3),'%2.0f') '-' num2str(t(2),'%2.0f') '-' num2str(t(1),'%4.0f') ', ' num2str(t(4),'%2.0f') ':' num2str(t(5),'%2.0f')];

save([dirLoc '/errors-lab-' labVal '.mat'],'errs')
