function varargout = [propDis](fcn,runVar,opts)

switch lower(deblank(fcn))
    case 'newpt'
        %<Generate pt1 and pdfBias based on a custom proposal distribution>
        varargout = {pt1,pdfBias};
        
    case 'adapt'
        %<Modify runVar.step. Information from opts.stepi is available>
 
        varargout = {runVar};

end

end