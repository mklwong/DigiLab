function [filename pathname namepath] = rmAnnot(fname)
%
% rmAnnot(sbmlModel)
%
%	This function removes all annotations from an SBML file. This is so any
%	parser will only need to work on the bare mimimum of information, and
%	checksum generators will only detect change in the actual guts of the
%	model, rather than changing the checksum after every little fit of
%	annotation change.

%Annotation "openers" to skip and their "closers"
toSkipHead = {'<annotation>','</annotation>';
	          '<notes>','</notes>'};

%Count number of lines in file for initialising cell size
fname = strrep(fname,'\','/');
sIndx = strfind(fname,'/');
if isempty(sIndx)
	pathname = [];
	filename = [fname(1:end-4) '_tmp.xml'];
else
	sIndx = sIndx(end);
	pathname = fname(1:(sIndx));
	filename = [fname((sIndx+1):end-4) '_tmp.xml'];
end
namepath = [pathname filename];
fid = fopen(fname, 'rt');
assert(fid ~= -1, [mfilename ':cannotReadFile'],['Could not read: ' fname]);
fileBin = fread( fid);  % read file in binary format
fclose(fid);
crLoc = find(fileBin==10);
crLoc = [0;crLoc];
fid = fopen(filename,'w');
annot = false;
for ii = 2:length(crLoc)
	lineReconst = char(fileBin((crLoc(ii-1)+1):crLoc(ii)-1))';
	toSkip = find(cellfun(@(x)strcmpi(lineReconst,x),toSkipHead(:,1)));
	if ~isempty(toSkip) && ~annot %check for annotation beginning
		annot = true;
		skipAnot = toSkip;
	elseif annot     %check for corresponding annotation ending if in annotation block
		if strcmpi(lineReconst,toSkipHead(skipAnot,2))
			annot = false;
		end
	elseif ~isempty(char(lineReconst)) %Not annotation and not empty: store.
		fprintf(fid,[char(lineReconst)]);
	end
end
fclose(fid);

end