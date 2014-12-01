function checkSum = sbmlCheckSum(fname)

file = rmAnnot(fname);
checkSum = md5(file);
delete(file);