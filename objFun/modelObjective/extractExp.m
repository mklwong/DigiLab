function [YExpMean YExpStd] = extractExp(expIndx,expMean,expStd)

YExpMean = expMean(:,expIndx);
YExpStd = expStd(:,expIndx);