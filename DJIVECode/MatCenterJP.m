function outMat = MatCenterJP(X, iColCent, iRowCent)

[d,n] = size(X) ;
outMat = X ;
if iColCent == 1
    outMat = outMat - repmat(mean(outMat,2),1,n) ;
end
if iRowCent == 1
    outMat = outMat - repmat(mean(outMat,1),d,1) ;
end

end