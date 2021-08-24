
load('toyDataThreeWay.mat')
out = DJIVEMainJP(datablock) ;
DJIVEAngleDiagnosticJP(datablock, dataname, out, 556, "Demo")