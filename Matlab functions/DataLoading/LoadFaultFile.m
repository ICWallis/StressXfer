function [P1,P2,P3,FaceNormalVector,X_fault,Y_fault,Z_fault,FaultMidPoint,FaultTriangles,FaultPoints] = LoadFaultFile(Faultfile)
%This function loads the fault file and parses it into the requistite
%data
    [ FaultPoints,FaultTriangles ] = GoCadAsciiReader( Faultfile );
    FaultPoints=[FaultPoints(:,1),FaultPoints(:,2),FaultPoints(:,3),(FaultPoints(:,4))];
  %%%PLOTS THE FAULT (UNSHIFTED)  
    [FaultMidPoint,FaceNormalVector] = MidPointCreate(FaultPoints,FaultTriangles);
    X_fault=FaultMidPoint(:,1);
    Y_fault=FaultMidPoint(:,2);
    Z_fault=FaultMidPoint(:,3);
    
    [P1,P2,P3] = CreateP1P2P3( FaultTriangles,FaultPoints );        
end

