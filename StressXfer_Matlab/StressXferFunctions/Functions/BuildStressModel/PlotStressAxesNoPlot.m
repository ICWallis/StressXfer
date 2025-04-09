function [Sxx3D,Syy3D,Szz3D,Sxy3D,Sxz3D,Syz3D] = PlotStressAxesNoPlot(Sv3D,Shmax3D,Shmin3D,DirS1,InputElevation,FaultFileString)

for j = 1
%%% THIS BLOCK PLOTS A STRESS TENSOR FOR A GIVEN DEPTH (i=2)
%%% CHECK AND MAKE SURE THE TENSOR LOOKS RIGHT
Test=Sv3D(j)>=Shmax3D(j);%Check to see if normal faulting stress regim; zero is false, one is true
Test2=Sv3D(j)<=Shmax3D(j) && Sv3D(j)>=Shmin3D(j);%Check to see if strikeslip faulting stress regime; zero is false, one is true

Sxx3D(1:(length(Sv3D)),1:1)=0;
Syy3D(1:(length(Sv3D)),1:1)=0;
Szz3D(1:(length(Sv3D)),1:1)=0;
Sxy3D(1:(length(Sv3D)),1:1)=0;
Sxz3D(1:(length(Sv3D)),1:1)=0;
Syz3D(1:(length(Sv3D)),1:1)=0;

if Test==1
   for i = 1:(length(Sv3D))
   [Sxx3D(i),Syy3D(i),Szz3D(i),Sxy3D(i),Sxz3D(i),Syz3D(i)]= TensorsFromPrincipal3d_normal(Sv3D(i),Shmax3D(i),Shmin3D(i),DirS1(i));
   end
    disp('Normal Faulting Stress Regime')
else
if Test2==1 && Test ~= 1
   for i = 1:(length(Sv3D))
   [Sxx3D(i),Syy3D(i),Szz3D(i),Sxy3D(i),Sxz3D(i),Syz3D(i)]= TensorsFromPrincipal3d_strikeslip(Shmax3D(i),Sv3D(i),Shmin3D(i),DirS1(i));
   end
    disp('Strike Slip Faulting Stress Regime')
else
    for i = 1:(length(Sv3D))
   [Sxx3D(i),Syy3D(i),Szz3D(i),Sxy3D(i),Sxz3D(i),Syz3D(i)]= TensorsFromPrincipal3d_thrust(Shmax3D(i),Shmin3D(i),Sv3D(i),DirS1(i));
   end
    disp('Reverse Faulting Stress Regime')
end
end

%This plots the Principal direction from the tensors, CHECK TO
%MAKE SURE IT LOOKS RIGHT
%DrawS1S2S3Directions([Sxx3D(j),Syy3D(j),Szz3D(j),Sxy3D(j),Sxz3D(j),Syz3D(j)],0,0,InputElevation,FaultFileString,'Scale', 5 )
end

