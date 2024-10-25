function [Sv3D,Shmax3D,Shmin3D,Sv3D_Obs,Shmax3D_Obs,Shmin3D_Obs,Pp,Pp3D] = Plot3DStressModel(InputElevation,Z_Obs,Sv1,Shmax1,Shmin1,Pp1,Z_Fault,MeanSurfaceElevation,FaultFileString)

%%%InputElevation           The elevation that the inital Sv (vertical
%                           stress) is calculated for in meters.

%%%Z_Obs                    The Z values of the obeservation points. This
%                           value comes out of the 'LoadData.m' function

%%%Sv1                      The vertical stress at the input elevation.
%                           This value comes out of the 'aphi_to_stressmags'
%                           function. Value in MPa

%%%Shmax1                   The maximum horiztontal stress at the input
%                           elevation. This value comes out of the
%                           'aphi_to_stressmags' function. Value in MPa

%%%Shmin1                   The minimum horiztontal stress at the input
%                           elevation. This value comes out of the
%                           'aphi_to_stressmags' function. Value in MPa

%%%Pp1                      The minimum horiztontal stress at the input
%                           elevation. This value comes out of the
%                           'BuildStressModel' function. Value in MPa

%%%Z_Fault                  The Z values of the fault facet midpoints. This
%                           value comes out of the 'LoadData.m' function

%%%MeanSurfaceElevation    The mean surface elevation in meters

%%%FaultFileString      A string used to append to the front of the file
%%%                     names of saved figures and exported data files

%%%SHIFTS THE FAULT TO DEPTH. The 'Adjust' sets the zero stress and Pp 
%%%to 25 m above the surface so there are no places on the faults with
%%%stess of zero, which screws up the calculation.
    Pp_depth = 0;
    Adjust = 25;
    ShiftAmount = MeanSurfaceElevation-InputElevation-Adjust;
%%%THIS BLOCK BUILDS A LINEARLY INCREASING WITH DEPTH STRESS MODEL BASED ON
%%%THE INPUT STRESSES IN THE FIRST BLOCK AND THE Z VALUES FROM THE FAULT
Z2 = (Z_Fault(:,1)-MeanSurfaceElevation-Adjust);
Sv3D(1:(length(Z_Fault(:,1))),1:1) = 0;
Shmax3D(1:(length(Z_Fault(:,1))),1:1) = 0;
Shmin3D(1:(length(Z_Fault(:,1))),1:1) = 0;
Pp3D(1:(length(Z_Fault(:,1))),1:1) = 0;

%%THE INPUT DEPTH IS THE DEPTH AT WHICH THE INPUT STRESSES ARE calculated

for i = 1:(length(Z_Fault(:,1)))
   Sv3D(i) = Sv1*(Z2(i)./(-1*(ShiftAmount+Adjust+Adjust)));
   Shmin3D(i) = Shmin1*(Z2(i)./(-1*(ShiftAmount+Adjust+Adjust)));
   Shmax3D(i) = Shmax1*(Z2(i)./(-1*(ShiftAmount+Adjust+Adjust)));
   Pp3D(i) = Pp1*(Z2(i)./(-1*(ShiftAmount+Adjust+Adjust)));
end

for j = 1:length(Z_Fault)
Test=Z_Fault(j) > MeanSurfaceElevation;
if Test == 0  
Sv3D(j) = Sv3D(j);
else
Sv3D(j) = 0;
end
end

for j = 1:length(Z_Fault)
Test=Z_Fault(j) > MeanSurfaceElevation;
if Test == 0  
Shmax3D(j) = Shmax3D(j);
else
Shmax3D(j) = 0;
end
end

for j = 1:length(Z_Fault)
Test=Z_Fault(j) > MeanSurfaceElevation;
if Test == 0  
Shmin3D(j) = Shmin3D(j);
else
Shmin3D(j) = 0;
end
end

Pp(1:(length(Z_Obs)),1:1) = 0;
Z3 = (Z_Obs-MeanSurfaceElevation-Adjust);
for i = 1:(length(Z3))
Pp(i) = Pp1*(Z3(i)/(-1*(ShiftAmount+Adjust+Adjust)));
end

lt500= lt(Pp_depth,abs(Z3));
Pp = Pp.*lt500;
dum = (min(Pp(Pp>0)));

for j = 1:length(Pp3D)
Test=Pp3D(j) == 0;
if Test == 0  
Pp3D(j) = Pp3D(j);
else
Pp3D(j) = 0;
end
end

for j = 1:length(Z_Fault)
Test=Z_Fault(j) > MeanSurfaceElevation;
if Test == 0  
Pp3D(j) = Pp3D(j);
else
Pp3D(j) = 0;
end
end

Z3 = (Z_Obs-MeanSurfaceElevation-Adjust);
Shmax3D_Obs(1:(length(Z_Obs)),1:1) = 0;
Shmin3D_Obs(1:(length(Z_Obs)),1:1) = 0;
Sv3D_Obs(1:(length(Z_Obs)),1:1) = 0;
Pp(1:(length(Z_Obs)),1:1) = 0;

for i = 1:(length(Z3))
   Sv3D_Obs(i) = Sv1*(Z3(i)/(-1*(ShiftAmount+Adjust+Adjust)));
   Shmin3D_Obs(i) = Shmin1*(Z3(i)/(-1*(ShiftAmount+Adjust+Adjust)));
   Shmax3D_Obs(i) = Shmax1*(Z3(i)/(-1*(ShiftAmount+Adjust+Adjust)));
   Pp(i) = Pp1*(Z3(i)/(-1*(ShiftAmount+Adjust+Adjust)));
end

for j = 1:length(Z_Obs)
Test=Z_Obs(j) > MeanSurfaceElevation;
if Test == 0  
Sv3D_Obs(j) = Sv3D_Obs(j);
else
Sv3D_Obs(j) = 0;
end
end

for j = 1:length(Z_Obs)
Test=Z_Obs(j) > MeanSurfaceElevation;
if Test == 0  
Shmax3D_Obs(j) = Shmax3D_Obs(j);
else
Shmax3D_Obs(j) = 0;
end
end

for j = 1:length(Z_Obs)
Test=Z_Obs(j) > MeanSurfaceElevation;
if Test == 0  
Shmin3D_Obs(j) = Shmin3D_Obs(j);
else
Shmin3D_Obs(j) = 0;
end
end



lt500= lt(Pp_depth,abs(Z3));
Pp = Pp.*lt500;
dum = (min(Pp(Pp>0)));

for j = 1:length(Pp)
Test=Pp(j) == 0;
if Test == 0  
Pp(j) = Pp(j);
else
Pp(j) = 0;
end
end

for j = 1:length(Pp)
Test=Z_Obs(j) > MeanSurfaceElevation;
if Test == 0  
Pp(j) = Pp(j);
else
Pp(j) = 0;
end
end
% This code plots and saves the vertical stress model. The plot shows both
% the input stress on the fault and the observations point. They should be
% the same. The figure is saved in the OutputFigures folder
VerticalStressModel = figure;
scatter(Sv3D,Z_Fault,'r');hold on; scatter(Shmax3D,Z_Fault,'g');hold on;
scatter(Shmin3D,Z_Fault,'b');hold on;scatter(Pp3D,Z_Fault,'k');hold off;
ax.YDir = 'reverse';
hold on
scatter(Sv3D_Obs,Z_Obs,200,'r','filled','square');hold on; 
scatter(Shmax3D_Obs,Z_Obs,200,'g','filled','square');hold on; 
scatter(Shmin3D_Obs,Z_Obs,200,'b','filled','square');hold on; 
scatter(Pp,Z_Obs,150,'k','filled','square');hold off
legend ('Sv fault', 'Shmax fault','Shmin fault','Sv obs',...
    'Shmax obs','Shmin obs', 'pore pressure');

hold on
scatter(Sv1,InputElevation,25,'k','d','filled');hold on; 
scatter(Shmax1,InputElevation,25,'k','d','filled');hold on; 
scatter(Shmin1,InputElevation,25,'k','d','filled');
hold on

Surface(1:100)= MeanSurfaceElevation;
SurfaceX = linspace(0,max(Sv3D));
scatter(SurfaceX,Surface,3,'c','o','filled')
legend ('Sv fault', 'Shmax fault','Shmin fault', 'PorePressure fault',...
    'Sv obs', 'Shmax obs','Shmin obs', 'pore pressure','Input Sv',...
    'Input Smax','Input Smin','Surface Elevation');
title ('Vertical Stress Model');
xlabel('Stress (MPa)');
ylabel('Elevation (m)');

filename = strcat(FaultFileString,'VerticalStressModel');
filename2 = strcat('OutputFigures/',filename);
saveas(VerticalStressModel,filename2, 'png')
end

