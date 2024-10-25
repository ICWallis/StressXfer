function [Dss,Dds,Dn,lambda,Mu,Cohesion,halfspace] = CalcSlip(FaultTriangles,FaultPoints,FaultMidPoint,P1,P2,P3,...
    FaceNormalVector,Sxx_Fault,Syy_Fault,Szz_Fault,Sxy_Fault,Sxz_Fault,Syz_Fault,...
    Z_Obs,Myu,Co,mu,nu,path,FaultFileString)

%This function does the slip calculation needed for the stress transfer
%modeling function 'RunAll.m'

%%%INPUTS%%%

%%%FaultTriangles           A three column vector containing the vertex
%%%                         number for each fault facet. This vector comes
%%%                         out of the 'LoadFaultFile' function.

%%%FaultPoints              A four column vector containing the vertices
%%%                         number and x,y,z coroodinates for each fault
%%%                         facet vertex. This vector comes out of the
%%%                         'LoadFaultFile' function.

%%%FaultMidPoint            A three column vector containing the x,y,z
%%%                         coordinate at the center of each fault facet
%%%                         This vector comes out of the 'LoadFaultFile'
%%%                         function.

%%%P1,P2,P3                 A three Column vector where each 'P' represents
%%%                         the different corner points of one of the
%%%                         triangles (XYZ). Not as efficient in terms of
%%%                         storage but easier to understand. These
%%%                         vectors comes out of the 'LoadFaultFile'
%%%                         function

%%%FaceNormalVector         A three column vector with the unit vector
%%%                         values for the normal to each fault facet. This
%%%                         vector comes out of the 'LoadFaultFile'
%%%                         function

%%%Sxx_Fault,Syy_Fault,     One column vectors containing the 
%%%Szz_Fault,Sxy_Fault,     magnitudes of the calculated Cartesian input
%%%Sxz_Fault,Syz_Fault      stress tensor for each fault facet. These
%%%                         vectors comes out of the 'BuildStressModel'
%%%                         function

%%%Z_Obs                    The Z values of the obeservation points. This
%%%                         value comes out of the 'LoadData.m' function

%%%mu                      Coefficient of friction. This contant is set in
%%%                        the 'BuildStressModel' function

%%%Co                       Cohesive strength. This contant is set in
%%%                         the 'BuildStressModel' function

%%%mu                       Shear modulus, relates shear stress to shear 
%%%                         strain. This contant is set in
%%%                         the 'BuildStressModel' function

%%nu                        Poisson's ratio. This contant is set in
%%%                         the 'BuildStressModel' function

%%%path                     The path to the folder containing 
%%%                         the 'RunAll.m' function

%%%FaultFileString      A string used to append to the front of the file
%%%                     names of saved figures and exported data files

%%%OUTPUTS%%%

%%%Dss                  A one column vector containing the calculated
%%%                     strike-slip displacement on every fault facet
%%%                     caused by the input stress. 

%%%Dds                  A one column vector containing the calculated
%%%                     dip-slip displacement on every fault facet
%%%                     caused by the input stress.  

%%%Dn                   A one column vector containing the calculated
%%%                     normal (opening) displacement on every fault facet
%%%                     caused by the input stress.  

%%lambda                   Lame's constant.
%%%  

%%%Mu                     Coefficient of friction. This contant is set in
%%%                        the 'BuildStressModel' function, here it is made
%%%                        into a vector with length = number of
%%%                        observation points.

%%%Cohesion                Cohesive strength. This contant is set in
%%%                        the 'BuildStressModel' function, here it is made
%%%                        into a vector with length = number of
%%%                        observation points.

%%%halfspace                This constant is 0

%The function relies hevily on functions written by Tim M. Davis 
%which are part of the 'Cut and Displace' code 
%(https://github.com/Timmmdavis/CutAndDisplace) 
%https://doi.org/10.5281/zenodo.3694164

addpath(genpath(path))
cmap = colormap_cpt('Ccool-warm');
cmap2 = colormap_cpt('Ccool-warm2');

%Checks and Creates Fdisp
chk=exist('Fdisp','var');
if chk==0 %exists in workspace
    Fdisp=[];
end

halfspace = 0; 
% Calculates some needed constants based on input nu and mu
[ K,E,lambda,nu,mu ] = ElasticConstantsCheck( mu,nu );
[Tn,Tss,Tds,Sf,strain,SecondSurface ] = CreateBlankVars;
Option='B'; %allow opening components on faults
%Option='D'; %no opening components on faults
Mu(1:(length(Z_Obs)),1:1)= Myu; 
Cohesion(1:(length(Z_Obs)),1:1)= Co;

Dss(1:(length(FaultTriangles)),1:1)=0;
Dds(1:(length(FaultTriangles)),1:1)=0;
Dn(1:(length(FaultTriangles)),1:1)=0;

for i = 1:length(FaultTriangles)
Dum =  FaultTriangles(i,:);
PointsForSlipCalculator = FaultPoints(Dum(:,:),:);
PointsForSlipCalculator = PointsForSlipCalculator(:,1:4);
[Dss(i),Dds(i),Dn(i)]...
    = SlipCalculator3d(FaultMidPoint(i,:),Sxx_Fault(i),Syy_Fault(i),Szz_Fault(i),...
    Sxy_Fault(i),Sxz_Fault(i),Syz_Fault(i),...
    Tn,Tss,Tds,mu,lambda,nu,P1(i,:),P2(i,:),P3(i,:),...
    halfspace,FaceNormalVector(i,:),Fdisp,strain,Myu,Sf,Option,...
    FaultTriangles(i,:),PointsForSlipCalculator);
end

DdsTable = array2table(Dds);
DssTable = array2table(Dss);
DnTable = array2table(Dn);
FaultMidPointTable = array2table(FaultMidPoint);
FaultSlipModel = [FaultMidPointTable,DnTable,DssTable,DdsTable];

filename = strcat(FaultFileString,'FaultSlipModel.csv');
filename2 = strcat('OutputData/',filename);


writetable(FaultSlipModel, filename2,'WriteVariableNames', true);

%%% saves the workspace so (if none of the input changes) don't need to
%%% recalculate the slip

% filename = strcat(FaultFileString,'CalcSlipWorkspaceVariables');
% filename2 = strcat('OutputData/',filename);
% save(filename2)

% Location = {'OutputFigures/'};
%    PlotSlipDistribution3d(FaultTriangles,FaultPoints,cmap2,Location,FaultFileString,Dds,Dss,Dn)
   
[rows, ~] = find(isnan(Dss));
Nans = unique(rows);
empty = isempty(Nans);
for empty = 0
i= 1:length(Nans);
Dss(Nans(i,:)) = zeros;
Dn(Nans(i,:)) = zeros;
Dds(Nans(i,:)) = zeros;
end
[rows, ~] = find(isnan(Dds));
Nans = unique(rows)

meen = mean(Dds);
stand = std(Dds);
Test=meen >= 0;
if Test == 1    
Range = meen+stand;
else
Range = meen-stand;
end


DipSlip = figure;
trisurf(FaultTriangles,FaultPoints(:,2),FaultPoints(:,3),FaultPoints(:,4),Dds);
    xlabel('x'); ylabel('y'); axis('equal');
    view(-45,25)
colormap
cmocean('curl')
colorbar('eastoutside')
clim([-abs(Range) abs(Range)]);
str = append("OutputFigures/",FaultFileString,"DipSlip");
title({'\fontsize{14}DsDisp','\fontsize{8}Positive (red) is normal faulting, negative (green) is reverse faulting'})
filename = append( str,'_displacement.png');
saveas(DipSlip,string(filename))

meen = mean(Dss);
stand = std(Dss);
Test=meen >= 0;
if Test == 1    
Range = meen+stand;
else
Range = meen-stand;
end

StrikeSlip = figure;
trisurf(FaultTriangles,FaultPoints(:,2),FaultPoints(:,3),FaultPoints(:,4),Dss);
    xlabel('x'); ylabel('y'); axis('equal');
    view(-45,25)
colormap
cmocean('delta')
colorbar('eastoutside')
clim([-abs(Range) abs(Range)]);
str = append("OutputFigures/",FaultFileString,"StrikeSlip");
title({'\fontsize{14}DsDisp','\fontsize{8}Positive (green) is right-lateral strike-slipe faulting, negative (blue) is left-lateral strike-slipe faulting'})
filename = append( str,'_displacement.png');
saveas(StrikeSlip,string(filename))

meen = mean(Dn);
stand = std(Dn);
Test=meen >= 0;
if Test == 1    
Range = meen+stand;
else
Range = meen-stand;
end

Dilatation = figure;
trisurf(FaultTriangles,FaultPoints(:,2),FaultPoints(:,3),FaultPoints(:,4),Dn);
    xlabel('x'); ylabel('y'); axis('equal');
    view(-45,25)
colormap
cmocean('-balance')
colorbar('eastoutside')
clim([-abs(Range) abs(Range)]);
str = append("OutputFigures/",FaultFileString,"Dilatation");
title({'\fontsize{14}DsDisp','\fontsize{8}Positive (blue) is opening, negative (red) shortening'})
filename = append( str,'_displacement.png');
saveas(Dilatation,string(filename))

end

