function [P1,P2,P3,FaceNormalVector,Z_Fault,FaultMidPoint,FaultTriangles,FaultPoints,X_Obs,Y_Obs,Z_Obs,X_cont,Y_cont,Z_cont] = LoadData_AutoGenObsPoints(Faults,X_ObsMax,Y_ObsMax,Z_ObsMax,X_ObsMin,Y_ObsMin,Z_ObsMin,path,FaultFileString,increment)
%This function loads the data needed for the stress transfer modeling
%function 'RunAll.m'

%%%INPUTS%%%

%%Faults                    .ts file containin the input faults. .ts files
%                           are a GOCAD format, Leapfrog can export them
%                           too. The fault should extend above the land
%                           surface, to work with the rest of the code. The
%                           below format is a fault made up of one
%                           triangular facet with vertices with coodiantes
%                           PVRTX1, PRVTX2, PRVTX3
                            
                            %GOCAD TSurf 1 
                            %HEADER {
                            %name:FaultBend
                            %}
                            %PROPERTIES 
                            %PROP_LEGAL_RANGES 
                            %NO_DATA_VALUES 
                                %PROPERTY_CLASSES 
                            %PROPERTY_KINDS 
                            %PROPERTY_SUBCLASSES QUANTITY 
                            %ESIZES 
                            %UNITS 
                            %TFACE
                            %PVRTX 1 433468 4397958 1372
                            %PVRTX 2 435226 4396927 234
                            %PVRTX 3 434987 4397251 221
                            %TRGL 1 2 3

%%%Observation points       .csv file containing the points around the
%                           faults for which stress change and strain will
%                           be calculated. Observations points should span
%                           the area of interest, but not extend past the
%                           fault. Use UTM so X,Y,Z are all in meters.
%                           Format: 
%                           X,Y,Z 
%                           439092.7,4373573.3,-5
%                           440092.7,4373573.3,-5

%%%path                     The path to the folder containing 
%                           the 'RunAll.m' function

%%%OUTPUTS%%%

%%%P1,P2,P3                 A three Column vector where each 'P' represents
%%%                         the different corner points of one of the
%%%                         triangles (XYZ). Not as efficient in terms of
%%%                         storage but easier to understand.

%%%FaceNormalVector         A three column vector with the unit vector
%%%                         values for the normal to each fault facet.

%%%Z_Fault                  The Z values of the fault facet midpoints. This
%%%                         value comes out of the 'LoadData.m' function

%%%FaultMidPoint            A three column vector containing the x,y,z
%%%                         coordinate at the center of each fault facet
                       
%%%FaultTriangles           A three column vector containing the vertex
%%%                         number for each fault facet. 

%%%FaultPoints              A four column vector containing the vertices
%%%                         number and x,y,z coroodinates for each fault
%%%                         facet vertex.

%%%X_Obs,Y_Obs,Z_Obs    One columns vectors containing the coordiates of
%%%                     the observation points.

%%%FaultFileString      A string used to append to the front of the file
%%%                     names of saved figures and exported data files

close all
addpath(genpath(path))


%The 'LoadFaultFile' function was written by Tim M. Davis as part of
%the 'Cut and Displace' code (https://github.com/Timmmdavis/CutAndDisplace)
%https://doi.org/10.5281/zenodo.3694164. It is slightly adapted for use in 
%this function

 [P1,P2,P3,FaceNormalVector,~,~,Z_Fault,FaultMidPoint,FaultTriangles,FaultPoints]...
 = LoadFaultFile(Faults);

%The 'LoadFaultFile' function was written by Tim M. Davis as part of
%the 'Cut and Displace' code (https://github.com/Timmmdavis/CutAndDisplace)
%https://doi.org/10.5281/zenodo.3694164. It is slightly adapted for use in 
%this function

X_cont = (X_ObsMin:increment:X_ObsMax);
Y_cont = (Y_ObsMin:increment:Y_ObsMax);
Z_cont = (Z_ObsMin:increment:Z_ObsMax);
sz = length(X_cont)*length(Y_cont);
[X_cont,Y_cont,Z_cont] = meshgrid(X_cont,Y_cont,Z_cont);
X_Obs = reshape(X_cont,sz,1);
Y_Obs = reshape(Y_cont,sz,1); 
Z_Obs = reshape(Z_cont,sz,1);


%This code prints the number of fault facets
NumberofFaultFacets = length(FaceNormalVector)
NumberofObservationPoints = length(Z_Obs)

%This code plots  a figure that shows the faults and the observation
%points and saves it in the 'OutputFigures' folder. The normals to the
%fault facets are plotted in red. Normals should point up, if not flip the
%fault(s)

FaultsandObsPoints = figure;
scatter3(X_Obs,Y_Obs,Z_Obs,10,'b','filled');hold on
trisurf(FaultTriangles,FaultPoints(:,2),...
   FaultPoints(:,3),FaultPoints(:,4),'FaceColor','k','FaceAlpha',(.1));
hold on
quiver3(FaultMidPoint(:,1),FaultMidPoint(:,2),FaultMidPoint(:,3), ...
FaceNormalVector(:,1),FaceNormalVector(:,2),FaceNormalVector(:,3),0.5, 'color','r');
title('Input faults and observation points')
axis equal
filename = strcat(FaultFileString,'FaultsandObsPoints');
filename2 = strcat('OutputFigures/',filename);
saveas(FaultsandObsPoints,filename2, 'png')



FaultPointData = array2table([FaultMidPoint(:,1),FaultMidPoint(:,2),FaultMidPoint(:,3)]);

filename = strcat(FaultFileString,'FaultPointData.csv');
filename2 = strcat('OutputData/',filename);

writetable(FaultPointData, filename2,'WriteVariableNames', true);
end


