function [DilatationData_NoOutliers] = PlotDilatation(FaultTriangles,...
    FaultPoints,X_Obs,Y_Obs,Z_Obs,Exx,Eyy,Ezz,FaultFileString)

%This function plots dilatation at the observation points as a result of
%the modeled slip. This is part of the stress transfer modeling function
%'RunAll.m'

%%%INPUTS%%%

%%%FaultTriangles           A three column vector containing the vertex
%%%                         number for each fault facet. These are
%%%                         calculated in the 'LoadData' function.

%%%FaultPoints              A four column vector containing the vertices
%%%                         number and x,y,z coroodinates for each fault
%%%                         facet vertex. These are calculated in the
%%%                         'LoadData' function.

%%%X_Obs,Y_Obs,Z_Obs    One columns vectors containing the coordiates of
%%%                     the observation points. These are calculated in the
%%%                     'LoadObservationPoints' function

%%%Exx,Eyy,Ezz          Strain change induced by movement of dislocations.
%%%                     These are calculated in the
%%%                     'StressDisplacementStran' function

%%%FaultFileString      A string used to append to the front of the file
%%%                     names of saved figures and exported data files

%%%OUTPUTS%%%

%%%DilatationData_NoOutliers    The dilataion at the observation points as
%%%                             a result of the modeled slip

Dilatation = Exx+Eyy+Ezz;

DilatationData  = [X_Obs Y_Obs Z_Obs Dilatation];

DilatationData_NoOutliers = filloutliers(Dilatation,'spline','mean');

DilatationFig = figure('Name','Dilatation on observation points',...
    'NumberTitle','off');
trisurf(FaultTriangles,FaultPoints(:,2),...
   FaultPoints(:,3),FaultPoints(:,4),'FaceColor','k','FaceAlpha',(.1));
hold on
scatter3(DilatationData(:,1),DilatationData(:,2),DilatationData(:,3),10,DilatationData_NoOutliers, 'filled')

hold on
xlabel('x'); ylabel('y'); axis('equal'); title('Dilatation');
hold on

meen = mean(DilatationData_NoOutliers);
stand = std(DilatationData_NoOutliers);
Test=meen >= 0;
if Test == 1    
Range = meen+stand/3;
else
Range = meen-stand/3;
end
caxis([-abs(Range) abs(Range)]);
hold on


colormap default;
colorbar('eastoutside')
view(-15,25)
filename = strcat(FaultFileString,'Dilatation');
filename2 = strcat('OutputFigures/',filename);

% DilatationGrid = reshape(DilatationData_NoOutliers,size(X_cont,1),size(X_cont,2));
% contourf(X_cont,Y_cont,DilatationGrid, 20, 'EdgeColor', 'none','ContourZLevel',Z_Obs(1,:)) ; 
% cmocean('thermal',20);

saveas(DilatationFig,filename2, 'png')
%     hold off
% histogram(DilatationData_NoOutliers)
end

