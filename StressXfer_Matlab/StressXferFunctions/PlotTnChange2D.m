function [PercentTnChange] = PlotTnChange2D(FaultTriangles,FaultPoints,...
    X_Obs,Y_Obs,Z_Obs,Tn_Change,Tn_In,X_cont,Y_cont,FaultFileString)

%This function plots the percent change in the Coulomb shear traction
%at the observation points as a result of the modeled slip. This is part 
%of the stress transfer modeling function 'RunAll.m'

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

%%%Tn_In,Tn_Change     The input and change in the 
%%%                    normal traction at the observation points 
%%%                    calculated for the surfaces defined by 
%%%                    TnFaceNormalVectorObs. These are calculated in the
%%%                    CalcStressChanges function.

%%%FaultFileString      A string used to append to the front of the file
%%%                     names of saved figures and exported data files

%%%OUTPUTS%%%

%%%PercentTnChange    The percent change in the normal traction relative
%                     to the traction induced by the input stress

TnchangeData  = [X_Obs Y_Obs Z_Obs Tn_Change];

TnchangeData_NoOutliers = filloutliers(TnchangeData(:,4),'spline','mean');

PercentTnChange = (TnchangeData_NoOutliers./Tn_In)*100;
PercentTnChange(isnan(PercentTnChange))=0;

NormalTraction = figure('Name','Normal Stress Change on observation points','NumberTitle','off');
trisurf(FaultTriangles,FaultPoints(:,2),...
   FaultPoints(:,3),FaultPoints(:,4),'FaceColor','k','FaceAlpha',(.1));
hold on
scatter3(TnchangeData(:,1),TnchangeData(:,2),TnchangeData(:,3),1,PercentTnChange, 'filled') 
hold on

meen = mean(PercentTnChange);
stand = std(PercentTnChange);
Test=meen >= 0;
if Test == 1    
Range = meen+stand;
else
Range = meen-stand;
end
caxis([-abs(Range) abs(Range)]);

colormap spring;
colorbar('eastoutside')
xlabel('x'); ylabel('y'); axis('equal'); title('Normal Traction Change');
view(-15,25)

filename = strcat(FaultFileString,'NormalTraction');
filename2 = strcat('OutputFigures/',filename);

TnChangeGrid = reshape(PercentTnChange,size(X_cont,1),size(X_cont,2));
contourf(X_cont,Y_cont,TnChangeGrid, 20, 'EdgeColor', 'none','ContourZLevel',Z_Obs(1,:)) ; 
cmocean('matter',20);

saveas(NormalTraction,filename2, 'png')
% hold off
% histogram(PercentTnChange)

end

