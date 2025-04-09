function [PercentTsChange] = PlotTsChange2D(FaultTriangles,FaultPoints,...
    X_Obs,Y_Obs,Z_Obs,Ts_Change,Ts_In,X_cont,Y_cont,FaultFileString)

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

%%%Ts_In,Ts_Change     The input and change in the 
%%%                    shear traction at the observation points 
%%%                    calculated for the surfaces defined by 
%%%                    TsFaceNormalVectorObs. These are calculated in the
%%%                    CalcStressChanges function.

%%%FaultFileString      A string used to append to the front of the file
%%%                     names of saved figures and exported data files

%%%OUTPUTS%%%

%%%PercentTsChange    The percent change in the  shear traction relative
%                     to the traction induced by the input stress
   
TschangeData  = [X_Obs Y_Obs Z_Obs Ts_Change];
TschangeData_NoOutliers = filloutliers(TschangeData(:,4),'spline','mean');
PercentTsChange = (TschangeData_NoOutliers./abs(Ts_In))*100;
PercentTsChange(isnan(PercentTsChange))=0;


ShearTraction = figure('Name','Shear tractions on observation points','NumberTitle','off');
trisurf(FaultTriangles,FaultPoints(:,2),...
   FaultPoints(:,3),FaultPoints(:,4),'FaceColor','k','FaceAlpha',(.1));hold on
scatter3(TschangeData(:,1),TschangeData(:,2),TschangeData(:,3),1,PercentTsChange, 'filled')
hold on

meen = median(PercentTsChange);
stand = std(PercentTsChange);
Test=meen >= 0;
if Test == 1    
Range = 50;
else
Range = 50;
end
caxis([0 Range]);

for i= 1:length(PercentTsChange)
    if PercentTsChange(i) > Range
  PercentTsChange(i) = Range;
    else
    if PercentTsChange(i) < -Range
  PercentTsChange(i) = -Range;
    else
  PercentTsChange(i) = PercentTsChange(i);
    end
    end
end


colormap(summer);
colorbar('eastoutside')
xlabel('x'); ylabel('y'); axis('equal'); title('Shear traction Change (increase is positive)');
view(-15,25)

filename = strcat(FaultFileString,'ShearTraction');
filename2 = strcat('OutputFigures/',filename);

TsChangeGrid = reshape(PercentTsChange,size(X_cont,1),size(X_cont,2));
contourf(X_cont,Y_cont,TsChangeGrid, 20, 'EdgeColor', 'none','ContourZLevel',Z_Obs(1,:)) ; 
cmocean('algae',20);
clim([0 abs(Range)]);

saveas(ShearTraction,filename2, 'png')

% hold off
% histogram(PercentTsChange)

end

