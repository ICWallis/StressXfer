function [PercentCssChange] = PlotCSSChange2D(FaultTriangles,FaultPoints,...
    X_Obs,Y_Obs,Z_Obs,CSS_Change,CSS_In,X_cont,Y_cont,FaultFileString)

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

%%%CSS_In,CSS_Change   The input and change in the Coulomb
%%%                    shear traction at the observation points 
%%%                    calculated for the surfaces defined by 
%%%                    CSSFaceNormalVectorObs. These are calculated in the
%%%                    CalcStressChanges function.

%%%FaultFileString      A string used to append to the front of the file
%%%                     names of saved figures and exported data files

%%%OUTPUTS%%%

%%%PercentCssChange    The percent change in the Coulomb shear traction
%                      relative to the traction induced by the input stress 

CSSChangeFig = figure;
trisurf(FaultTriangles,FaultPoints(:,2),...
   FaultPoints(:,3),FaultPoints(:,4),'FaceColor','k','FaceAlpha',(.1));
hold on

%Use this to look at a particular map slice
CSSChangeData  = [X_Obs Y_Obs Z_Obs CSS_Change];

%This code move outlier values back toward the population (some really high
%and really low values occur when an observation point is very very close
%to a fault facet, so these need to be smooth off)
CSSChangeData_NoOutliers = filloutliers(CSS_Change,'spline','mean');
PercentCssChange = (CSSChangeData_NoOutliers./CSS_In)*-100;

meen = mean(PercentCssChange);
stand = std(PercentCssChange);

Test=meen >= 0;
if Test == 1    
Range = meen+stand;
else
Range = meen-stand;
end
caxis([-abs(Range) abs(Range)]);

cmap2 = colormap_cpt('Ccool-warm2');
colormap(cmap2);
colorbar('eastoutside');
scatter3(CSSChangeData(:,1),CSSChangeData(:,2),CSSChangeData(:,3),1,PercentCssChange, 'filled') 

xlabel('x'); ylabel('y'); axis('equal'); title('Coulomb Shear Traction Change');
view(-15,25)

filename = strcat(FaultFileString,'CSS_change');
filename2 = strcat('OutputFigures/',filename);

CssChangeGrid = reshape(PercentCssChange,size(X_cont,1),size(X_cont,2));
contourf(X_cont,Y_cont,CssChangeGrid, 20, 'EdgeColor', 'none','ContourZLevel',Z_Obs(1,:)) ; 
cmocean('-deep',20);

saveas(CSSChangeFig,filename2, 'png')
hold on

% %plots a histogram of the CSSChange data
% hold off
% histogram(PercentCssChange)
% %Plots the CSSChange data with respect to Z
end

