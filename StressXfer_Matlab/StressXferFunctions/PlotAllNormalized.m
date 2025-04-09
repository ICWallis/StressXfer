function [AllSum] = PlotAllNormalized(FaultTriangles,FaultPoints,...
    X_Obs,Y_Obs,Z_Obs,Ts_Change,Ts_In,Tn_Change,Tn_In,CSS_Change,CSS_In ...
    ,Exx,Eyy,Ezz,FaultFileString)

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

Dilatation = Exx+Eyy+Ezz;
   
TschangeData  = [X_Obs Y_Obs Z_Obs Ts_Change];
TschangeData_NoOutliers = filloutliers(TschangeData(:,4),'spline','mean');
PercentTsChange = (TschangeData_NoOutliers./abs(Ts_In))*100;
PercentTsChange(isnan(PercentTsChange))=0;

TnchangeData  = [X_Obs Y_Obs Z_Obs Tn_Change];
TnchangeData_NoOutliers = filloutliers(TnchangeData(:,4),'spline','mean');
PercentTnChange = (TnchangeData_NoOutliers./Tn_In)*100;
PercentTnChange(isnan(PercentTnChange))=0;

CSSChangeData  = [X_Obs Y_Obs Z_Obs CSS_Change];
CSSChangeData_NoOutliers = filloutliers(CSS_Change,'spline','mean');
PercentCssChange = (CSSChangeData_NoOutliers./CSS_In)*-100;
PercentCssChange(isnan(PercentCssChange))=0;

PercentCssChangeShift = PercentCssChange+min(PercentCssChange);
CSSNorm = (PercentCssChangeShift- min(PercentCssChangeShift))/(max(PercentCssChangeShift)-min(PercentCssChangeShift));

PercentTsChangeShift = PercentTsChange+min(PercentTsChange);
TsNorm = (PercentTsChangeShift- min(PercentTsChangeShift))/(max(PercentTsChangeShift)-min(PercentTsChangeShift));

PercentTnChangeShift = -1*PercentTnChange+min(PercentTnChange);
TnNorm = (PercentTnChangeShift- min(PercentTnChangeShift))/(max(PercentTnChangeShift)-min(PercentTnChangeShift));

DilitationChangeShift = -1*Dilatation+min(Dilatation);
DilitationNorm = (DilitationChangeShift- min(DilitationChangeShift))/(max(DilitationChangeShift)-min(DilitationChangeShift));

AllSum=(CSSNorm+TnNorm+TsNorm+DilitationNorm);

AllSumFig = figure('Name','Summed Stress Change and Strain','NumberTitle','off');
trisurf(FaultTriangles,FaultPoints(:,2),...
   FaultPoints(:,3),FaultPoints(:,4),'FaceColor','k','FaceAlpha',(.1));hold on
scatter3(TschangeData(:,1),TschangeData(:,2),TschangeData(:,3),10,AllSum, 'filled')
hold on

meen = mean(AllSum);
stand = std(AllSum);
Test=meen >= 0;
if Test == 1    
Range = meen+stand;
else
Range = meen-stand;
end
caxis([0 abs(Range)]);

colormap(jet);
colorbar('eastoutside')
xlabel('x'); ylabel('y'); axis('equal'); title('All Sum');
view(-15,25)

filename = strcat(FaultFileString,'AllSumFig');
filename2 = strcat('OutputFigures/',filename);

% AllSumGrid = reshape(AllSum,size(X_cont,1),size(X_cont,2));
% contourf(X_cont,Y_cont,AllSumGrid, 20, 'EdgeColor', 'none','ContourZLevel',Z_Obs(1,:)) ; 
% cmocean('-haline',20);

saveas(AllSumFig,filename2, 'png')

% hold off
% histogram(PercentTsChange)

end

