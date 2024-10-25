function [S1In,S2In,S3In,S1dirIn,S2dirIn,S3dirIn,S1Out,S2Out,S3Out,S1dirOut,S2dirOut,S3dirOut,S1,S2,S3,S1dir,S2dir,S3dir] = PlotOutputStresses(Sxx_Obs_Out,Syy_Obs_Out,Szz_Obs_Out,...
    Sxy_Obs_Out,Sxz_Obs_Out,Syz_Obs_Out,Sxx_Obs_In,Syy_Obs_In,Szz_Obs_In,...
    Sxy_Obs_In,Sxz_Obs_In,Syz_Obs_In,X_Obs, Y_Obs, Z_Obs,FaultPoints,FaultTriangles,FaultFileString)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
StressTensorOut = [Sxx_Obs_Out,Syy_Obs_Out,Szz_Obs_Out,...
    Sxy_Obs_Out,Sxz_Obs_Out,Syz_Obs_Out];
StressTensorIn = [Sxx_Obs_In,Syy_Obs_In,Szz_Obs_In,...
    Sxy_Obs_In,Sxz_Obs_In,Syz_Obs_In];

[S1In,S2In,S3In,S1dirIn,S2dirIn,S3dirIn] = EigCalc3d(Sxx_Obs_In(:),Syy_Obs_In(:),Szz_Obs_In(:),Sxy_Obs_In(:),Sxz_Obs_In(:),Syz_Obs_In(:));
[S1Out,S2Out,S3Out,S1dirOut,S2dirOut,S3dirOut] = EigCalc3d(Sxx_Obs_Out(:),Syy_Obs_Out(:),Szz_Obs_Out(:),Sxy_Obs_Out(:),Sxz_Obs_Out(:),Syz_Obs_Out(:));
%DrawOutputS1S2S3Directions(StressTensorOut, X_Obs, Y_Obs, Z_Obs,FaultFileString,'Scale', 20, 'Points', FaultPoints, 'Triangles', FaultTriangles);
[S1,S2,S3,S1dir,S2dir,S3dir] = EigCalc3d(Sxx_Obs_Out(:)-Sxx_Obs_In(:),Syy_Obs_Out(:)-Syy_Obs_In(:),Szz_Obs_Out(:)-Szz_Obs_In(:),Sxy_Obs_Out(:)-Sxy_Obs_In(:),Sxz_Obs_Out(:)-Sxz_Obs_In(:),Syz_Obs_Out(:)-Syz_Obs_In(:));
S1Mag_pchange = S1./S1In.*100;
S2Mag_pchange = S2./S2In.*100;
S3Mag_pchange = S3./S3In.*100;
Sdiffchange = S1-S3;
SdiffIn = S1In-S3In;

Sdiff_pchange = Sdiffchange./SdiffIn.*100;

StressChangeModel = array2table([X_Obs,Y_Obs,Z_Obs,S1,S2,S3,S1dir,S2dir,S3dir]);
StressChangeModel = renamevars(StressChangeModel,["Var1","Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10","Var11","Var12","Var13","Var14","Var15"],["X_Obs", "Y_Obs","Z_Obs","S1_mag","S2_mag", "S3_mag", "S1DirX", "S1DirY", "S1DirZ", "S2DirX", "S2DirY", "S2DirZ", "S3DirX", "S3DirY", "S3DirZ"]);
filename = strcat(FaultFileString,'ChangeFromInitialStress.csv');
filename2 = strcat('OutputData/',filename);
writetable(StressChangeModel, filename2,'WriteVariableNames', true);
hold off


DiffStressHistogram = figure('Name','Percent Change in Stress Magnitudes','NumberTitle','off');
StdS = std(vertcat(S1Mag_pchange,S2Mag_pchange,S3Mag_pchange));
StdR = 1.5*StdS;
edges = linspace((0-StdR),(0+StdR),30);
h1 = histogram(S1Mag_pchange,'BinEdges',edges);
h1.FaceColor = 'r';
hold on
h2 = histogram(S2Mag_pchange,'BinEdges',edges);
h2.FaceColor = 'g';
hold on
h3 = histogram(S3Mag_pchange,'BinEdges',edges);
h3.FaceColor = 'b';
h4 = histogram(Sdiff_pchange,'BinEdges',edges);
h4.FaceColor = 'y';

xlim ([(-1*StdR) (StdR)])
title('Percent Change in Stress Magnitudes');
xlabel('Change in Stress Magnitudes (%)');
legend('S1','S2', 'S3','SDiff', 'Location','northeast')
filename = strcat(FaultFileString,'StressChangeHistogram');
filename2 = strcat('OutputFigures/',filename);
saveas(DiffStressHistogram,filename2, 'png')


DiffStressChange = figure('Name','Percent Change in Differential Stress','NumberTitle','off');
trisurf(FaultTriangles,FaultPoints(:,2),...
   FaultPoints(:,3),FaultPoints(:,4),'FaceColor','k','FaceAlpha',(.1));
hold on
scatter3(X_Obs(:),Y_Obs(:),Z_Obs(:),10,Sdiff_pchange(:) ,'filled') 
hold on
meen = mean(Sdiff_pchange);
Range = meen;
clim([-abs(Range) abs(Range)]);
hold on
% xlim([min(X_Obs) max(X_Obs)])
% hold on
% ylim([min(Y_Obs) max(Y_Obs)])
colormap jet;
colorbar('eastoutside')
xlabel('x'); ylabel('y');zlabel('z'); axis('equal'); title('Percent Change in Differential Stress'); subtitle('circle size and color correspond to change')
view(-15,25)
filename = strcat(FaultFileString,'DiffStressChange');
filename2 = strcat('OutputFigures/',filename);

% DiffStressChangeGrid = reshape(Sdiff_pchange,size(X_cont,1),size(X_cont,2),size(X_cont,3));
% contourf(X_cont,Y_cont,DiffStressChangeGrid, 20, 'EdgeColor', 'none','ContourZLevel',Z_Obs(1,:)) ; 
% cmocean('-speed',20);

saveas(DiffStressChange,filename2, 'png')

end