function [TsFaceNormalVectorObs,TnFaceNormalVectorObs,CSSFaceNormalVectorObs,...
    TsRake,TnRake,CSSRake] = SetObsPointAttitudes(FaceNormalVector,Z_Obs...
    ,Dss,Dds,Sxx_Change,Syy_Change,Szz_Change,Sxy_Change,Sxz_Change,...
    Syz_Change,Mu,Cohesion,spacing,path,FaultFileString)

%This function calculates the geometry of the surfaces at the observation
%points that will be examined to evaluate stress change.These are needed
%for the stress transfer modeling function 'RunAll.m'

%%%INPUTS%%%

%%%FaceNormalVector         A three column vector with the unit vector
%%%                         values for the normal to each fault facet. This
%%%                         vector comes out of the 'LoadFaultFile'
%%%                         function

%%%Z_Obs                    The Z values of the obeservation points. This
%%%                         value comes out of the 'LoadData.m' function

%%%Dss                  A one column vector containing the calculated
%%%                     strike-slip displacement on every fault facet
%%%                     caused by the input stress. This vector is
%%%                     calculated in the 'CalcSlip' function

%%%Dds                  A one column vector containing the calculated
%%%                     dip-slip displacement on every fault facet
%%%                     caused by the input stress.  This vector is
%%%                     calculated in the 'CalcSlip' function

%%%Dn                   A one column vector containing the calculated
%%%                     normal (opening) displacement on every fault facet
%%%                     caused by the input stress.  This vector is
%%%                     calculated in the 'CalcSlip' function

%%%Sxx_Obs_Change,Syy_Obs_Change,     One column vectors containing the 
%%%Szz_Obs_Change,Sxy_Obs_Change,     magnitudes of the calculated Cartesian
%%%Sxz_Obs_Change,Syz_Obs_Change      stress tensor change, that is the stress 
%%%                                   difference between the input and output
%%%                                   stress,for each observation point.
%%%                                   These are calculated in the 'Stress
%%%                                   DisplacementStrain' function

%%%Mu                      Coefficient of friction. This contant is set in
%%%                        the 'BuildStressModel' function, here it is a
%%%                        vector with length = number of observation
%%%                        points.

%%%Cohesion                Cohesive strength. This contant is set in
%%%                        the 'BuildStressModel' function, here it is a
%%%                        vector with length = number of observation
%%%                        points.

%%% spacing                The spacing for the equally spaced trends and
%%%                        plunges that will be evaluated to see which is
%%%                        the most optimally oriented surfaces. This
%%%                        constant is set in the 'RunAll' function.

%%%path                     The path to the folder containing 
%%%                         the 'RunAll.m' function

%%%FaultFileString      A string used to append to the front of the file
%%%                     names of saved figures and exported data files

%%%OUTPUTS%%%
%%%TsFaceNormalVectorObs    The normals to the surfaces that result in the
%%%TnFaceNormalVectorObs    maximum shear traction increase, normal
%%%CSSFaceNormalVectorObs   traction decrease, and Coulomb shear traction
%%%                         increase at each observation point

%%%TsRake,TnRake,CSSRake    The rakes that result in the
%%%                         maximum shear traction increase, normal
%%%                         traction decrease, and Coulomb shear traction
%%%                         increase at each observation point

addpath(genpath(path))
[PlungeFault,TrendFault] = PlungeAndTrend(FaceNormalVector);
TrendFault = TrendFault+180;
for i = 1:length(TrendFault)
    Test  = TrendFault(i) > 360;
        if Test == 1
        TrendFault(i) = TrendFault(i) - 360;
        else
        TrendFault(i) = TrendFault(i);
        end
end
TrendFault = filloutliers(TrendFault,'spline','mean');
PlungeFault = filloutliers(PlungeFault,'spline','mean');
MinPlunge = min(PlungeFault); MaxPlunge = max(PlungeFault); stdPlunge = std(PlungeFault); 
MinTrend = min(TrendFault); MaxTrend = max(TrendFault); stdTrend = std(TrendFault); 

Test = TrendFault < 180;
Test2 = TrendFault >= 180;
        MinTrendEast = min(TrendFault(Test));
        MaxTrendEast = max(TrendFault(Test));
        stdTrendEast = std(TrendFault(Test));
        MinTrendWest = min(TrendFault(Test2));
        stdTrendWest = std(TrendFault(Test2));
        MaxTrendWest = max(TrendFault(Test2));

AbsDss =abs(Dss);
AbsDds =abs(Dds);
Slip(1:(length(Dds)),1:1)=0;

for i = 1:length(AbsDss)
%for i = 1:10
    Index = AbsDds(i) == 0;
    Index2 = AbsDds(i)>= AbsDss(i);
        if Index == 1
            Slip(i) = 180;
        else
            if Index2 == 1
                Slip(i) = 180+((AbsDss(i)/AbsDds(i))*90);
            else
                Slip(i) = 90+cos(AbsDds(i)/AbsDss(i))*90;
            end
            
        end
end
RakeFault = filloutliers(Slip,'spline','mean');
MinRake = min(RakeFault); MaxRake = max(RakeFault); stdRake = std(RakeFault); 

PlngSpc= spacing; %Set the increment of the plunge vector 
TrndSpc1= spacing/2; %Set the increment of the trend vector 
TrndSpc= spacing; %Set the increment of the trend vector 
RakeSpc = spacing/2; %Set the increment of the rake vector 

P = linspace(MinPlunge-stdPlunge,MaxPlunge+stdPlunge,PlngSpc)';P(P > 90) = 90;P(P < 0) = 0;
empty = isempty( MinTrendEast );
empty2 = isempty( MinTrendWest ); 
if empty == 0 && empty2 ==0
T1 = linspace(MinTrendWest-stdTrendWest,MaxTrendWest+stdTrendWest,TrndSpc1)'-180;
T2 = linspace(MinTrendEast-stdTrendEast,MaxTrendEast+stdTrendEast,TrndSpc1)'-180;
  
T = vertcat(T1,T2);
else
T = linspace(MinTrend-stdTrend,MaxTrend+stdTrend,TrndSpc)'-180;
end
for i = 1:length(T)
    Test  = T(i) > 360;
        if Test == 1
        T(i) = T(i) - 360;
        else
        T(i) = T(i);
        end
end
for i = 1:length(T)
    Test  = T(i) < 0;
        if Test == 1
        T(i) = T(i) + 360;
        else
        T(i) = T(i);
        end
end

R = linspace(MinRake-stdRake,MaxRake+stdRake,RakeSpc)'; R(R > 270) = 270;R(R < 90) = 90;

Trend = sort(repmat(T',1,PlngSpc*RakeSpc))';
Plunge = repmat(P',1,TrndSpc*RakeSpc)';
Rke = sort(repmat(R',1,PlngSpc))';
Rake = repmat(Rke',1,TrndSpc)';

TPR = [Trend Plunge Rake]; %should be one unique set of TPR values at the increments above

%calculate the components of the vector from the T and P values
x =  cos(deg2rad(Plunge)).*cos(deg2rad(-Trend+90));
y =  cos(deg2rad(Plunge)).*sin(deg2rad(-Trend+90));
z =  sin(deg2rad(Plunge));

FaceNormalVectorObs = [x y z]; 
%vector = FNVNear;
 Ts_Chnge = zeros(length(FaceNormalVectorObs(:,1)),1); Tn_Chnge = zeros(length(FaceNormalVectorObs(:,1)),1); CSS_Chnge = zeros(length(FaceNormalVectorObs(:,1)),1); ...
 Tss_Chnge = zeros(length(FaceNormalVectorObs(:,1)),1); Tds_Chnge = zeros(length(FaceNormalVectorObs(:,1)),1);
 Ts_Change = zeros(length(Z_Obs),1); Tn_Change = zeros(length(Z_Obs),1); CSS_Change = zeros(length(Z_Obs),1);...
      Tss_Change = zeros(length(Z_Obs),1); Tds_Change = zeros(length(Z_Obs),1);
    Index_Ts = zeros(length(Z_Obs),1);Index_Tn = zeros(length(Z_Obs),1);Index_CSS = zeros(length(Z_Obs),1);...
        Index_Tss = zeros(length(Z_Obs),1);Index_Tds = zeros(length(Z_Obs),1);
    
%use a loop to calculate the CSS, Tn, Ts, Tss, and Tds for each possible observation point geometry    
%for i = 1:1
for i = 1:length(Z_Obs) 
    %for j = 1:1
    for j = 1:length(FaceNormalVectorObs(:,1))
   [ Ts_Chnge(j), Tn_Chnge(j), CSS_Chnge(j), Tss_Chnge(j), Tds_Chnge(j)  ] = CalculateCoulombStressOnPlane...
    (FaceNormalVectorObs(j,:),Sxx_Change(i),Syy_Change(i),Szz_Change(i),...
Sxy_Change(i),Sxz_Change(i),Syz_Change(i),Mu(i),Cohesion(i),Rake(j));    
    [~,Index_Ts(i)] = max(Ts_Chnge); %take the max Ts Change
    [~,Index_Tn(i)] = min(Tn_Chnge); %take the min Tn Change
    [~,Index_CSS(i)] = max(CSS_Chnge); %take the max CSS Change
    [~,Index_Tss(i)] = max(Tss_Chnge); %take the max Tss Change
    [~,Index_Tds(i)] = max(Tds_Chnge); %take the max Tds Change
    end
end
MinTnNormals = FaceNormalVectorObs(Index_Tn,:); % get the observation point geometries that result in the minimum Tn change
% get the rake of the observation points that results in the highest CSS change
TsRake = Rake(Index_Ts,:); 
TnRake = Rake(Index_Tn,:); 
CSSRake = Rake(Index_CSS,:); 

TsFaceNormalVectorObs = FaceNormalVectorObs(Index_Ts,:);
TnFaceNormalVectorObs = FaceNormalVectorObs(Index_Tn,:);
CSSFaceNormalVectorObs = FaceNormalVectorObs(Index_CSS,:);

[TsPlungeObs,TsTrendObs] = PlungeAndTrend(TsFaceNormalVectorObs);
TP = {'Trend', 'Plunge'};
TsTrendObs = TsTrendObs+180;
for i = 1:length(TsTrendObs)
    Test  = TsTrendObs(i) > 360;
        if Test == 1
        TsTrendObs(i) = TsTrendObs(i) - 360;
        else
        TsTrendObs(i) = TsTrendObs(i);
        end
end

ObsPoles1 = table(TsTrendObs,TsPlungeObs, 'VariableNames',TP);

filename = strcat(FaultFileString,'TsObservationPointPoles.txt');
filename2 = strcat('OutputData/',filename);

writetable(ObsPoles1, filename2,'WriteVariableNames', true);

[TnPlungeObs,TnTrendObs] = PlungeAndTrend(TnFaceNormalVectorObs);
TnTrendObs = TnTrendObs+180;
for i = 1:length(TnTrendObs)
    Test  = TnTrendObs(i) > 360;
        if Test == 1
        TnTrendObs(i) = TnTrendObs(i) - 360;
        else
        TnTrendObs(i) = TnTrendObs(i);
        end
end

ObsPoles2 = table(TnTrendObs,TnPlungeObs, 'VariableNames',TP);

filename = strcat(FaultFileString,'TnObservationPointPoles.txt');
filename2 = strcat('OutputData/',filename);

writetable(ObsPoles2, filename2,'WriteVariableNames', true);


[CSSPlungeObs,CSSTrendObs] = PlungeAndTrend(CSSFaceNormalVectorObs);
CSSTrendObs = CSSTrendObs+180;
for i = 1:length(CSSTrendObs)
    Test  = CSSTrendObs(i) > 360;
        if Test == 1
        CSSTrendObs(i) = CSSTrendObs(i) - 360;
        else
        CSSTrendObs(i) = CSSTrendObs(i);
        end
end

ObsPoles3 = table(CSSTrendObs,CSSPlungeObs, 'VariableNames',TP);

filename = strcat(FaultFileString,'CSSObservationPointPoles.txt');
filename2 = strcat('OutputData/',filename);

writetable(ObsPoles3, filename2,'WriteVariableNames', true);

TrendObs = vertcat(CSSTrendObs,TsTrendObs,TnTrendObs);
PlungeObs = vertcat(CSSPlungeObs,TsPlungeObs,TnPlungeObs);

ObsPoles = table(TrendObs,PlungeObs, 'VariableNames',TP);
filename = strcat(FaultFileString,'AllObservationPoles.txt');
filename2 = strcat('OutputData/',filename);
writetable(ObsPoles, filename2,'WriteVariableNames', true);

[FaultPlunge,FaultTrend] = PlungeAndTrend(FaceNormalVector);
FaultTrend = FaultTrend+180;
for i = 1:length(FaultTrend)
    Test  = FaultTrend(i) > 360;
        if Test == 1
        FaultTrend(i) = FaultTrend(i) - 360;
        else
        FaultTrend(i) = FaultTrend(i);
        end
end
FaultTrend = filloutliers(FaultTrend,'spline','mean');
FaultPlunge = filloutliers(FaultPlunge,'spline','mean');
ObsPoles = table(FaultTrend,FaultPlunge, 'VariableNames',TP);

filename = strcat(FaultFileString,'InputFaultPoles.txt');
filename2 = strcat('OutputData/',filename);
writetable(ObsPoles, filename2,'WriteVariableNames', true);

% filename = strcat(FaultFileString,'SetObsPointsWorkspaceVariables');
% filename2 = strcat('OutputData/',filename);
% save(filename2)

end

