function [] = RunAll(~)
%This Matlab function calculates the slip on input faults as a result of an
%applied stress and the resultant stress changes and strain in the area
%surrounding the faults. This function was written by Drew Siler, 2022. The
%function calls many functions written by Tim M. Davis which are part of
%the 'Cut and Displace' code (https://github.com/Timmmdavis/CutAndDisplace)
%https://doi.org/10.5281/zenodo.3694164

%Define the folder path, fault (.ts file), and observation points (x,y,z
%points. Fault files with ~2600 facets and ~300 observation points run in
%about 3 minutes. Larger files will take longer.
%******STEP 1: SPECIFY INPUT DATA****
path = 'C:\Users\dsile\Geologica Dropbox\Drew Siler\VOI_workshop\StressXfer\'; %%%The path to the folder containing
addpath(genpath(path));                                       %%%the 'RunAll.m' function
%addpath(genpath('C:\Users\dsile\Geologica Dropbox\Drew Siler\DynamicStressModeling\StressTransferModeling'))

Fault = 'SimpleFault.ts';      %%%The fault file, see inclosed .ts 
                                       %%%file as an example of the format
                                       
ObservationPoints = 'Block Model_1000mcell_0melev.csv'; %%% The observation 
                                                %%%points around the fault
                                                %%%see enclosed .csv file 0
                                                %%%as an example of the
                                                %%%format
FaultFileString = 'SimpleFault';          %%% A string to append to 
                                                %%% to the front of the
                                                %%% output file, leave
                                                %%% blank ('') if you don't want
increment = 1000;                                                %%% to use it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                
%******STEP 2: LOAD DATA****

%The 'LoadData.m' function loads the fault file and observation points file
%and uses those data to define several values needed for the calculation
[P1,P2,P3,FaceNormalVector,Z_Fault,FaultMidPoint...
    ,FaultTriangles,FaultPoints,X_Obs,Y_Obs,Z_Obs,X_cont,Y_cont]...
    = LoadData_InputObsPoints(Fault, ObservationPoints,path,FaultFileString,increment);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%******STEP 3: BUILD STRESS MODEL****

%Define S3 (minimum horizontal stress) direction. In this code we assume
%one principal stress is vertical and the other two are horizontal.

%Define Aphi, a numeric measure of the stress regime. 0.5 is normal
%faulting, 1.5 is strike-slip faulting, 2.5 is reverse faulting. See
%(Simpson, 1997) for more on Aphi. This block calls the function
%'aphi_to_stressmags' written by Jens Lund Snee. The formalae can be found
%here; https://doi.org/10.1306/08102120151

S3_Direction = 101; %%%degrees azimuth clockwise, i.e. 90 equals due east  
                   %%%or an east-west S3.
                   
Aphi = 0.5;        %%%0.5 is normal faulting, 1.5 is strike-slip faulting, 
                   %%%2.5 is reverse faulting. See (Simpson, 1997) 
                   %%%for more on Aphi.
                   
Elevation = 750; %%%elevation that the vertical stress is calculated for in
               %%%meters
               
MeanSurfaceElevation = 1750; %%%%mean surface elevation in meters

%The 'BuildStressModel.m' function calculates S1,S2,S3 at the input
%elevation from the inputs. Then a vertical stress model is built assumeing
%all stresses are zero at the surface and increase linearly with depth. The
%calculated stresses are projected to all input fault facets and
%observation points for later use.


[Sxx_Fault,Syy_Fault,Szz_Fault,Sxy_Fault,...
    Sxz_Fault,Syz_Fault,Sxx_Obs_In,Syy_Obs_In,...
    Szz_Obs_In,Sxy_Obs_In,Sxz_Obs_In,Syz_Obs_In,mu,nu,Myu,Co]...
    = BuildStressModel(S3_Direction,Aphi, Z_Obs, Z_Fault, Elevation,...
    MeanSurfaceElevation,path,FaultFileString);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%******STEP 4: CALCULATE SLIP****

%The CalcSlip.m function calculates the slip (strike-slip, dip-slip, and
%normal (opening) that occurs as a result of the input stress model

%The function relies hevily on functions written by Tim M. Davis 
%which are part of the 'Cut and Displace' code 
%(https://github.com/Timmmdavis/CutAndDisplace) 
%https://doi.org/10.5281/zenodo.3694164

%This code outputs 'FaultSlipModel.csv' is the which is thecalculated slip
%on the input faults as a result of the input stress. In this file the mid
%point of each fault facet (x,y,z) and the Dss, Dds, and Dn are reported.

[Dss,Dds,Dn,lambda,Mu,Cohesion,halfspace] = CalcSlip(FaultTriangles,...
    FaultPoints,FaultMidPoint,P1,P2,P3,FaceNormalVector,Sxx_Fault,...
    Syy_Fault,Szz_Fault,Sxy_Fault,Sxz_Fault,Syz_Fault,Z_Obs,Myu,...
    Co,mu,nu,path,FaultFileString);

disp('Done with slip calculation')

filename = strcat(FaultFileString,'STEP_4_Variables.mat');
filename2 = strcat('OutputData/',filename);
save(filename2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%******STEP 5: CALCULATE STRESS CHANGE AND STRAIN****

%The 'StressDisplacementsStrain.m' function calculates the output stress
%and stress changes that occur on the observation points as a result of the
%fault slip

%This function relies hevily on functions written by Tim M. Davis 
%which are part of the 'Cut and Displace' code 
%(https://github.com/Timmmdavis/CutAndDisplace) 
%https://doi.org/10.5281/zenodo.3694164

[Sxx_Change,Syy_Change,Szz_Change,Sxy_Change,Sxz_Change,Syz_Change,...
    Sxx_Obs_Out,Syy_Obs_Out,Szz_Obs_Out,Sxy_Obs_Out,Sxz_Obs_Out,...
    Syz_Obs_Out,Exx,Eyy,Ezz] = StressDisplacementStrain(Dss,Dds,Dn,mu,...
    X_Obs,Y_Obs,Z_Obs,Sxx_Obs_In,Syy_Obs_In,Szz_Obs_In,Sxy_Obs_In,...
    Sxz_Obs_In,Syz_Obs_In,P1,P2,P3,halfspace,nu,lambda,path, ...
    FaultFileString,FaultTriangles, FaultPoints);

filename = strcat(FaultFileString,'STEP_5_Variables.mat');
filename2 = strcat('OutputData/',filename);
save(filename2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%******STEP 6: CALCULATE OBSERVATION POINT GEOMETRIES****

%The 'SetObsAttitudes.m' function calculates the shear traction (Ts),
%normal (Tn), and Coulomb Shear Traction (CSS) for an array (length defined
%below) of fracture orientations within 1 standard deviation of the trend,
%plunge, of the input faults, and the rake of the initial slip). The code
%calculates the plane and rake that will result in the maximum CSS
%increase, maximum Ts increase and maximum Tn decrease. Essentially, it
%calculates the most optimally oriented fracutres with similar geometry to
%the input faults.

spacing = 20; %%%this is the spacing of trend, plunge, and rake 
%(within 1 standard deviation) points the input fault that that code will
%evaluate for the maximum stress changes. 20 is a good number here, much
%more than that and the code slows way down. Note that for input fault data
%with oppositely dipping faults, this code will break the search for
%optimate orientation into eastern and western hemisphere. I recommend
%taking the 'InputFaultPoles.txt', 'CSSObservationPointPoles.txt',
%'TsObservationPointPoles.txt', 'TsObservationPointPoles.txt,' which are
%outputs from this function, and plotting them on stereonet (try Stereonet
%http://www.geo.cornell.edu/geology/faculty/RWA/programs/stereonet.html) to
%make sure this part is working the way you want it to.

[TsFaceNormalVectorObs,TnFaceNormalVectorObs,CSSFaceNormalVectorObs,...
    TsRake,TnRake,CSSRake] = SetObsPointAttitudes(FaceNormalVector,...
    Z_Obs,Dss,Dds,Sxx_Change,Syy_Change,Szz_Change,Sxy_Change,...
    Sxz_Change,Syz_Change,Mu,Cohesion,spacing,path,FaultFileString);

filename = strcat(FaultFileString,'STEP_6_Variables.mat');
filename2 = strcat('OutputData/',filename);
save(filename2)

disp('Done calculating observation points')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The 'CalcStressChanges.m' function calculates the tractions (Ts, Tn, and
%CSS) on the planes/rakes defined by the SetObsPointAttitudes function from
%the input and output (stress resulting from the modeled slip) stresses. It
%calculates the stress changes too.
%******STEP 7: CALCULATE STRESS CHANGES ON PLANES****

[Ts_In,Ts_Out,Ts_Change,Tn_In,Tn_Out,Tn_Change,CSS_In,CSS_Out,CSS_Change]... 
     = CalcStressChanges(TsFaceNormalVectorObs,TnFaceNormalVectorObs,...
    CSSFaceNormalVectorObs,Sxx_Obs_In,Syy_Obs_In,Szz_Obs_In,Sxy_Obs_In,...
    Sxz_Obs_In,Syz_Obs_In,Sxx_Obs_Out,Syy_Obs_Out,Szz_Obs_Out,...
    Sxy_Obs_Out,Sxz_Obs_Out,Syz_Obs_Out,Sxx_Change,Syy_Change,...
    Szz_Change,Sxy_Change,Sxz_Change,Syz_Change,Mu,Cohesion,TsRake,...
    TnRake,CSSRake,FaultFileString);

filename = strcat(FaultFileString,'STEP_7_Variables.mat');
filename2 = strcat('OutputData/',filename);
save(filename2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The 'PlotOutputStresses.m' function plots the new stresses (the stress
% that results from fault slip) at each observation point.  

filename = strcat(FaultFileString,'WorkspaceVariables');
filename2 = strcat('OutputData/',filename);
save(filename2)

PlotOutputStresses(Sxx_Obs_Out,Syy_Obs_Out,Szz_Obs_Out,...
    Sxy_Obs_Out,Sxz_Obs_Out,Syz_Obs_Out,Sxx_Obs_In,Syy_Obs_In,Szz_Obs_In,...
    Sxy_Obs_In,Sxz_Obs_In,Syz_Obs_In,X_Obs, Y_Obs, Z_Obs,FaultPoints,...
    FaultTriangles,FaultFileString);
% 
filename = strcat(FaultFileString,'STEP_7_Variables.mat');
filename2 = strcat('OutputData/',filename);
load (filename2)

StressChangePlaying
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The four functions below plot the CSSchange, Tn change, Ts change, and
%dilatation that occur at the observation points as a result of the slip.
%Results are in percent change for the stresses and percent volume change
%for dilatation.
%% 

[PercentCSSChange] = PlotCSSChange2D(FaultTriangles,FaultPoints,X_Obs,...
    Y_Obs,Z_Obs,CSS_Change,CSS_In,X_cont,Y_cont,FaultFileString);
[PercentTnChange] = PlotTnChange2D(FaultTriangles,FaultPoints,X_Obs,...
    Y_Obs,Z_Obs,Tn_Change,Tn_In,X_cont,Y_cont,FaultFileString);
% [PercentTsChange] = PlotTsChange2D(FaultTriangles,FaultPoints,X_Obs,...
%     Y_Obs,Z_Obs,Ts_Change,Ts_In,X_cont,Y_cont,FaultFileString);
[Dilatation] = PlotDilatation2D(FaultTriangles,FaultPoints,X_Obs,...
    Y_Obs,Z_Obs,Exx,Eyy,Ezz,X_cont,Y_cont,FaultFileString);
[AllSum] = PlotAllNormalized2D(FaultTriangles,FaultPoints,...
    X_Obs,Y_Obs,Z_Obs,Ts_Change,Ts_In,Tn_Change,Tn_In,CSS_Change, ...
    CSS_In,Exx,Eyy,Ezz,X_cont,Y_cont,FaultFileString);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The 'WriteResults.m' function write two tables. The
%'StressChangeModel.csv' reports the x,y,z of the obseration points in the
%percent CSS, Tn, and Ts that occur as a result of the modeled slip.
%'DiliationModel.csv' reports the x,y,z and dilatation at the observation
%points.

WriteResults(PercentCSSChange,PercentTnChange,Dilatation,AllSum,...
   X_Obs,Y_Obs,Z_Obs,FaultFileString);

%load handel  %%% plays the chorus from Handel's 'Messiah' so that 
             %%%you know it is done ;-)
%sound(y,Fs)

end

