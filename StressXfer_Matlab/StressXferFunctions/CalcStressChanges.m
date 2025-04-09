function [Ts_In,Ts_Out,Ts_Change,Tn_In,Tn_Out,Tn_Change,CSS_In,CSS_Out,CSS_Change] =...
    CalcStressChanges(TsFaceNormalVectorObs,TnFaceNormalVectorObs,CSSFaceNormalVectorObs...
    ,Sxx_Obs_In,Syy_Obs_In,Szz_Obs_In,Sxy_Obs_In,Sxz_Obs_In,Syz_Obs_In,...
    Sxx_Obs_Out,Syy_Obs_Out,Szz_Obs_Out,Sxy_Obs_Out,Sxz_Obs_Out,Syz_Obs_Out,...
    Sxx_Change,Syy_Change,Szz_Change,Sxy_Change,Sxz_Change,Syz_Change,...
    Mu,Cohesion,TsRake,TnRake,CSSRake,FaultFileString)

%%%This function calculates the input,output,and change in the shear
%%%traction, normal traction and Coulomb shear traction at each observation
%%%point. This is part of the stress transfer modeling function 'RunAll.m'

%%%INPUTS%%%

%%%TsFaceNormalVectorObs    The normals to the surfaces that result in the
%%%TnFaceNormalVectorObs    maximum shear traction increase, normal
%%%CSSFaceNormalVectorObs   traction decrease, and Coulomb shear traction
%%%                         increase at each observation point

%%%Sxx_Obs_Change,Syy_Obs_Change,     One column vectors containing the 
%%%Szz_Obs_Change,Sxy_Obs_Change,     magnitudes of the calculated Cartesian
%%%Sxz_Obs_Change,Syz_Obs_Change      stress tensor change, that is the stress 
%%%                                   difference between the input and output
%%%                                   stress,for each observation point.
%%%                                   These are calculated in the 'Stress
%%%                                   DisplacementStrain' function

%%%Sxx_Obs_Out,Syy_Obs_Out,     One column vectors containing the 
%%%Szz_Obs_Out,Sxy_Obs_Out,     magnitudes of the calculated Cartesian
%%%Sxz_Obs_Out,Syz_Obs_Out      output stress tensor, that is the stress 
%%%                             that results from the slip,for each
%%%                             observation point.

%%%Sxx_Obs_In,Syy_Obs_In,     One column vectors containing the 
%%%Szz_Obs_In,Sxy_Obs_In,     magnitudes of the calculated Cartesian input
%%%Sxz_Obs_In,Syz_Obs_In      stress tensor for each observation point. 
%%%                           These vectors comes out of the 
%%%                           'BuildStressModel' function

%%%Mu                      Coefficient of friction. This contant is set in
%%%                        the 'BuildStressModel' function, here it is a
%%%                        vector with length = number of observation
%%%                        points.

%%%Cohesion                Cohesive strength. This contant is set in
%%%                        the 'BuildStressModel' function, here it is a
%%%                        vector with length = number of observation
%%%                        points.

%%%TsRake,TnRake,CSSRake    The rakes that result in the
%%%                         maximum shear traction increase, normal
%%%                         traction decrease, and Coulomb shear traction
%%%                         increase at each observation point

%%%OUTPUTS%%%
%%%Ts_In,Ts_Out,Ts_Change      The input, output, and change in the shear
%%%                            traction at the observation points 
%%%                            calculated for the surfaces defined by 
%%%                            TsFaceNormalVectorObs

%%%Tn_In,Tn_Out,Tn_Change      The input, output, and change in the normal
%%%                            traction at the observation points 
%%%                            calculated for the surfaces defined by 
%%%                            TnFaceNormalVectorObs

%%%CSS_In,CSS_Out,CSS_Change   The input, output, and change in the Coulomb
%%%                            shear traction at the observation points 
%%%                            calculated for the surfaces defined by 
%%%                            CSSFaceNormalVectorObs

%The function relies hevily on functions written by Tim M. Davis 
%which are part of the 'Cut and Displace' code 
%(https://github.com/Timmmdavis/CutAndDisplace) 
%https://doi.org/10.5281/zenodo.3694164

% % This code calculates the shear, normal, Coulomb shear, stirke-slip shear,
% % and dip-slip shear stress  on the observation points in initial stress state
[Ts_In, ~, ~, ~, ~ ] = CalculateCoulombStressOnPlane...
(TsFaceNormalVectorObs,Sxx_Obs_In,Syy_Obs_In,Szz_Obs_In,...
Sxy_Obs_In,Sxz_Obs_In,Syz_Obs_In,Mu,Cohesion,TsRake);

% This code calculates the shear, normal, Coulomb shear, stirke-slip shear,
% and dip-slip shear stress  on the observation points in initial stress state
[Ts_Out,~, ~, ~, ~ ] = CalculateCoulombStressOnPlane...
(TsFaceNormalVectorObs,Sxx_Obs_Out,Syy_Obs_Out,Szz_Obs_Out,...
Sxy_Obs_Out,Sxz_Obs_Out,Syz_Obs_Out,Mu,Cohesion,TsRake);

[ Ts_Change, ~,~, ~, ~ ] = CalculateCoulombStressOnPlane...
(TsFaceNormalVectorObs,Sxx_Change,Syy_Change,Szz_Change,...
Sxy_Change,Sxz_Change,Syz_Change,Mu,Cohesion,TsRake);

% % This code calculates the shear, normal, Coulomb shear, stirke-slip shear,
% % and dip-slip shear stress  on the observation points in initial stress state
[~, Tn_In, ~, ~, ~ ] = CalculateCoulombStressOnPlane...
(TnFaceNormalVectorObs,Sxx_Obs_In,Syy_Obs_In,Szz_Obs_In,...
Sxy_Obs_In,Sxz_Obs_In,Syz_Obs_In,Mu,Cohesion,TnRake);

% This code calculates the shear, normal, Coulomb shear, stirke-slip shear,
% and dip-slip shear stress  on the observation points in initial stress state
[ ~, Tn_Out, ~, ~, ~ ] = CalculateCoulombStressOnPlane...
(TnFaceNormalVectorObs,Sxx_Obs_Out,Syy_Obs_Out,Szz_Obs_Out,...
Sxy_Obs_Out,Sxz_Obs_Out,Syz_Obs_Out,Mu,Cohesion,TnRake);

[ ~, Tn_Change, ~, ~, ~ ] = CalculateCoulombStressOnPlane...
(TnFaceNormalVectorObs,Sxx_Change,Syy_Change,Szz_Change,...
Sxy_Change,Sxz_Change,Syz_Change,Mu,Cohesion,TnRake);
% % This code calculates the shear, normal, Coulomb shear, stirke-slip shear,
% % and dip-slip shear stress  on the observation points in initial stress state
[ ~, ~, CSS_In, ~, ~ ] = CalculateCoulombStressOnPlane...
(CSSFaceNormalVectorObs,Sxx_Obs_In,Syy_Obs_In,Szz_Obs_In,...
Sxy_Obs_In,Sxz_Obs_In,Syz_Obs_In,Mu,Cohesion,CSSRake);

% This code calculates the shear, normal, Coulomb shear, stirke-slip shear,
% and dip-slip shear stress  on the observation points in initial stress state
[ ~, ~,CSS_Out, ~, ~ ] = CalculateCoulombStressOnPlane...
(CSSFaceNormalVectorObs,Sxx_Obs_Out,Syy_Obs_Out,Szz_Obs_Out,...
Sxy_Obs_Out,Sxz_Obs_Out,Syz_Obs_Out,Mu,Cohesion,CSSRake);

[ ~, ~,CSS_Change, ~, ~ ] = CalculateCoulombStressOnPlane...
(CSSFaceNormalVectorObs,Sxx_Change,Syy_Change,Szz_Change,...
Sxy_Change,Sxz_Change,Syz_Change,Mu,Cohesion,CSSRake);

% filename = strcat(FaultFileString,'CalcStressChangesWorkspaceVariables');
% filename2 = strcat('OutputData/',filename);
% save(filename2)

end

