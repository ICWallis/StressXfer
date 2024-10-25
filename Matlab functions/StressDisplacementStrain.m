function [Sxx_Change,Syy_Change,Szz_Change,Sxy_Change,Sxz_Change,Syz_Change,Sxx_Obs_Out,Syy_Obs_Out,...
    Szz_Obs_Out,Sxy_Obs_Out,Sxz_Obs_Out,Syz_Obs_Out,Exx,Eyy,Ezz] = StressDisplacementStrain(Dss,Dds,Dn,mu,X_Obs,Y_Obs,Z_Obs,...
Sxx_Obs_In,Syy_Obs_In,Szz_Obs_In,Sxy_Obs_In,Sxz_Obs_In,Syz_Obs_In,P1,P2,P3,halfspace,nu,lambda,path,FaultFileString,FaultTriangles, FaultPoints)

%This function calculates the stress changes and strain that result from
%the modelled slip. These are needed for the stress transfer modeling
%function 'RunAll.m'

%%%INPUTS%%%

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

%%%mu                   Shear modulus, relates shear stress to shear 
%%%                     strain. This contant is set in
%%%                     the 'BuildStressModel' function

%%%X_Obs,Y_Obs,Z_Obs    One columns vectors containing the coordiates of
%%%                     the observation points. These are calculated in the
%%%                     'LoadObservationPoints' function

%%%Sxx_Obs_In,Syy_Obs_In,     One column vectors containing the 
%%%Szz_Obs_In,Sxy_Obs_In,     magnitudes of the calculated Cartesian input
%%%Sxz_Obs_In,Syz_Obs_In      stress tensor for each observation point. 
%%%                           These vectors comes out of the 
%%%                           'BuildStressModel' function

%%%P1,P2,P3                 A three Column vector where each 'P' represents
%%%                         the different corner points of one of the
%%%                         triangles (XYZ). Not as efficient in terms of
%%%                         storage but easier to understand. These
%%%                         vectors comes out of the 'LoadFaultFile'
%%%                         function

%%%halfspace                This constant is 0, set in the CalcSlip
%%%                         function

%%nu                        Poisson's ratio. This contant is set in
%%%                         the 'BuildStressModel' function

%%lambda                    Lame's constant. This contant is calculated in
%%%                         the 'BuildStressModel' function

%%%path                     The path to the folder containing 
%%%                         the 'RunAll.m' function

%%%FaultFileString      A string used to append to the front of the file
%%%                     names of saved figures and exported data files

%%%OUTPUTS%%%

%%%Sxx_Obs_Change,Syy_Obs_Change,     One column vectors containing the 
%%%Szz_Obs_Change,Sxy_Obs_Change,     magnitudes of the calculated Cartesian
%%%Sxz_Obs_Change,Syz_Obs_Change      stress tensor change, that is the stress 
%%%                                   difference between the input and output
%%%                                   stress,for each observation point.

%%%Sxx_Obs_Out,Syy_Obs_Out,     One column vectors containing the 
%%%Szz_Obs_Out,Sxy_Obs_Out,     magnitudes of the calculated Cartesian
%%%Sxz_Obs_Out,Syz_Obs_Out      output stress tensor, that is the stress 
%%%                             that results from the slip,for each
%%%                             observation point.

%%%Exx,Eyy,Ezz                  Strain change induced by movement of 
%%%                             dislocations.

addpath(genpath(path))

%This function relies hevily on functions written by Tim M. Davis 
%which are part of the 'Cut and Displace' code 
%(https://github.com/Timmmdavis/CutAndDisplace) 
%https://doi.org/10.5281/zenodo.3694164

%%% Calculating the stress changes and strains that occur at each
%%% observation point based on the slip and the initial stresses
[StressTTotal,StrainTTotal,StressTChg,StrainTChg,StressTReg,StrainTReg]=...
CalculateStressOnSurroundingPoints3d(Dss,Dds,Dn,mu,lambda,X_Obs,Y_Obs,Z_Obs,...
Sxx_Obs_In,Syy_Obs_In,Szz_Obs_In,Sxy_Obs_In,Sxz_Obs_In,Syz_Obs_In,P1,P2,P3,halfspace,nu);
close all

[Exx,Eyy,Ezz,Exy,Exz,Eyz ] = ExtractCols( StrainTChg );

%% Calculating the displacement that occurs on each of the observation
%%% points as a result of the slip
 [Ux,Uy,Uz] = CalculateDisplacementOnSurroundingPoints3d...
 (Dss,Dds,Dn,nu,X_Obs,Y_Obs,Z_Obs,P1,P2,P3,halfspace);

% Calculating the resultant stress vectors that are imparted on each
% observation points as a result of the slip
[Sxx_Obs_Out,Syy_Obs_Out,Szz_Obs_Out,Sxy_Obs_Out,Sxz_Obs_Out,Syz_Obs_Out]...
    = ExtractCols( StressTTotal );

% Calculating the change between the initial stress vecotrs and the output
% stress vectors at each observation points
[Sxx_Change,Syy_Change,Szz_Change,Sxy_Change,Sxz_Change,Syz_Change]...
    = ExtractCols( StressTChg );

%Converts output stress tensors on observation points to principal 
% directions and magnitudes
[S1_Obs_Out,S2_Obs_Out,S3_Obs_Out,S1dir_Obs_Out,S2dir_Obs_Out,S3dir_Obs_Out]...
    = EigCalc3d(Sxx_Obs_Out,Syy_Obs_Out,Szz_Obs_Out,...
    Sxy_Obs_Out,Sxz_Obs_Out,Syz_Obs_Out);
% 
%Converts input stress tensors on observation points to principal
% directions and magnitudes
[S1_Obs_In,S2_Obs_In,S3_Obs_In,S1dir_Obs_In,S2dir_Obs_In,S3dir_Obs_In]...
    = EigCalc3d(Sxx_Obs_In,Syy_Obs_In,Szz_Obs_In,...
    Sxy_Obs_In,Sxz_Obs_In,Syz_Obs_In);
% a check
[ Uniform ] = UniformGridCheck3d( X_Obs,Y_Obs,Z_Obs );

% filename = strcat(FaultFileString,'StressDisplacementStrainWorkspaceVariables');
% filename2 = strcat('OutputData/',filename);
% save(filename2)
end

