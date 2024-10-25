function [Sxx_Fault,Syy_Fault,Szz_Fault,Sxy_Fault,Sxz_Fault,Syz_Fault,Sxx_Obs_In,Syy_Obs_In,Szz_Obs_In,...
    Sxy_Obs_In,Sxz_Obs_In,Syz_Obs_In,mu,nu,Myu,Co] = BuildStressModel(S3_Azimuth,Aphi, Z_Obs, Z_Fault,InputElevation,MeanSurfaceElevation, path,FaultFileString)
%This function builds the stress model data needed for the stress transfer
%modeling function 'RunAll.m'

%%%INPUTS%%%

%%%S3_Azimuth               The azimuth (0-360 clockwise) of the minimum
%                           horizontal stress in Mpa

%%%Aphi                     A numerica indicator of the faulting regime.0.5
%                           is normal faulting, 1.5 is strikeslip faulting,
%                           2.5 is reverse faulting. See (Simpson, 1997) 
%                           for more on Aphi.

%%%Z_Obs                    The Z values of the obeservation points. This
%                           value comes out of the 'LoadData.m' function

%%%Z_Fault                  The Z values of the fault facet midpoints. This
%                           value comes out of the 'LoadData.m' function

%%%InputElevation           The elevation that the inital Sv (vertical
%                           stress) is calculated for in meters.

%%%MeanSurfaceElevation    The mean surface elevation in meters

%%%path                     The path to the folder containing 
%                           the 'RunAll.m' function

%%%OUTPUTS%%%

%%%Sxx_Fault,Syy_Fault,     One column vectors containing the 
%%%Szz_Fault,Sxy_Fault,     magnitudes of the calculated Cartesian input
%%%Sxz_Fault,Syz_Fault      stress tensor for each fault facet.

%%%Sxx_Obs_In,Syy_Obs_In,     One column vectors containing the 
%%%Szz_Obs_In,Sxy_Obs_In,     magnitudes of the calculated Cartesian input
%%%Sxz_Obs_In,Syz_Obs_In      stress tensor for each observation point.

%%%Myu                      Coefficient of friction.

%%%Co                       Cohesive strength. 

%%%mu                       Shear modulus, relates shear stress to shear 
%%%                         strain. 

%%nu                        Poisson's ratio. 

%%%FaultFileString      A string used to append to the front of the file
%%%                     names of saved figures and exported data files


addpath(genpath(path))
Depth = MeanSurfaceElevation-InputElevation; %Depth of interest in meters

%Define a number of needed parameters
ym = 75000 ;        %Young's modulus Mpa (E) (~10-100 Gpa for ingenous 
                    %rocks, so 75,000 Mpa is reasonable)
                    
nu = 0.25;          %Poisson's ratio, Nu or V. Rubber 0.5, Cork 0,
                    %Rock 0.1-0.3;
                    
mu = ym/(2*(1+nu)); %Shear modulus, relates shear stress to shear 
                    %strain. 
                    
Myu = 0.60;         %Coefficient of friction

Co = 0.0;           %Cohesive strength


AphiAll = linspace(0,3,1001)'; %% Initiates a matrix of values for Aphi for
                               %%  plotting

%This block calls the function 'aphi_to_stressmags' written by Jens Lund 
%Snee (Lundstern). The formalae can be found here; 
%https://doi.org/10.1306/08102120151

%This code calculates the stresses for all values of aphi at a given depth
[Shmax_MagAll, Shmin_MagAll, Sv_MagAll, ~] = aphi_to_stressmags(AphiAll, Depth/1000, Myu);
%% first value is aphi, second is depth of interest in km, third is coefficient of friction

L = length(AphiAll);
Sv_MagAll = ones(L,1)*Sv_MagAll;

%This code calculates the stress for the input aphi at a given depth
[S2_Magnitude, S3_Magnitude, S1_Magnitude, ~] = aphi_to_stressmags(Aphi, Depth/1000, Myu);

s1 = num2str( S1_Magnitude );
s2 = num2str( S2_Magnitude );
s3 = num2str( S3_Magnitude );
%Pp1 = 0.45*(S1_Magnitude); %pore pressure in Mpa 
Pp1 = (0.433*3280)/145.038; %pore pressure in Mpa 
PpAll = ones(L,1)*Pp1;

%This code plots and saves a figure of all Aphi values and the Aphi
%values for this the input depth. The figure is saved in the 
%'OutputFigures' folder.

NormalFaulting_InitialStressValues = figure;
scatter(AphiAll,Sv_MagAll,2,'r');
hold on 
scatter(AphiAll,Shmax_MagAll,2,'g');
hold on 
scatter(AphiAll,Shmin_MagAll,2,'b');
hold on
scatter(AphiAll,PpAll,2,'k');
hold on
scatter(Aphi,S1_Magnitude,25,'r');
hold on 
scatter(Aphi,S2_Magnitude,25,'g');
hold on 
scatter(Aphi,S3_Magnitude,25,'b');
hold on
scatter(Aphi,Pp1,25,'k');
legend on
title(['stress at ',num2str(Depth),' m depth; normal faulting']);
legend('S1', 'S2', 'S3','pore pressure', 'Location','southeast');
ylabel('Stress (MPa)');
xlabel('A_{φ}');
xlim([0 3])
ylim([0 85])
daspect([1,65,1]);
str = ['Sv = ' num2str( round(S1_Magnitude,2) )];
text(.02,80,str,'Color','red');
str = ['SHmax = ' num2str( round(S2_Magnitude,2) )];
text(.55,80,str,'Color','green')
str = ['Shmin = ' num2str( round(S3_Magnitude,2) )];
text(1.25,80,str,'Color','blue')
str = ['Pp = ' num2str( round(Pp1,2) )];
text(1.9,80,str,'Color','black')

filename = strcat(FaultFileString,'InitialStressValues');
filename2 = strcat('OutputFigures/',filename);
saveas(NormalFaulting_InitialStressValues,filename2, 'png')

%This code builds, plots and saves it a figure of the vertical stress
%model. The figure is saved in the 'OutputFigures' folder


[Sv_Fault,Shmax_Fault,Shmin_Fault,Sv_Obs,Shmax_Obs,Shmin_Obs,PpObs,Pp_Fault]...
    = Plot3DStressModel(InputElevation,Z_Obs,S1_Magnitude,S2_Magnitude,S3_Magnitude,Pp1,Z_Fault,MeanSurfaceElevation,FaultFileString);


%This saves the vertical stress model in 'OutputData' folder
SvFaultTable = array2table(Sv_Fault);
Shmax_FaultTable = array2table(Shmax_Fault);
Shmin_FaultTable = array2table(Shmin_Fault);
Z_FaultTable = array2table(Z_Fault);
Pp_FaultTable = array2table(Pp_Fault);

StressModel = [Z_FaultTable, SvFaultTable, Shmax_FaultTable, Shmin_FaultTable,Pp_FaultTable];

filename = strcat(FaultFileString,'StressModel.csv');
filename2 = strcat('OutputData/',filename);

writetable(StressModel, filename2,'WriteVariableNames', true);

%The 'PlotStressAxes' function was written by Tim M. Davis as part of
%the 'Cut and Displace' code (https://github.com/Timmmdavis/CutAndDisplace)
%https://doi.org/10.5281/zenodo.3694164. It is slightly adapted for use in 
%this function

% This code plots the axis of the three principal stresses, one for the
% observation points and one for the faults. Only one is saved in the
% OutputFigures folder
[Sxx_Obs_In,Syy_Obs_In,Szz_Obs_In,Sxy_Obs_In,Sxz_Obs_In,Syz_Obs_In] ...
    = PlotStressAxes(Sv_Obs,Shmax_Obs,Shmin_Obs,S3_Azimuth,InputElevation,FaultFileString);

[Sxx_Fault,Syy_Fault,Szz_Fault,Sxy_Fault,Sxz_Fault,Syz_Fault] ...
    = PlotStressAxes(Sv_Fault,Shmax_Fault,Shmin_Fault,S3_Azimuth,InputElevation,FaultFileString);

S1_TP = [0 0];
S2_TP = [0 0];
S3_TP = [0 0];

for i = 1:length(Aphi)
if Aphi(i) <= 1
    S1_TP = array2table([0 -90]);
    S2_TP = array2table([S3_Azimuth-90 0]);
    S3_TP = array2table([S3_Azimuth 0]);
    if S3_Azimuth-90 <0 
        S2_TP = array2table([S3_Azimuth+270 0]);
    end
elseif Aphi(i) > 1 && Aphi(i)<=2
    S2_TP = array2table([0 -90]);
    S1_TP = array2table([S3_Azimuth-90 0]);
    S3_TP = array2table([S3_Azimuth 0]);
    if S3_Azimuth-90 <0 
        S1_TP = array2table([S3_Azimuth+270 0]);
    end
else
    S3_TP = array2table([0 -90]);
    S1_TP = array2table([S3_Azimuth-90 0]);
    S2_TP = array2table([S3_Azimuth 0]);
    if S3_Azimuth-90 <0 
        S1_TP = array2table([S3_Azimuth+270 0]);
    end
end
end


S1_TP = renamevars(S1_TP,["Var1","Var2"],["Trend", "Plunge"]);
S2_TP = renamevars(S2_TP,["Var1","Var2"],["Trend", "Plunge"]);
S3_TP = renamevars(S3_TP,["Var1","Var2"],["Trend", "Plunge"]);
filename = strcat(FaultFileString,'S1_TP_Input.csv');
filename2 = strcat('OutputData/',filename);
writetable(S1_TP, filename2,'WriteVariableNames', true);

filename = strcat(FaultFileString,'S2_TP_Input.csv');
filename2 = strcat('OutputData/',filename);
writetable(S2_TP, filename2,'WriteVariableNames', true);

filename = strcat(FaultFileString,'S3_TP_Input.csv');
filename2 = strcat('OutputData/',filename);
writetable(S3_TP, filename2,'WriteVariableNames', true);
end


