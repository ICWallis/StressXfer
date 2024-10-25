function [] = WriteResults(PercentCSSChange,PercentTnChange,Dilatation,AllSum,...
    X_Obs,Y_Obs,Z_Obs,FaultFileString)

%This function writes the final results from stress transfer modeling
%function 'RunAll.m'

%%%INPUTS%%%

%%%PercentCSSChange   The percent change in the Coulomb shear traction relative
%                     to the traction induced by the input stress

%%%PercentTsChange    The percent change in the  shear traction relative
%                     to the traction induced by the input stress

%%%PercentTnChange    The percent change in the normal traction relative
%                     to the traction induced by the input stress

%%%Dilatation         The dilataion at the observation points as
%%%                   a result of the modeled slip

%%%X_Obs,Y_Obs,Z_Obs    One columns vectors containing the coordiates of
%%%                     the observation points. These are calculated in the
%%%                     'LoadObservationPoints' function.

%%%FaultFileString      A string used to append to the front of the file
%%%                     names of saved figures and exported data files

    CSS_ChangeTable = array2table(PercentCSSChange);
   % ShearTractionTable = array2table(PercentTsChange);
    NormalTractionTable = array2table(PercentTnChange);
    AllSumTable = array2table(AllSum);
    X_ObsTable = array2table(X_Obs);
    Y_ObsTable = array2table(Y_Obs);
    Z_ObsTable = array2table(Z_Obs);
    StressChangeModel = [X_ObsTable,Y_ObsTable,Z_ObsTable,CSS_ChangeTable,NormalTractionTable,AllSumTable];
    
    filename = strcat(FaultFileString,'StressChangeModel.csv');
    filename2 = strcat('OutputData/',filename);
    writetable(StressChangeModel, filename2,'WriteVariableNames', true);
    
  
    DilatationTable = array2table(Dilatation);
    DilatationModel = [X_ObsTable,Y_ObsTable,Z_ObsTable,DilatationTable];
    
    filename = strcat(FaultFileString,'DilatationModel.csv');
    filename2 = strcat('OutputData/',filename);
    writetable(DilatationModel, filename2,'WriteVariableNames', true);
    
% %   filename = strcat(FaultFileString,'WorkspaceVariables');
% %   filename2 = strcat('OutputData/',filename);
% %   save(filename2)
end

