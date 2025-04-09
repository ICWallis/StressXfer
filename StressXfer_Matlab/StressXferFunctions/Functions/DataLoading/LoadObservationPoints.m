function [X_Obs,Y_Obs,Z_Obs,ObsPoints2] = LoadObservationPoints(Obsfile,Z_Fault)
%%%This function loads the observation point data
%%% DATA SHOULD BE X,Y,Z IN UTM, Z IN ELEVATION
%%% 'Block Model_72 Pts.csv' is good example
% %Import Observation Points from Brady
  ObsPoints = readtable(Obsfile);
  ObsPoints2 = table2array(ObsPoints);
  rowsToDelete = ObsPoints2(:,3) > max(Z_Fault)+std(Z_Fault);
  
  X_Obs = (ObsPoints2(:,1));
  X_Obs(rowsToDelete) = [];
  Y_Obs = (ObsPoints2(:,2));
  Y_Obs(rowsToDelete) = [];
  Z_Obs = (ObsPoints2(:,3));
  Z_Obs(rowsToDelete) = [];
  
  X_Obs = round(X_Obs,2);
  Y_Obs = round(Y_Obs,2);
  Z_Obs = round(Z_Obs,2); 
end

