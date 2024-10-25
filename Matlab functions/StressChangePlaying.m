
% filename = strcat(FaultFileString,'WorkSpaceVariables');
%load('OutputData\CarversCalcStressChangesWorkspaceVariables')

[S1In,S2In,S3In,S1dirIn,S2dirIn,S3dirIn,S1Out,S2Out,S3Out,S1dirOut,S2dirOut,S3dirOut,S1,S2,S3,S1dir,S2dir,S3dir] = PlotOutputStresses(Sxx_Obs_Out,Syy_Obs_Out,Szz_Obs_Out,...
    Sxy_Obs_Out,Sxz_Obs_Out,Syz_Obs_Out,Sxx_Obs_In,Syy_Obs_In,Szz_Obs_In,...
    Sxy_Obs_In,Sxz_Obs_In,Syz_Obs_In,X_Obs, Y_Obs, Z_Obs,FaultPoints,...
    FaultTriangles,FaultFileString);

S1Plunge = -1*asind(S1dirOut(:,3));
S1Trend(1:(length(S1dirOut)),1:1)=0;
for i = 1:length(S1dirOut)
    if S1dirOut(i,1) == 0 && S1dirOut(i,2) == 0  
        S1Trend(i) = 0; 
    else
        if S1dirOut(i,1) == 0
           S1Trend(i) = asind(0);
        elseif S1dirOut(i,2) == 0
           S1Trend(i) = asind(1);
        else
            if S1dirOut(i,2) > 0
            S1Trend(i) = atand(S1dirOut(i,1)./S1dirOut(i,2));
            else
            S1Trend(i) = atand(S1dirOut(i,1)./S1dirOut(i,2))+180;  
            end
        end
    end
end
S2Plunge = -1*asind(S2dirOut(:,3));
S2Trend(1:(length(S2dirOut)),1:1)=0;
for i = 1:length(S2dirOut)
    if S2dirOut(i,1) == 0 && S2dirOut(i,2) == 0  
        S2Trend(i) = 0; 
    else
        if S2dirOut(i,1) == 0
           S2Trend(i) = asind(0);
        elseif S2dirOut(i,2) == 0
           S2Trend(i) = asind(1);
        else
            if S2dirOut(i,2) > 0
            S2Trend(i) = atand(S2dirOut(i,1)./S2dirOut(i,2));
            else
            S2Trend(i) = atand(S2dirOut(i,1)./S2dirOut(i,2))+180;  
            end
        end
    end
end
S3Plunge = -1*asind(S3dirOut(:,3));
S3Trend(1:(length(S3dirOut)),1:1)=0;
for i = 1:length(S3dirOut)
    if S3dirOut(i,1) == 0 && S3dirOut(i,2) == 0  
        S3Trend(i) = 0; 
    else
        if S3dirOut(i,1) == 0
           S3Trend(i) = asind(0);
        elseif S3dirOut(i,2) == 0
           S3Trend(i) = asind(1);
        else
            if S3dirOut(i,2) > 0
            S3Trend(i) = atand(S3dirOut(i,1)./S3dirOut(i,2));
            else
            S3Trend(i) = atand(S3dirOut(i,1)./S3dirOut(i,2))+180;  
            end
        end
    end
end
S1TP = [S1Trend, S1Plunge];
S2TP = [S2Trend, S2Plunge];
S3TP = [S3Trend, S3Plunge];
S1TP = array2table(S1TP); S2TP = array2table(S2TP); S3TP = array2table(S3TP); 
S1TP = renamevars(S1TP,["S1TP1","S1TP2"],["Trend", "Plunge"]);
S2TP = renamevars(S2TP,["S2TP1","S2TP2"],["Trend", "Plunge"]);
S3TP = renamevars(S3TP,["S3TP1","S3TP2"],["Trend", "Plunge"]);
filename = strcat(FaultFileString,'PostSlip_S1_TP.csv');
filename2 = strcat('OutputData/',filename);
writetable(S1TP, filename2,'WriteVariableNames', true);
filename = strcat(FaultFileString,'PostSlip_S2_TP.csv');
filename2 = strcat('OutputData/',filename);
writetable(S2TP, filename2,'WriteVariableNames', true);
filename = strcat(FaultFileString,'PostSlip_S3_TP.csv');
filename2 = strcat('OutputData/',filename);
writetable(S3TP, filename2,'WriteVariableNames', true);
S1_Trend = 0;
S1_Plunge = -90;
S2_Trend = 0;
S2_Plunge = 0;
S3_Trend = 90;
S3_Plunge = 0;