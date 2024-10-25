function [] = StressChangeStereonet(S1_Trend,S1_Plunge,S2_Trend,S2_Plunge,S3_Trend,S3_Plunge,S1dirOut,S2dirOut,S3dirOut,FaultFileString)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
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
filename = strcat(FaultFileString,'S1TP.csv');
filename2 = strcat('OutputData/',filename);
writetable(S1TP, filename2,'WriteVariableNames', true);
filename = strcat(FaultFileString,'S2TP.csv');
filename2 = strcat('OutputData/',filename);
writetable(S2TP, filename2,'WriteVariableNames', true);
filename = strcat(FaultFileString,'S3TP.csv');
filename2 = strcat('OutputData/',filename);
writetable(S3TP, filename2,'WriteVariableNames', true);


trend1 = table2array(S1TP(:,1));
plunge1 = table2array(S1TP(:,2));
num = height(S1TP(:,1));
R = 1;
trendr1 = trend1*pi/180;
plunger1 = plunge1(:,1)*pi/180;
rho1 = R.*tan(pi/4 - ((plunger1)/2));
figure(1)
for i=1:num
% polarb plots ccl from 3:00, so convert to cl from 12:00
polarplot(pi/2-trendr1(i),rho1(i),'o')
hold on
end