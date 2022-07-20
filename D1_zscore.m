% ---------------------------
% Selection algorithm: z-score
% Match patients and healthy participants speeds
% Take all discrete gait mtrics as input
% Mathilde 20.07.2022
% ---------------------------

clear all
close all
clc

% Define parameters/ file characteristics for patients
day = {'190913','','190924','191003','191004','191004','191022','191023'};
speeds = {'FST0.15' 'PGV0.1' ''; '' '' '';'FST0.7' 'PGV0.4' 'SLW0.4'; 'FST0.15' 'PGV0.05' '';'FST0.75' 'PGV0.55' 'SLW0.3'; 'PGV0.25' 'PGV0.5' 'PGV0.8'; 'FST0.35' 'PGV0.20' 'SLW0.07crop'; 'FST0.3' 'PGV0.17' 'SLW0.08'};
velP = [0.15 0.1 0; 0 0 0; 0.7 0.4 0.4; 0.15 0.05 0; 0.75 0.55 0.3; 0.25 0.5 0.8; 0.35 0.20 0.07; 0.3 0.17 0.08];
assist ={'_NH_NS' '_NRH_NS1' '';'' '' ''; '_BH_NS' '_BH_NS' '_BH_NS'; '_BH_NS' '_BH_NS' '_BH_NS'; '_BH_NS' '_BH_NS' '_BH_NS'; '_BH_NS' '_BH_NS' '_BH_NS'; '_BH_NS' '_BH_NS' '_BH_NS'; '_BH_NS' '_BH_NS' '_BH_NS'};
suffixe = '_MoCgapfilled';

velH = [0.11 0.19 0.28 0.36 0.42 0.53 0.61 0.69 0.78 0.86 0.94 1.03 1.11]; % velocities healthy participants

subject = 1;

% Load data
% Mean healthy
fileHealthy = 'D:\StimuLOOP\DataGait\NM_Reference\statHealthy.mat';
load(fileHealthy);
% Patient
if subject == 4
    subjectN = [day{subject},'_S0',num2str(subject),'_T2'];
else
    subjectN = [day{subject},'_S0',num2str(subject)];
end
folder = ['D:\StimuLOOP\DataGait\NM_GaitSegmentation\',subjectN,'\04_Visual3D\'];
filePatient = [folder,'MatlabData\',subjectN,'_parameters'];
load(filePatient);
for speedN = 1:size(speeds,2)
    if velP(subject,speedN) == 0
        % do nothing
    else
        % find the corresponding velocity of healthy participants (closest one)
        [minDistance,indOfMin] = min(abs(velP(subject,speedN)-velH));
        % compute z-score for all metrics
        fields = fieldnames(param{speedN});
        for k = 1:size(fields,1)
            s = size(meanSubject{speedN}.(fields{k}),2);
            zsc(k,1:s) = (meanSubject{speedN}.(fields{k})-meanHealthy{indOfMin}.(fields{k}))./stdHealthy{indOfMin}.(fields{k});
            zsc(k,s:3) = NaN; % put nan when there is only one dimension
        end
    end
end
