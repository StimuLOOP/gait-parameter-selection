% ---------------------------
% Selection algorithm: z-score
% Match patients and healthy participants speeds
% Take all discrete gait metrics as input
% Mathilde 20.07.2022
% ---------------------------

clear all
close all
clc

type = 'patient';

if strcmp(type,'patient')
    % Define parameters/ file characteristics for patients
    day = {'190913','','190924','191003','191004','191004','191022','191023'};
    speeds = {'FST0.15' 'PGV0.1' ''; '' '' '';'FST0.7' 'PGV0.4' 'SLW0.4'; 'FST0.15' 'PGV0.05' '';'FST0.75' 'PGV0.55' 'SLW0.3'; 'PGV0.25' 'PGV0.5' 'PGV0.8'; 'FST0.35' 'PGV0.20' 'SLW0.07crop'; 'FST0.3' 'PGV0.17' 'SLW0.08'};
    velP = [0.15 0.1 0; 0 0 0; 0.7 0.4 0.4; 0.15 0.05 0; 0.75 0.55 0.3; 0.25 0.5 0.8; 0.35 0.20 0.07; 0.3 0.17 0.08];
    assist ={'_NH_NS' '_NRH_NS1' '';'' '' ''; '_BH_NS' '_BH_NS' '_BH_NS'; '_BH_NS' '_BH_NS' '_BH_NS'; '_BH_NS' '_BH_NS' '_BH_NS'; '_BH_NS' '_BH_NS' '_BH_NS'; '_BH_NS' '_BH_NS' '_BH_NS'; '_BH_NS' '_BH_NS' '_BH_NS'};
    suffixe = '_MoCgapfilled';

elseif strcmp(type,'healthy')
    day = {'200113';'200302';'200306';'200306';'200306';'200306';'200617'};
    speeds = {'0.11'; '0.19';'0.28';'0.36';'0.42';'0.53';'0.61';'0.69';'0.78';'0.86';'0.94';'1.03';'1.11'};
end
velH = [0.11 0.19 0.28 0.36 0.42 0.53 0.61 0.69 0.78 0.86 0.94 1.03 1.11]; % velocities healthy participants

maxM = 0;
for subject = [1 3:8]
    %% Load data
    % Mean healthy
    fileHealthy = 'C:\Users\mlestoille\Documents\DataGait\NM_Reference\statHealthy.mat';
    load(fileHealthy);

    if strcmp(type,'patient')
        % Patient
        if subject == 4
            subjectN = [day{subject},'_S0',num2str(subject),'_T2'];
        else
            subjectN = [day{subject},'_S0',num2str(subject)];
        end
        folder = ['C:\Users\mlestoille\Documents\DataGait\NM_GaitSegmentation\',subjectN,'\04_Visual3D\'];
        filePatient = [folder,'MatlabData\',subjectN,'_parameters'];
        load(filePatient);
        speedSize = size(speeds,2);

    elseif strcmp(type,'healthy')
        velP(subject,:) = velH;
        subjectN = ['REF0',num2str(subject)];
        folder = ['C:\Users\mlestoille\Documents\DataGait\NM_Reference\ReferenceData\Data\',subjectN,'\3_C3D_Files\'];
        fileHealthy = [folder,'MatlabData\',day{subject},'_',subjectN,'_parameters'];
        load(fileHealthy);
        speedSize = size(speeds,1);
    end
    dimName = {' /X','Y','Z'};

    %% Zscore and selection
    for speedN = 1:speedSize
        if velP(subject,speedN) == 0
            % do nothing
        else
            % find the corresponding velocity of healthy participants (closest one)
            if strcmp(type,'patient')
                [minDistance,indOfMin] = min(abs(velP(subject,speedN)-velH));
            elseif strcmp(type,'healthy')
                indOfMin = speedN;
            end
            % compute z-score for all metrics
            fields = fieldnames(param{speedN});
            for k = 1:size(fields,1)
                s = size(meanSubject{speedN}.(fields{k}),2);
                zsc(k,1:s) = (meanSubject{speedN}.(fields{k})-meanHealthy{indOfMin}.(fields{k}))./stdHealthy{indOfMin}.(fields{k});
                if s<3
                    zsc(k,s+1:3) = NaN; % put nan when there is only one dimension
                end
            end
            % Find 4 max values
            % find in each column
            [M{speedN,subject},indOfMax] = maxk(zsc,4);
            % gather in one column to find the 4 max of the whole zscore matrix
            M2 = []; ind = [];
            for i = 1:3
                M2 = [M2;M{speedN,subject}(:,i)];
                ind = [ind; indOfMax(:,i)];
            end
            [maxZscore,indOfMaxBis] = maxk(M2,4);
            numImpParam = ind(indOfMaxBis);
            dim = fix((indOfMaxBis-1)/4)+1;
            impParam{subject}(:,2*(speedN-1)+1:2+2*(speedN-1)) = [fields(numImpParam) dimName(dim)'];

            % find max zscore over all healthy participants
            if strcmp(type,'healthy')
                if maxM < max(M2)
                    maxM = max(M2);
                end
            end
            
            % store all values in single array
            zscore{speedN,subject} = zsc; 
        end
    end
end