% -------------------
% 1. Take a c3d file (not necessarily processed by visual3D)
% 2. Use Chris' codes to detect events based on Zeni's method
% 3. Check the events
% Mathilde, 16.03.2022
% ------------------

clear all
close all
clc
addpath('btk');
addpath('ChrisFunctions');
% Choose between healthy and patient (different file characteristics)
type = 'healthy'; 

plotFig = 0;  % choose or not to plot the figures
saveFig = 0; % choose or not to save the figures

% Define parameters/ file characteristics
k = 0.75; % threshold for detecting too short time between similar events
if strcmp(type,'healthy')
    day = {'200113';'200302';'200306';'200306';'200306';'200306';'200617'};
    speeds = {'0.11'; '0.19';'0.28';'0.36';'0.42';'0.53';'0.61';'0.69';'0.78';'0.86';'0.94';'1.03';'1.11'};
    suffixe = {'_MoC'; '_MoC';'';'';'_gapfilled';'';''};
elseif strcmp(type,'patient')
    day = {'190913','','190924','191003','191004','191004','191022','191023'};
    speeds = {'FST0.15' 'PGV0.1' ''; '' '' '';'FST0.7' 'PGV0.4' 'SLW0.4'; 'FST0.15' 'PGV0.05' '';'FST0.75' 'PGV0.55' 'SLW0.3'; 'PGV0.5' 'PGV0.8' 'PGV0.25'; 'FST0.35' 'PGV0.20' 'SLW0.07crop'; 'FST0.3' 'PGV0.17' 'SLW0.08'};
    assist ={'_NH_NS' '_NRH_NS1' '';'' '' ''; '_BH_NS' '_BH_NS' '_BH_NS'; '_BH_NS' '_BH_NS' '_BH_NS'; '_BH_NS' '_BH_NS' '_BH_NS'; '_BH_NS' '_BH_NS' '_BH_NS'; '_BH_NS' '_BH_NS' '_BH_NS'; '_BH_NS' '_BH_NS' '_BH_NS'};
    suffixe = '_MoCgapfilled';
end

for subject = [1 3:7]%[1 3:8]
    if strcmp(type,'healthy')
        if subject < 10
            subjectN = ['REF0', num2str(subject)];
        else
            subjectN = ['REF', num2str(subject)];
        end
        folder = ['D:\StimuLOOP\DataGait\NM_Reference\ReferenceData\Data\',subjectN,'\3_C3D_Files\'];
        % create a folder to store the Matlab figures
        if not(isfolder([folder,'MatlabData']))
            mkdir(fullfile(folder, 'MatlabData'))
        end
        fileMatlab = [folder,'\MatlabData\',day{subject},'_',subjectN,'_TM_RawEvents'];
        %     folder = ['C:\Users\Mathilde\Documents\StimuLOOP\DataGaitSamples\',subjectN,'\'];
    elseif strcmp(type,'patient')
        if subject == 4
            subjectN = [day{subject},'_S0',num2str(subject),'_T2'];
        else
            subjectN = [day{subject},'_S0',num2str(subject)];
        end
        folder = ['D:\StimuLOOP\DataGait\NM_GaitSegmentation\',subjectN,'\04_Visual3D\'];
        % create a folder to store the Matlab figures
        if not(isfolder([folder,'MatlabData']))
            mkdir(fullfile(folder, 'MatlabData'))
        end
        fileMatlab = [folder,'\MatlabData\',subjectN,'_TM_RawEvents'];
    end
    close all

    for speedN = 1:size(speeds,1)
        clearvars -except k folder subjectN speeds speedN assist ...
            day fileMatlab timepb subject suffixe plotFig saveFig events type
        if strcmp(type,'healthy')
            filename = [folder,day{subject},'_',subjectN,'_TM_',speeds{speedN},'_NH_NS',suffixe{subject}];
            figureMatlabh = [folder,'\MatlabData\',day{subject},'_',subjectN,'_',speeds{speedN},'_events.fig'];
        elseif strcmp(type,'patient')
            if isempty(speeds{subject, speedN})
                break
            end
            filename = [folder,subjectN,'_TM_',speeds{subject,speedN},assist{subject,speedN},suffixe];
            figureMatlabh = [folder,'\MatlabData\',day{subject},'_',subjectN,'_TM_',speeds{subject,speedN},assist{subject,speedN},'_events.fig'];
        end
        c3d = btkReadAcquisition([filename,'.c3d']);
        % Add subjects Meta data
        info.format = 'Char';
        info.values = {subjectN};
        btkAppendMetaData(c3d,'SUBJECTS','NAMES',info);
        info.format = 'Char';
        if strcmp(type,'healthy')
            info.values = speeds(speedN);
        elseif strcmp(type,'patient')
            info.values = speeds(subject,speedN);
        end
        btkAppendMetaData(c3d,'SUBJECTS','SPEED',info);
        % Save metadata
        btkWriteAcquisition(c3d, [filename,'.c3d']);
        % Filtering data
        c3d = filteringMarkersData(c3d);

        % Gait event detection. Output = structure with four matrices
        [success,events{speedN},h(1)] = setKinematicEvents(c3d, 50, 100, 2, 2, plotFig);

        % Check gait events
        RFS = events{speedN}.RStrike;
        RFO = events{speedN}.ROff;
        LFS = events{speedN}.LStrike;
        LFO = events{speedN}.LOff;
        Labels = [repmat({'RFO'},size(RFO,1),1);repmat({'RFS'},size(RFS,1),1);...
            repmat({'LFO'},size(LFO,1),1);repmat({'LFS'},size(LFS,1),1)];
        GE_table = table([RFO; RFS; LFO; LFS],Labels);
        GE_table_sorted = sortrows(GE_table);

        z1 = [RFS RFS];  z2 = [RFO RFO];
        v1 = [LFS LFS];  v2 = [LFO LFO];
        a = [0; 1]; b = [1 2];

        if plotFig
            h(2) = figure(2);
            hold on;
            title(['Gait events, ', char(subjectN),' Speed ', char(speeds(speedN))])
            plot(z1.',a, 'b');  plot(z2.',a, 'r');
            plot(v1.',b, 'g');  plot(v2.',b, 'c'); % LOFF/LTO
        end

        % save figures for later manual correction
        if plotFig && saveFig
            savefig(h,figureMatlabh);
        end

        % Close btk handle
        btkCloseAcquisition(c3d);

    end
    % save data in matlab file for manual correction
    save([fileMatlab,'.mat'],'events');
    close all
    clear events
end