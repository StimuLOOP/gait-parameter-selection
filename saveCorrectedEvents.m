% -------------------
% Gather parameters + gait events in one matlab file
% Mathilde, 10.05.2022
% ------------------

clear all
close all
clc
addpath('btk');
addpath('ChrisFunctions');
% Choose between healthy and patient (different file characteristics)
type = 'patient'; 

% Define parameters/ file characteristics
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

for subject = [1 3:8]  %[1 3:7]
    if strcmp(type,'healthy')
        if subject < 10
            subjectN = ['REF0', num2str(subject)];
        else
            subjectN = ['REF', num2str(subject)];
        end
        folder = ['D:\StimuLOOP\DataGait\NM_Reference\ReferenceData\Data\',subjectN,'\3_C3D_Files\'];

         % open Matlab files with gait events detected in B1
        fileMatlab = [folder,'MatlabData\',day{subject},'_',subjectN,'_TM_RawEvents.mat'];
        load(fileMatlab);
        fileMatlabAll = [folder,'MatlabData\',day{subject},'_',subjectN,'_parameters.mat'];
    elseif strcmp(type,'patient')
        if subject == 4
            subjectN = [day{subject},'_S0',num2str(subject),'_T2'];
        else
            subjectN = [day{subject},'_S0',num2str(subject)];
        end
        folder = ['D:\StimuLOOP\DataGait\NM_GaitSegmentation\',subjectN,'\04_Visual3D\'];
         % open Matlab files with gait events detected in B1
        fileMatlab = [folder,'MatlabData\',subjectN,'_TM_RawEvents.mat'];
        load(fileMatlab);
        fileMatlabAll = [folder,'MatlabData\',subjectN,'_parameters.mat'];
    end

        if strcmp(type,'healthy')
             sizeSpeed = size(speeds,1);
        elseif strcmp(type,'patient')
            sizeSpeed = size(speeds,2);
        end

        % Save events before manual correction
        rawEvents = events';
        clear events

    for speedN = 1:sizeSpeed
        % Load c3d file with corrected gait events
         if strcmp(type,'healthy')
            filename = [folder,day{subject},'_',subjectN,'_TM_',speeds{speedN},'_NH_NS',suffixe{subject},'_EventsAdded'];
        elseif strcmp(type,'patient')
            if isempty(speeds{subject, speedN})
                break
            end
            filename = [folder,subjectN,'_TM_',speeds{subject,speedN},assist{subject,speedN},suffixe,'_EventsAdded'];
         end

        c3d = btkReadAcquisition([filename,'.c3d']);
        evts = btkGetEvents(c3d);
        freq = btkGetPointFrequency(c3d);
        ff = btkGetFirstFrame(c3d);
        p = btkGetMetaData(c3d, 'SUBJECTS', 'NAMES');
        subjectName = p.info.values;

        % Transform from events in s to events in frame number
        events{speedN,1}.LOff = round(evts.Left_Foot_Off*freq - ff);
        events{speedN,1}.LStrike = round(evts.Left_Foot_Strike*freq - ff);
        events{speedN,1}.ROff = round(evts.Right_Foot_Off*freq - ff);
        events{speedN,1}.RStrike = round(evts.Right_Foot_Strike*freq - ff);

        % Close btk handle
        btkCloseAcquisition(c3d);
    end
    save(fileMatlabAll,'events','rawEvents','-append');
end