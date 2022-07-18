% -------------------
% 1. Correct manually the events detected in B1 (with Zeni's method),
% through a GUI
% 3. Check the corrected events
% Mathilde, 30.03.2022
% ------------------

clear all
close all
clc
addpath('btk');
addpath('ChrisFunctions');
% Choose between healthy and patient (different file characteristics)
type = 'patient'; 

insert = @(a, x, n)cat(2,  x(1:n), a, x(n+1:end)); % define insert functions
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

for subject = 6%3:7  %[1 3:7]
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
        %     folder = ['C:\Users\Mathilde\Documents\StimuLOOP\DataGaitSamples\',subjectN,'\'];

        % open Matlab files with gait events detected in B1
        fileMatlab = [folder,'MatlabData\',day{subject},'_',subjectN,'_TM_rawEvents.mat'];
        newfileMatlab = [folder,'MatlabData\',day{subject},'_',subjectN,'_TM_Events.mat'];
        load(fileMatlab);
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

        % open Matlab files with gait events detected in B1
        fileMatlab = [folder,'MatlabData\',subjectN,'_TM_rawEvents.mat'];
        newfileMatlab = [folder,'MatlabData\',subjectN,'_TM_Events.mat'];
        load(fileMatlab);
    end


        % Open the app
        app = correctEventApp;
        if strcmp(type,'healthy')
             sizeSpeed = size(speeds,1);
        elseif strcmp(type,'patient')
            sizeSpeed = size(speeds,2);
        end
    for speedN = 1:sizeSpeed
        % Load c3d file
         if strcmp(type,'healthy')
            filename = [folder,day{subject},'_',subjectN,'_TM_',speeds{speedN},'_NH_NS',suffixe{subject}];
        elseif strcmp(type,'patient')
            if isempty(speeds{subject, speedN})
                break
            end
            filename = [folder,subjectN,'_TM_',speeds{subject,speedN},assist{subject,speedN},suffixe];
            end
        newfilename = [filename,'_EventsAdded.c3d'];

        c3d = btkReadAcquisition([filename,'.c3d']);
        evts = btkGetEvents(c3d);
        freq = btkGetPointFrequency(c3d);
        ff = btkGetFirstFrame(c3d);
        p = btkGetMetaData(c3d, 'SUBJECTS', 'NAMES');
        subjectName = p.info.values;

        % Transform from events in s to events in frame number (easier to
        % work with with the figures)
        Left_Off = events{speedN}.LOff;
        Left_Strike = events{speedN}.LStrike;
        Right_Off = events{speedN}.ROff;
        Right_Strike = events{speedN}.RStrike;

        % Plot gait events
        figure(1); hold on;
        subplot(1,2,1)
        hold on;
        title([subjectName,'Left, speed ',speeds{speedN}]);
        plot(events{speedN}.signalPosLeft);
        s1 = scatter(Left_Strike, events{speedN}.signalPosLeft(Left_Strike),'d','m');
        s2 = scatter(Left_Off, events{speedN}.signalPosLeft(Left_Off),'^','m');
        legend([s1 s2], {'FootStrike','FootOff'});

        figure(2);
        subplot(1,2,1)
        hold on;
        title([subjectName,'Right, speed ',speeds{speedN}]);
        plot(events{speedN}.signalPosRight);
        s1 = scatter(Right_Strike, events{speedN}.signalPosRight(Right_Strike),'d','m');
        s2 = scatter(Right_Off, events{speedN}.signalPosRight(Right_Off),'^','m');
        legend([s1 s2], {'FootStrike','FootOff'});
        

        % Correct event(s) through a small app
        text = 'When you are done with the app, press enter to continue';
        input(text);

        % check the new events. Replot only the new elements
        figure(1); hold on;
        subplot(1,2,2);
        title('Corrected events');
        hold on;
        plot(events{speedN}.signalPosLeft);
        scatter(Left_Strike, events{speedN}.signalPosLeft(Left_Strike),'d','b');
        scatter(Left_Off, events{speedN}.signalPosLeft(Left_Off),'^','b');

        figure(2); hold on;
        subplot(1,2,2);
        title('Corrected events');
        hold on;
        plot(events{speedN}.signalPosRight);
        scatter(Right_Strike, events{speedN}.signalPosRight(Right_Strike),'d','b');
        scatter(Right_Off, events{speedN}.signalPosRight(Right_Off),'^','b');

        % Ask if a new correction is required
        prompt = 'Do you need more corrections? y/n [n]: ';
        str = input(prompt,'s');
        if isempty(str)
            str = 'n';
        end
        if strcmp(str,'y')
            text = 'Correct again and press enter when you are done';
            input(text);
        else
            % do nothing
        end

        % Save the new events vectors
        for i_evt = 1:length(Left_Strike)
            btkAppendEvent(c3d, 'Foot_Strike', (Left_Strike(i_evt) + ff) / freq, 'Left', char(subjectName), 'positionKinematic', 1);
        end
        for i_evt = 1:length(Left_Off)
            btkAppendEvent(c3d, 'Foot_Off', (Left_Off(i_evt) + ff) / freq, 'Left', char(subjectName), 'positionKinematic', 2);
        end

        for i_evt = 1:length(Right_Strike)
            btkAppendEvent(c3d, 'Foot_Strike', (Right_Strike(i_evt) + ff) / freq, 'Right', char(subjectName), 'positionKinematic', 1);
        end
        for i_evt = 1:length(Right_Off)
            btkAppendEvent(c3d, 'Foot_Off', (Right_Off(i_evt) + ff) / freq, 'Right', char(subjectName), 'positionKinematic', 2);
        end

        text = 'Press enter to close the figures and continue to the next speed';
        input(text);
        close all
        % save new eventsin c3d
        btkWriteAcquisition(c3d, newfilename);
        % in matlab file
        events{speedN}.LOff = Left_Off;
        events{speedN}.LStrike= Left_Strike;
        events{speedN}.ROff = Right_Off;
        events{speedN}.RStrike = Right_Strike;
        % Close btk handle
        btkCloseAcquisition(c3d);
        save(newfileMatlab,'events','-append');
    end
    save(newfileMatlab,'events','-append');
    app.delete;
end