% -------------------
% Correct individual gait events
% Mathilde, 04.07.2022
% ------------------

clear all
close all
clc
addpath('btk');
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

subject = 8;
if strcmp(type,'healthy')
    if subject < 10
        subjectN = ['REF0', num2str(subject)];
    else
        subjectN = ['REF', num2str(subject)];
    end
    folder = ['D:\StimuLOOP\DataGait\NM_Reference\ReferenceData\Data\',subjectN,'\3_C3D_Files\'];

    % open Matlab files with gait events manually corrected
    fileMatlabAll = [folder,'MatlabData\',day{subject},'_',subjectN,'_parameters.mat'];
    load(fileMatlabAll);
elseif strcmp(type,'patient')
    if subject == 4
        subjectN = [day{subject},'_S0',num2str(subject),'_T2'];
    else
        subjectN = [day{subject},'_S0',num2str(subject)];
    end
    folder = ['D:\StimuLOOP\DataGait\NM_GaitSegmentation\',subjectN,'\04_Visual3D\'];
    % open Matlab files
    fileMatlabAll = [folder,'MatlabData\',subjectN,'_parameters.mat'];
    load(fileMatlabAll);
end

if strcmp(type,'healthy')
    sizeSpeed = size(speeds,1);
elseif strcmp(type,'patient')
    sizeSpeed = size(speeds,2);
end


speedN = 3;
app = correctEventApp;
% Load c3d file with corrected gait events
if strcmp(type,'healthy')
    filename = [folder,day{subject},'_',subjectN,'_TM_',speeds{speedN},'_NH_NS',suffixe{subject},'_EventsAdded'];
elseif strcmp(type,'patient')
    filename = [folder,subjectN,'_TM_',speeds{subject,speedN},assist{subject,speedN},suffixe,'_EventsAdded'];
end

c3d = btkReadAcquisition([filename,'.c3d']);
evts = btkGetEvents(c3d);
freq = btkGetPointFrequency(c3d);
ff = btkGetFirstFrame(c3d);
p = btkGetMetaData(c3d, 'SUBJECTS', 'NAMES');
subjectName = p.info.values;

Left_Off = events{speedN}.LOff;
Left_Strike = events{speedN}.LStrike;
Right_Off = events{speedN}.ROff;
Right_Strike = events{speedN}.RStrike;

% Plot gait events
figure(1); hold on;
subplot(1,2,1)
hold on;
title([subjectName,'Left, speed ',speeds{speedN}]);
plot(rawEvents{speedN}.signalPosLeft);
s1 = scatter(Left_Strike, rawEvents{speedN}.signalPosLeft(Left_Strike),'d','m');
s2 = scatter(Left_Off, rawEvents{speedN}.signalPosLeft(Left_Off),'^','m');
legend([s1 s2], {'FootStrike','FootOff'});

figure(2);
subplot(1,2,1)
hold on;
title([subjectName,'Right, speed ',speeds{speedN}]);
plot(rawEvents{speedN}.signalPosRight);
s1 = scatter(Right_Strike, rawEvents{speedN}.signalPosRight(Right_Strike),'d','m');
s2 = scatter(Right_Off, rawEvents{speedN}.signalPosRight(Right_Off),'^','m');
legend([s1 s2], {'FootStrike','FootOff'});

% Correct event(s) through a small app
text = 'When you are done with the app, press enter to continue';
input(text);

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

close all
% save new events in c3d
btkWriteAcquisition(c3d, [filename,'.c3d']);
% in matlab file
events{speedN}.LOff = Left_Off;
events{speedN}.LStrike= Left_Strike;
events{speedN}.ROff = Right_Off;
events{speedN}.RStrike = Right_Strike;

% Close btk handle
btkCloseAcquisition(c3d);

save(fileMatlabAll,'events','-append');
