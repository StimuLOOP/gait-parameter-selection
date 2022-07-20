%--------------------------------
% Compute mean for healthy participants
% 1. either a global mean taken over all participants (take individual step-by-step values from all
% healthy participants and compute mean from them)
% 2. or mean(meanSubject, with meanSubject being the mean over gait cycles
% for each participant
% Possibility to draw some barplots to visualize the metrics
% For discrete parameters
% Mathilde 18.07.2022
% ------------------------------

clear all
close all
clc
addpath('btk');
addpath('ChrisFunctions');

type = 'healthy'; % no patient used in this script
day = {'200113';'200302';'200306';'200306';'200306';'200306';'200617'};
speeds = {'0.11'; '0.19';'0.28';'0.36';'0.42';'0.53';'0.61';'0.69';'0.78';'0.86';'0.94';'1.03';'1.11'};
fileMatlab = 'D:\StimuLOOP\DataGait\NM_Reference\statHealthy.mat';
load(fileMatlab);

saveFile = 1;
plotFigures =0;
varToPlot = 'RStanceTime';

for speedN = 1:size(speeds,1)
   
    % Gather healthy values
    % Start with REF01
    subject = 1;
    healthy = paramAll{speedN,subject};
    folder = ['D:\StimuLOOP\DataGait\NM_Reference\ReferenceData\Data\REF01\3_C3D_Files\'];
    fileMatlab = [folder,'MatlabData\',day{subject},'_REF01_parameters'];
    load(fileMatlab);
    storeMean{1} = meanSubject;
    meanAll = meanSubject{speedN};

    for subject = [3:7]
        if speedN == 1
            if subject < 10
                subjectN = ['REF0', num2str(subject)];
            else
                subjectN = ['REF', num2str(subject)];
            end
            folder = ['D:\StimuLOOP\DataGait\NM_Reference\ReferenceData\Data\',subjectN,'\3_C3D_Files\'];
            fileMatlab = [folder,'MatlabData\',day{subject},'_',subjectN,'_parameters'];
            load(fileMatlab);
            storeMean{subject} = meanSubject;
        end
        fields = fieldnames(paramAll{speedN,subject});
        for i = 1:size(fields,1)
            % For global mean
            healthy.(fields{i}) = cat(1,healthy.(fields{i}),paramAll{speedN,subject}.(fields{i}));

            % for mean of mean
            meanAll.(fields{i}) = cat(1,meanAll.(fields{i}),storeMean{subject}{speedN}.(fields{i}));
        end
    end

    %% Global
    clear meanAllHealthy stdAllHealthy varAllHealthy
    % Compute mean, std and coeff of variation
    meanGlobalHealthy{speedN,1} = structfun(@mean,healthy,'UniformOutput',false);
    stdGlobalHealthy{speedN,1} = structfun(@(x)(std(x,0,1)),healthy,'UniformOutput',false);
    varGlobalHealthy{speedN,1} = struct(fields{:});
    for k = 1: size(fields,1)
        varGlobalHealthy{speedN,1}.(fields{k}) = stdGlobalHealthy{speedN,1}.(fields{k})./meanGlobalHealthy{speedN,1}.(fields{k});
    end

    %% Mean of participants' mean
    meanHealthy{speedN,1} = structfun(@mean,meanAll,'UniformOutput',false);
    stdHealthy{speedN,1} = structfun(@(x)(std(x,0,1)),meanAll,'UniformOutput',false);

    %% Bar plots
    if plotFigures
        names = {'Max';'Min';'ROM';'SwMax';'SwMin';'SwROM';'StMax';'StMin';'StROM'};
        
        a = [varToPlot,'Max'];
        try  meanGlobalHealthy{speedN}.(a); % see if it is a joint param (max,min,rom + sw and stance) or a spatiotemporal one (only one value)
            figure(speedN); hold on;
            for i = 1:9
                a = [varToPlot,names{i}];
                subplot(3,3,i); 
                barweb(meanGlobalHealthy{speedN}.(a)', stdGlobalHealthy{speedN}.(a)');
                axis([0.5 1.5 -inf inf]);
                title(a);
                xticklabels('');
            end
        catch
            figure(1); hold on;
            subplot(2,7,speedN);
            barweb(meanGlobalHealthy{speedN}.(varToPlot)', stdGlobalHealthy{speedN}.(varToPlot)');
            axis([0.5 1.5 -inf inf]);
                title(a);
                xticklabels('');
        end
    end
end

if saveFile
    save(fileMatlab,'meanGlobalHealthy','stdGlobalHealthy','varGlobalHealthy','meanHealthy','stdHealthy','-append');
end