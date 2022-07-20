%--------------------------------
% Compute mean for all healthy participants (take individual time-normalized values from all
% healthy participants and compute the global mean)
% Possibility to draw some plots to visualize the parameters
% For continuous parameters (joint angles trajectories)
% Mathilde 20.07.2022
% ------------------------------

clear all
close all
clc

type = 'healthy'; % no patient used in this script
speeds = {'0.11'; '0.19';'0.28';'0.36';'0.42';'0.53';'0.61';'0.69';'0.78';'0.86';'0.94';'1.03';'1.11'};
day = {'200113';'200302';'200306';'200306';'200306';'200306';'200617'};

saveFile = 0;
plotFigures =1;
varToPlot = 'RHip';

for speedN = 1:size(speeds,1)
    %% Gather all healthy values in one matrix for the different parameters
    
    for subject = [1 3:7]
        if subject < 10
            subjectN = ['REF0', num2str(subject)];
        else
            subjectN = ['REF', num2str(subject)];
        end
        folder = ['D:\StimuLOOP\DataGait\NM_Reference\ReferenceData\Data\',subjectN,'\3_C3D_Files\'];
        fileMatlab = [folder,'MatlabData\',day{subject},'_',subjectN,'_parameters'];
        load(fileMatlab,'kinNorm');
        if subject == 1
            dataAll{speedN} = kinNorm{speedN};
        else
            fields = fieldnames(kinNorm{speedN});
            for i = 1:size(fields,1)
                dataAll{speedN}.(fields{i}) = cat(2,dataAll{speedN}.(fields{i}),kinNorm{speedN}.(fields{i}));
            end
        end
    end

    %% Compute mean, std and coeff of variation
    meanAllTraj{speedN,1} = structfun(@(x)(mean(x,2,'omitnan')),dataAll{speedN},'UniformOutput',false);
    stdAllTraj{speedN,1} = structfun(@(x)(std(x,0,2,'omitnan')),dataAll{speedN},'UniformOutput',false);
    
    %% Plots
    if plotFigures 
        colorm = [[255 0 0]; [240 195 0]; [128 255 0]; [0 153 0];[0 255 255];[127 0 255];...
        [0 0 0]; [102 51 0];[255 128 0];[171 0 171];[153 0 76];[32 32 32];[102 102 0]]/255;
        figure(1); hold on;
        if strcmp(varToPlot,'RFPA') || strcmp(varToPlot,'LFPA')
            nbTrials = size(dataAll{speedN}.(varToPlot),2);
            plot([0:100],meanAllTraj{speedN}.(varToPlot),'Color',colorm(speedN,:));
            plot([0:100], meanAllTraj{speedN}.(varToPlot) + 1.96/sqrt(nbTrials)*stdAllTraj{speedN}.(varToPlot),'LineStyle',':','Color',colorm(speedN,:));
            plot([0:100], meanAllTraj{speedN}.(varToPlot) - 1.96/sqrt(nbTrials)*stdAllTraj{speedN}.(varToPlot),'LineStyle',':','Color',colorm(speedN,:));
            fillhandle = jbfill([0:100],meanAllTraj{speedN}.(varToPlot)',meanAllTraj{speedN}.(varToPlot)' + 1.96/sqrt(nbTrials)*stdAllTraj{speedN}.(varToPlot)',colorm(speedN,:),colorm(speedN,:),[],0.1);
            fillhandle = jbfill([0:100],meanAllTraj{speedN}.(varToPlot)',meanAllTraj{speedN}.(varToPlot)' - 1.96/sqrt(nbTrials)*stdAllTraj{speedN}.(varToPlot)',colorm(speedN,:),colorm(speedN,:),[],0.1);
            title(varToPlot);
        else
            ax = {'X';'Y';'Z'};
            for i = 1:3
                figure(1); hold on;
                varToPlota = [varToPlot,ax{i}];
                nbTrials = size(dataAll{speedN}.(varToPlota),2);
                subplot(2,2,i); hold on;
                plot([0:100], meanAllTraj{speedN}.(varToPlota),'Color',colorm(speedN,:));
                plot([0:100], meanAllTraj{speedN}.(varToPlota) + 1.96/sqrt(nbTrials)*stdAllTraj{speedN}.(varToPlota),'LineStyle',':','Color',colorm(speedN,:));
                plot([0:100], meanAllTraj{speedN}.(varToPlota) - 1.96/sqrt(nbTrials)*stdAllTraj{speedN}.(varToPlota),'LineStyle',':','Color',colorm(speedN,:));
                fillhandle = jbfill([0:100],meanAllTraj{speedN}.(varToPlota)',meanAllTraj{speedN}.(varToPlota)' + 1.96/sqrt(nbTrials)*stdAllTraj{speedN}.(varToPlota)',colorm(speedN,:),colorm(speedN,:),[],0.1);
                fillhandle = jbfill([0:100],meanAllTraj{speedN}.(varToPlota)',meanAllTraj{speedN}.(varToPlota)' - 1.96/sqrt(nbTrials)*stdAllTraj{speedN}.(varToPlota)',colorm(speedN,:),colorm(speedN,:),[],0.1);
                title(varToPlota);
                
                figure(100); hold on;
                subplot(2,2,i); hold on;
                plot([0:100], meanAllTraj{speedN}.(varToPlota),'Color',colorm(speedN,:));
                plot([0:100], meanAllTraj{speedN}.(varToPlota) + stdAllTraj{speedN}.(varToPlota),'LineStyle',':','Color',colorm(speedN,:));
                plot([0:100], meanAllTraj{speedN}.(varToPlota) -stdAllTraj{speedN}.(varToPlota),'LineStyle',':','Color',colorm(speedN,:));
                fillhandle = jbfill([0:100],meanAllTraj{speedN}.(varToPlota)',meanAllTraj{speedN}.(varToPlota)' + stdAllTraj{speedN}.(varToPlota)',colorm(speedN,:),colorm(speedN,:),[],0.1);
                fillhandle = jbfill([0:100],meanAllTraj{speedN}.(varToPlota)',meanAllTraj{speedN}.(varToPlota)' - stdAllTraj{speedN}.(varToPlota)',colorm(speedN,:),colorm(speedN,:),[],0.1);
                title(varToPlota);
            end
        end     
    end
end

if saveFile
    save('D:\StimuLOOP\DataGait\NM_Reference\statHealthy.mat','meanAllTraj','stdAllTraj','dataAll','-append');
end