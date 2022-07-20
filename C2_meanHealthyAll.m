%--------------------------------
% Compute mean for all healthy participants (take individual step-by-step values from all
% healthy participants and compute stat from them)
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
speeds = {'0.11'; '0.19';'0.28';'0.36';'0.42';'0.53';'0.61';'0.69';'0.78';'0.86';'0.94';'1.03';'1.11'};
fileMatlab = 'D:\StimuLOOP\DataGait\NM_Reference\statHealthy.mat';
load(fileMatlab);

saveFile = 0;
plotFigures =1;
varToPlot = 'RStanceTime';

for speedN = 1:size(speeds,1)
    %% Gather all healthy values in one matrix for the different parameters
    healthy = paramAll{speedN,1};
    for subject = [3:7]
        fields = fieldnames(paramAll{speedN,subject});
        for i = 1:size(fields,1)
        healthy.(fields{i}) = cat(1,healthy.(fields{i}),paramAll{speedN,subject}.(fields{i}));
        end
    end

    %% Compute mean, std and coeff of variation
    meanAllHealthy{speedN,1} = structfun(@mean,healthy,'UniformOutput',false);
    stdAllHealthy{speedN,1} = structfun(@(x)(std(x,0,1)),healthy,'UniformOutput',false);
    varAllHealthy{speedN,1} = struct(fields{:});
    for k = 1: size(fields,1)
        varAllHealthy{speedN,1}.(fields{k}) = stdAllHealthy{speedN,1}.(fields{k})./meanAllHealthy{speedN,1}.(fields{k});
    end

    %% Bar plots
    if plotFigures
        names = {'Max';'Min';'ROM';'SwMax';'SwMin';'SwROM';'StMax';'StMin';'StROM'};
        
        a = [varToPlot,'Max'];
        try  meanAllHealthy{speedN}.(a); % see if it is a joint param (max,min,rom + sw and stance) or a spatiotemporal one (only one value)
            figure(speedN); hold on;
            for i = 1:9
                a = [varToPlot,names{i}];
                subplot(3,3,i); 
                barweb(meanAllHealthy{speedN}.(a)', stdAllHealthy{speedN}.(a)');
                axis([0.5 1.5 -inf inf]);
                title(a);
                xticklabels('');
            end
        catch
            figure(1); hold on;
            subplot(2,7,speedN);
            barweb(meanAllHealthy{speedN}.(varToPlot)', stdAllHealthy{speedN}.(varToPlot)');
            axis([0.5 1.5 -inf inf]);
                title(a);
                xticklabels('');
        end
    end
end

if saveFile
    save(fileMatlab,'meanAllHealthy','stdAllHealthy','varAllHealthy','-append');
end