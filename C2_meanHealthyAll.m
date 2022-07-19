%--------------------------------
% Compute mean for all healthy participants (take individual step-by-step values from all
% healthy participants and compute stat from them)
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
    meanAllHealthy = structfun(@mean,healthy,'UniformOutput',false);
    stdAllHealthy{speedN,1} = structfun(@(x)(std(x,0,1)),healthy,'UniformOutput',false);
    varAllHealthy{speedN,1} = struct(fields{:});
    for i = 1: size(fields,1)
        varAllHealthy{speedN,1}.(fields{k}) = stdAllHealthy{speedN,1}.(fields{k})./meanAllHealthy{speedN,1}.(fields{k});
    end

end