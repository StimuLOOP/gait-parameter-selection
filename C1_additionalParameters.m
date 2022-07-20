%--------------------------------
% After visual3D analysis pipeline, 
% 1. compute additional gait parameters
% 2. segment with gait cycle
% 3. compute statistics
% Mathilde 10.05.2022
% ------------------------------

clear all
close all
clc
addpath('btk');
addpath('ChrisFunctions');
% Choose between healthy and patient (different file characteristics)
type = 'healthy'; 

% Define parameters/ file characteristics
if strcmp(type,'healthy')
    day = {'200113';'200302';'200306';'200306';'200306';'200306';'200617'};
    speeds = {'0.11'; '0.19';'0.28';'0.36';'0.42';'0.53';'0.61';'0.69';'0.78';'0.86';'0.94';'1.03';'1.11'};
    suffixe = {'_MoC'; '_MoC';'';'';'_gapfilled';'';''};
elseif strcmp(type,'patient')
    day = {'190913','','190924','191003','191004','191004','191022','191023'};
    speeds = {'FST0.15' 'PGV0.1' ''; '' '' '';'FST0.7' 'PGV0.4' 'SLW0.4'; 'FST0.15' 'PGV0.05' '';'FST0.75' 'PGV0.55' 'SLW0.3'; 'PGV0.25' 'PGV0.5' 'PGV0.8'; 'FST0.35' 'PGV0.20' 'SLW0.07crop'; 'FST0.3' 'PGV0.17' 'SLW0.08'};
    assist ={'_NH_NS' '_NRH_NS1' '';'' '' ''; '_BH_NS' '_BH_NS' '_BH_NS'; '_BH_NS' '_BH_NS' '_BH_NS'; '_BH_NS' '_BH_NS' '_BH_NS'; '_BH_NS' '_BH_NS' '_BH_NS'; '_BH_NS' '_BH_NS' '_BH_NS'; '_BH_NS' '_BH_NS' '_BH_NS'};
    suffixe = '_MoCgapfilled';
end

for subject = [1 3:7]%[1 3:8]
    clearvars -except type day speeds suffixe subject RKinTable RSTTable LKinTable LSTTable globalSTTable assist paramAll
    % Load data
    if strcmp(type,'healthy')
        if subject < 10
            subjectN = ['REF0', num2str(subject)];
        else
            subjectN = ['REF', num2str(subject)];
        end
        folder = ['D:\StimuLOOP\DataGait\NM_Reference\ReferenceData\Data\',subjectN,'\3_C3D_Files\'];
        fileMatlab = [folder,'MatlabData\',day{subject},'_',subjectN,'_parameters'];
    elseif strcmp(type,'patient')
        if subject == 4
            subjectN = [day{subject},'_S0',num2str(subject),'_T2'];
        else
            subjectN = [day{subject},'_S0',num2str(subject)];
        end
        folder = ['D:\StimuLOOP\DataGait\NM_GaitSegmentation\',subjectN,'\04_Visual3D\'];
        fileMatlab = [folder,'MatlabData\',subjectN,'_parameters'];
    end
    load(fileMatlab);
    
    %% Compute additional parameters
    if strcmp(type,'healthy')
             sizeSpeed = size(speeds,1);
        elseif strcmp(type,'patient')
            sizeSpeed = size(speeds,2);
    end
    X = [1 0 0];
    Y = [0 1 0];
    Z = [0 0 1];

    for speedN = 1:sizeSpeed
        if strcmp(type,'healthy')
            filec3d = [folder, day{subject},'_',subjectN,'_TM_',speeds{speedN},'_NH_NS',suffixe{subject}];
        elseif strcmp(type,'patient')
            if isempty(speeds{subject, speedN})
                break
            end
            filec3d = [folder,subjectN,'_TM_',speeds{subject,speedN},assist{subject,speedN},suffixe];
        end
        c3d = btkReadAcquisition([filec3d,'.c3d']);
        markers = btkGetMarkers(c3d);
        freq = btkGetPointFrequency(c3d);

        % Rename markers if names not conform to the ones used here
        try
            lpsi = markers.LPSI;
        catch
            markers = renameMarkers(markers);
        end

         % correct scale for S01
        if (strcmp(subject,'190913_S01') && strcmp(speed,'PGV0.1'))
            mult = @(x)(x*1000);
            markers = structfun(mult,markers,'UniformOutput',false);
        end
        % Rename markers if names not conform to the ones used here
        try
            lpsi = markers.LPSI;
        catch
            markers = renameMarkers(markers);
        end
        
        nbSamples = size(markers.LASI,1);

        % Trunk motions (not available for stroke patients)
        if strcmp(type,'healthy')
            % define trunk frame
            Ztrunk = (markers.CLAV - markers.TV10);
            Ztrunk = Ztrunk./vecnorm(Ztrunk,2,2);
            shoulVec = markers.RSHO - markers.LSHO;
            shoulVec = shoulVec./vecnorm(shoulVec,2,2);
            Xtrunk = cross(shoulVec,Ztrunk);
            Ytrunk = cross(Ztrunk, Xtrunk);

            % Compute rotation matrix
            for i=1:nbSamples
                R(1,1) = dot(Xtrunk(i,:),X); R(1,2) = dot(Xtrunk(i,:),Y);  R(1,3) = dot(Xtrunk(i,:),Z);
                R(2,1) = dot(Ytrunk(i,:),X); R(2,2) = dot(Ytrunk(i,:),Y);  R(2,3) = dot(Ytrunk(i,:),Z);
                R(3,1) = dot(Ztrunk(i,:),X); R(3,2) = dot(Ztrunk(i,:),Y);  R(3,3) = dot(Ztrunk(i,:),Z);

                % Compute Euler angles between the lab and the trunk frames
                eulerT{speedN,1}(i,:) = rotm2eul(R)*180/pi;
            end

            % Arm swing
            centerPoint = (markers.RASI + markers.LASI)./2;
            LArmSw{speedN} = ((markers.LRSP + markers.LUSP)./2) - centerPoint;
            RArmSw{speedN} = ((markers.RRSP + markers.RUSP)./2) - centerPoint;
        else
            eulerT{speedN,1} = [];
            LArmSw{speedN} = [];
            RArmSw{speedN} = [];
        end

        
    
        % Foot progression angles 
        Rfootvec = markers.RFM2 - markers.RFCC;
        for sampleNo = 1:size(Rfootvec,1)
            RFPA{speedN}(sampleNo,:) = atan2d(dot(cross(Rfootvec(sampleNo,:),-X),Z),dot(-X,Rfootvec(sampleNo,:)));
        end
        Lfootvec = markers.LFM2 - markers.LFCC;
        for sampleNo = 1:size(Lfootvec,1)
            LFPA{speedN}(sampleNo,:) = atan2d(dot(cross(-X,Lfootvec(sampleNo,:)),Z),dot(-X,Lfootvec(sampleNo,:)));
        end

        % Filter FPA
%         fc = 5;
%         fs = freq;
%         [b,a] = butter(2,fc/(fs/2));
%         LFPA_filter = filter(b,a,LFPA{speedN},[],1);
%         LFPA_filter(1:10,:) = LFPA{speedN}(1:10); % remove the peak at the beginning due to the filtering
% %         figure;plot(correctLFP_filter)
%         LFPA{speedN,1} = LFPA_filter;

        %% Segment kinetic and kinematic parameters with gait events + time normalization + max and mean
        for e = 1:min([size(events{speedN}.LStrike,2) size(events{speedN}.LOff,2) size(events{speedN}.RStrike,2) size(events{speedN}.ROff,2)])-1
            % Segmentation - Step level
            kinSeg{speedN,1}.RAnkle{e,1} = RAnkle{speedN}(events{speedN}.RStrike(e):events{speedN}.RStrike(e+1),:);
            kinSeg{speedN,1}.RHip{e,1} = RHip{speedN}(events{speedN}.RStrike(e):events{speedN}.RStrike(e+1),:);
            kinSeg{speedN,1}.RKnee{e,1} = RKnee{speedN}(events{speedN}.RStrike(e):events{speedN}.RStrike(e+1),:);
            kinSeg{speedN,1}.RKnee{e,1}(:,1) = - kinSeg{speedN,1}.RKnee{e,1}(:,1); % take the opposite so that knee flexion is positive

            kinSeg{speedN,1}.LAnkle{e,1} = LAnkle{speedN}(events{speedN}.LStrike(e):events{speedN}.LStrike(e+1),:);
            kinSeg{speedN,1}.LHip{e,1} = LHip{speedN}(events{speedN}.LStrike(e):events{speedN}.LStrike(e+1),:);
            kinSeg{speedN,1}.LKnee{e,1} = LKnee{speedN}(events{speedN}.LStrike(e):events{speedN}.LStrike(e+1),:);
            kinSeg{speedN,1}.LKnee{e,1}(:,1) = - kinSeg{speedN,1}.LKnee{e,1}(:,1); % take the opposite so that knee flexion is positive
            
            if strcmp(type,"healthy")
            kinSeg{speedN,1}.RArmSw{e,1} = RArmSw{speedN}(events{speedN}.RStrike(e):events{speedN}.RStrike(e+1),:);
            kinSeg{speedN,1}.LArmSw{e,1} = LArmSw{speedN}(events{speedN}.LStrike(e):events{speedN}.LStrike(e+1),:);
            end

            kinSeg{speedN,1}.RFPA{e,1} = RFPA{speedN}(events{speedN}.RStrike(e):events{speedN}.RStrike(e+1),:);
            kinSeg{speedN,1}.LFPA{e,1} = LFPA{speedN}(events{speedN}.LStrike(e):events{speedN}.LStrike(e+1),:);
            
                % Sample rate of force plates = 600Hz-> rescaling required
            ff = 600; % forceplate frequency
            fmc = 100; % motion capture frequency
            kinSeg{speedN,1}.R_FP{e,1} = FP1{speedN}(events{speedN}.RStrike(e)*ff/fmc:events{speedN}.RStrike(e+1)*ff/fmc,:);
            kinSeg{speedN,1}.L_FP{e,1} = FP2{speedN}(events{speedN}.LStrike(e)*ff/fmc:events{speedN}.LStrike(e+1)*ff/fmc,:);

            if strcmp(type,'healthy')
                kinSeg{speedN,1}.ReulerT{e,1} = eulerT{speedN}(events{speedN}.RStrike(e):events{speedN}.RStrike(e+1),:);
                kinSeg{speedN,1}.LeulerT{e,1} = eulerT{speedN}(events{speedN}.LStrike(e):events{speedN}.LStrike(e+1),:);
            end

             % Segmentation - Stance phase
             % Adapt the index of the off event, depending on which event
             % is detected first
             if events{speedN}.RStrike(e) < events{speedN}.ROff(e)
                 indexROff = e;
             else
                 indexROff = e+1;
             end
             if events{speedN}.LStrike(e) < events{speedN}.LOff(e)
                 indexLOff = e;
             else
                 indexLOff = e+1;
             end
            kinStance{speedN,1}.RAnkle{e,1} = RAnkle{speedN}(events{speedN}.RStrike(e):events{speedN}.ROff(indexROff),:);
            kinStance{speedN,1}.RHip{e,1} = RHip{speedN}(events{speedN}.RStrike(e):events{speedN}.ROff(indexROff),:);
            kinStance{speedN,1}.RKnee{e,1} = RKnee{speedN}(events{speedN}.RStrike(e):events{speedN}.ROff(indexROff),:);
            kinStance{speedN,1}.RKnee{e,1}(:,1) = - kinStance{speedN,1}.RKnee{e,1}(:,1); % take the opposite so that knee flexion is positive

            kinStance{speedN,1}.LAnkle{e,1} = LAnkle{speedN}(events{speedN}.LStrike(e):events{speedN}.LOff(indexLOff),:);
            kinStance{speedN,1}.LHip{e,1} = LHip{speedN}(events{speedN}.LStrike(e):events{speedN}.LOff(indexLOff),:);
            kinStance{speedN,1}.LKnee{e,1} = LKnee{speedN}(events{speedN}.LStrike(e):events{speedN}.LOff(indexLOff),:);
            kinStance{speedN,1}.LKnee{e,1}(:,1) = - kinStance{speedN,1}.LKnee{e,1}(:,1); % take the opposite so that knee flexion is positive
            
            if strcmp(type,"healthy")
            kinStance{speedN,1}.RArmSw{e,1} = RArmSw{speedN}(events{speedN}.RStrike(e):events{speedN}.ROff(indexROff),:);
            kinStance{speedN,1}.LArmSw{e,1} = LArmSw{speedN}(events{speedN}.LStrike(e):events{speedN}.LOff(indexLOff),:);
            end

            kinStance{speedN,1}.RFPA{e,1} = RFPA{speedN}(events{speedN}.RStrike(e):events{speedN}.ROff(indexROff),:);
            kinStance{speedN,1}.LFPA{e,1} = LFPA{speedN}(events{speedN}.LStrike(e):events{speedN}.LOff(indexLOff),:);
            
            ff = 600; % forceplate frequency
            fmc = 100; % motion capture frequency
            kinStance{speedN,1}.R_FP{e,1} = FP1{speedN}(events{speedN}.RStrike(e)*ff/fmc:events{speedN}.ROff(indexROff)*ff/fmc,:);
            kinStance{speedN,1}.L_FP{e,1} = FP2{speedN}(events{speedN}.LStrike(e)*ff/fmc:events{speedN}.LOff(indexLOff)*ff/fmc,:);

            if strcmp(type,'healthy')
                kinStance{speedN,1}.ReulerT{e,1} = eulerT{speedN}(events{speedN}.RStrike(e):events{speedN}.ROff(indexROff),:);
                kinStance{speedN,1}.LeulerT{e,1} = eulerT{speedN}(events{speedN}.LStrike(e):events{speedN}.LOff(indexLOff),:);
            end

            % Segmentation - Swing phase
             % Adapt the index of the off event, depending on which event
             % is detected first
             if events{speedN}.ROff(e) < events{speedN}.RStrike(e)
                 indexRStrike = e;
             else
                 indexRStrike = e+1;
             end
             if events{speedN}.LOff(e) < events{speedN}.LStrike(e)
                 indexLStrike = e;
             else
                 indexLStrike = e+1;
             end
            kinSwing{speedN,1}.RAnkle{e,1} = RAnkle{speedN}(events{speedN}.ROff(e):events{speedN}.RStrike(indexRStrike),:);
            kinSwing{speedN,1}.RHip{e,1} = RHip{speedN}(events{speedN}.ROff(e):events{speedN}.RStrike(indexRStrike),:);
            kinSwing{speedN,1}.RKnee{e,1} = RKnee{speedN}(events{speedN}.ROff(e):events{speedN}.RStrike(indexRStrike),:);
            kinSwing{speedN,1}.RKnee{e,1}(:,1) = - kinSwing{speedN,1}.RKnee{e,1}(:,1); % take the opposite so that knee flexion is positive

            kinSwing{speedN,1}.LAnkle{e,1} = LAnkle{speedN}(events{speedN}.LOff(e):events{speedN}.LStrike(indexLStrike),:);
            kinSwing{speedN,1}.LHip{e,1} = LHip{speedN}(events{speedN}.LOff(e):events{speedN}.LStrike(indexLStrike),:);
            kinSwing{speedN,1}.LKnee{e,1} = LKnee{speedN}(events{speedN}.LOff(e):events{speedN}.LStrike(indexLStrike),:);
            kinSwing{speedN,1}.LKnee{e,1}(:,1) = - kinSwing{speedN,1}.LKnee{e,1}(:,1); % take the opposite so that knee flexion is positive
            
            if strcmp(type,"healthy")
            kinSwing{speedN,1}.RArmSw{e,1} = RArmSw{speedN}(events{speedN}.ROff(e):events{speedN}.RStrike(indexRStrike),:);
            kinSwing{speedN,1}.LArmSw{e,1} = LArmSw{speedN}(events{speedN}.LOff(e):events{speedN}.LStrike(indexLStrike),:);
            end

            kinSwing{speedN,1}.RFPA{e,1} = RFPA{speedN}(events{speedN}.ROff(e):events{speedN}.RStrike(indexRStrike),:);
            kinSwing{speedN,1}.LFPA{e,1} = LFPA{speedN}(events{speedN}.LOff(e):events{speedN}.LStrike(indexLStrike),:);
            
            ff = 600; % forceplate frequency
            fmc = 100; % motion capture frequency
            kinSwing{speedN,1}.R_FP{e,1} = FP1{speedN}(events{speedN}.ROff(e)*ff/fmc:events{speedN}.RStrike(indexRStrike)*ff/fmc,:);
            kinSwing{speedN,1}.L_FP{e,1} = FP2{speedN}(events{speedN}.LOff(e)*ff/fmc:events{speedN}.LStrike(indexLStrike)*ff/fmc,:);

            if strcmp(type,'healthy')
                kinSwing{speedN,1}.ReulerT{e,1} = eulerT{speedN}(events{speedN}.ROff(e):events{speedN}.RStrike(indexRStrike),:);
                kinSwing{speedN,1}.LeulerT{e,1} = eulerT{speedN}(events{speedN}.LOff(e):events{speedN}.LStrike(indexLStrike),:);
            end

            % Time normalization
            % Right ankle
            sizeVar = size(kinSeg{speedN}.RAnkle{e}(:,1),1);
            kinNorm{speedN,1}.RAnkleX(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.RAnkle{e}(:,1),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.RAnkleY(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.RAnkle{e}(:,2),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.RAnkleZ(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.RAnkle{e}(:,3),linspace(1, sizeVar, 101)))';
            % Right hip
            kinNorm{speedN,1}.RHipX(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.RHip{e}(:,1),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.RHipY(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.RHip{e}(:,2),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.RHipZ(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.RHip{e}(:,3),linspace(1, sizeVar, 101)))';
            % Right knee
            kinNorm{speedN,1}.RKneeX(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.RKnee{e}(:,1),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.RKneeY(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.RKnee{e}(:,2),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.RKneeZ(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.RKnee{e}(:,3),linspace(1, sizeVar, 101)))';
            % Right arm swing
            if strcmp(type,"healthy")
            kinNorm{speedN,1}.RArmSwX(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.RArmSw{e}(:,1),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.RArmSwY(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.RArmSw{e}(:,2),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.RArmSwZ(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.RArmSw{e}(:,3),linspace(1, sizeVar, 101)))';
            end
            % Right foot progression angle
            kinNorm{speedN,1}.RFPA(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.RFPA{e}(:,1),linspace(1, sizeVar, 101)))';
            % Trunk angles on right cycles
            if strcmp(type,"healthy")
            kinNorm{speedN,1}.ReulerTX(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.ReulerT{e}(:,1),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.ReulerTY(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.ReulerT{e}(:,2),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.ReulerTZ(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.ReulerT{e}(:,3),linspace(1, sizeVar, 101)))';
            end
            % Right forceplate
            sizeForce = size(kinSeg{speedN,1}.R_FP{e}(:,1),1);
            kinNorm{speedN,1}.R_FPX(:,e) = (interp1([1:sizeForce], kinSeg{speedN,1}.R_FP{e}(:,1),linspace(1, sizeForce, 101)))';
            kinNorm{speedN,1}.R_FPY(:,e) = (interp1([1:sizeForce], kinSeg{speedN,1}.R_FP{e}(:,2),linspace(1, sizeForce, 101)))';
            kinNorm{speedN,1}.R_FPZ(:,e) = (interp1([1:sizeForce], kinSeg{speedN,1}.R_FP{e}(:,3),linspace(1, sizeForce, 101)))';
            
            % Left ankle
            sizeVar = size(kinSeg{speedN,1}.LAnkle{e}(:,1),1);
            kinNorm{speedN,1}.LAnkleX(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.LAnkle{e}(:,1),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.LAnkleY(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.LAnkle{e}(:,2),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.LAnkleZ(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.LAnkle{e}(:,3),linspace(1, sizeVar, 101)))';
            % Left hip
            kinNorm{speedN,1}.LHipX(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.LHip{e}(:,1),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.LHipY(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.LHip{e}(:,2),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.LHipZ(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.LHip{e}(:,3),linspace(1, sizeVar, 101)))';
            % Left knee
            kinNorm{speedN,1}.LKneeX(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.LKnee{e}(:,1),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.LKneeY(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.LKnee{e}(:,2),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.LKneeZ(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.LKnee{e}(:,3),linspace(1, sizeVar, 101)))';
            % Left arm swing
            if strcmp(type,"healthy")
            kinNorm{speedN,1}.LArmSwX(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.LArmSw{e}(:,1),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.LArmSwY(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.LArmSw{e}(:,2),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.LArmSwZ(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.LArmSw{e}(:,3),linspace(1, sizeVar, 101)))';
            end
            % Left foot progression angle
            kinNorm{speedN,1}.LFPA(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.LFPA{e}(:,1),linspace(1, sizeVar, 101)))';
            % Trunk angles on left cycles
            if strcmp(type,"healthy")
            kinNorm{speedN,1}.LeulerTX(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.LeulerT{e}(:,1),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.LeulerTY(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.LeulerT{e}(:,2),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.LeulerTZ(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.LeulerT{e}(:,3),linspace(1, sizeVar, 101)))';
            end
            % Left forceplate
            sizeForce = size(kinSeg{speedN,1}.L_FP{e}(:,1),1);
            kinNorm{speedN,1}.L_FPX(:,e) = (interp1([1:sizeForce], kinSeg{speedN,1}.L_FP{e}(:,1),linspace(1, sizeForce, 101)))';
            kinNorm{speedN,1}.L_FPY(:,e) = (interp1([1:sizeForce], kinSeg{speedN,1}.L_FP{e}(:,2),linspace(1, sizeForce, 101)))';
            kinNorm{speedN,1}.L_FPZ(:,e) = (interp1([1:sizeForce], kinSeg{speedN,1}.L_FP{e}(:,3),linspace(1, sizeForce, 101)))';
            
            % Compute extrema and ROM within gait cycles/swing phase/stance phase
            param{speedN,1}.RAnkleMax(e,:) = max(kinSeg{speedN,1}.RAnkle{e,1});
            param{speedN,1}.RAnkleMin(e,:) = min(kinSeg{speedN,1}.RAnkle{e,1});
            param{speedN,1}.RAnkleSwMax(e,:) = max(kinSwing{speedN,1}.RAnkle{e,1});
            param{speedN,1}.RAnkleSwMin(e,:) = min(kinSwing{speedN,1}.RAnkle{e,1});
            param{speedN,1}.RAnkleStMax(e,:) = max(kinStance{speedN,1}.RAnkle{e,1});
            param{speedN,1}.RAnkleStMin(e,:) = min(kinStance{speedN,1}.RAnkle{e,1});
            param{speedN,1}.RAnkleROM(e,:) = param{speedN,1}.RAnkleMax(e,:) - param{speedN,1}.RAnkleMin(e,:);
            param{speedN,1}.RAnkleSwROM(e,:) = param{speedN,1}.RAnkleSwMax(e,:) - param{speedN,1}.RAnkleSwMin(e,:);
            param{speedN,1}.RAnkleStROM(e,:) = param{speedN,1}.RAnkleStMax(e,:) - param{speedN,1}.RAnkleStMin(e,:);

            param{speedN,1}.RHipMax(e,:) = max(kinSeg{speedN,1}.RHip{e,1});
            param{speedN,1}.RHipMin(e,:) = min(kinSeg{speedN,1}.RHip{e,1});
            param{speedN,1}.RHipSwMax(e,:) = max(kinSwing{speedN,1}.RHip{e,1});
            param{speedN,1}.RHipSwMin(e,:) = min(kinSwing{speedN,1}.RHip{e,1});
            param{speedN,1}.RHipStMax(e,:) = max(kinStance{speedN,1}.RHip{e,1});
            param{speedN,1}.RHipStMin(e,:) = min(kinStance{speedN,1}.RHip{e,1});
            param{speedN,1}.RHipROM(e,:) = param{speedN,1}.RHipMax(e,:) - param{speedN,1}.RHipMin(e,:);
            param{speedN,1}.RHipSwROM(e,:) = param{speedN,1}.RHipSwMax(e,:) - param{speedN,1}.RHipSwMin(e,:);
            param{speedN,1}.RHipStROM(e,:) = param{speedN,1}.RHipStMax(e,:) - param{speedN,1}.RHipStMin(e,:);

            param{speedN,1}.RKneeMax(e,:) = max(kinSeg{speedN,1}.RKnee{e,1});
            param{speedN,1}.RKneeMin(e,:) = min(kinSeg{speedN,1}.RKnee{e,1});
            param{speedN,1}.RKneeSwMax(e,:) = max(kinSwing{speedN,1}.RKnee{e,1});
            param{speedN,1}.RKneeSwMin(e,:) = min(kinSwing{speedN,1}.RKnee{e,1});
            param{speedN,1}.RKneeStMax(e,:) = max(kinStance{speedN,1}.RKnee{e,1});
            param{speedN,1}.RKneeStMin(e,:) = min(kinStance{speedN,1}.RKnee{e,1});
            param{speedN,1}.RKneeROM(e,:) = param{speedN,1}.RKneeMax(e,:) - param{speedN,1}.RKneeMin(e,:);
            param{speedN,1}.RKneeSwROM(e,:) = param{speedN,1}.RKneeSwMax(e,:) - param{speedN,1}.RKneeSwMin(e,:);
            param{speedN,1}.RKneeStROM(e,:) = param{speedN,1}.RKneeStMax(e,:) - param{speedN,1}.RKneeStMin(e,:);

            if strcmp(type,"healthy")
            param{speedN,1}.RArmSwMax(e,:) = max(kinSeg{speedN,1}.RArmSw{e,1});
            param{speedN,1}.RArmSwMin(e,:) = min(kinSeg{speedN,1}.RArmSw{e,1});
            param{speedN,1}.RArmSwSwMax(e,:) = max(kinSwing{speedN,1}.RArmSw{e,1});
            param{speedN,1}.RArmSwSwMin(e,:) = min(kinSwing{speedN,1}.RArmSw{e,1});
            param{speedN,1}.RArmSwStMax(e,:) = max(kinStance{speedN,1}.RArmSw{e,1});
            param{speedN,1}.RArmSwStMin(e,:) = min(kinStance{speedN,1}.RArmSw{e,1});
            param{speedN,1}.RArmSwROM(e,:) = param{speedN,1}.RArmSwMax(e,:) - param{speedN,1}.RArmSwMin(e,:);
            param{speedN,1}.RArmSwSwROM(e,:) = param{speedN,1}.RArmSwSwMax(e,:) - param{speedN,1}.RArmSwSwMin(e,:);
            param{speedN,1}.RArmSwStROM(e,:) = param{speedN,1}.RArmSwStMax(e,:) - param{speedN,1}.RArmSwStMin(e,:);
            
            param{speedN,1}.ReulerTMax(e,:) = max(kinSeg{speedN,1}.ReulerT{e,1});
            param{speedN,1}.ReulerTMin(e,:) = min(kinSeg{speedN,1}.ReulerT{e,1});
            param{speedN,1}.ReulerTSwMax(e,:) = max(kinSwing{speedN,1}.ReulerT{e,1});
            param{speedN,1}.ReulerTSwMin(e,:) = min(kinSwing{speedN,1}.ReulerT{e,1});
            param{speedN,1}.ReulerTStMax(e,:) = max(kinStance{speedN,1}.ReulerT{e,1});
            param{speedN,1}.ReulerTStMin(e,:) = min(kinStance{speedN,1}.ReulerT{e,1});
            param{speedN,1}.ReulerTROM(e,:) = param{speedN,1}.ReulerTMax(e,:) - param{speedN,1}.ReulerTMin(e,:);
            param{speedN,1}.ReulerTSwROM(e,:) = param{speedN,1}.ReulerTSwMax(e,:) - param{speedN,1}.ReulerTSwMin(e,:);
            param{speedN,1}.ReulerTStROM(e,:) = param{speedN,1}.ReulerTStMax(e,:) - param{speedN,1}.ReulerTStMin(e,:);
            end

            param{speedN,1}.R_FPMax(e,:) = max(kinSeg{speedN,1}.R_FP{e,1});
            param{speedN,1}.R_FPMin(e,:) = min(kinSeg{speedN,1}.R_FP{e,1});
            param{speedN,1}.R_FPSwMax(e,:) = max(kinSwing{speedN,1}.R_FP{e,1});
            param{speedN,1}.R_FPSwMin(e,:) = min(kinSwing{speedN,1}.R_FP{e,1});
            param{speedN,1}.R_FPStMax(e,:) = max(kinStance{speedN,1}.R_FP{e,1});
            param{speedN,1}.R_FPStMin(e,:) = min(kinStance{speedN,1}.R_FP{e,1});
            param{speedN,1}.R_FPROM(e,:) = param{speedN,1}.R_FPMax(e,:) - param{speedN,1}.R_FPMin(e,:);
            param{speedN,1}.R_FPSwROM(e,:) = param{speedN,1}.R_FPSwMax(e,:) - param{speedN,1}.R_FPSwMin(e,:);
            param{speedN,1}.R_FPStROM(e,:) = param{speedN,1}.R_FPStMax(e,:) - param{speedN,1}.R_FPStMin(e,:);

            param{speedN,1}.RFPAMax(e,:) = max(kinSeg{speedN,1}.RFPA{e,1});
            param{speedN,1}.RFPAMin(e,:) = min(kinSeg{speedN,1}.RFPA{e,1});
            param{speedN,1}.RFPASwMax(e,:) = max(kinSwing{speedN,1}.RFPA{e,1});
            param{speedN,1}.RFPASwMin(e,:) = min(kinSwing{speedN,1}.RFPA{e,1});
            param{speedN,1}.RFPAStMax(e,:) = max(kinStance{speedN,1}.RFPA{e,1});
            param{speedN,1}.RFPAStMin(e,:) = min(kinStance{speedN,1}.RFPA{e,1});
            param{speedN,1}.RFPAROM(e,:) = param{speedN,1}.RFPAMax(e,:) - param{speedN,1}.RFPAMin(e,:);
            param{speedN,1}.RFPASwROM(e,:) = param{speedN,1}.RFPASwMax(e,:) - param{speedN,1}.RFPASwMin(e,:);
            param{speedN,1}.RFPAStROM(e,:) = param{speedN,1}.RFPAStMax(e,:) - param{speedN,1}.RFPAStMin(e,:);

            param{speedN,1}.LAnkleMax(e,:) = max(kinSeg{speedN,1}.LAnkle{e,1});
            param{speedN,1}.LAnkleMin(e,:) = min(kinSeg{speedN,1}.LAnkle{e,1});
            param{speedN,1}.LAnkleSwMax(e,:) = max(kinSwing{speedN,1}.LAnkle{e,1});
            param{speedN,1}.LAnkleSwMin(e,:) = min(kinSwing{speedN,1}.LAnkle{e,1});
            param{speedN,1}.LAnkleStMax(e,:) = max(kinStance{speedN,1}.LAnkle{e,1});
            param{speedN,1}.LAnkleStMin(e,:) = min(kinStance{speedN,1}.LAnkle{e,1});
            param{speedN,1}.LAnkleROM(e,:) = param{speedN,1}.LAnkleMax(e,:) - param{speedN,1}.LAnkleMin(e,:);
            param{speedN,1}.LAnkleSwROM(e,:) = param{speedN,1}.LAnkleSwMax(e,:) - param{speedN,1}.LAnkleSwMin(e,:);
            param{speedN,1}.LAnkleStROM(e,:) = param{speedN,1}.LAnkleStMax(e,:) - param{speedN,1}.LAnkleStMin(e,:);

            param{speedN,1}.LHipMax(e,:) = max(kinSeg{speedN,1}.LHip{e,1});
            param{speedN,1}.LHipMin(e,:) = min(kinSeg{speedN,1}.LHip{e,1});
             param{speedN,1}.LHipSwMax(e,:) = max(kinSwing{speedN,1}.LHip{e,1});
            param{speedN,1}.LHipSwMin(e,:) = min(kinSwing{speedN,1}.LHip{e,1});
             param{speedN,1}.LHipStMax(e,:) = max(kinStance{speedN,1}.LHip{e,1});
            param{speedN,1}.LHipStMin(e,:) = min(kinStance{speedN,1}.LHip{e,1});
            param{speedN,1}.LHipROM(e,:) = param{speedN,1}.LHipMax(e,:) - param{speedN,1}.LHipMin(e,:);
            param{speedN,1}.LHipSwROM(e,:) = param{speedN,1}.LHipSwMax(e,:) - param{speedN,1}.LHipSwMin(e,:);
            param{speedN,1}.LHipStROM(e,:) = param{speedN,1}.LHipStMax(e,:) - param{speedN,1}.LHipStMin(e,:);

            param{speedN,1}.LKneeMax(e,:) = max(kinSeg{speedN,1}.LKnee{e,1});
            param{speedN,1}.LKneeMin(e,:) = min(kinSeg{speedN,1}.LKnee{e,1});
            param{speedN,1}.LKneeSwMax(e,:) = max(kinSwing{speedN,1}.LKnee{e,1});
            param{speedN,1}.LKneeSwMin(e,:) = min(kinSwing{speedN,1}.LKnee{e,1});
            param{speedN,1}.LKneeStMax(e,:) = max(kinStance{speedN,1}.LKnee{e,1});
            param{speedN,1}.LKneeStMin(e,:) = min(kinStance{speedN,1}.LKnee{e,1});
            param{speedN,1}.LKneeROM(e,:) = param{speedN,1}.LKneeMax(e,:) - param{speedN,1}.LKneeMin(e,:);
            param{speedN,1}.LKneeSwROM(e,:) = param{speedN,1}.LKneeSwMax(e,:) - param{speedN,1}.LKneeSwMin(e,:);
            param{speedN,1}.LKneeStROM(e,:) = param{speedN,1}.LKneeStMax(e,:) - param{speedN,1}.LKneeStMin(e,:);
            
            if strcmp(type,"healthy")
            param{speedN,1}.LArmSwMax(e,:) = max(kinSeg{speedN,1}.LArmSw{e,1});
            param{speedN,1}.LArmSwMin(e,:) = min(kinSeg{speedN,1}.LArmSw{e,1});
            param{speedN,1}.LArmSwSwMax(e,:) = max(kinSwing{speedN,1}.LArmSw{e,1});
            param{speedN,1}.LArmSwSwMin(e,:) = min(kinSwing{speedN,1}.LArmSw{e,1});
            param{speedN,1}.LArmSwStMax(e,:) = max(kinStance{speedN,1}.LArmSw{e,1});
            param{speedN,1}.LArmSwStMin(e,:) = min(kinStance{speedN,1}.LArmSw{e,1});
            param{speedN,1}.LArmSwROM(e,:) = param{speedN,1}.LArmSwMax(e,:) - param{speedN,1}.LArmSwMin(e,:);
            param{speedN,1}.LArmSwSwROM(e,:) = param{speedN,1}.LArmSwSwMax(e,:) - param{speedN,1}.LArmSwSwMin(e,:);
            param{speedN,1}.LArmSwStROM(e,:) = param{speedN,1}.LArmSwStMax(e,:) - param{speedN,1}.LArmSwStMin(e,:);

            param{speedN,1}.LeulerTMax(e,:) = max(kinSeg{speedN,1}.LeulerT{e,1});
            param{speedN,1}.LeulerTMin(e,:) = min(kinSeg{speedN,1}.LeulerT{e,1});
            param{speedN,1}.LeulerTSwMax(e,:) = max(kinSwing{speedN,1}.LeulerT{e,1});
            param{speedN,1}.LeulerTSwMin(e,:) = min(kinSwing{speedN,1}.LeulerT{e,1});
            param{speedN,1}.LeulerTStMax(e,:) = max(kinStance{speedN,1}.LeulerT{e,1});
            param{speedN,1}.LeulerTStMin(e,:) = min(kinStance{speedN,1}.LeulerT{e,1});
            param{speedN,1}.LeulerTROM(e,:) = param{speedN,1}.LeulerTMax(e,:) - param{speedN,1}.LeulerTMin(e,:);
            param{speedN,1}.LeulerTSwROM(e,:) = param{speedN,1}.LeulerTSwMax(e,:) - param{speedN,1}.LeulerTSwMin(e,:);
            param{speedN,1}.LeulerTStROM(e,:) = param{speedN,1}.LeulerTStMax(e,:) - param{speedN,1}.LeulerTStMin(e,:);
            end

            param{speedN,1}.L_FPMax(e,:) = max(kinSeg{speedN,1}.L_FP{e,1});
            param{speedN,1}.L_FPMin(e,:) = min(kinSeg{speedN,1}.L_FP{e,1});
            param{speedN,1}.L_FPSwMax(e,:) = max(kinSwing{speedN,1}.L_FP{e,1});
            param{speedN,1}.L_FPSwMin(e,:) = min(kinSwing{speedN,1}.L_FP{e,1});
            param{speedN,1}.L_FPStMax(e,:) = max(kinStance{speedN,1}.L_FP{e,1});
            param{speedN,1}.L_FPStMin(e,:) = min(kinStance{speedN,1}.L_FP{e,1});
            param{speedN,1}.L_FPROM(e,:) = param{speedN,1}.L_FPMax(e,:) - param{speedN,1}.L_FPMin(e,:);
            param{speedN,1}.L_FPSwROM(e,:) = param{speedN,1}.L_FPSwMax(e,:) - param{speedN,1}.L_FPSwMin(e,:);
            param{speedN,1}.L_FPStROM(e,:) = param{speedN,1}.L_FPStMax(e,:) - param{speedN,1}.L_FPStMin(e,:);

            param{speedN,1}.LFPAMax(e,:) = max(kinSeg{speedN,1}.LFPA{e,1});
            param{speedN,1}.LFPAMin(e,:) = min(kinSeg{speedN,1}.LFPA{e,1});
            param{speedN,1}.LFPASwMax(e,:) = max(kinSwing{speedN,1}.LFPA{e,1});
            param{speedN,1}.LFPASwMin(e,:) = min(kinSwing{speedN,1}.LFPA{e,1});
            param{speedN,1}.LFPAStMax(e,:) = max(kinStance{speedN,1}.LFPA{e,1});
            param{speedN,1}.LFPAStMin(e,:) = min(kinStance{speedN,1}.LFPA{e,1});
            param{speedN,1}.LFPAROM(e,:) = param{speedN,1}.LFPAMax(e,:) - param{speedN,1}.LFPAMin(e,:);
            param{speedN,1}.LFPASwROM(e,:) = param{speedN,1}.LFPASwMax(e,:) - param{speedN,1}.LFPASwMin(e,:);
            param{speedN,1}.LFPAStROM(e,:) = param{speedN,1}.LFPAStMax(e,:) - param{speedN,1}.LFPAStMin(e,:);
        
        end

        % Add spatiotemporal parameters in the param structure
        param{speedN,1}.RStanceTime = RStanceTime{speedN};
        param{speedN,1}.RStepLength = RStepLength{speedN};
        param{speedN,1}.RStepTime = RStepTime{speedN};
        param{speedN,1}.RStrideLength = RStrideLength{speedN};
        param{speedN,1}.RSwingTime = RSwingTime{speedN};

        param{speedN,1}.LStanceTime = LStanceTime{speedN};
        param{speedN,1}.LStepLength = LStepLength{speedN};
        param{speedN,1}.LStepTime = LStepTime{speedN};
        param{speedN,1}.LStrideLength = LStrideLength{speedN};
        param{speedN,1}.LSwingTime = LSwingTime{speedN};

        param{speedN,1}.CycleTime = CycleTime{speedN};
        param{speedN,1}.StrideLength = StrideLength{speedN};
        param{speedN,1}.StrideWidth = StrideWidth{speedN};
        param{speedN,1}.DSupportMean = DSupportMean{speedN};
        
        % Gather parameters from all healthy subjects in one structure
        paramAll{speedN,subject} = param{speedN,1};
        % Plot time-normalized data to check gait events
%             figure(subject); hold on;
%             subplot(2,7,speedN);
%             plot(kinNorm{speedN,1}.LKneeX);

        %% Compute paramistics (mean and std of spatiotemporal obtained with Visual3D) 
        % Mean over gait cycles + std + coeff of variation
        fields = fieldnames(param{speedN});
        meanAllKin{speedN,1} = structfun(@mean,param{speedN},'UniformOutput',false);
        stdAllKin{speedN,1} = structfun(@(x)(std(x,0,1)),param{speedN},'UniformOutput',false);
        varAllKin{speedN,1} = struct(fields{:});
        for k = 1: size(fields,1)
            varAllKin{speedN,1}.(fields{k}) = stdAllKin{speedN,1}.(fields{k})./meanAllKin{speedN,1}.(fields{k});
        end

        varAllST{speedN,1}.RStanceTime = RStanceTimeStd{speedN}/RStanceTimeMean{speedN};
        varAllST{speedN,1}.RStepLength = RStepLengthStd{speedN}/RStepLengthMean{speedN};
        varAllST{speedN,1}.RStepTime = RStepTimeStd{speedN}/RStepTimeMean{speedN};
        varAllST{speedN,1}.RStrideLength = RStrideLengthStd{speedN}/RStrideLengthMean{speedN};
        varAllST{speedN,1}.RSwingTime = RSwingTimeStd{speedN}/RSwingTimeMean{speedN};

        varAllST{speedN,1}.LStanceTime = LStanceTimeStd{speedN}/LStanceTimeMean{speedN};
        varAllST{speedN,1}.LStepLength = LStepLengthStd{speedN}/LStepLengthMean{speedN};
        varAllST{speedN,1}.LStepTime = LStepTimeStd{speedN}/LStepTimeMean{speedN};
        varAllST{speedN,1}.LStrideLength = LStrideLengthStd{speedN}/LStrideLengthMean{speedN};
        varAllST{speedN,1}.LSwingTime = LSwingTimeStd{speedN}/LSwingTimeMean{speedN};

        varAllST{speedN,1}.CycleTime = CycleTimeStd{speedN}/CycleTimeMean{speedN};
        varAllST{speedN,1}.StrideLength = StrideLengthStd{speedN}/StrideLengthMean{speedN};
        varAllST{speedN,1}.StrideWidth = StrideWidthStd{speedN}/StrideWidthMean{speedN};
        varAllST{speedN,1}.DSupport = DSupportStd{speedN}/DSupportMean{speedN};
        
        %% Save in tables
        namesRow = {'Mean','Std','Coeff var'};

        % Kinetics/kinematics
        RightKinM = [];
        LeftKinM = [];
        counterR = 1;
        counterL = 1;
        for k = 1:size(fields,1)
            dim = size(stdAllKin{speedN}.(fields{k}),2);
            if fields{k}(1) == 'R' % Right side
                RightKinM = [RightKinM [meanAllKin{speedN}.(fields{k}); stdAllKin{speedN}.(fields{k}); varAllKin{speedN}.(fields{k})]];
                if dim == 1
                    RightKinN{counterR} = fields{k};
                elseif dim == 3
                    RightKinN(1,counterR:counterR+2) = {[fields{k},' X'],[fields{k},' Y'],[fields{k},' Z']};
                end
                counterR = counterR + dim;
            elseif fields{k}(1) == 'L' % Left side
                LeftKinM = [LeftKinM [meanAllKin{speedN}.(fields{k}); stdAllKin{speedN}.(fields{k}); varAllKin{speedN}.(fields{k})]];
                if dim == 1
                    LeftKinN{counterL} = fields{k};
                elseif dim == 3
                    LeftKinN(1,counterL:counterL+2) = {[fields{k},' X'],[fields{k},' Y'],[fields{k},' Z']};
                end
                counterL = counterL + dim;
            end
        end

        RKinTable{subject,speedN} = array2table(RightKinM,'RowNames',namesRow,'VariableNames',RightKinN);
        LKinTable{subject,speedN} = array2table(LeftKinM,'RowNames',namesRow,'VariableNames',LeftKinN);
        
        % Spatio Temporal
        % Right
        RightSTempM = [RStanceTimeMean{speedN} RStepLengthMean{speedN} RStepTimeMean{speedN} RStrideLengthMean{speedN} RSwingTimeMean{speedN};...
                        RStanceTimeStd{speedN} RStepLengthStd{speedN} RStepTimeStd{speedN} RStrideLengthStd{speedN} RSwingTimeStd{speedN};...
                        varAllST{speedN,1}.RStanceTime varAllST{speedN,1}.RStepLength varAllST{speedN,1}.RStepTime varAllST{speedN,1}.RStrideLength varAllST{speedN,1}.RSwingTime];
        RightSTempN = {'RStanceTime', 'RStepLength', 'RStepTime', 'RStrideLength', 'RSwingTime'};

        
        RSTTable{subject,speedN} = array2table(RightSTempM,'RowNames',namesRow,'VariableNames',RightSTempN);

        %Left
        LeftSTempM = [LStanceTimeMean{speedN} LStepLengthMean{speedN} LStepTimeMean{speedN} LStrideLengthMean{speedN} LSwingTimeMean{speedN};...
                        LStanceTimeStd{speedN} LStepLengthStd{speedN} LStepTimeStd{speedN} LStrideLengthStd{speedN} LSwingTimeStd{speedN};...
                        varAllST{speedN,1}.LStanceTime varAllST{speedN,1}.LStepLength varAllST{speedN,1}.LStepTime varAllST{speedN,1}.LStrideLength varAllST{speedN,1}.LSwingTime];
        LeftSTempN = {'LStanceTime', 'LStepLength', 'LStepTime', 'LStrideLength', 'LSwingTime'};

        LSTTable{subject,speedN} = array2table(LeftSTempM,'RowNames',namesRow,'VariableNames',LeftSTempN);
        
        % Global
        globalM = [CycleTimeMean{speedN} StrideLengthMean{speedN} StrideWidthMean{speedN} DSupportMean{speedN};...
                CycleTimeStd{speedN} StrideLengthStd{speedN} StrideWidthStd{speedN} DSupportStd{speedN};...
                varAllST{speedN,1}.CycleTime varAllST{speedN,1}.StrideLength varAllST{speedN,1}.StrideWidth varAllST{speedN,1}.DSupport];
        globalN ={'CycleTime', 'StrideLength', 'StrideWidth', 'DSupport'};

        globalSTTable{subject,speedN} = array2table(globalM,'RowNames',namesRow,'VariableNames',globalN);

%         close all
    end
    % save individual data
    save(fileMatlab,'kinSeg','kinSwing','kinStance','kinNorm','eulerT','RArmSw','LArmSw','RFPA','LFPA','param','-append');
end

% Save file with data from all subjects (1 file healthy, 1 file stroke)
if strcmp(type,'healthy')
    save('D:\StimuLOOP\DataGait\NM_Reference\statHealthy.mat','RKinTable','RSTTable','LKinTable','LSTTable','globalSTTable','paramAll');
elseif strcmp (type,'patient')
    save('D:\StimuLOOP\DataGait\NM_GaitSegmentation\statPatient.mat','RKinTable','RSTTable','LKinTable','LSTTable','globalSTTable');
end