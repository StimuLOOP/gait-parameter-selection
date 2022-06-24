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
    speeds = {'FST0.15' 'PGV0.1' ''; '' '' '';'FST0.7' 'PGV0.4' 'SLW0.4'; 'FST0.15' 'PGV0.05' '';'FST0.75' 'PGV0.55' 'SLW0.3'; 'PGV0.5' 'PGV0.8' 'PGV0.25'; 'FST0.35' 'PGV0.20' 'SLW0.07crop'; 'FST0.3' 'PGV0.17' 'SLW0.08'};
    assist ={'_NH_NS' '_NRH_NS1' '';'' '' ''; '_BH_NS' '_BH_NS' '_BH_NS'; '_BH_NS' '_BH_NS' '_BH_NS'; '_BH_NS' '_BH_NS' '_BH_NS'; '_BH_NS' '_BH_NS' '_BH_NS'; '_BH_NS' '_BH_NS' '_BH_NS'; '_BH_NS' '_BH_NS' '_BH_NS'};
    suffixe = '_MoCgapfilled';
end

for subject = [1 3:7]%[1 3:8]
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
            Ztrunk = markers.CLAV - markers.TV10;
            shoulVec = markers.RSHO - markers.LSHO;
            Xtrunk = cross(shoulVec,Ztrunk);
            Ytrunk = cross(Ztrunk, Xtrunk);

            % Compute rotation matrix
            for i=1:nbSamples
                R(1,1) = dot(Xtrunk(i,:),X); R(1,2) = dot(Xtrunk(i,:),Y);  R(1,3) = dot(Xtrunk(i,:),Z);
                R(2,1) = dot(Ytrunk(i,:),X); R(2,2) = dot(Ytrunk(i,:),Y);  R(2,3) = dot(Ytrunk(i,:),Z);
                R(3,1) = dot(Ztrunk(i,:),X); R(3,2) = dot(Ztrunk(i,:),Y);  R(3,3) = dot(Ztrunk(i,:),Z);
            end

            % Compute Euler angles between the lab and the trunk frames
            eulerT{speedN,1}(i,:) = rotm2eul(R);
        else
            eulerT{speedN,1} = [];
        end

        % Arm swing
        centerPoint = (markers.RASI + markers.LASI)./2;
        LArmSw{speedN} = ((markers.LRSP + markers.LUSP)./2) - centerPoint;
        RArmSw{speedN} = ((markers.RRSP + markers.RUSP)./2) - centerPoint;
    
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
        for e = 1:min(size(events{speedN}.LStrike,2),size(events{speedN}.RStrike,2))-1
            % Segmentation - Step level
            kinSeg{speedN,1}.RAnkle{e,1} = RAnkle{speedN}(events{speedN}.RStrike(e):events{speedN}.RStrike(e+1),:);
            kinSeg{speedN,1}.RHip{e,1} = RHip{speedN}(events{speedN}.RStrike(e):events{speedN}.RStrike(e+1),:);
            kinSeg{speedN,1}.RKnee{e,1} = RKnee{speedN}(events{speedN}.RStrike(e):events{speedN}.RStrike(e+1),:);
            kinSeg{speedN,1}.RKnee{e,1}(:,1) = - kinSeg{speedN,1}.RKnee{e,1}(:,1); % take the opposite so that knee flexion is positive

            kinSeg{speedN,1}.LAnkle{e,1} = LAnkle{speedN}(events{speedN}.LStrike(e):events{speedN}.LStrike(e+1),:);
            kinSeg{speedN,1}.LHip{e,1} = LHip{speedN}(events{speedN}.LStrike(e):events{speedN}.LStrike(e+1),:);
            kinSeg{speedN,1}.LKnee{e,1} = LKnee{speedN}(events{speedN}.LStrike(e):events{speedN}.LStrike(e+1),:);
            kinSeg{speedN,1}.LKnee{e,1}(:,1) = - kinSeg{speedN,1}.LKnee{e,1}(:,1); % take the opposite so that knee flexion is positive
            
            kinSeg{speedN,1}.RArmSw{e,1} = RArmSw{speedN}(events{speedN}.RStrike(e):events{speedN}.RStrike(e+1),:);
            kinSeg{speedN,1}.LArmSw{e,1} = LArmSw{speedN}(events{speedN}.LStrike(e):events{speedN}.LStrike(e+1),:);

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
            
            kinStance{speedN,1}.RArmSw{e,1} = RArmSw{speedN}(events{speedN}.RStrike(e):events{speedN}.ROff(indexROff),:);
            kinStance{speedN,1}.LArmSw{e,1} = LArmSw{speedN}(events{speedN}.LStrike(e):events{speedN}.LOff(indexLOff),:);

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
            
            kinSwing{speedN,1}.RArmSw{e,1} = RArmSw{speedN}(events{speedN}.ROff(e):events{speedN}.RStrike(indexRStrike),:);
            kinSwing{speedN,1}.LArmSw{e,1} = LArmSw{speedN}(events{speedN}.LOff(e):events{speedN}.LStrike(indexLStrike),:);

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
            kinNorm{speedN,1}.RArmSwX(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.RArmSw{e}(:,1),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.RArmSwY(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.RArmSw{e}(:,2),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.RArmSwZ(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.RArmSw{e}(:,3),linspace(1, sizeVar, 101)))';
            % Right foot progression angle
            kinNorm{speedN,1}.RFPA(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.RFPA{e}(:,1),linspace(1, sizeVar, 101)))';
            % Trunk angles on right cycles
            kinNorm{speedN,1}.ReulerTX(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.ReulerT{e}(:,1),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.ReulerTY(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.ReulerT{e}(:,2),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.ReulerTZ(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.ReulerT{e}(:,3),linspace(1, sizeVar, 101)))';
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
            kinNorm{speedN,1}.LArmSwX(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.LArmSw{e}(:,1),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.LArmSwY(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.LArmSw{e}(:,2),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.LArmSwZ(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.LArmSw{e}(:,3),linspace(1, sizeVar, 101)))';
            % Left foot progression angle
            kinNorm{speedN,1}.LFPA(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.LFPA{e}(:,1),linspace(1, sizeVar, 101)))';
            % Trunk angles on left cycles
            kinNorm{speedN,1}.LeulerTX(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.LeulerT{e}(:,1),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.LeulerTY(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.LeulerT{e}(:,2),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.LeulerTZ(:,e) = (interp1([1:sizeVar], kinSeg{speedN,1}.LeulerT{e}(:,3),linspace(1, sizeVar, 101)))';
            % Left forceplate
            sizeForce = size(kinSeg{speedN,1}.L_FP{e}(:,1),1);
            kinNorm{speedN,1}.L_FPX(:,e) = (interp1([1:sizeForce], kinSeg{speedN,1}.L_FP{e}(:,1),linspace(1, sizeForce, 101)))';
            kinNorm{speedN,1}.L_FPY(:,e) = (interp1([1:sizeForce], kinSeg{speedN,1}.L_FP{e}(:,2),linspace(1, sizeForce, 101)))';
            kinNorm{speedN,1}.L_FPZ(:,e) = (interp1([1:sizeForce], kinSeg{speedN,1}.L_FP{e}(:,3),linspace(1, sizeForce, 101)))';
            
            % Compute extrema and ROM within gait cycles/swing phase/stance phase
            stat{speedN,1}.RAnkleMax(e,:) = max(kinSeg{speedN,1}.RAnkle{e,1});
            stat{speedN,1}.RAnkleMin(e,:) = min(kinSeg{speedN,1}.RAnkle{e,1});
            stat{speedN,1}.RAnkleSwMax(e,:) = max(kinSwing{speedN,1}.RAnkle{e,1});
            stat{speedN,1}.RAnkleSwMin(e,:) = min(kinSwing{speedN,1}.RAnkle{e,1});
            stat{speedN,1}.RAnkleStMax(e,:) = max(kinStance{speedN,1}.RAnkle{e,1});
            stat{speedN,1}.RAnkleStMin(e,:) = min(kinStance{speedN,1}.RAnkle{e,1});
            stat{speedN,1}.RAnkleROM(e,:) = stat{speedN,1}.RAnkleMax(e,:) - stat{speedN,1}.RAnkleMin(e,:);
            stat{speedN,1}.RAnkleSwROM(e,:) = stat{speedN,1}.RAnkleSwMax(e,:) - stat{speedN,1}.RAnkleSwMin(e,:);
            stat{speedN,1}.RAnkleStROM(e,:) = stat{speedN,1}.RAnkleStMax(e,:) - stat{speedN,1}.RAnkleStMin(e,:);

            stat{speedN,1}.RHipMax(e,:) = max(kinSeg{speedN,1}.RHip{e,1});
            stat{speedN,1}.RHipMin(e,:) = min(kinSeg{speedN,1}.RHip{e,1});
            stat{speedN,1}.RHipSwMax(e,:) = max(kinSwing{speedN,1}.RHip{e,1});
            stat{speedN,1}.RHipSwMin(e,:) = min(kinSwing{speedN,1}.RHip{e,1});
            stat{speedN,1}.RHipStMax(e,:) = max(kinStance{speedN,1}.RHip{e,1});
            stat{speedN,1}.RHipStMin(e,:) = min(kinStance{speedN,1}.RHip{e,1});
            stat{speedN,1}.RHipROM(e,:) = stat{speedN,1}.RHipMax(e,:) - stat{speedN,1}.RHipMin(e,:);
            stat{speedN,1}.RHipSwROM(e,:) = stat{speedN,1}.RHipSwMax(e,:) - stat{speedN,1}.RHipSwMin(e,:);
            stat{speedN,1}.RHipStROM(e,:) = stat{speedN,1}.RHipStMax(e,:) - stat{speedN,1}.RHipStMin(e,:);

            stat{speedN,1}.RKneeMax(e,:) = max(kinSeg{speedN,1}.RKnee{e,1});
            stat{speedN,1}.RKneeMin(e,:) = min(kinSeg{speedN,1}.RKnee{e,1});
            stat{speedN,1}.RKneeSwMax(e,:) = max(kinSwing{speedN,1}.RKnee{e,1});
            stat{speedN,1}.RKneeSwMin(e,:) = min(kinSwing{speedN,1}.RKnee{e,1});
            stat{speedN,1}.RKneeStMax(e,:) = max(kinStance{speedN,1}.RKnee{e,1});
            stat{speedN,1}.RKneeStMin(e,:) = min(kinStance{speedN,1}.RKnee{e,1});
            stat{speedN,1}.RKneeROM(e,:) = stat{speedN,1}.RKneeMax(e,:) - stat{speedN,1}.RKneeMin(e,:);
            stat{speedN,1}.RKneeSwROM(e,:) = stat{speedN,1}.RKneeSwMax(e,:) - stat{speedN,1}.RKneeSwMin(e,:);
            stat{speedN,1}.RKneeStROM(e,:) = stat{speedN,1}.RKneeStMax(e,:) - stat{speedN,1}.RKneeStMin(e,:);

            stat{speedN,1}.RArmSwMax(e,:) = max(kinSeg{speedN,1}.RArmSw{e,1});
            stat{speedN,1}.RArmSwMin(e,:) = min(kinSeg{speedN,1}.RArmSw{e,1});
            stat{speedN,1}.RArmSwSwMax(e,:) = max(kinSwing{speedN,1}.RArmSw{e,1});
            stat{speedN,1}.RArmSwSwMin(e,:) = min(kinSwing{speedN,1}.RArmSw{e,1});
            stat{speedN,1}.RArmSwStMax(e,:) = max(kinStance{speedN,1}.RArmSw{e,1});
            stat{speedN,1}.RArmSwStMin(e,:) = min(kinStance{speedN,1}.RArmSw{e,1});
            stat{speedN,1}.RArmSwROM(e,:) = stat{speedN,1}.RArmSwMax(e,:) - stat{speedN,1}.RArmSwMin(e,:);
            stat{speedN,1}.RArmSwSwROM(e,:) = stat{speedN,1}.RArmSwSwMax(e,:) - stat{speedN,1}.RArmSwSwMin(e,:);
            stat{speedN,1}.RArmSwStROM(e,:) = stat{speedN,1}.RArmSwStMax(e,:) - stat{speedN,1}.RArmSwStMin(e,:);

            stat{speedN,1}.ReulerTMax(e,:) = max(kinSeg{speedN,1}.ReulerT{e,1});
            stat{speedN,1}.ReulerTMin(e,:) = min(kinSeg{speedN,1}.ReulerT{e,1});
            stat{speedN,1}.ReulerTSwMax(e,:) = max(kinSwing{speedN,1}.ReulerT{e,1});
            stat{speedN,1}.ReulerTSwMin(e,:) = min(kinSwing{speedN,1}.ReulerT{e,1});
            stat{speedN,1}.ReulerTStMax(e,:) = max(kinStance{speedN,1}.ReulerT{e,1});
            stat{speedN,1}.ReulerTStMin(e,:) = min(kinStance{speedN,1}.ReulerT{e,1});
            stat{speedN,1}.ReulerTROM(e,:) = stat{speedN,1}.ReulerTMax(e,:) - stat{speedN,1}.ReulerTMin(e,:);
            stat{speedN,1}.ReulerTSwROM(e,:) = stat{speedN,1}.ReulerTSwMax(e,:) - stat{speedN,1}.ReulerTSwMin(e,:);
            stat{speedN,1}.ReulerTStROM(e,:) = stat{speedN,1}.ReulerTStMax(e,:) - stat{speedN,1}.ReulerTStMin(e,:);

            stat{speedN,1}.R_FPMax(e,:) = max(kinSeg{speedN,1}.R_FP{e,1});
            stat{speedN,1}.R_FPMin(e,:) = min(kinSeg{speedN,1}.R_FP{e,1});
            stat{speedN,1}.R_FPSwMax(e,:) = max(kinSwing{speedN,1}.R_FP{e,1});
            stat{speedN,1}.R_FPSwMin(e,:) = min(kinSwing{speedN,1}.R_FP{e,1});
            stat{speedN,1}.R_FPStMax(e,:) = max(kinStance{speedN,1}.R_FP{e,1});
            stat{speedN,1}.R_FPStMin(e,:) = min(kinStance{speedN,1}.R_FP{e,1});
            stat{speedN,1}.R_FPROM(e,:) = stat{speedN,1}.R_FPMax(e,:) - stat{speedN,1}.R_FPMin(e,:);
            stat{speedN,1}.R_FPSwROM(e,:) = stat{speedN,1}.R_FPSwMax(e,:) - stat{speedN,1}.R_FPSwMin(e,:);
            stat{speedN,1}.R_FPStROM(e,:) = stat{speedN,1}.R_FPStMax(e,:) - stat{speedN,1}.R_FPStMin(e,:);

            stat{speedN,1}.RFPAMax(e,:) = max(kinSeg{speedN,1}.RFPA{e,1});
            stat{speedN,1}.RFPAMin(e,:) = min(kinSeg{speedN,1}.RFPA{e,1});
            stat{speedN,1}.RFPASwMax(e,:) = max(kinSwing{speedN,1}.RFPA{e,1});
            stat{speedN,1}.RFPASwMin(e,:) = min(kinSwing{speedN,1}.RFPA{e,1});
            stat{speedN,1}.RFPAStMax(e,:) = max(kinStance{speedN,1}.RFPA{e,1});
            stat{speedN,1}.RFPAStMin(e,:) = min(kinStance{speedN,1}.RFPA{e,1});
            stat{speedN,1}.RFPAROM(e,:) = stat{speedN,1}.RFPAMax(e,:) - stat{speedN,1}.RFPAMin(e,:);
            stat{speedN,1}.RFPASwROM(e,:) = stat{speedN,1}.RFPASwMax(e,:) - stat{speedN,1}.RFPASwMin(e,:);
            stat{speedN,1}.RFPAStROM(e,:) = stat{speedN,1}.RFPAStMax(e,:) - stat{speedN,1}.RFPAStMin(e,:);

            stat{speedN,1}.LAnkleMax(e,:) = max(kinSeg{speedN,1}.LAnkle{e,1});
            stat{speedN,1}.LAnkleMin(e,:) = min(kinSeg{speedN,1}.LAnkle{e,1});
            stat{speedN,1}.LAnkleSwMax(e,:) = max(kinSwing{speedN,1}.LAnkle{e,1});
            stat{speedN,1}.LAnkleSwMin(e,:) = min(kinSwing{speedN,1}.LAnkle{e,1});
            stat{speedN,1}.LAnkleStMax(e,:) = max(kinStance{speedN,1}.LAnkle{e,1});
            stat{speedN,1}.LAnkleStMin(e,:) = min(kinStance{speedN,1}.LAnkle{e,1});
            stat{speedN,1}.LAnkleROM(e,:) = stat{speedN,1}.LAnkleMax(e,:) - stat{speedN,1}.LAnkleMin(e,:);
            stat{speedN,1}.LAnkleSwROM(e,:) = stat{speedN,1}.LAnkleSwMax(e,:) - stat{speedN,1}.LAnkleSwMin(e,:);
            stat{speedN,1}.LAnkleStROM(e,:) = stat{speedN,1}.LAnkleStMax(e,:) - stat{speedN,1}.LAnkleStMin(e,:);

            stat{speedN,1}.LHipMax(e,:) = max(kinSeg{speedN,1}.LHip{e,1});
            stat{speedN,1}.LHipMin(e,:) = min(kinSeg{speedN,1}.LHip{e,1});
             stat{speedN,1}.LHipSwMax(e,:) = max(kinSwing{speedN,1}.LHip{e,1});
            stat{speedN,1}.LHipSwMin(e,:) = min(kinSwing{speedN,1}.LHip{e,1});
             stat{speedN,1}.LHipStMax(e,:) = max(kinStance{speedN,1}.LHip{e,1});
            stat{speedN,1}.LHipStMin(e,:) = min(kinStance{speedN,1}.LHip{e,1});
            stat{speedN,1}.LHipROM(e,:) = stat{speedN,1}.LHipMax(e,:) - stat{speedN,1}.LHipMin(e,:);
            stat{speedN,1}.LHipSwROM(e,:) = stat{speedN,1}.LHipSwMax(e,:) - stat{speedN,1}.LHipSwMin(e,:);
            stat{speedN,1}.LHipStROM(e,:) = stat{speedN,1}.LHipStMax(e,:) - stat{speedN,1}.LHipStMin(e,:);

            stat{speedN,1}.LKneeMax(e,:) = max(kinSeg{speedN,1}.LKnee{e,1});
            stat{speedN,1}.LKneeMin(e,:) = min(kinSeg{speedN,1}.LKnee{e,1});
            stat{speedN,1}.LKneeSwMax(e,:) = max(kinSwing{speedN,1}.LKnee{e,1});
            stat{speedN,1}.LKneeSwMin(e,:) = min(kinSwing{speedN,1}.LKnee{e,1});
            stat{speedN,1}.LKneeStMax(e,:) = max(kinStance{speedN,1}.LKnee{e,1});
            stat{speedN,1}.LKneeStMin(e,:) = min(kinStance{speedN,1}.LKnee{e,1});
            stat{speedN,1}.LKneeROM(e,:) = stat{speedN,1}.LKneeMax(e,:) - stat{speedN,1}.LKneeMin(e,:);
            stat{speedN,1}.LKneeSwROM(e,:) = stat{speedN,1}.LKneeSwMax(e,:) - stat{speedN,1}.LKneeSwMin(e,:);
            stat{speedN,1}.LKneeStROM(e,:) = stat{speedN,1}.LKneeStMax(e,:) - stat{speedN,1}.LKneeStMin(e,:);

            stat{speedN,1}.LArmSwMax(e,:) = max(kinSeg{speedN,1}.LArmSw{e,1});
            stat{speedN,1}.LArmSwMin(e,:) = min(kinSeg{speedN,1}.LArmSw{e,1});
            stat{speedN,1}.LArmSwSwMax(e,:) = max(kinSwing{speedN,1}.LArmSw{e,1});
            stat{speedN,1}.LArmSwSwMin(e,:) = min(kinSwing{speedN,1}.LArmSw{e,1});
            stat{speedN,1}.LArmSwStMax(e,:) = max(kinStance{speedN,1}.LArmSw{e,1});
            stat{speedN,1}.LArmSwStMin(e,:) = min(kinStance{speedN,1}.LArmSw{e,1});
            stat{speedN,1}.LArmSwROM(e,:) = stat{speedN,1}.LArmSwMax(e,:) - stat{speedN,1}.LArmSwMin(e,:);
            stat{speedN,1}.LArmSwSwROM(e,:) = stat{speedN,1}.LArmSwSwMax(e,:) - stat{speedN,1}.LArmSwSwMin(e,:);
            stat{speedN,1}.LArmSwStROM(e,:) = stat{speedN,1}.LArmSwStMax(e,:) - stat{speedN,1}.LArmSwStMin(e,:);

            stat{speedN,1}.LeulerTMax(e,:) = max(kinSeg{speedN,1}.LeulerT{e,1});
            stat{speedN,1}.LeulerTMin(e,:) = min(kinSeg{speedN,1}.LeulerT{e,1});
            stat{speedN,1}.LeulerTSwMax(e,:) = max(kinSwing{speedN,1}.LeulerT{e,1});
            stat{speedN,1}.LeulerTSwMin(e,:) = min(kinSwing{speedN,1}.LeulerT{e,1});
            stat{speedN,1}.LeulerTStMax(e,:) = max(kinStance{speedN,1}.LeulerT{e,1});
            stat{speedN,1}.LeulerTStMin(e,:) = min(kinStance{speedN,1}.LeulerT{e,1});
            stat{speedN,1}.LeulerTROM(e,:) = stat{speedN,1}.LeulerTMax(e,:) - stat{speedN,1}.LeulerTMin(e,:);
            stat{speedN,1}.LeulerTSwROM(e,:) = stat{speedN,1}.LeulerTSwMax(e,:) - stat{speedN,1}.LeulerTSwMin(e,:);
            stat{speedN,1}.LeulerTStROM(e,:) = stat{speedN,1}.LeulerTStMax(e,:) - stat{speedN,1}.LeulerTStMin(e,:);

            stat{speedN,1}.L_FPMax(e,:) = max(kinSeg{speedN,1}.L_FP{e,1});
            stat{speedN,1}.L_FPMin(e,:) = min(kinSeg{speedN,1}.L_FP{e,1});
            stat{speedN,1}.L_FPSwMax(e,:) = max(kinSwing{speedN,1}.L_FP{e,1});
            stat{speedN,1}.L_FPSwMin(e,:) = min(kinSwing{speedN,1}.L_FP{e,1});
            stat{speedN,1}.L_FPStMax(e,:) = max(kinStance{speedN,1}.L_FP{e,1});
            stat{speedN,1}.L_FPStMin(e,:) = min(kinStance{speedN,1}.L_FP{e,1});
            stat{speedN,1}.L_FPROM(e,:) = stat{speedN,1}.L_FPMax(e,:) - stat{speedN,1}.L_FPMin(e,:);
            stat{speedN,1}.L_FPSwROM(e,:) = stat{speedN,1}.L_FPSwMax(e,:) - stat{speedN,1}.L_FPSwMin(e,:);
            stat{speedN,1}.L_FPStROM(e,:) = stat{speedN,1}.L_FPStMax(e,:) - stat{speedN,1}.L_FPStMin(e,:);

            stat{speedN,1}.LFPAMax(e,:) = max(kinSeg{speedN,1}.LFPA{e,1});
            stat{speedN,1}.LFPAMin(e,:) = min(kinSeg{speedN,1}.LFPA{e,1});
            stat{speedN,1}.LFPASwMax(e,:) = max(kinSwing{speedN,1}.LFPA{e,1});
            stat{speedN,1}.LFPASwMin(e,:) = min(kinSwing{speedN,1}.LFPA{e,1});
            stat{speedN,1}.LFPAStMax(e,:) = max(kinStance{speedN,1}.LFPA{e,1});
            stat{speedN,1}.LFPAStMin(e,:) = min(kinStance{speedN,1}.LFPA{e,1});
            stat{speedN,1}.LFPAROM(e,:) = stat{speedN,1}.LFPAMax(e,:) - stat{speedN,1}.LFPAMin(e,:);
            stat{speedN,1}.LFPASwROM(e,:) = stat{speedN,1}.LFPASwMax(e,:) - stat{speedN,1}.LFPASwMin(e,:);
            stat{speedN,1}.LFPAStROM(e,:) = stat{speedN,1}.LFPAStMax(e,:) - stat{speedN,1}.LFPAStMin(e,:);
        
        end

        %% Compute statistics (mean and std of spatiotemporal obtained with Visual3D) 
        % Mean over gait cycles + std + coeff of variation
        fields = fieldnames(stat{speedN});
        meanAllKin{speedN,1} = structfun(@mean,stat{speedN},'UniformOutput',false);
        stdAllKin{speedN,1} = structfun(@(x)(std(x,0,1)),stat{speedN},'UniformOutput',false);
        varAllKin{speedN,1} = struct(fields{:});
        for k = 1: size(fields,1)
            varAllKin{speedN,1}.(fields{k}) = stdAllKin{speedN,1}.(fields{k})./meanAllKin{speedN,1}.(fields{k});
        end

        varAllST{speedN,1}.RStanceTime = RStanceTimeMean{speedN}/RStanceTimeStd{speedN};
        varAllST{speedN,1}.RStepLength = RStepLengthMean{speedN}/RStepLengthStd{speedN};
        varAllST{speedN,1}.RStepTime = RStepTimeMean{speedN}/RStepTimeStd{speedN};
        varAllST{speedN,1}.RStrideLength = RStrideLengthMean{speedN}/RStrideLengthStd{speedN};
        varAllST{speedN,1}.RSwingTime = RSwingTimeMean{speedN}/RSwingTimeStd{speedN};

        varAllST{speedN,1}.LStanceTime = LStanceTimeMean{speedN}/LStanceTimeStd{speedN};
        varAllST{speedN,1}.LStepLength = LStepLengthMean{speedN}/LStepLengthStd{speedN};
        varAllST{speedN,1}.LStepTime = LStepTimeMean{speedN}/LStepTimeStd{speedN};
        varAllST{speedN,1}.LStrideLength = LStrideLengthMean{speedN}/LStrideLengthStd{speedN};
        varAllST{speedN,1}.LSwingTime = LSwingTimeMean{speedN}/LSwingTimeStd{speedN};

        varAllST{speedN,1}.CycleTime = CycleTimeMean{speedN}/CycleTimeStd{speedN};
        varAllST{speedN,1}.StrideLength = StrideLengthMean{speedN}/StrideLengthStd{speedN};
        varAllST{speedN,1}.StrideWidth = StrideWidthMean{speedN}/StrideWidthStd{speedN};
        varAllST{speedN,1}.DSupport = DSupportMean{speedN}/DSupportStd{speedN};
        
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

        close all
    end
    % save individual data
    save(fileMatlab,'kinSeg','kinSwing','kinStance','kinNorm','eulerT','RArmSw','LArmSw','RFPA','LFPA','stat','-append');
end

% Save file with data from all subjects (1 file healthy, 1 file stroke)
if strcmp(type,'healthy')
    save('D:\StimuLOOP\DataGait\statHealthy.mat','RKinTable','RSTTable','LKinTable','LSTTable','globalSTTable');
elseif strcmp (type,'patient')
    save('D:\StimuLOOP\DataGait\statPatient.mat','RKinTable','RSTTable','LKinTable','LSTTable','globalSTTable');
end