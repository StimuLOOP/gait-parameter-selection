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

for subject = 1%3:7  %[1 3:7]
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

    for speedN = 8%1:sizeSpeed
        if strcmp(type,'healthy')
            filec3d = [folder, day{subject},'_',subjectN,'_TM_',speeds{speedN},'_NH_NS',suffixe{subject}];
        elseif strcmp(type,'patient')
            filec3d = [folder,subjectN,'_TM_',speeds{subject,speedN},assist{subject,speedN},suffixe];
        end
        c3d = btkReadAcquisition([filec3d,'.c3d']);
        markers = btkGetMarkers(c3d);

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
        end

        % Arm swing
        centerPoint = (markers.RASI + markers.LASI)./2;
        LArmSw{speedN} = ((markers.LRSP + markers.LUSP)./2) - centerPoint;
        RArmSw{speedN} = ((markers.RRSP + markers.RUSP)./2) - centerPoint;
    
        %% Correct foot progression angles --> still on process. 
        correctLFP = LFootProg{speedN}(:,2);
        % interpolation of the NaN space
        indexNan = find(isnan(correctLFP));
        if find(diff(indexNan)
        a = find(diff(indexNan)==1);
        indexNan = [indexNan; size(correctLFP,1)];
        for j = 1:2:length(indexNan)-1
            if ismember(j,a) % consecutive NaNs
                
            else
            if abs(correctLFP(indexNan(j)+1,:) - correctLFP(indexNan(j)-1,:)) > 50 
                correctLFP(indexNan(j)+1:indexNan(j+1),:) = - correctLFP(indexNan(j)+1:indexNan(j+1),:);
            end
            end
        end
        % correction sign
        indexJump = find(abs(diff(LFootProg{speedN}(:,2)))>100);
        for j=1:2:length(indexJump)-1
            if correctLFP(indexJump(j)+1) + correctLFP(indexJump(j)) < 1
            %change = LFootProg{speedN}(indexJump(j)+1,3) - correctLFP(indexJump(j),:);
            correctLFP(indexJump(j)+1:indexJump(j+1)) = -LFootProg{speedN}(indexJump(j)+1:indexJump(j+1),2);   
            end
        end
        

        %% Segment kinetic and kinematic parameters with gait events + time normalization + max and mean
        for e = 1:min(size(events{speedN}.LStrike,2),size(events{speedN}.RStrike,2))-1
            % Segmentation
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
            
                % Sample rate of force plates = 600Hz-> rescaling required
            ff = 600; % forceplate frequency
            fmc = 100; % motion capture frequency
            kinSeg{speedN,1}.R_FP{e,1} = FP1{speedN}(events{speedN}.RStrike(e)*ff/fmc:events{speedN}.RStrike(e+1)*ff/fmc,:);
            kinSeg{speedN,1}.L_FP{e,1} = FP2{speedN}(events{speedN}.LStrike(e)*ff/fmc:events{speedN}.LStrike(e+1)*ff/fmc,:);

            if strcmp(type,'healthy')
                kinSeg{speedN,1}.ReulerT{e,1} = eulerT{speedN}(events{speedN}.RStrike(e):events{speedN}.RStrike(e+1),:);
                kinSeg{speedN,1}.LeulerT{e,1} = eulerT{speedN}(events{speedN}.LStrike(e):events{speedN}.LStrike(e+1),:);
            end

            % Time normalization
            % Right ankle
            sizeVar = size(kinSeg{speedN}.RAnkle{e}(:,1),1);
            kinNorm{speedN,1}.RAnkleX(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.RAnkle{e}(:,1),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.RAnkleY(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.RAnkle{e}(:,2),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.RAnkleZ(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.RAnkle{e}(:,3),linspace(1, sizeVar, 101)))';
            % Right hip
            kinNorm{speedN,1}.RHipX(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.RHip{e}(:,1),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.RHipY(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.RHip{e}(:,2),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.RHipZ(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.RHip{e}(:,3),linspace(1, sizeVar, 101)))';
            % Right knee
            kinNorm{speedN,1}.RKneeX(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.RKnee{e}(:,1),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.RKneeY(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.RKnee{e}(:,2),linspace(1, sizeVar, 101)))';
            kinNorm{speedN,1}.RKneeZ(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.RKnee{e}(:,3),linspace(1, sizeVar, 101)))';
            % Right arm swing
            kinNorm{speedN}.RArmSwX(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.RArmSw{e}(:,1),linspace(1, sizeVar, 101)))';
            kinNorm{speedN}.RArmSwY(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.RArmSw{e}(:,2),linspace(1, sizeVar, 101)))';
            kinNorm{speedN}.RArmSwZ(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.RArmSw{e}(:,3),linspace(1, sizeVar, 101)))';
            % Trunk angles on right cycles
            kinNorm{speedN}.ReulerTX(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.ReulerT{e}(:,1),linspace(1, sizeVar, 101)))';
            kinNorm{speedN}.ReulerTY(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.ReulerT{e}(:,2),linspace(1, sizeVar, 101)))';
            kinNorm{speedN}.ReulerTZ(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.ReulerT{e}(:,3),linspace(1, sizeVar, 101)))';
            % Right forceplate
            sizeForce = size(kinSeg{speedN}.R_FP{e}(:,1),1);
            kinNorm{speedN}.R_FPX(:,e) = (interp1([1:sizeForce], kinSeg{speedN}.R_FP{e}(:,1),linspace(1, sizeForce, 101)))';
            kinNorm{speedN}.R_FPY(:,e) = (interp1([1:sizeForce], kinSeg{speedN}.R_FP{e}(:,2),linspace(1, sizeForce, 101)))';
            kinNorm{speedN}.R_FPZ(:,e) = (interp1([1:sizeForce], kinSeg{speedN}.R_FP{e}(:,3),linspace(1, sizeForce, 101)))';
            
            % Left ankle
            sizeVar = size(kinSeg{speedN}.LAnkle{e}(:,1),1);
            kinNorm{speedN}.LAnkleX(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.LAnkle{e}(:,1),linspace(1, sizeVar, 101)))';
            kinNorm{speedN}.LAnkleY(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.LAnkle{e}(:,2),linspace(1, sizeVar, 101)))';
            kinNorm{speedN}.LAnkleZ(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.LAnkle{e}(:,3),linspace(1, sizeVar, 101)))';
            % Left hip
            kinNorm{speedN}.LHipX(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.LHip{e}(:,1),linspace(1, sizeVar, 101)))';
            kinNorm{speedN}.LHipY(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.LHip{e}(:,2),linspace(1, sizeVar, 101)))';
            kinNorm{speedN}.LHipZ(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.LHip{e}(:,3),linspace(1, sizeVar, 101)))';
            % Left knee
            kinNorm{speedN}.LKneeX(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.LKnee{e}(:,1),linspace(1, sizeVar, 101)))';
            kinNorm{speedN}.LKneeY(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.LKnee{e}(:,2),linspace(1, sizeVar, 101)))';
            kinNorm{speedN}.LKneeZ(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.LKnee{e}(:,3),linspace(1, sizeVar, 101)))';
            % Left arm swing
            kinNorm{speedN}.LArmSwX(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.LArmSw{e}(:,1),linspace(1, sizeVar, 101)))';
            kinNorm{speedN}.LArmSwY(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.LArmSw{e}(:,2),linspace(1, sizeVar, 101)))';
            kinNorm{speedN}.LArmSwZ(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.LArmSw{e}(:,3),linspace(1, sizeVar, 101)))';
            % Trunk angles on left cycles
            kinNorm{speedN}.LeulerTX(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.LeulerT{e}(:,1),linspace(1, sizeVar, 101)))';
            kinNorm{speedN}.LeulerTY(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.LeulerT{e}(:,2),linspace(1, sizeVar, 101)))';
            kinNorm{speedN}.LeulerTZ(:,e) = (interp1([1:sizeVar], kinSeg{speedN}.LeulerT{e}(:,3),linspace(1, sizeVar, 101)))';
            % Left forceplate
            sizeForce = size(kinSeg{speedN}.L_FP{e}(:,1),1);
            kinNorm{speedN}.L_FPX(:,e) = (interp1([1:sizeForce], kinSeg{speedN}.L_FP{e}(:,1),linspace(1, sizeForce, 101)))';
            kinNorm{speedN}.L_FPY(:,e) = (interp1([1:sizeForce], kinSeg{speedN}.L_FP{e}(:,2),linspace(1, sizeForce, 101)))';
            kinNorm{speedN}.L_FPZ(:,e) = (interp1([1:sizeForce], kinSeg{speedN}.L_FP{e}(:,3),linspace(1, sizeForce, 101)))';
            
            % Compute extrema
            stat{speedN,1}.maxRAnkle(e,:) = max(kinSeg{speedN,1}.RAnkle{e,1});
            stat{speedN,1}.minRAnkle(e,:) = min(kinSeg{speedN,1}.RAnkle{e,1});

            stat{speedN,1}.maxRHip(e,:) = max(kinSeg{speedN,1}.RHip{e,1});
            stat{speedN,1}.minRHip(e,:) = min(kinSeg{speedN,1}.RHip{e,1});

            stat{speedN,1}.maxRKnee(e,:) = max(kinSeg{speedN,1}.RKnee{e,1});
            stat{speedN,1}.minRKnee(e,:) = min(kinSeg{speedN,1}.RKnee{e,1});

            stat{speedN,1}.maxRArmSw(e,:) = max(kinSeg{speedN,1}.RArmSw{e,1});
            stat{speedN,1}.minRArmSw(e,:) = min(kinSeg{speedN,1}.RArmSw{e,1});

            stat{speedN,1}.maxReulerT(e,:) = max(kinSeg{speedN,1}.ReulerT{e,1});
            stat{speedN,1}.minReulerT(e,:) = min(kinSeg{speedN,1}.ReulerT{e,1});

            stat{speedN,1}.maxR_FP(e,:) = max(kinSeg{speedN,1}.R_FP{e,1});
            stat{speedN,1}.minR_FP(e,:) = min(kinSeg{speedN,1}.R_FP{e,1});

            stat{speedN,1}.maxLAnkle(e,:) = max(kinSeg{speedN,1}.LAnkle{e,1});
            stat{speedN,1}.minLAnkle(e,:) = min(kinSeg{speedN,1}.LAnkle{e,1});

            stat{speedN,1}.maxLHip(e,:) = max(kinSeg{speedN,1}.LHip{e,1});
            stat{speedN,1}.minLHip(e,:) = min(kinSeg{speedN,1}.LHip{e,1});

            stat{speedN,1}.maxLKnee(e,:) = max(kinSeg{speedN,1}.LKnee{e,1});
            stat{speedN,1}.minLKnee(e,:) = min(kinSeg{speedN,1}.LKnee{e,1});

            stat{speedN,1}.maxLArmSw(e,:) = max(kinSeg{speedN,1}.LArmSw{e,1});
            stat{speedN,1}.minLArmSw(e,:) = min(kinSeg{speedN,1}.LArmSw{e,1});

            stat{speedN,1}.maxLeulerT(e,:) = max(kinSeg{speedN,1}.LeulerT{e,1});
            stat{speedN,1}.minLeulerT(e,:) = min(kinSeg{speedN,1}.LeulerT{e,1});

            stat{speedN,1}.maxL_FP(e,:) = max(kinSeg{speedN,1}.L_FP{e,1});
            stat{speedN,1}.minL_FP(e,:) = min(kinSeg{speedN,1}.L_FP{e,1});
        
        end

        %% Compute statistics (mean and std of spatiotemporal obtained with Visual3D) 
        % Mean over gait cycles + std
        meanMaxRAnkle = mean(stat{speedN,1}.maxRAnkle,1);
        stdMaxRAnkle = std(stat{speedN,1}.maxRAnkle,0,1);
        meanMinRAnkle = mean(stat{speedN,1}.minRAnkle,1);
        stdMinRAnkle = std(stat{speedN,1}.minRAnkle,0,1);

        meanMaxRHip = mean(stat{speedN,1}.maxRHip,1);
        stdMaxRHip = std(stat{speedN,1}.maxRHip,0,1);
        meanMinRHip = mean(stat{speedN,1}.minRHip,1);
        stdMinRHip = std(stat{speedN,1}.minRHip,0,1);

        meanMaxRKnee = mean(stat{speedN,1}.maxRKnee,1);
        stdMaxRKnee = std(stat{speedN,1}.maxRKnee,0,1);
        meanMinRKnee = mean(stat{speedN,1}.minRKnee,1);
        stdMinRKnee = std(stat{speedN,1}.minRKnee,0,1);

        meanMaxRArmSw = mean(stat{speedN,1}.maxRArmSw,1);
        stdMaxRArmSw = std(stat{speedN,1}.maxRArmSw,0,1);
        meanMinRArmSw = mean(stat{speedN,1}.minRArmSw,1);
        stdMinRArmSw = std(stat{speedN,1}.minRArmSw,0,1);

        meanMaxReulerT = mean(stat{speedN,1}.maxReulerT,1);
        stdMaxReulerT = std(stat{speedN,1}.maxReulerT,0,1);
        meanMinReulerT = mean(stat{speedN,1}.minReulerT,1);
        stdMinReulerT = std(stat{speedN,1}.minReulerT,0,1);

        meanMaxR_FP = mean(stat{speedN,1}.maxR_FP,1);
        stdMaxR_FP = std(stat{speedN,1}.maxR_FP,0,1);
        meanMinR_FP = mean(stat{speedN,1}.minR_FP,1);
        stdMinR_FP = std(stat{speedN,1}.minR_FP,0,1);

        meanMaxLAnkle = mean(stat{speedN,1}.maxLAnkle,1);
        stdMaxLAnkle = std(stat{speedN,1}.maxLAnkle,0,1);
        meanMinLAnkle = mean(stat{speedN,1}.minLAnkle,1);
        stdMinLAnkle = std(stat{speedN,1}.minLAnkle,0,1);

        meanMaxLHip = mean(stat{speedN,1}.maxLHip,1);
        stdMaxLHip = std(stat{speedN,1}.maxLHip,0,1);
        meanMinLHip = mean(stat{speedN,1}.minLHip,1);
        stdMinLHip = std(stat{speedN,1}.minLHip,0,1);

        meanMaxLKnee = mean(stat{speedN,1}.maxLKnee,1);
        stdMaxLKnee = std(stat{speedN,1}.maxLKnee,0,1);
        meanMinLKnee = mean(stat{speedN,1}.minLKnee,1);
        stdMinLKnee = std(stat{speedN,1}.minLKnee,0,1);

        meanMaxLArmSw = mean(stat{speedN,1}.maxLArmSw,1);
        stdMaxLArmSw = std(stat{speedN,1}.maxLArmSw,0,1);
        meanMinLArmSw = mean(stat{speedN,1}.minLArmSw,1);
        stdMinLArmSw = std(stat{speedN,1}.minLArmSw,0,1);

        meanMaxLeulerT = mean(stat{speedN,1}.maxLeulerT,1);
        stdMaxLeulerT = std(stat{speedN,1}.maxLeulerT,0,1);
        meanMinLeulerT = mean(stat{speedN,1}.minLeulerT,1);
        stdMinLeulerT = std(stat{speedN,1}.minLeulerT,0,1);

        meanMaxL_FP = mean(stat{speedN,1}.maxL_FP,1);
        stdMaxL_FP = std(stat{speedN,1}.maxL_FP,0,1);
        meanMinL_FP = mean(stat{speedN,1}.minL_FP,1);
        stdMinL_FP = std(stat{speedN,1}.minL_FP,0,1);

        % Coefficients of variation
        varMaxRAnkle = stdMaxRAnkle./meanMaxRAnkle;
        varMinRAnkle = stdMinRAnkle./meanMinRAnkle;

        varMaxRHip = stdMaxRHip./meanMaxRHip;
        varMinRHip = stdMinRHip./meanMinRHip;

        varMaxRKnee = stdMaxRKnee./meanMaxRKnee;
        varMinRKnee = stdMinRKnee./meanMinRKnee;

        varMaxRArmSw = stdMaxRArmSw./meanMaxRArmSw;
        varMinRArmSw = stdMinRArmSw./meanMinRArmSw;

        varMaxReulerT = stdMaxReulerT./meanMaxReulerT;
        varMinReulerT = stdMinReulerT./meanMinReulerT;

        varMaxR_FP = stdMaxR_FP./meanMaxR_FP;
        varMinR_FP = stdMinR_FP./meanMinR_FP;

        varMaxLAnkle = stdMaxLAnkle./meanMaxLAnkle;
        varMinLAnkle = stdMinLAnkle./meanMinLAnkle;

        varMaxLHip = stdMaxLHip./meanMaxLHip;
        varMinLHip = stdMinLHip./meanMinLHip;

        varMaxLKnee = stdMaxLKnee./meanMaxLKnee;
        varMinLKnee = stdMinLKnee./meanMinLKnee;

        varMaxLArmSw = stdMaxLArmSw./meanMaxLArmSw;
        varMinLArmSw = stdMinLArmSw./meanMinLArmSw;

        varMaxLeulerT = stdMaxLeulerT./meanMaxLeulerT;
        varMinLeulerT = stdMinLeulerT./meanMinLeulerT;

        varMaxL_FP = stdMaxL_FP./meanMaxL_FP;
        varMinL_FP = stdMinL_FP./meanMinL_FP;

        varRStanceTime = RStanceTimeMean{speedN}/RStanceTimeStd{speedN};
        varRStepLength = RStepLengthMean{speedN}/RStepLengthStd{speedN};
        varRStepTime = RStepTimeMean{speedN}/RStepTimeStd{speedN};
        varRStrideLength = RStrideLengthMean{speedN}/RStrideLengthStd{speedN};
        varRSwingTime = RSwingTimeMean{speedN}/RSwingTimeStd{speedN};

        varLStanceTime = LStanceTimeMean{speedN}/LStanceTimeStd{speedN};
        varLStepLength = LStepLengthMean{speedN}/LStepLengthStd{speedN};
        varLStepTime = LStepTimeMean{speedN}/LStepTimeStd{speedN};
        varLStrideLength = LStrideLengthMean{speedN}/LStrideLengthStd{speedN};
        varLSwingTime = LSwingTimeMean{speedN}/LSwingTimeStd{speedN};

        varCycleTime = CycleTimeMean{speedN}/CycleTimeStd{speedN};
        varStrideLength = StrideLengthMean{speedN}/StrideLengthStd{speedN};
        varStrideWidth = StrideWidthMean{speedN}/StrideWidthStd{speedN};
        varDSupport = DSupportMean{speedN}/DSupportStd{speedN};
        
        %% Save in tables
        namesRow = {'Mean','Std','Coeff var'};
        % Right side
        RightKinM = [meanMaxRAnkle meanMinRAnkle meanMaxRHip meanMinRHip meanMaxRKnee meanMinRKnee meanMaxRArmSw meanMinRArmSw meanMaxReulerT meanMinReulerT meanMaxR_FP meanMinR_FP;...
                        stdMaxRAnkle stdMinRAnkle stdMaxRHip stdMinRHip stdMaxRKnee stdMinRKnee stdMaxRArmSw stdMinRArmSw stdMaxReulerT stdMinReulerT stdMaxR_FP stdMinR_FP;...
                        varMaxRAnkle varMinRAnkle varMaxRHip varMinRHip varMaxRKnee varMinRKnee varMaxRArmSw varMinRArmSw varMaxReulerT varMinReulerT varMaxR_FP varMinR_FP];
        RightKinN = {'MaxRAnkle X', 'MaxRAnkle Y', 'MaxRAnkle Z','MinRAnkle X', 'MinRAnkle Y', 'MinRAnkle Z', ...
                'MaxRHip X', 'MaxRHip Y', 'MaxRHip Z', 'MinRHip X', 'MinRHip Y', 'MinRHip Z', ...
                'MaxRKnee X', 'MaxRKnee Y', 'MaxRKnee Z', 'MinRKNee X', 'MinRKNee Y', 'MinRKNee Z',...
                'MaxRArmSw X', 'MaxRArmSw Y', 'MaxRArmSw Z', 'MinRArmSw X', 'MinRArmSw Y', 'MinRArmSw Z', ...
                'MaxReulerT X', 'MaxReulerT Y', 'MaxReulerT Z', 'MinReulerT X', 'MinReulerT Y','MinReulerT Z',...
                'MaxR_FP X', 'MaxR_FP Y', 'MaxR_FP Z', 'MinR_FP X', 'MinR_FP Y','MinR_FP Z'};
        
        RightSTempM = [RStanceTimeMean{speedN} RStepLengthMean{speedN} RStepTimeMean{speedN} RStrideLengthMean{speedN} RSwingTimeMean{speedN};...
                        RStanceTimeStd{speedN} RStepLengthStd{speedN} RStepTimeStd{speedN} RStrideLengthStd{speedN} RSwingTimeStd{speedN};...
                        varRStanceTime varRStepLength varRStepTime varRStrideLength varRSwingTime];
        RightSTempN = {'RStanceTime', 'RStepLength', 'RStepTime', 'RStrideLength', 'RSwingTime'};

        RKinTable{speedN} = array2table(RightKinM,'RowNames',namesRow,'VariableNames',RightKinN);
        RSTTable{speedN} = array2table(RightSTempM,'RowNames',namesRow,'VariableNames',RightSTempN);

        % Left side
        LeftKinM = [meanMaxLAnkle meanMinLAnkle meanMaxLHip meanMinLHip meanMaxLKnee meanMinLKnee meanMaxLArmSw meanMinLArmSw meanMaxLeulerT meanMinLeulerT meanMaxL_FP meanMinL_FP;...
                        stdMaxLAnkle stdMinLAnkle stdMaxLHip stdMinLHip stdMaxLKnee stdMinLKnee stdMaxLArmSw stdMinLArmSw stdMaxLeulerT stdMinLeulerT stdMaxL_FP stdMinL_FP;...
                        varMaxLAnkle varMinLAnkle varMaxLHip varMinLHip varMaxLKnee varMinLKnee varMaxLArmSw varMinLArmSw varMaxLeulerT varMinLeulerT varMaxL_FP varMinL_FP];
        LeftKinN = {'MaxLAnkle X', 'MaxLAnkle Y', 'MaxLAnkle Z','MinLAnkle X', 'MinLAnkle Y', 'MinLAnkle Z', ...
                'MaxLHip X', 'MaxLHip Y', 'MaxLHip Z', 'MinLHip X', 'MinLHip Y', 'MinLHip Z', ...
                'MaxLKnee X', 'MaxLKnee Y', 'MaxLKnee Z', 'MinLKNee X', 'MinLKNee Y', 'MinLKNee Z',...
                'MaxLArmSw X', 'MaxLArmSw Y', 'MaxLArmSw Z', 'MinLArmSw X', 'MinLArmSw Y', 'MinLArmSw Z', ...
                'MaxLeulerT X', 'MaxLeulerT Y', 'MaxLeulerT Z', 'MinLeulerT X', 'MinLeulerT Y','MinLeulerT Z',...
                'MaxL_FP X', 'MaxL_FP Y', 'MaxL_FP Z', 'MinL_FP X', 'MinL_FP Y','MinL_FP Z'};
        
        LeftSTempM = [LStanceTimeMean{speedN} LStepLengthMean{speedN} LStepTimeMean{speedN} LStrideLengthMean{speedN} LSwingTimeMean{speedN};...
                        LStanceTimeStd{speedN} LStepLengthStd{speedN} LStepTimeStd{speedN} LStrideLengthStd{speedN} LSwingTimeStd{speedN};...
                        varLStanceTime varLStepLength varLStepTime varLStrideLength varLSwingTime];
        LeftSTempN = {'LStanceTime', 'LStepLength', 'LStepTime', 'LStrideLength', 'LSwingTime'};

        LKinTable{speedN} = array2table(LeftKinM,'RowNames',namesRow,'VariableNames',LeftKinN);
        LSTTable{speedN} = array2table(LeftSTempM,'RowNames',namesRow,'VariableNames',LeftSTempN);
        
        % Global
        globalM = [CycleTimeMean{speedN} StrideLengthMean{speedN} StrideWidthMean{speedN} DSupportMean{speedN};...
                CycleTimeStd{speedN} StrideLengthStd{speedN} StrideWidthStd{speedN} DSupportStd{speedN};...
                varCycleTime varStrideLength varStrideWidth varDSupport];
        globalN ={'CycleTime', 'StrideLength', 'StrideWidth', 'DSupport'};

        globalSTTable{speedN} = array2table(globalM,'RowNames',namesRow,'VariableNames',globalN);

    end
end