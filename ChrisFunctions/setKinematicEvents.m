function [success,events,h] = setKinematicEvents(c3d, saliance, depolarization, terrain, type, plotFlag)

% If there are no function arguments, let the user designate a folder in
% which all trials should be evented with the below specifications
if nargin <1
    saliance = 50; %100
    depolarization = 100; %350 for 1kmh,
    terrain = 2; % 1:overground, 2:treadmill, 3:GRAIL (negative is forwards)
    type = 2; % 1: velocity, 2: position
    plotFlag = 1;

    %     %% Select folder with trials to process
    %     [path] = uigetdir('Choose folder to process trials');
    %     cd(path)
    %
    %     a = dir ('*.c3d');
    %     for i = 1:length(a)
    %         files(i) = cellstr(a(i).name);
    %     end
end


%% Main loop through all files
for j = 1 %1:length(files)
    try

        btkClearEvents(c3d);

        % Extract some shit
        p = btkGetMetaData(c3d, 'SUBJECTS', 'NAMES');
        subject = p.info.values;
        p = btkGetMetaData(c3d, 'SUBJECTS', 'SPEED');
        speed = p.info.values;
        freq = btkGetPointFrequency(c3d);
        ff = btkGetFirstFrame(c3d);
        evts = btkGetEvents(c3d);
        evtsTrue=fields(evts);
        for i = 1:length(evtsTrue)
            eval(sprintf('evts.%s = evts.%s*freq;', evtsTrue{i}, evtsTrue{i}))
            eval(sprintf('evts.%s = evts.%s-ff;', evtsTrue{i}, evtsTrue{i}))
        end



        % Extract marker positions
        markers = btkGetMarkers(c3d);
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

        %         %% if the markers have a subject prefix, erase it!
        %         clear longest test n m skip
        %         test = fieldnames(markers);
        %         for i = 1:length(subject)
        %             longest(i) = length (subject{i});
        %         end
        %         [n,m] = sort(longest, 'descend');
        %
        %
        %         % if ANY of the subjects are found in the markers definition
        %         for i = 1:length(test)
        %             skip = 0;
        %             clear newfield
        %             for ii = 1:length(m)
        %                 if regexp(test{i}, subject{m(ii)}) == 1
        %                     if skip == 0
        %                         newfield = test{i}; newfield(1:length(subject{m(ii)})+1) = [];
        %                         [markers.(newfield)] = markers.(test{i});
        %                         markers = rmfield(markers,test{i});
        %                         skip = 1;
        %                     end
        %                 end
        %             end
        %         end
        %
        %         % change the subject to reflect the shortes denominator
        %         subject = subject{m(end)};



        %% LEFT SIDE
        heel(:,1) = nanmean([markers.LFAL(:,1)  markers.LTAM(:,1)  markers.LFCC(:,1)],2);
        heel(:,2) = nanmean([markers.LFAL(:,2)  markers.LTAM(:,2)  markers.LFCC(:,2)],2);
        heel(:,3) = nanmean([markers.LFAL(:,3)  markers.LTAM(:,3)  markers.LFCC(:,3)],2);
        toe(:,1) = nanmean([markers.LFM1(:,1)  markers.LFM2(:,1)  markers.LFM5(:,1)],2);
        toe(:,2) = nanmean([markers.LFM1(:,2)  markers.LFM2(:,2)  markers.LFM5(:,2)],2);
        toe(:,3) = nanmean([markers.LFM1(:,3)  markers.LFM2(:,3)  markers.LFM5(:,3)],2);
        try
            sacr = markers.SACR;
        catch
            if ismember('LPSI', fieldnames(markers))
                sacr = (markers.LPSI(:,:) + markers.RPSI(:,:))/2;
                if isnan(sacr)
                    sacr = 10000 * ones(length(heel),3);
                end
            else
                sacr = 10000 * ones(length(heel),3);
            end
        end

        % Determine walking direction
        if terrain == 1
            if sacr(1,2) > sacr(end,2)
                inverted = 0;
            else
                inverted = 1;
            end
        elseif terrain == 3
            inverted = -1;
        else
            inverted = -1;
        end

        signalPos = sacr(:,1) - heel(:,1);
        signalVel = sacr(:,1) - toe(:,1);
        signalVel = diff(signalVel);

        % Use position algorithm
        if inverted == 1
            [val, posStrike] = findpeaks(signalPos*-1, 'MINPEAKDISTANCE', depolarization, 'MinPeakProminence', saliance);
            [val, posOff] = findpeaks(signalPos, 'MINPEAKDISTANCE', depolarization, 'MinPeakProminence', saliance);
        else
            % automatic tuning of depolarization and salience parameters
            [val, peaks] = findpeaks(signalPos);
            meanDep = mean(diff(peaks));
            depolarization = meanDep - 100; % for patients % depolarization = meanDep - 10; % for healthy
            meanSal = mean(val);
            saliance = meanSal - 150; % for patients % saliance = meanSal - 50; % for healthy
            if saliance < 0
                saliance = 50;
            end
            [val, posStrike] = findpeaks(signalPos, 'MINPEAKDISTANCE', depolarization, 'MinPeakProminence', saliance);

            [val, posOff] = findpeaks(signalPos*-1, 'MINPEAKDISTANCE', depolarization, 'MinPeakProminence', saliance);
        end

        % use Velocity algorithm
        tmp = crossing(signalVel);

        % remove close crosings
        [a,b] = find (diff(tmp) < 50);
        tmp(b) = [];

        % Sort strikes from Offs
        velStrike = []; velOff = [];
        if inverted == 1
            for ii = 1 : length(tmp)
                if signalVel(tmp(ii)) < signalVel(tmp(ii)+1)
                    velStrike = [velStrike tmp(ii)];
                else
                    velOff = [velOff; tmp(ii)];
                end
            end
        else
            for ii = 1 : length(tmp)
                if signalVel(tmp(ii)) > signalVel(tmp(ii)+1)
                    velStrike = [velStrike tmp(ii)];
                else
                    velOff = [velOff; tmp(ii)];
                end
            end
        end

        % Remove extra strikes
        [a,b] = find (diff(velStrike) < 200);
        velStrike(b+1) = [];

        if plotFlag ==1
            % Double check
            h = figure(1);
            subplot 221
            plot (signalPos);
            title('Position left side')
            a = get(gca, 'ylim');
            lim = a(2) - a(1);
            hold on
            s1 = scatter(posStrike, signalPos(posStrike), 'd', 'm');
            s2 = scatter(posOff, signalPos(posOff), '^', 'm');
            legend([s1 s2], {'FootStrike','FootOff'});
            if length (evtsTrue) == 4
                scatter(round(abs(evts.Left_Foot_Strike)), signalPos(abs(round(evts.Left_Foot_Strike))), 'd', 'g')
                scatter(abs(round(evts.Left_Foot_Off)), signalPos(abs(round(evts.Left_Foot_Off))), '^', 'g')
                title('LEFT Position based events')
                try
                    cut = min(size(posStrike,1) , size(evts.Left_Foot_Strike',1));
                    offset = evts.Left_Foot_Strike(1:cut)'- posStrike(1:cut);
                    writePos = a(1) + 4*(lim-lim/8)/4;
                    val = mean(offset);
                    if abs(val) < 20
                        text(10, writePos, ['LFS: ' num2str(val) ], 'BackgroundColor',[.7 .9 .7]);
                    else text(10, writePos, ['LFS: ' num2str(val) ], 'BackgroundColor',[.8 .4 .4]);
                    end
                    cut = min(size(posOff,1) , size(evts.Left_Foot_Off',1));
                    offset = evts.Left_Foot_Off(1:cut)'- posOff(1:cut);
                    writePos = a(1) + 3*(lim-lim/8)/4;
                    val = mean(offset);
                    if abs(val) < 20
                        text(10, writePos, ['LTO: ' num2str(val) ], 'BackgroundColor',[.7 .9 .7]);
                    else text(10, writePos, ['LTO: ' num2str(val) ], 'BackgroundColor',[.8 .4 .4]);
                    end
                catch msgbox('Unequal LFS nr detected...')
                end
            end
            hold off
            clear offset cut

            subplot 223
            plot (signalVel)
            title('Velocity left side')
            a = get(gca, 'ylim');
            lim = a(2) - a(1);
            hold on
            scatter(velStrike, signalVel(velStrike), 'd', 'm')
            scatter(velOff, signalVel(velOff), '^', 'm')
            if length (evtsTrue) == 4
                scatter(abs(round(evts.Left_Foot_Strike)), signalVel(abs(round(evts.Left_Foot_Strike))), 'd', 'g')
                scatter(abs(round(evts.Left_Foot_Off)), signalVel(abs(round(evts.Left_Foot_Off))), '^', 'g')
                title('LEFT Velocity based events')
                try
                    cut = min(size(velStrike,1) , size(evts.Left_Foot_Strike',1));
                    offset = evts.Left_Foot_Strike(1:cut)'- velStrike(1:cut);
                    writePos = a(1) + 4*(lim-lim/8)/4;
                    val = mean(offset);
                    if abs(val) < 20
                        text(10, writePos, ['LFS: ' num2str(val) ], 'BackgroundColor',[.7 .9 .7]);
                    else text(10, writePos, ['LFS: ' num2str(val) ], 'BackgroundColor',[.8 .4 .4]);
                    end
                    cut = min(size(velOff,1) , size(evts.Left_Foot_Off',1));
                    offset = evts.Left_Foot_Off(1:cut)'- velOff(1:cut);
                    writePos = a(1) + 3*(lim-lim/8)/4;
                    val = mean(offset);
                    if abs(val) < 20
                        text(10, writePos, ['LTO: ' num2str(val) ], 'BackgroundColor',[.7 .9 .7]);
                    else text(10, writePos, ['LTO: ' num2str(val) ], 'BackgroundColor',[.8 .4 .4]);
                    end
                catch msgbox('Unequal LTO nr detected...')
                end
            end
            hold off
        else
            h = 0;
        end

        % Set Events, defined as a frame number (type 1 == velocity)
        if type == 1
            events.LStrike = velStrike+ff;
            events.LOff = velOff+ff;
        else
            events.LStrike = posStrike+ff;
            events.LOff = posOff+ff;
        end
        events.signalPosLeft = signalPos;
        events.signalVelLeft = signalVel;
        % Do some cleanup
        clear posStrike posOff heel toe signalPos tmp


        %% RIGHT SIDE
        heel(:,1) = nanmean([markers.RFAL(:,1)  markers.RTAM(:,1)  markers.RFCC(:,1)],2);
        heel(:,2) = nanmean([markers.RFAL(:,2)  markers.RTAM(:,2)  markers.RFCC(:,2)],2);
        heel(:,3) = nanmean([markers.RFAL(:,3)  markers.RTAM(:,3)  markers.RFCC(:,3)],2);
        toe(:,1) = nanmean([markers.RFM1(:,1)  markers.RFM2(:,1)  markers.RFM5(:,1)],2);
        toe(:,2) = nanmean([markers.RFM1(:,2)  markers.RFM2(:,2)  markers.RFM5(:,2)],2);
        toe(:,3) = nanmean([markers.RFM1(:,3)  markers.RFM2(:,3)  markers.RFM5(:,3)],2);
        clear tmp a b
        try
            sacr = markers.SACR;
        catch
            if ismember('LPSI', fieldnames(markers))
                sacr = (markers.LPSI(:,:) + markers.RPSI(:,:))./2;
                if isnan(sacr)
                    sacr = 10000 * ones(length(heel),3);
                end
            else
                sacr = 10000 * ones(length(heel),3);
            end
        end



        signalPos = sacr(:,1) - heel(:,1);
        signalVel = sacr(:,1) - toe(:,1);
        signalVel = diff(signalVel);

        % Use position algorithm
        if inverted == 1
            [val, posStrike] = findpeaks(signalPos*-1, 'MINPEAKDISTANCE', depolarization, 'MinPeakProminence', saliance);
            [val, posOff] = findpeaks(signalPos, 'MINPEAKDISTANCE', depolarization, 'MinPeakProminence', saliance);
        else
            [val, posStrike] = findpeaks(signalPos, 'MINPEAKDISTANCE', depolarization, 'MinPeakProminence', saliance);
            [val, posOff] = findpeaks(signalPos*-1, 'MINPEAKDISTANCE', depolarization, 'MinPeakProminence', saliance);
        end

        % use Velocity algorithm
        tmp = crossing(signalVel);

        % remove close crosings
        [a,b] = find (diff(tmp) < 50);
        tmp(b) = [];

        % Sort strikes from Offs
        velStrike = []; velOff = [];

        if inverted == 1
            for ii = 1 : length(tmp)
                if signalVel(tmp(ii)) < signalVel(tmp(ii)+1)
                    velStrike = [velStrike tmp(ii)];
                else
                    velOff = [velOff tmp(ii)];
                end
            end
        else
            for ii = 1 : length(tmp)
                if signalVel(tmp(ii)) > signalVel(tmp(ii)+1)
                    velStrike = [velStrike tmp(ii)];
                else
                    velOff = [velOff tmp(ii)];
                end
            end
        end

        % Remove extra strikes
        clear a b
        [a,b] = find (diff(velStrike) < 200);
        velStrike(b+1) = [];

        if plotFlag == 1
            % Double Check
            subplot 222
            plot (signalPos)
            title('Position right side')
            a = get(gca, 'ylim');
            lim = a(2) - a(1);
            hold on
            s1 = scatter(posStrike, signalPos(posStrike), 'd', 'm');
            s2 = scatter(posOff, signalPos(posOff), '^', 'm');
            legend([s1 s2], {'FootStrike','FootOff'});
            linkdata on
            brush on
            if length (evtsTrue) == 4
                scatter(abs(round(evts.Right_Foot_Strike)), signalPos(abs(round(evts.Right_Foot_Strike))), 'd', 'g')
                scatter(abs(round(evts.Right_Foot_Off)), signalPos(abs(round(evts.Right_Foot_Off))), '^', 'g')
                title('RIGHT Position based events')
                try
                    cut = min(size(posStrike,1) , size(evts.Right_Foot_Strike',1));
                    offset = evts.Right_Foot_Strike(1:cut)'- posStrike(1:cut);
                    writePos = a(1) + 4*(lim-lim/8)/4;
                    val = mean(offset);
                    if abs(val) < 20
                        text(10, writePos, ['RFS: ' num2str(val) ], 'BackgroundColor',[.7 .9 .7]);
                    else text(10, writePos, ['RFS: ' num2str(val) ], 'BackgroundColor',[.8 .4 .4]);
                    end
                    cut = min(size(posOff,1) , size(evts.Right_Foot_Off',1));
                    offset = evts.Right_Foot_Off(1:cut)'- posOff(1:cut);
                    writePos = a(1) + 3*(lim-lim/8)/4;
                    val = mean(offset);
                    if abs(val) < 20
                        text(10, writePos, ['RTO: ' num2str(val) ], 'BackgroundColor',[.7 .9 .7]);
                    else text(10, writePos, ['RTO: ' num2str(val) ], 'BackgroundColor',[.8 .4 .4]);
                    end
                    text(10, writePos, ['RTO: ' num2str(mean(offset))])
                catch msgbox('Unequal R event nr detected...')
                end
            end
            hold off

            subplot 224
            plot (signalVel)
            title('Velocity right side')
            a = get(gca, 'ylim');
            lim = a(2) - a(1);
            hold on
            scatter(velStrike, signalVel(velStrike), 'd', 'm')
            scatter(velOff, signalVel(velOff), '^', 'm')
            if length (evtsTrue) == 4
                scatter(abs(round(evts.Right_Foot_Strike)), signalVel(abs(round(evts.Right_Foot_Strike))), 'd', 'g')
                scatter(abs(round(evts.Right_Foot_Off)), signalVel(abs(round(evts.Right_Foot_Off))), '^', 'g')
                title('RIGHT Velocity based events')
                try
                    cut = min(size(velStrike,1) , size(evts.Right_Foot_Strike',1));
                    offset = evts.Right_Foot_Strike(1:cut)'- velStrike(1:cut);
                    writePos = a(1) + 4*(lim-lim/8)/4;
                    val = mean(offset);
                    if abs(val) < 20
                        text(10, writePos, ['RFS: ' num2str(val) ], 'BackgroundColor',[.7 .9 .7]);
                    else text(10, writePos, ['RFS: ' num2str(val) ], 'BackgroundColor',[.8 .4 .4]);
                    end
                    cut = min(size(velOff,1) , size(evts.Right_Foot_Off',1));
                    offset = evts.Right_Foot_Off(1:cut)'- velOff(1:cut);
                    writePos = a(1) + 3*(lim-lim/8)/4;
                    val = mean(offset);
                    if abs(val) < 20
                        text(10, writePos, ['RTS: ' num2str(val) ], 'BackgroundColor',[.7 .9 .7]);
                    else text(10, writePos, ['RTS: ' num2str(val) ], 'BackgroundColor',[.8 .4 .4]);
                    end
                catch msgbox('Unequal R event nr detected...')
                end
            end
            hold off
            suptitleBGA([char(subject),' Speed ', char(speed)]);
            %             suptitleBGA(filename(end-33:end));
            %     waitforbuttonpress()
        end

        % Set events
        if type == 1
            events.RStrike = velStrike + ff;
            events.ROff = velOff + ff;
        else
            events.RStrike = posStrike + ff;
            events.ROff = posOff + ff;
        end

        events.signalPosRight = signalPos;
        events.signalVelRight = signalVel;
        % Do some cleanup
        clear posStrike posOff heel toe signalPos tmp

        success(j) = 1;
    catch
        success(j) = 0;
    end

end


end