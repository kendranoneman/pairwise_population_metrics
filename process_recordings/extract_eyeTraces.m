function extract_eyeTraces(datafolder,datafiles,kernel,colors,window_s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for extracting and preprocessing eye traces from .ns5 files
%%%%%%%%%% INPUTS %%%%%%%%%%%
% datafolder = string/char of full path where data is stored
% datafiles = cell array of datafile names
% kernel = size of gaussian kernel for smoothing (bigger = more smoothing)
% colors = 1x2 cell array (rgb color for horizontal and vertical traces)
% window_s = window (in sec) that you want to plot

%%%%%%%%%% OUTPUTS %%%%%%%%%%%
% eyeTraces --> table with eye traces
%    - length_s = length of recording session in seconds
%    - time_s = time points (in sec) to coincide with each datapoint
%    - pos_raw = [HE VE] eye position, after downsampling Fs to 1000 Hz
%    - vel_raw = [HE VE] eye velocity, after downsampling Fs to 1000 Hz
%    - vel_raw = [HE VE] pupil diameter, after downsampling Fs to 1000 Hz
%    - pos_filtered = [HE VE] eye position, after "removing" blinks/artifacts
%                     (found points with velocity > 600 deg/s & interpolated)
%          ...
%    - pos_smoothed = [HE VE] eye position, after "smoothing" filtered data
%                     (filtering data with a Gaussian window, based on kernel input)

%%%%%%%%%% EXAMPLE %%%%%%%%%%%
% e.g. datafolder = '/Users/kendranoneman/Projects/mayo/n64/pb18n64_K_trials';
%      datafiles = {'pb18n64_K_trial50011';'pb18n64_K_trial650012';'pb18n64_K_trial7650013'};
%      colors = {[255,0,0]./255; [0,0,255]./255};
%      window_s = [0 1]
% "extract_eyeTraces(dataFolder,datafiles,10,colors)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varnames = ["session_name","subject_name","trial_number","length_s","time_s","pos_raw","vel_raw","pup_raw","pos_filtered","vel_filtered","pup_filtered","pos_smoothed","vel_smoothed","pup_smoothed"];
tt = cell(length(datafiles),length(varnames));
for df = 1:length(datafiles) % for each datafile
    thisdataset = datafiles{df};

    % Extract data from file
    dat = read_nsx(sprintf('%s/%s.ns5',datafolder,thisdataset),'chanindx',1:3, 'readdata',true);
    seps = strfind(datafiles{df},'_');
    if isequal(thisdataset(seps(2)-1),'K')
        subject = 'kendra';
    elseif isequal(thisdataset(seps(2)-1),'J')
        subject = 'patrick';
    end

    % Filter/smooth data and calculate velocity
    Fs_old = double(dat.hdr.Fs);
    Fs = 1000;    % desired sampling rate
    downsampleFactor = Fs_old / Fs;
    
    pos = downsample(dat.data',downsampleFactor);
    x = (0:1/Fs:(size(pos, 1)-1)/Fs);

    vel = zeros(size(pos));
    for p = 1:size(pos,2)
        vel(:,p) = (gradient(pos(:,p)) ./ gradient(x(:))) ./ 1000;
    end
    
    % removing blinks
    peakIndices = find(abs(vel(:,1))>600 | abs(vel(:,2))>600);
    groups = {};
    currentGroup = {peakIndices(1)};
    for i = 2:numel(peakIndices)
        if peakIndices(i) == peakIndices(i-1) + 1
            currentGroup{end+1} = peakIndices(i);
        else
            groups{end+1} = currentGroup;
            currentGroup = {peakIndices(i)};
        end
    end
    groups{end+1} = currentGroup;
    groups(cellfun(@length, groups, 'UniformOutput', true) == 1) = [];

    filteredVel = vel;
    filteredPos = pos;
    for b = 1:length(groups)
        blinkStart = groups{b}{1}-1;
        blinkEnd = groups{b}{end}+1;

        xq = x(blinkStart-200:blinkEnd+200);
        ipt = findchangepts(pos(blinkStart-200:blinkEnd+200,1),MaxNumChanges=2);
        ipt(1) = ipt(1); ipt(2) = ipt(2);
        for i = 1:size(filteredVel,2)
            vq = vel(blinkStart-200:blinkEnd+200,i)';
            vq1 = interp1([xq(1:ipt(1)-5),xq(ipt(2)+5:end)],[vq(1:ipt(1)-5),vq(ipt(2)+5:end)],xq(ipt(1)-4:ipt(2)+4));
            vq(ipt(1)-4:ipt(2)+4) = vq1;
            filteredVel(blinkStart-200:blinkEnd+200,i) = vq';

            pq = pos(blinkStart-200:blinkEnd+200,i)';
            pq1 = interp1([xq(1:ipt(1)-5),xq(ipt(2)+5:end)],[pq(1:ipt(1)-5),pq(ipt(2)+5:end)],xq(ipt(1)-4:ipt(2)+4));
            pq(ipt(1)-4:ipt(2)+4) = pq1;
            filteredPos(blinkStart-200:blinkEnd+200,i) = pq';
        end
    end
    
    vel_filtered = medfilt1(filteredVel,2,[],1);
    pos_filtered = medfilt1(filteredPos,2,[],1);
    
    pos_smoothed = smoothdata(pos_filtered,1,"gaussian",kernel);
    vel_smoothed = smoothdata(vel_filtered,1,"gaussian",kernel);
    
    tt(df,:) = [{thisdataset(1:seps(1)-1)} {subject} {str2double(thisdataset(seps(2)+6:end))}  {x(end)} {x} {pos(:,1:2)} {vel(:,1:2)} {pos(:,3)} {pos_filtered(:,1:2)} {vel_filtered(:,1:2)} {pos_filtered(:,3)} {pos_smoothed(:,1:2)} {vel_smoothed(:,1:2)} {pos_smoothed(:,3)}];

end

eyeTraces = cell2table(tt,'VariableNames',varnames);
eyeTraces.session_name = categorical(string(eyeTraces.session_name)); eyeTraces.subject_name = categorical(string(eyeTraces.subject_name)); %trialTbl.subject_name = categorical(string(trialTbl.subject_name));
    
save(sprintf('%s/%s_processed.mat',datafolder,thisdataset(1:seps(1)-1)),'eyeTraces','-v7');

%% Plot eye traces
if nargin > 3
    f1 = figure;
    f1.Position = [100 100 900 900];
    tl = tiledlayout(height(eyeTraces),2);
    tl.TileSpacing = 'loose'; 
    tl.Padding = 'compact';
    
    for t = 1:height(eyeTraces)
        nexttile
        x = eyeTraces.time_s{t};
        %xline(x(abs(trialTbl.vel_smoothed{t}(:,1))>100),'k','alpha',0.2)

        plot(eyeTraces.time_s{t}, eyeTraces.pos_smoothed{t}(:,1),'Color',colors{1});
        hold on
        plot(eyeTraces.time_s{t}, eyeTraces.pos_smoothed{t}(:,2),'Color',colors{2});
        xlim(window_s)
        if t==1
            title('position')
        end

        nexttile
        plot(eyeTraces.time_s{t}, eyeTraces.vel_smoothed{t}(:,1),'Color',colors{1});
        hold on
        plot(eyeTraces.time_s{t}, eyeTraces.vel_smoothed{t}(:,2),'Color',colors{2});
        xlim(window_s)
         if t==1
            title('velocity')
        end
    end
    
    xlabel(tl,'time (s)');
    ylabel(tl,'amplitude');
end
end