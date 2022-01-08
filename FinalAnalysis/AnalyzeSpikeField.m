%% Plot figures for spike-field synchrony in each state
% Phase-locking value

%% Get data from excel log
clear; close all;

user = getenv('username');

metafile = ['C:\Users\',user,'\Dropbox\Fetz Lab\_Brain States\Experiments.xlsx'];
opts = detectImportOptions(metafile);
metadata = readtable(metafile,opts);

% Combined data
SFC = {};
PLV = [];
Counts = {};
AmpBins = {};
PhaseBins = -pi:pi/180:pi;
M = [];
C = [];

%% For each experiment
for m = 1:size(metadata,1)
    
    %% Path logistics
    animal = metadata.Animal{m};
    exp = metadata.Experiment{m};
    
    fprintf('%s - %s\n', animal, exp);
    
    filepath = fullfile('R:\Yun',animal,'Neurochip');
    
    files = dir(fullfile(filepath,exp,'*.mat'));
    filenames = extractfield(files,'name');
    
    % Load metadata
    [fpath,fname,Channels,fs,session_time] = getNCData(filepath,exp);
    
    parampath = fullfile(fpath,'SpikeParams');
    load(fullfile(parampath,'Channels.mat'))
    
    spectpath = fullfile(fpath,'Spectra');
    phaseamppath = fullfile(fpath,'PhaseAmp');
    
    %% Load states
    % idx, smoothidx, bins
    bins = []; ts = []; smoothidx = [];
    load(fullfile(fpath,'SortedIdx'));
    
    % Determine spikes to skip
    spikepath = fullfile(fpath,'Spikes');
    load(fullfile(spikepath,'CoD'));
    
    coef = cellfun(@mean,coef);
    coefcont = cellfun(@mean,coefcont);
    skip = find(coef < 0.8 & coefcont > coef);
    
    %% Loop through channels
    for c = 1:length(spikechannels)
        
        fprintf('%d\n',spikechannels(c));
        
        if(any(c==skip)), continue; end
        
        spikefile = [num2str(spikechannels(c)),'.mat'];
        
        fprintf('Calculating...'); tic;
        
        %% Spike Field Coherence
        % trig,range,lfpfs,STA,f,STAspect,Trialspect
        load(fullfile(spectpath,spikefile));
        SFC(end+1,:) = cellfun(@(x,y) x./y, STAspect, Trialspect, 'UniformOutput', false);
        
        %% Phase and amplitudes at spike timing
        % trig,bands,names,phase,amp
        load(fullfile(phaseamppath,spikefile));
        bw = [1,0.5,0.5,0.3,0.2,0.1]; 
        
        ts = trig./lfpfs;
        % Find out which bin each spike is in
        state = zeros(1,length(ts));
        parfor t = 1:length(ts)
            ind = find(ts(t) > bins(:,1) & ts(t) < bins(:,2));
            if(~isempty(ind))
                state(t) = smoothidx(ind);
            end
        end
        
        ind = size(PLV,1);
        
        for s = 1:4
            
            statespike = state==s;
            
            for b = 1:size(bands,1)
                % Phase locking value
                PLV(ind+1,s,b) = getPLV(phase{b}(statespike));
                
                % 2D histogram
                lim = prctile(amp{b}(statespike),99);
                AmpBins{ind+1,s,b} = 0:bw(b):lim;
                Counts{ind+1,s,b} = histcounts2(phase{b}(statespike),amp{b}(statespike),PhaseBins,AmpBins{ind+1,s,b});
            end
            
        end
                
        M(end+1) = m;
        C(end+1) = c;
        
        t = toc;
        fprintf('%f s\n',t);
        
    end
end
temppath = ['C:\Users\',user,'\Dropbox\Fetz Lab\_Brain States'];
save(fullfile(temppath,'SpikeField.mat'),'f','bands','names','bw','SFC','PLV','Counts','AmpBins','PhaseBins','M','C','-v7.3');

%% STA

%% PLV
% Phase histograms of PLVs
figure;
for b = 1:6
    subplot(4,3,b);
    for s = 1:4
        polarhistogram(angle(PLV(:,s,b)),-pi:pi/8:pi,...
            'displaystyle','stairs','linewidth',1.5)
        set(gca,'Thetaaxisunits','radians');
        set(gca,'Gridlinestyle','--');
        thetaticks(0:pi/4:2*pi);
        set(gca,'FontSize',10); 
        hold on;
    end
    title(names{b});
end

% Plot distribution of locked phases per state per frequency band 
% Plot distribution of magnitude of PLV per state per frequency band
direction = [];
magnitude = [];
bad = [];
for b = 1:6
    phases = angle(PLV(:,:,b));
    phases(bad,:) = [];
    direction(b,:) = circMean(phases);
    mags = abs(PLV(:,:,b));
    mags(bad,:) = [];
    magnitude(b,:) = mean(mags);
    subplot(4,3,b+6);
    colors = get(gca,'colororder');
    for s = 1:4
        polarplot([0,direction(b,s)],[0,magnitude(b,s)],'--','color',colors(s,:));
        polarscatter(direction(b,s),magnitude(b,s),50,colors(s,:),'x','linewidth',2);
        set(gca,'Thetaaxisunits','radians');
        set(gca,'Gridlinestyle','--');
        thetaticks(0:pi/4:2*pi);
        set(gca,'FontSize',10);
        hold on;
    end
    title(names{b});
end

set(gcf,'renderer','painters');

%% Counts
%% Generate packets to parse through
for b = 1:6
    packet = ['C:\Users\Richy Yun\Dropbox\Fetz Lab\_Brain States\Packets\SpikeField_',names{b},'.ps'];
    for i = 1:size(Counts,1)
        disp(i);
        figure('visible','off');
        for s = 1:4
            subplot(2,2,s)
            temp = squeeze(Counts{i,s,b});
            temp = (temp)./sum(temp,'omitnan');
            temp = [temp;temp(end,:)];
            temp(isnan(temp)) = 0;
            temp = imgaussfilt(temp,2);
            polarPcolor(AmpBins{i,s,b}(1:end-1),rad2deg(PhaseBins),imgaussfilt(temp,2)')
            if(s==1)
                title([num2str(M(i)),'-',num2str(C(i))]);
            end
        end
        
        print(packet,'-append','-dpsc2','-fillpage');
        close(gcf);
    end
    callps2pdf(packet);
end

%% Normalize the counts for each spike then combine them per state per
% frequency band and plot 
bad = [remove,maybe];
figure;
for b = 1:6
    
    clim = [inf,-inf];
    for s = 1:4
        subplot(6,4,4*(b-1)+s);
        tempcounts = Counts(:,s,b);
        tempcounts(bad) = [];
        sz = max(cellfun(@(x) size(x,2), tempcounts));
        total = nan(length(tempcounts),361,sz);
        for i = 1:length(tempcounts)
            temp = squeeze(tempcounts{i});
            temp = [temp;temp(end,:)];
            total(i,:,1:size(temp,2)) = temp;
        end
        
        temp = sum(total,1,'omitnan'); 
        temp = sum(temp,2,'omitnan');
        temp = squeeze(temp);
        temp = cumsum(temp);
        keep = find(temp>0.90*temp(end),1);

        total = total(:,:,1:keep);
        
        norm = sum(total,2,'omitnan');
        total = total./norm;
        
        total = squeeze(nanmean(total,1));
        
        ampbin = 0:bw(b):(size(total,2).*bw(b));
        polarPcolor(ampbin(1:end-1),rad2deg(PhaseBins),imgaussfilt(total,2)')
        
        temp = caxis;
        clim(1) = min([clim(1),temp(1)]);
        clim(2) = max([clim(2),temp(2)]);
        
        if(s==1), title(names{b}); end
    end
    
    for s = 1:4
        subplot(6,4,4*(b-1)+s);
        caxis(clim);
    end
    
end
set(gcf,'renderer','painters');









