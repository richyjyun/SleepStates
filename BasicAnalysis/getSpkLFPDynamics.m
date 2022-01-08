function getSpkLFPDynamics(fpath,fname,session_time)
% Calculates spike field coherence, instantaneous amplitude and phase of 
% spike timing, and cross-frequency coupling for each spike channel

%% Load channels
parampath = fullfile(fpath,'SpikeParams');
load(fullfile(parampath,'Channels.mat'))

%% Load states
% idx, smoothidx, bins
bins = []; ts = []; smoothidx = [];
load(fullfile(fpath,'SortedIdx'));

% Find data indices for each state
stateinds = cell(1,4); lfpfs = 1000;
for s = 1:4
    statebins = bins(smoothidx==s,:);
    statebins = round(statebins.*lfpfs);
    
    inds = arrayfun(@(x) (statebins(x,1)+1):statebins(x,2), 1:size(statebins,1), 'UniformOutput', false);
    stateinds{s} = cell2mat(inds);
end

%% Paths
spectpath = fullfile(fpath,'Spectra');
if(~exist(spectpath)), mkdir(spectpath); end

phaseamppath = fullfile(fpath,'PhaseAmp');
if(~exist(phaseamppath)), mkdir(phaseamppath); end

crossfreqpath = fullfile(fpath,'CrossFreq');
if(~exist(crossfreqpath)), mkdir(crossfreqpath); end

%% Loop through channels
spikepath = fullfile(fpath,'Spikes');
f = 0:0.1:200;
for c = 1:length(spikechannels)
    
    %% Logistics
    fprintf('Channel %d\n', spikechannels(c));
    spikefile = [num2str(spikechannels(c)),'.mat'];
%     if(exist(fullfile(spectpath,spikefile))), continue; end
    
    %% Load data
    fprintf('Loading...');
    tic;
    % Load spike times
    load(fullfile(spikepath,spikefile));
    
    % Load whole channel
    [data, ~, ~] = nc3data(spikechannels(c), 0, round(session_time), lfpfs, [0.1,200], fname);

    % Find out which bin each spike is in
    state = zeros(1,length(ts));
    parfor t = 1:length(ts)
        ind = find(ts(t) > bins(:,1) & ts(t) < bins(:,2));
        if(~isempty(ind))
            state(t) = smoothidx(ind);
        end
    end
    
    t = toc;
    fprintf('%f s\n',t);
    
    %% Spike field coherence
    fprintf('Spike field coherence...');
    tic;

    % Define spiketimes and windows for STAs
    trig = round(ts*lfpfs);
    range = round(-1*lfpfs):(1*lfpfs);
    
    %% Get spike triggered trials and find spike field coherence
    % Get STA with the trials then the spectra
    % Get spectra of each trial and the average
    STA = cell(1,4);
    STAspect = cell(1,4);
    Trialspect = cell(1,4);
    
    % Loop through all states
    for s = 1:4
        
        fprintf('%d...',s);
        
        % Spikes that occured in that state
        statetrig = trig(state==s);
        
        % Split up the spikes to avoid running out of memory
        bw = 1e4; %10k spikes seems to be reasonable for speed
        lims = [0:bw:length(statetrig),length(statetrig)];
        
        % Find spike triggered averages
        sta = cell(1,length(lims)-1);
        tspect = cell(1,length(lims)-1);        
        parfor b = 1:(length(lims)-1)
                       
            % Get windows around each spike
            tempTrig = statetrig((lims(b)+1):lims(b+1));
            trials = getTrialinds(tempTrig, range, length(data));
            trials = data(trials);
            
            % Get STA
            sta{b} = sum(trials,2);
            
            % Get spectra
            [tempSpect, ~] = pwelch(trials,[],[],f,lfpfs);
            tspect{b} = sum(tempSpect,2);

        end
        
        % Sum and divide for averaging
        STA{s} = sum(cell2mat(sta),2)./length(statetrig);
        Trialspect{s} = sum(cell2mat(tspect),2)./length(statetrig);
        
        % Get spectra of STA
        temp = pwelch(STA{s},[],[],f,lfpfs);
        STAspect{s} = temp(:);
    end
    
    save(fullfile(spectpath,spikefile), 'trig','range','lfpfs','STA','f','STAspect','Trialspect','-v7.3');
    
    % Delete variables to free up memory
    clearvars trials temp STA STAspect Trialspect
    
    t = toc;
    fprintf('%f s\n',t);
    
    %% Phases and amplitudes
    fprintf('Oscillation phases and amplitudes...');
    tic;
    bands = [0.5,4; 4,8; 8,12; 15,30; 30,70; 70,150];
    names = {'Delta','Theta','Alpha','Beta','Low Gamma','High Gamma'};
    bw = [1,0.5,0.5,0.3,0.2,0.1];
    
    % Save hilbert transform of each band
    % AmpBins are for cross frequency coupling
    AmpBins = cell(1,size(bands,1));
    PhaseBins = -pi:pi/180:pi;
    h = cell(1,size(bands,1));
    for b = 1:size(bands,1)
        filt = bpfilt(data, bands(b,:), lfpfs, 2);
        h{b} = hilbert(filt);
        
        % assign bins
        amp = abs(h{b});
%         [~,edges] = histcounts(amp);
        lim = prctile(amp,99);
        AmpBins{b} = 0:bw(b):lim;
    end
    
    % Just in case
    trig(trig>length(data)) = [];
    
    % Get all phases and amplitudes of each frequency band at the time of
    % spiking
    phase = cell(1,size(bands,1));
    amp = cell(1,size(bands,1));
    for b = 1:size(bands,1)
        temp = angle(h{b}); phase{b} = temp(trig);
        temp = abs(h{b}); amp{b} = temp(trig);
    end
    
    % Save
    save(fullfile(phaseamppath,spikefile),'trig','bands','names','phase','amp','-v7.3');
    t = toc;
    fprintf('%f s\n',t);
    
%     % For plotting later
%     b = 3;
%     temp = histcounts2(phase{b},amp{b},PhaseBins,AmpBins{b});
%     temp = zscore(temp); temp = [temp;temp(end,:)]; temp = imgaussfilt(temp,2); 
%     figure; polarPcolor(AmpBins{b}(1:end-1),rad2deg(PhaseBins),temp)
    
    % Delete variables to free up memory
    clearvars phase amp
    
    %% Cross frequency coupling - using mean vector length (Canolty et al. (2006))
    fprintf('Cross Frequency Coupling...');
    tic;
    
    clearvars filt amp edges
    
    % MVL using amplitude from low freq, phase from high freq, and vice
    % versa
    MVL_PhasetoAmp = zeros(4,size(bands,1),size(bands,1));
    Counts_PhasetoAmp = cell(4,size(bands,1),size(bands,1));
    % Go through each pairwise comparison
    for b1 = 1:(size(bands,1)-1)
        for b2 = (b1+1):size(bands,1)          
            % Loop through each state
            for s = 1:4       
                temp = stateinds{s};
                temp(temp>length(data)) = [];
                
                b1amp = abs(h{b1}(temp));
                b1phase = angle(h{b1}(temp));
                b2amp = abs(h{b2}(temp));
                b2phase = angle(h{b2}(temp));
                
                % Calculate mean vector length
                % Indices are: state, phase band, amplitude band
                MVL_PhasetoAmp(s,b1,b2) = mean(b2amp.*exp(1i.*b1phase));
                MVL_PhasetoAmp(s,b2,b1) = mean(b1amp.*exp(1i.*b2phase));

                % Calculate counts in each bin
                % Indices are: state, phase band, amplitude band
                Counts_PhasetoAmp{s,b1,b2} = histcounts2(b1phase,b2amp,PhaseBins,AmpBins{b2});
                Counts_PhasetoAmp{s,b2,b1} = histcounts2(b2phase,b1amp,PhaseBins,AmpBins{b1});

%                 % For plotting phase amplitude plots later
%                 temp = Counts_PhasetoAmp{s,b1,b2}; temp = zscore(temp); 
%                 temp = [temp;temp(end,:)]; temp = imgaussfilt(temp,2); 
%                 figure; polarPcolor(AmpBins{b2}(1:end-1),rad2deg(PhaseBins),temp)
            end
            
        end
    end
    
    save(fullfile(crossfreqpath,spikefile),'MVL_PhasetoAmp','Counts_PhasetoAmp',...
        'AmpBins','PhaseBins','names','bands','bw','-v7.3');
    
    clearvars h MVL_PhasetoAmp Counts_PhasetoAmp b1amp b1phase b1phase b2phase
    
    t = toc;
    fprintf('%f s\n',t);
end
end