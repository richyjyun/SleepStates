function getAllSpikes(fpath,fname,session_time,fs,stim)
%% Saves all spikes 
% Load in all of a single channel, sort the spikes, save, and plot 

% Path all the spikes are saved into
savepath = fullfile(fpath,'Spikes');
if(~exist(savepath)), mkdir(savepath); end

% Packet to print data into
user = getenv('username');
temp = strsplit(fpath,'\');
packet = ['C:\Users\',user,'\Dropbox\Fetz Lab\_Brain States\Packets\Spikes\',temp{end},'.ps'];

%% Load bins saved from saveSpectra
load(fullfile(fpath,'Bins.mat'))

%% Load channels
parampath = fullfile(fpath,'SpikeParams');
load(fullfile(parampath,'Channels.mat'))
spikerange = round(-0.0005*fs):1:round(0.0015*fs); % window for detecting spikes

%% Loop through each spike channel
for c = 1:length(spikechannels)
    
    fprintf('%d - Loading data...',spikechannels(c));
    tic;
    
    % File name
    spikefile = [num2str(spikechannels(c)),'.mat'];
    
    %% Load data
    % Load sorting parameters
    load(fullfile(parampath,spikefile));
       
    % Load data and filter
    [data, ~, ~] = nc3data(spikechannels(c), 0, session_time, fs, [], fname);
    if(isempty(data))
        [data, ~, ~] = nc3data(spikechannels(c), 0, round(session_time), fs, [], fname);
    end
    filt = HPF(data,fs,1000);
    filt = LPF(filt,fs,2000);
    t = toc;
    fprintf('%f seconds\n',t);
    
    %% Find all threshold crossings and ones that cross both windows
    fprintf('%d - Sorting data... ',spikechannels(c));
    tic;
    [cross, ind] = GetSnips(filt, p.thresh, spikerange, 1000, stim,  0.0012, fs);
    detected = TwoWindowDiscrim(cross, p);
    
    %% Remove traces that go through the windows but are very off
    template = mean(snips,2);
    stddev = 2*std(snips,[],2);
    traceerr = sum(abs(cross - template));
    keep = traceerr < sum(stddev);
    ts = ind(detected & keep)/fs;
    traces = cross(:,detected & keep);
               
    %% Save
    savefile = fullfile(savepath,spikefile);
    save(savefile,'ts','traces','-v7.3');
    
    %% Plot for debugging
    figure('visible','off');
    x = spikerange/fs*1000;
    
    % From original sorting
    subplot(2,2,1);
    temp1 = mean(snips,2); temp2 = std(snips,[],2);
    patch([x ,fliplr(x)],...
        [temp1+temp2;flipud(temp1-temp2)],'k','facealpha',0.4,'edgealpha',0)
    hold on; plot(x ,temp1,'k','linewidth',2);
    title(spikechannels(c));
    
    % From new sorted spikes
    subplot(2,2,2);
    temp1 = mean(traces,2); temp2 = std(traces,[],2);
    patch([x ,fliplr(x)],...
        [temp1+temp2;flipud(temp1-temp2)],'k','facealpha',0.4,'edgealpha',0)
    hold on; plot(x ,temp1,'k','linewidth',2);
    title(spikechannels(c));
    
    % Firing rate 
    subplot(2,1,2); 
    bins = 0:60:max(ts); counts = histcounts(ts,bins);
    plot(bins(1:end-1)/60/60,counts/60); ylabel('Hz'); xlabel('h');
    xlim([0,bins(end)/60/60]);
    
    % Print to file
    print(packet,'-append','-dpsc2','-fillpage');
    t = toc;
    fprintf('%f seconds\n',t);
    
end

% Convert ps to pdf file
callps2pdf(packet,0,1);


