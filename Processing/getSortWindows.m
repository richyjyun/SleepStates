function getSortWindows(fpath,fname,Channels,fs,stim)
%% Loops through each channel to manually determine which channels have spikes
% Then saves all two window discrimination of the spikes
% Assumes 1 spike per channel

% Path all the parameters are saved into
savepath = fullfile(fpath,'SpikeParams');
if(~exist(savepath)), mkdir(savepath); end

% If already sorted, skip
channelsave = fullfile(savepath,'Channels.mat');
if(exist(channelsave)), return; end

% Loading parameters
skip = 6*60*60; % check middle of the night (6 hours into recording) for stable units
dur = 10*60; % check for 10 minutes to sort with

%% Go through each channel, determine channels with spikes, get sorting parameters
spikechannels = [];
for c = 1:length(Channels)
    
    fprintf('Channel %d\n',Channels(c));
        
    % Load data
    data = nc3data(Channels(c),skip, dur, fs, [], fname);
    
    % Filter
    data = HPF(data,fs,1000);
    data = LPF(data,fs,2000);
    
    [spikes, snips, p] = SortSpikes(data,fs,stim);
    % If there's a spike, save channel number and sorting parameters
    if ~isempty(spikes)
        spikechannels(end+1) = Channels(c);
    end
    
    % Save
    paramfile = fullfile(savepath,[num2str(Channels(c)),'.mat']);
    save(paramfile,'spikes','snips','p');
    close all;
end

% Save channels
save(channelsave,'spikechannels');

end