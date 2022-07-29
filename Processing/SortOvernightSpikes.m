%% Saves all sorting parameters for each dataset then sorts the spikes

%% Get data from excel log
clear; close all;

user = getenv('username');

metafile = 'Experiments.xlsx';
opts = detectImportOptions(metafile);
metadata = readtable(metafile,opts);

%% For each experiment
for m = 1:size(metadata,1)
    % Path logistics
    animal = metadata.Animal{m};
    exp = metadata.Experiment{m};
    
    fprintf('\n\n%s - %s\n', animal, exp);
    
    filepath = fullfile('R:\Yun',animal,'Neurochip');
    
    % Load metadata
    [fpath,fname,Channels,fs,session_time] = getNCData(filepath,exp);
    events = nc3events(fname);
    
    % Was there stimulation
    stim = ~isnan(metadata.Stim(m));
    if(stim) % events are sometimes saved even if there's no stim
        stim = events.stim{1};
    else
        stim = [];
    end
    
    % Data was bad after 17.6 hours on this day
    if(strcmp(exp,'Kronk_20191118_01'))
        session_time = 17.6*60*60;
    end

    % Last 2 channels reserved for EOGs in these experiments
    if(~isnan(metadata.EOG(m))), Channels = Channels(1:14); end

    % Get sorting parameters
    fprintf('Getting sorting parameters...\n');
    channelsave = fullfile(savepath,'Channels.mat');
    if(~exist(channelsave))
        getSortWindows(fpath,fname,Channels,fs,stim);
    end
    fprintf('Done\n');
    
    % Sort spikes throughout the night
    fprintf('Sorting all spikes...\n');
    getAllSpikes(fpath,fname,session_time,fs,stim);
    fprintf('Done\n');
    
end

