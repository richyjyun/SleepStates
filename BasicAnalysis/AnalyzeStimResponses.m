%% Analysis of excitatory and inhibitory stimulus response

%% Get data from excel log
clear; close all;

user = getenv('username');

metafile = 'Experiments.xlsx';
opts = detectImportOptions(metafile);
metadata = readtable(metafile,opts);

%% For each experiment
for m = 1:size(metadata,1)
    
    % Skip if no stim
    if(isnan(metadata.Stim(m))), continue; end
    
    % Path logistics
    animal = metadata.Animal{m};
    exp = metadata.Experiment{m};
    
    fprintf('%s - %s\n', animal, exp);
    
    filepath = fullfile('R:\Yun',animal,'Neurochip');
    
    files = dir(fullfile(filepath,exp,'*.mat'));
    filenames = extractfield(files,'name');    
    
    %% Load data
    [fpath,fname,Channels,fs,session_time] = getNCData(filepath,exp);
    
    % Load channels
    % spikechannels
    parampath = fullfile(fpath,'SpikeParams');
    load(fullfile(parampath,'Channels.mat'))
    
    % Load CoD to ignore spikes that are bad
    % coef, coefcont, pval
    spikepath = fullfile(fpath,'Spikes');
    load(fullfile(spikepath,'CoD'));
    
    % Determine channels to skip
    coef = cellfun(@mean,coef); 
    coefcont = cellfun(@mean,coefcont);
    skip = find(coef < 0.8 & coefcont > coef);
    
    % Load sorted states
    % bins, smoothidx
    load(fullfile(fpath,'SortedIdx'));
    
    % Get stimtimes 
    events = nc3events(fname);
    stim = cellfun(@(x) x',events.stim,'uniformoutput',0);
    stim = sort(cell2mat(stim));
    
    % Determine which stimtimes belong to which state
    ind = discretize(stim,[bins(:,1);bins(end,2)]);
    stim(isnan(ind)) = [];
    ind(isnan(ind)) = [];
    state = smoothidx(ind);
    
    %% Loop through each channel and get evoked spikes and inhibition during each state
    fprintf('Going through channels...'); tic;
    dt = cell(4,length(spikechannels));
    prc = cell(4,length(spikechannels));
    norm = cell(4,length(spikechannels));
    inhib = cell(4,length(spikechannels));    
    for c = 1:length(spikechannels)
        
        fprintf('%d...',spikechannels(c));
        
        if(any(c==skip)), continue; end
        
        % Load spike times
        % ts, traces
        spikefile = [num2str(spikechannels(c)),'.mat'];
        load(fullfile(spikepath,spikefile),'ts');
        
        % Loop through each state and find evoked spikes and inhibition
        for s = 1:4
            [~, dt{s,c}, prc{s,c}, norm{s,c}] = getESpikes(ts, stim(state==s));
            if ~isnan(prc{s,c})
                inhib{s,c} = getInhib(dt{s,c},ts, stim(state==s));
            else
                inhib{s,c} = nan;
            end
        end
    end   
    
    save(fullfile(spikepath,'StimResponse.mat'),'dt','prc','norm','inhib');
    
    t = toc;
    fprintf('%f s\n',t);
    
end











