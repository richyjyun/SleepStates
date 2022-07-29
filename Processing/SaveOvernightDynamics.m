%% Saves spectra, average pairwise coherence, and EOG metrics for each experiment

%% Get data from excel log
clear; close all;

user = getenv('username');

metafile = 'Experiments.xlsx';
opts = detectImportOptions(metafile);
metadata = readtable(metafile,opts);

%% For each experiment
for i = 1:size(metadata,1)
    
    % Path logistics
    animal = metadata.Animal{i};
    exp = metadata.Experiment{i};
    
    fprintf('%s - %s\n', animal, exp);
    
    filepath = fullfile('R:\Yun',animal,'Neurochip');
    
    files = dir(fullfile(filepath,exp,'*.mat'));
    filenames = extractfield(files,'name');    
    
    % Was there stimulation
    stim = ~isnan(metadata.Stim(i));
    
    % Load metadata
    [fpath,fname,Channels,fs,session_time] = getNCData(filepath,exp);

    % Data was bad after certain times for some days
    if(~isnan(metadata.Finish(i)))
        session_time = metadata.Finish(i)*60*60;
    end

    % Last 2 channels reserved for EOGs in these experiments
    EOG = [];
    if ~isnan(metadata.EOG(i)) 
        EOG = Channels(15:16); 
        Channels = Channels(1:14);
    end

    % Sampling rate at 1kHz for LFP data
    dwnfs = 1000;
    
    % Save spectra in each epoch
    if ~any(cellfun(@(x) strcmp(x,'Spectra.mat'), filenames))
        fprintf('Saving Spectra...'); tic;
        getSpectra(fpath,fname,Channels,session_time,dwnfs,stim);
        t = toc;
        fprintf('%f s\n',t);
    end
    
    % Save coherence in each epoch
    if ~any(cellfun(@(x) strcmp(x,'Coherence.mat'), filenames))
        fprintf('Saving Coherence...');  tic;
        getCoherence(fpath,fname,Channels,session_time,dwnfs,stim);
        t = toc;
        fprintf('%f s\n',t);
    end
    
    % Save EOG data for days with EOG recordings
    if ~isnan(metadata.EOG(i)) && ~any(cellfun(@(x) strcmp(x,'EOG.mat'), filenames))
        fprintf('Saving EOG...');  tic;
        getEOG(fpath,fname,EOG,session_time,dwnfs,stim);
        t = toc;
        fprintf('%f s\n',t);
    end
    
end


