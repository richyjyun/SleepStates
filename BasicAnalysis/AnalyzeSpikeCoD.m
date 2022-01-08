%% Calculates the coefficient of determination 
% Between each spike's first and last 1000 instances and between each
% spike's first 1000 and a random spike's last 1000 as a control

%% Get data from excel log
clear; close all;

user = getenv('username');

metafile = ['C:\Users\',user,'\Dropbox\Fetz Lab\_Brain States\Experiments.xlsx'];
opts = detectImportOptions(metafile);
metadata = readtable(metafile,opts);

packet = ['C:\Users\',user,'\Dropbox\Fetz Lab\_Brain States\Packets\CoD.ps'];

%% For each experiment
for m = 1:size(metadata,1)
    
    % Path logistics
    animal = metadata.Animal{m};
    exp = metadata.Experiment{m};
    
    fprintf('%s - %s\n', animal, exp);
    
    filepath = fullfile('R:\Yun',animal,'Neurochip');
    
    files = dir(fullfile(filepath,exp,'*.mat'));
    filenames = extractfield(files,'name');    
    
    %% Load data
    fprintf('Loading Data...');
    tic;
    [fpath,fname,Channels,fs,session_time] = getNCData(filepath,exp);
    
    % Load channels
    parampath = fullfile(fpath,'SpikeParams');
    load(fullfile(parampath,'Channels.mat'))
    
    % Load first and last 1000 spike traces for each spike
    spikepath = fullfile(fpath,'Spikes');
    allTraces = cell(2,length(spikechannels));
    for c = 1:length(spikechannels)
        spikefile = [num2str(spikechannels(c)),'.mat'];
        load(fullfile(spikepath,spikefile));
        allTraces{1,c} = traces(:,1:1000);
        allTraces{2,c} = traces(:,end-999:end);
    end
    
    % Delete unused variables
    clearvars ts traces
    t = toc;
    fprintf('%f s\n',t);
    
    %% Compare spike waveform from first 1000 to last 1000
    % Get the average of a random set of 100 spikes within the 1000 to do
    % coefficient of determination. Repeat it 1000 times to get a
    % distribution.
    % Compare it with random 100 spikes from another spike
    fprintf('Calculating CoD...');
    tic;
    rep = 1000;
    coef = {};
    coefcont = {};
    for c = 1:size(allTraces,2)
        cod = zeros(1,rep);
        codcont = zeros(1,rep);
        for r = 1:rep
            temp = randperm(1000);
            pre = mean(allTraces{1,c}(:,temp(1:100)),2);
            post = mean(allTraces{2,c}(:,temp(end-99:end)),2);
            
            compareind = 1+mod(c,size(allTraces,2));
            comppost = mean(allTraces{2,compareind}(:,temp(end-99:end)),2);
            
            cod(r) = getCoD(pre,post);
            codcont(r) = getCoD(pre,comppost);
        end
        coef{c} = cod;
        coefcont{c} = codcont;
    end

    % Statistics
    pval = zeros(1,length(coef));
    for c = 1:length(coef)
        if(mean(coef{c}) < mean(coefcont{c}))
            pval(c) = 1;
        else
            [~,pval(c)] = ttest(coef{c},coefcont{c});
        end
    end
    
    % Save
    codfile = fullfile(spikepath,'CoD');
    save(codfile,'coef','coefcont','pval');
    
    %% Plot for sanity check
    figure('visible','off');
    subplot(8,2,1); text(0,0.5,exp, 'Interpreter', 'none'); axis off;
    for c = 1:length(coef)
       subplot(8,2,c+1);
       histogram(coef{c},'edgealpha',0); 
       hold on; histogram(coefcont{c},'edgealpha',0);
       title(spikechannels(c));
    end
    print(packet,'-append','-dpsc2','-fillpage');
    close(gcf);
    
    t = toc;
    fprintf('%f s\n',t);
    
end

callps2pdf(packet);

















