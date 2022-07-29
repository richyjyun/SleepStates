%% Combine all spectra and pairwise-coherence for each state
% Separated by animal

%% Get data from excel log
clear; close all;

user = getenv('username');

metafile = 'Experiments.xlsx';
opts = detectImportOptions(metafile);
metadata = readtable(metafile,opts);

Coh = {};
stddev = {};
N = [];
Spectra = {};
notsaved = [];
packet = 'Coherence.ps';
for m = 1:size(metadata,1)
    
    %% Path logistics
    animal = metadata.Animal{m};
    exp = metadata.Experiment{m};
    
    fprintf('%s - %s\n', animal, exp);
    
    filepath = fullfile('R:\Yun',animal,'Neurochip');
    
    files = dir(fullfile(filepath,exp,'*.mat'));
    filenames = extractfield(files,'name');
    
    % Load data
    [fpath,fname,Channels,fs,session_time] = getNCData(filepath,exp);
    
    load(fullfile(fpath,'Coherence.mat'));
    cohbins = bins; cohf = f;
    
    load(fullfile(fpath,'SortedIdx.mat'));
    sortbins = bins;
    
    % Find corresponding indices
    inds = size(sortbins,1);
    for i = 1:size(sortbins,1)
        temp = find(cohbins(:,2) == sortbins(i,2));
        if(isempty(temp))
            temp = find(cohbins(:,1) == sortbins(i,1));
        end
        if(isempty(temp))
            inds(i) = nan;
        else
            inds(i) = temp;
        end
    end
    
    % Go through each state
    figure('visible','off');
    colors = get(gca,'colororder');
    for s = 1:4
        temp = inds(smoothidx==s);
        temp(isnan(temp)) = [];
        temp = coherence_raw(:,temp);
        Coh{s,m} = nanmean(temp,2);
        stddev{s,m} = nanstd(temp,[],2);
        N(s,m) = size(temp,2);
        
        hold on; 
        plot(f,mean(temp,2),'color',colors(s,:),'linewidth',2);
    end
    
    print(packet,'-append','-dpsc2','-fillpage');
    close(gcf);
    
    % Spectra
    load(fullfile(fpath,'SpectraLong.mat'));
    spectbins = bins; spectf = f;
    
    if(size(spectbins,1) ~= size(spectra,1))
        spectbins = spectbins(1:size(spectra,1),:);
    end
    
    % Find corresponding indices
    inds = size(sortbins,1);
    for i = 1:size(sortbins,1)
        temp = find(cohbins(:,2) == sortbins(i,2));
        if(isempty(temp))
            temp = find(cohbins(:,1) == sortbins(i,1));
        end
        if(isempty(temp))
            inds(i) = nan;
        else
            inds(i) = temp;
        end
    end
    
    % Go through each state
    for s = 1:4
        temp = inds(smoothidx==s);
        temp(isnan(temp)) = [];
        temp = spectra(temp,:);
        Spectra{s,m} = nanmean(temp);
    end
        
end

callps2pdf(packet);

temppath = ['C:\Users\',user,'\Dropbox\Fetz Lab\_Brain States'];
save(fullfile(temppath,'Coherence.mat'),'Coh','stddev','N','cohf');
save(fullfile(temppath,'SpectraAvg.mat'),'Spectra','spectf');

%% Average spectra per animal
f = spectf; 
figure; 
colors = get(gca,'colororder');
for s = 4:-1:1
    subplot(2,2,1);
    temp = cell2mat(Spectra(s,1:16)');
    hold on; patch([f,fliplr(f)],[mean(temp)+std(temp)./4,...
        fliplr(mean(temp)-std(temp)./4)],colors(s,:),...
        'edgealpha',0,'facealpha',0.2);
    plot(f,mean(temp),'color',colors(s,:),'linewidth',1.5);
    xlabel('Frequency (Hz)'); ylabel('dB');
    title('Monkey K');
    xlim([0,120]); ylim([-4,8]);
    set(gca,'FontSize',12);

    subplot(2,2,2);
    temp = cell2mat(Spectra(s,1:17:end)');
           hold on; patch([f,fliplr(f)],[mean(temp)+std(temp)./3,...
       fliplr(mean(temp)-std(temp)./3)],colors(s,:),...
       'edgealpha',0,'facealpha',0.2);
    plot(f,mean(temp),'color',colors(s,:),'linewidth',1.5);
    xlabel('Frequency (Hz)'); ylabel('dB');
    title('Monkey J');
    xlim([0,120]); ylim([-4,8]);
    set(gca,'FontSize',12);
end

%% Average coherence per animal
f = cohf;
colors = get(gca,'colororder');
for s = 4:-1:1
    subplot(2,2,3);
    temp = cell2mat(Coh(s,1:16));
       hold on; patch([f,fliplr(f)],[mean(temp,2)+std(temp,[],2)./4;...
       flipud(mean(temp,2)-std(temp,[],2)./4)],colors(s,:),...
       'edgealpha',0,'facealpha',0.2);
    plot(f,mean(temp,2),'color',colors(s,:),'linewidth',2.5);
    ylim([0.2,0.9])
    xlabel('Frequency (Hz)'); ylabel('Coherence');
    xlim([0,120]);
    set(gca,'FontSize',12);

    subplot(2,2,4);
    temp = cell2mat(Coh(s,1:17:end));
           hold on; patch([f,fliplr(f)],[mean(temp,2)+std(temp,[],2)./3;...
       flipud(mean(temp,2)-std(temp,[],2)./3)],colors(s,:),...
       'edgealpha',0,'facealpha',0.2);
    plot(f,mean(temp,2),'color',colors(s,:),'linewidth',2.5);
    ylim([0.2,0.9])
    xlabel('Frequency (Hz)'); ylabel('Coherence');
    xlim([0,120]);
    set(gca,'FontSize',12);
end











