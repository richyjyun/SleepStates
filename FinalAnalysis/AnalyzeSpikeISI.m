%% Get data from excel log
clear; close all;

metafile = 'Experiments.xlsx';
opts = detectImportOptions(metafile);
metadata = readtable(metafile,opts);

AutoCor = {};
ISI = {};
ind = 1;
bw = 0.0005;
Traces = {};
tracefs = 20000;
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
    
    parampath = fullfile(fpath,'SpikeParams');
    load(fullfile(parampath,'Channels.mat'));
    
    % Determine spikes to skip
    spikepath = fullfile(fpath,'Spikes');
    load(fullfile(spikepath,'CoD'));
    
    coef = cellfun(@mean,coef);
    coefcont = cellfun(@mean,coefcont);
    skip = find(coef < 0.8 & coefcont > coef);
    
    % Load state
    % bins, smoothidx
    load(fullfile(fpath,'SortedIdx'));
        
    for c = 1:length(spikechannels)
        
        fprintf('%d\n',spikechannels(c));
        tic;
        
        if(any(c==skip)), continue; end
        
        spikefile = [num2str(spikechannels(c)),'.mat'];
       
        % 'ts', 'traces'
        load(fullfile(spikepath,spikefile),'ts','traces');
        
        Traces{ind} = mean(traces,2);
        
        ts(ts < bins(1,1)) = [];
        binid = arrayfun(@(x) find(bins(:,1) <= ts(x), 1, 'last'),1:length(ts));
        
        states = smoothidx(binid);
        
        for s = 1:4
            [AutoCor{ind,s},lags] = CrossCorr(ts(states==s),'ts2',ts(states==s),'binsize', 0.0005,'lag',[-0.2,0.2],'suppress_plot',1);
            temp = diff(ts(states==s));
            ISI{ind,s} = histcounts(temp,0:bw:2);
        end
        
        ind = ind+1;
        
        t = toc;
        fprintf('%f seconds\n',t);
        
    end
    
end
temppath = ['_Brain States'];
save(fullfile(temppath,'ISI'),'ISI','bw','AutoCor','Traces','lags','tracefs','-v7.3');

%% Just plotting 
packet = 'ISI.ps';
for i = 1:length(ISI)
    disp(i)
    figure('visible','off');
    for j = 1:4
        temp = ISI{i,j}./sum(ISI{i,j});
        plot(bins*1000,smooth(temp),'linewidth',2); hold on;
    end
    xlim([0,100]);
    print('-painters',packet,'-append','-dpsc','-fillpage');
end

% Plot example ISIs
spks = [43, 45, 26];
figure;
for i = 1:3
    subplot(3,3,i);
    for j = 1:4
        temp = ISI{spks(i),j}./sum(ISI{spks(i),j})*100;
        plot(bins*1000,smooth(temp),'linewidth',2); hold on;
    end
    xlim([0,100]); 
    xlabel('ms'); ylabel('ISI (%)');
    set(gca,'FontSize',10); box off;
    if(i==1)
        legend({'Move','Rest','REM','NREM'},'box','off');
    end
end

%% Classify into fast or regular spiking w/ ISI distribution
% Classify with ISI
bins = 0:bw:(length(ISI{1,1})-1)*bw;
rlim = find(bins>=0.1,1);
type = [];
bad = 13:16;
peaks = [];
for i = 1:length(ISI)
    if(any(i==bad))
        type(i) = 0;
        continue;
    end
    temp = smoothdata(ISI{i,1},'gaussian',20); 
    temp = temp(1:rlim);
%     [~,pks] = findpeaks(temp);
    [~,pks] = max(temp);
    peaks(i) = bins(pks);
    if(length(pks)>1)
        type(i) = 3;
    else
        if(bins(pks) < 0.01)
            type(i) = 2;
        else
            type(i) = 1;
        end
    end
end
figure; 
subplot(3,2,5);
histogram(peaks(peaks<=0.01)*1000,(0:0.002:0.05)*1000,'edgecolor',[0.4,0.4,0.4]);
hold on; histogram(peaks(peaks>0.01)*1000,(0:0.002:0.05)*1000,'edgecolor',[0.4,0.4,0.4]);
xlim([0,50]); box off; ylabel('Spikes'); xlabel('ISI peak (ms)');
set(gca,'fontsize',10); legend({'FS','RS'},'box','off');

% Classify with spike shape
bad = 13:16;
[~,maxind] = cellfun(@max,Traces);
[~,minind] = cellfun(@min,Traces);
width = (maxind-minind)./tracefs*1000;
type(width<=0.4) = 2; type(width>0.4) = 1;
type(bad) = 0;


%% Plot ISIs and their differences for classified spikes
% Average ISI distribution of each
figure; 
titles = {'Regular spiking','Fast spiking','Combined'};
xl = [0,150];
for t = 1:3
    if(t~=3)
        ind = type==t;
    else
        ind = 1:length(ISI);
        ind(bad) = [];
    end
    subplot(3,3,t);
    for s = 1:4
        temp = ISI(ind,s);
        temp = cell2mat(temp);
        temp = temp./sum(temp,2);
        temp = smoothdata(temp','gaussian',5)';
        plot(bins*1000,mean(temp)*100,'linewidth',2); hold on;
    end
    xlim(xl); title(titles{t}); ylabel('%');
    xlabel('ms'); set(gca,'fontsize',10); box off;
end

% Average differences of each
colors = get(gca,'colororder');
xl1 = [0,150];
xl2 = [0,200];
for t = 1:3
    if(t~=3)
        ind = type==t;
    else
        ind = 1:length(ISI);
        ind(bad) = [];
    end
    base = ISI(ind,1);
    base = cell2mat(base);
    base = base./sum(base,2);
    base = smoothdata(base','gaussian',5)';
    base = mean(base);
    
    subplot(3,3,3+t);
    plot(xl1,[0,0],'--','color',colors(1,:),'linewidth',2); hold on;
    subplot(3,3,6+t);
    plot(xl2,[0,0],'--','color',colors(1,:),'linewidth',2); hold on;
    
    for s = 2:4
        temp = ISI(ind,s);
        temp = cell2mat(temp);
        temp = temp./sum(temp,2);
        temp = smoothdata(temp','gaussian',5)';
        temp = temp-base;
        
        subplot(3,3,3+t);
        plot(bins*1000,mean(temp)*100,'linewidth',2,'color',colors(s,:)); hold on;
        xlim(xl1); ylabel('% difference'); xlabel('ms'); set(gca,'fontsize',10); box off;

        subplot(3,3,6+t);
        plot(1./bins,mean(temp)*100,'linewidth',2,'color',colors(s,:)); hold on;
        xlim(xl2); ylabel('% difference'); xlabel('Hz'); set(gca,'fontsize',10); box off;
    end
end


%% Compare ISIs by subtracting?
bins = 0:bw:(length(ISI{1,1})-1)*bw;
bad = 13:16; % weird ISIs

% Plot average ISI
subplot(3,2,3);
for j = 1:4
    temp = ISI(:,j);
    temp(bad) = [];
    temp = cell2mat(temp);
    temp = temp./sum(temp,2);
    temp = smoothdata(temp','gaussian',5)';
    plot(bins*1000,mean(temp)*100,'linewidth',2); hold on;
end
xlim([0,100]); xlabel('ms'); ylabel('%'); box off;
ylim([0,0.8]); set(gca,'FontSize',10);

subplot(3,2,4);
for j = 1:4
    temp = ISI(:,j);
    temp(bad) = [];
    temp = cell2mat(temp);
    temp = temp./sum(temp,2);
    temp = smoothdata(temp','gaussian',5)';
    plot(1./bins,mean(temp)*100,'linewidth',2); hold on;
end
xlim([0,200]); xlabel('Hz'); ylabel('%'); box off; ylim([0,0.8]); 
set(gca,'FontSize',10);

% Plot difference in ISIs
Sub = {};
for i = 1:length(ISI)
    for j = 2:4
        Sub{i,j} = (ISI{i,j}./sum(ISI{i,j}))-(ISI{i,1}./sum(ISI{i,1})); 
    end
end
Sub(bad,:) = [];

subplot(3,2,5);
colors = get(gca,'colororder');
plot([0,150],[0,0],'--','Color',colors(1,:),'linewidth',2); hold on;
for i = 2:4
    temp = Sub(:,i);
    temp = cell2mat(temp);
    plot(bins*1000,smooth(mean(temp))*100,'Color',colors(i,:),'linewidth',2);
%     plot(lags,mean(temp,2),'Color',colors(i,:));
    hold on;
end
ylabel('% Difference'); xlim([0,150]); xlabel('ms'); ylim([-.3,0.101]); box off;
set(gca,'FontSize',10);

subplot(3,2,6);
plot([0,200],[0,0],'--','Color',colors(1,:),'linewidth',2); hold on;
for i = 2:4
    temp = Sub(:,i);
    temp = cell2mat(temp);
    plot(1./bins,smooth(mean(temp))*100,'Color',colors(i,:),'linewidth',2);
%     plot(lags,mean(temp,2),'Color',colors(i,:));
    hold on; xl = xlim; 
end
ylabel('% Difference'); xlim([0,200]); xlabel('Hz'); ylim([-.3,0.101]); box off;
set(gca,'FontSize',10);












