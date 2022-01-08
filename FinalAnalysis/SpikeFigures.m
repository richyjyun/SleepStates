%% Spike example and its firing rate over time
path = 'R:\Yun\Kronk\Neurochip';
day = 'Kronk_20191128_01';
[fpath,fname,Channels,fs,session_time] = getNCData(path,day);

parampath = fullfile(fpath,'SpikeParams');
load(fullfile(parampath,'Channels.mat'));

figure;
%% 7 and 8 are good waveform examples
c = 8;
spikefile = [num2str(spikechannels(c)),'.mat'];
load(fullfile(spikepath,spikefile));

subplot(2,2,1); 
inds = randperm(length(traces));
inds = inds(1:1000);
x = (0:size(traces,1)-1)/fs*1000;
stddev = std(traces,[],2);
plot(x,traces(:,inds)','color',[0.7,0.7,0.7]);
hold on; plot(x,mean(traces,2),'k','linewidth',2);
hold on; plot(x,mean(traces,2)+stddev,'r--');
hold on; plot(x,mean(traces,2)-stddev,'r--');
xlabel('Time (ms)'); ylabel('\muV'); box off;
set(gca,'FontSize',12);

%% CoD
spikepath = fullfile(fpath,'Spikes');
codfile = fullfile(spikepath,'CoD');
load(codfile);

c = 8;

subplot(2,2,2);
histogram(coef{c},'edgealpha',0); 
hold on; histogram(coefcont{c},'edgealpha',0);
box off;
xlabel('CoD');
ylabel('Pairs');
legend({'Same','Different'},'box','off');
set(gca,'FontSize',12);

%% Spikerate 
load(fullfile(fpath,'SortedIdx'));
subplot(2,4,5:7); 

tic
chns = [2,4,6,9];
rates = {};
for c = 1:length(chns)
    spikefile = [num2str(spikechannels(chns(c))),'.mat'];
    load(fullfile(spikepath,spikefile),'ts');
    count = arrayfun(@(x) sum(ts >= bins(x,1) & ts <= bins(x,2)), 1:size(bins,1));
    rates{c} = count./diff(bins');
    rates{c} = smooth(rates{c},8);
    if(c~=1)
        rates{c} = rates{c} + max(rates{c-1}) - min(rates{c});
    end
    hold on; plot(bins(:,1)./3600,rates{c},'k');
end
toc; ylim([0,130]);

colors = get(gca,'colororder');
edges = [bins(:,1);bins(end,2)]./3600;
yl = ylim; ylim(yl);
for i = 1:max(smoothidx)
    hold on;
    ind = find(smoothidx==i);
    left = edges(ind)'; right = edges(ind+1)';
    bottom = yl(1)*ones(1,length(ind));
    top = yl(2)*ones(1,length(ind));
    patch([left;left;right;right],[bottom;top;top;bottom],colors(i,:),...
        'edgealpha',0,'facealpha',0.75);
end
xlim([0,max(edges)]);
xlabel('Time (h)');
ylabel('Firing rate (Hz)');
box off;

for c = 1:length(chns)
    hold on; plot(bins(:,1)./3600,rates{c},'k');
end
set(gca,'FontSize',12);

%% Combined
subplot(2,4,8);

allratefile = fullfile(path,'SpikeRateStates.mat');
if(exist(allratefile))
    load(allratefile)
else
    user = getenv('username');
    
    metafile = ['C:\Users\',user,'\Dropbox\Fetz Lab\_Brain States\Experiments.xlsx'];
    opts = detectImportOptions(metafile);
    metadata = readtable(metafile,opts);
    % For each experiment
    staterate = []; sz = 1;
    for m = 1:size(metadata,1)
        
        % Path logistics
        animal = metadata.Animal{m};
        exp = metadata.Experiment{m};
        
        fprintf('%s - %s\n', animal, exp);
        
        filepath = fullfile('R:\Yun',animal,'Neurochip');
        
        files = dir(fullfile(filepath,exp,'*.mat'));
        filenames = extractfield(files,'name');
        
        % Load data
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
        
        % Loop through each channel and get FR for each state
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
                edges = bins(smoothidx==s,:);
                counts = arrayfun(@(x) sum(ts >= edges(x,1) & ts <= edges(x,2)),1:size(edges,1));
                staterate(sz,s) = sum(counts)./sum(diff(edges'));
            end
            sz = sz+1;
        end
        
    end
    
    save(allratefile,'staterate');
    
end

boxplot((staterate-mean(staterate,2))./mean(staterate,2)*100,'notch','on','symbol','w');
xticklabels({'Move','Rest','REM','NREM'});
xl = xlim; hold on; plot(xl,[0,0],'k--');
ylabel('% difference');
box off;
set(gca,'FontSize',12);
ylim([-100,100]);
set(gcf,'renderer','painters');





















