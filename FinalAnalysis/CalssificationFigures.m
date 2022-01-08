%% Classification example
% diagram showing process of classification
% Example of LFP traces during different states

path = 'R:\Yun\Kronk\Neurochip';
day = 'Kronk_20191128_01';
[fpath,fname,Channels,fs,session_time] = getNCData(path,day);

spectrafile = fullfile(fpath,'Spectra.mat');
load(spectrafile);
load(fullfile(fpath,'SortedIdx'));

spectra = spectra_raw;
spectra = pow2db(spectra);
spectra = trimmean(spectra,20,3);

% All states and spectral power 
figure; subplot(20,1,1); 
colors = get(gca,'colororder');
edges = [bins(:,1);bins(end,2)]./3600;
yl = [0,1];
for i = 1:max(smoothidx)
    hold on;
    ind = find(smoothidx==i);
    left = edges(ind)'; right = edges(ind+1)';
    bottom = yl(1)*ones(1,length(ind));
    top = yl(2)*ones(1,length(ind));
    patch([left;left;right;right],[bottom;top;top;bottom],colors(i,:),...
        'edgealpha',0,'facealpha',1);
end
ylim(yl); xlim([bins(1,1),bins(end,2)]/3600);
yticks([]); xticks([]);

keep = find(f>=35,1);
temp = subplot(20,1,2:10);
temppos = temp.Position;
imagesc(bins(:,1)/3600,f(1:keep),imgaussfilt(spectra(:,1:keep)',5)); box off;
ylim([0,35]);
colormap turbo
ylabel('Frequency (Hz)'); xlabel('Time (h)');
c = colorbar; ylabel(c, 'dB');
set(temp, 'Position', temppos);
set(gca,'FontSize',12);


% Average spectra
subplot(2,3,4); 
colors = get(gca,'colororder'); p = [];
for i = 1:max(smoothidx)
    x = f';
    y = mean(spectra(smoothidx==i,:))';
    dy = std(spectra(smoothidx==i,:))';
    h = fill([x;flipud(x)],[y-dy;flipud(y+dy)],colors(i,:),'edgealpha',0);
    set(h,'facealpha',0.15);
    hold on;
    p(i) = plot(x,mean(spectra(smoothidx==i,:)),'Color',colors(i,:),'linewidth',2);
end
legend(p,{'AM','AR','REM','NREM'},'box','off');
xlabel('Frequency (Hz)'); box off; 
ylabel('dB'); xlim([0,50]); ylim([-5,40])
set(gca,'FontSize',12);

% Example LFP traces during different states
fs = 500;
AM = nc3data(17,bins(5755,1),diff(bins(5755,:)),fs,[],fname);
AR = nc3data(17,bins(6477,1),diff(bins(6477,:)),fs,[],fname);
REM = nc3data(17,bins(1913,1),diff(bins(1913,:)),fs,[],fname);
NREM = nc3data(17,bins(1335,1),diff(bins(1335,:)),fs,[],fname);

subplot(2,3,5:6);
plot((1:length(AM))/fs,AM);
hold on; plot((1:length(AR))/fs,AR-400);
hold on; plot((1:length(REM))/fs,REM-800);
hold on; plot((1:length(NREM))/fs,NREM-1600);

yticks(fliplr([mean(AM),mean(AR)-400,mean(REM)-800,mean(NREM)-1600])); 
yticklabels({'NREM','REM','AR','AM'});
xlabel('Time (s)'); box off;
set(gca,'FontSize',12);

set(gcf,'renderer','painters');

%% Validation
%% EOG figures
% Figure of EOG trace with denoted eye movements
path = 'R:\Yun\Jafar\Neurochip';
day = 'Jafar_20210113_01';

[fpath,fname,Channels,fs,session_time] = getNCData(path,day);

fs = 500;
data = nc3data([31,32],22*60*60+47*60+36,6*60,fs,[],fname);

EOG = data(:,1)-data(:,2);
EOG = bpfilter(EOG,[1,20],fs,2);
en = envelope(abs(EOG),fs/2,'rms');
en = smooth(abs(EOG),fs*2);

figure; subplot(4,1,1); plot((1:length(EOG))/fs-5,EOG,'k','linewidth',1.5);
xlim([0,15]); xlabel('Time (s)'); ylabel('mV'); ylim([-300,500]); box off;
set(gca,'FontSize',12);

% EOG overnight (overlaid on classficiation example?)
load(fullfile(fpath,'EOG'));
eogbins = bins; 
load(fullfile(fpath,'SortedIdx'));

subplot(4,1,2);
colors = get(gca,'colororder'); 
eogbins = eogbins-bins(1,1);
plot(eogbins(:,1)./3600,smooth(EOGstd,4));
edges = [bins(:,1);bins(end,2)]./3600;
edges = edges-edges(1);
ylim([0,500]);
yl = ylim;
for i = 1:max(smoothidx)
    hold on;
    ind = find(smoothidx==i);
    left = edges(ind)'; right = edges(ind+1)';
    bottom = yl(1)*ones(1,length(ind));
    top = yl(2)*ones(1,length(ind));
    patch([left;left;right;right],[bottom;top;top;bottom],colors(i,:),...
        'edgealpha',0,'facealpha',0.75);
end
xlim([edges(1),edges(end)]); 

% subplot(40,1,12:20);
hold on; plot(eogbins(:,1)./3600,smooth(EOGstd),'k');
xlim([edges(1),edges(end)]); ylim([0,300]);

xlabel('Time(h)');

% xticklabels({'6pm','6am'}); 
ylabel('EOG stddev (mV)'); box off; set(gca,'FontSize',12);

%% Kinect figures
% Kinect movement vs accelerometer

% Load Kinect data
fname = 'R:\Yun\Jafar\OvernightTracking\20210322_Overnight.bin';
FID = fopen(fname,'r');

fields = {'double','single','single','single','single','int32','int32','single','single'};
fieldsizes = [8, 4, 4, 4, 4, 4, 4, 4, 4];
skip = @(n) sum(fieldsizes) - fieldsizes(n); 
offset = @(n) sum(fieldsizes(1:n)); 

data = [];
for i = 1:length(fieldsizes)
    data(:,i) = fread(FID, inf, fields{i}, skip(i));
    fseek(FID, offset(i), -1);
end

time = data(:,1);

posWeight = data(:,6);
negWeight = data(:,7);

weight = (posWeight+negWeight)/2;
smoothweight = smooth(weight);

% Load accelerometer data
filepath = fullfile('R:\Yun\Jafar\Neurochip');
exp = 'Jafar_20210323_01';

[fpath,fname,Channels,fs,session_time] = getNCData(filepath,exp);

load(fullfile(fpath,'SortedIdx'));
sortbins = bins;

X = nc3chan(0, ceil(session_time), 100, [], 100, fullfile(fpath,[exp,'_AccelX.i16']));
Y = nc3chan(0, ceil(session_time), 100, [], 100, fullfile(fpath,[exp,'_AccelY.i16']));
Z = nc3chan(0, ceil(session_time), 100, [], 100, fullfile(fpath,[exp,'_AccelZ.i16']));

dX = diff(X);
dY = diff(Y);
dZ = diff(Z);

mfs = 100;

movement = (dX.^2+dY.^2+dZ.^2);

% Get variance of movement
move = zeros(1,length(bins));

for i = 1:length(bins)
    left = round(bins(i,1)*mfs);
    right = round(bins(i,2)*mfs);
    if(left<=0)
        left = 1;
    end
    if(right > length(movement))
        right = length(movement);
    end
    move(i) = var(movement(left:right));
end

move = zscore(log(move));


bins = [sortbins(:,1);sortbins(end,2)];
% 1.0374 hour difference in recordings
x = (time-time(1))+1.0374*3600;
inds = discretize(x,bins);
binnedkinect = arrayfun(@(x) var(smoothweight(inds==x)), min(inds):max(inds));
binnedkinect = log(binnedkinect);

subplot(4,1,3); 
bins = bins - sortbins(1,1);
sortbins = sortbins-sortbins(1,1);
edges = [sortbins(:,1);sortbins(end,2)]./3600;
yl = [0,1];
for i = 1:max(smoothidx)
    hold on;
    ind = find(smoothidx==i);
    left = edges(ind)'; right = edges(ind+1)';
    bottom = yl(1)*ones(1,length(ind));
    top = yl(2)*ones(1,length(ind));
    patch([left;left;right;right],[bottom;top;top;bottom],colors(i,:),...
        'edgealpha',0,'facealpha',0.75);
end

tempmove = (move-min(move))./(max(move)-min(move));
L(1) = plot(bins(1:end-1)/3600, smooth(tempmove),'k');
tempkinect = (binnedkinect-min(binnedkinect))./(max(binnedkinect)-min(binnedkinect));
hold on; L(2) = plot(bins(1:length(binnedkinect))/3600, smooth(tempkinect),'r');
xlim([0,bins(end)/3600]);

xlabel('Time (h)'); ylabel('a.u.');
yticks([]); box off; set(gca,'FontSize',12);
legend(L,'Accelerometer','Kinect','box','off'); 

%% Cross validation figures
% Test and traininig error in terms of percent but also in terms of total
% time

if(exist(fullfile(fpath,'kfoldavg.mat')))
    load(fullfile(fpath,'kfoldavg.mat'))
else
    
    user = getenv('username');
    
    metafile = ['C:\Users\',user,'\Dropbox\Fetz Lab\_Brain States\Experiments.xlsx'];
    opts = detectImportOptions(metafile);
    metadata = readtable(metafile,opts);
    
    TestErr = {};
    TrainErr = {};
    TotalTime = {};
    
    for m = 1:size(metadata,1)
        % Path logistics
        animal = metadata.Animal{m};
        exp = metadata.Experiment{m};
        
        fprintf('%s - %s\n', animal, exp);
        
        filepath = fullfile('R:\Yun',animal,'Neurochip');
        
        files = dir(fullfile(filepath,exp,'*.mat'));
        filenames = extractfield(files,'name');
        
        [fpath,fname,Channels,fs,session_time] = getNCData(filepath,exp);
        
        load(fullfile(fpath,'k-fold.mat'))
        
        load(fullfile(fpath,'SortedIdx'),'bins');
        
        TestErr{m} = testerror;
        TrainErr{m} = trainerror;
        
        TotalTime{m} = bins(end,2)-bins(1,1);
        
    end
    
    avgTrainErr = cellfun(@(x) median(median(x)), TrainErr);
    avgTestErr = cellfun(@(x) median(median(x)), TestErr);
    
    avgTrainErrTime = cellfun(@(x,y) median(median(x.*y)), TrainErr, TotalTime);
    avgTestErrTime = cellfun(@(x,y) median(median(x.*y)), TestErr, TotalTime);
    
    save(fullfile(fpath,'kfoldavg.mat'),'avgTrainErr','avgTestErr','avgTrainErrTime','avgTestErrTime');
    
end

% Typically less than 4% test and train error, average of about 3.5% which
% translates to roughly 40 minutes of error. Fairly robust, most of the
% error is likely due to the majority filter but we need it to make sure
% it's closer to the "true" classification.
a = subplot(4,2,7); boxplot([avgTestErr,avgTrainErr].*100,[zeros(1,25),ones(1,25)],'notch','on','symbol','w');
xticklabels({'Test Error','TrainError'}); box off;
ylabel('Percent'); set(gca,'FontSize',12);
b = subplot(4,2,8); boxplot([avgTestErrTime,avgTrainErrTime]./60,[zeros(1,25),ones(1,25)],'notch','on');
xticklabels({'Test Error','TrainError'});
ylabel('Minutes'); box off;
set(gca,'FontSize',12);

a.Position(2:4) = b.Position(2:4);

set(gcf,'renderer','painters');



















