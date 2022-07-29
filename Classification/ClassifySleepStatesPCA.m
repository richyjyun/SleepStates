%% Classifies sleep states for each experiment using PCA and saves them
% Additionally plots the sorting into a pdf for debugging

%% Get data from excel log
clear; close all;

user = getenv('username');

metafile = 'Experiments.xlsx';
opts = detectImportOptions(metafile);
metadata = readtable(metafile,opts);

packet = 'ClassifiedStates.ps';

%% For each experiment
for m = 1:size(metadata,1)
    
    %% Path logistics
    animal = metadata.Animal{m};
    exp = metadata.Experiment{m};
    
    fprintf('%s - %s\n', animal, exp);
    
    filepath = fullfile('R:\Yun',animal,'Neurochip');
    
    %% Load meta data
    [fpath,fname,Channels,fs,session_time] = getNCData(filepath,exp);
    
    %% Load Spectra
    spectrafile = fullfile(fpath,'Spectra.mat');
    if(exist(spectrafile))
        fprintf('Loading spectra\n');
        load(spectrafile); 
    else
        fprintf('No spectra saved\n');
        continue;
    end
    
    % Some days had sedation, need to remove first x hours
    if(~isnan(metadata.Start(m)))
        temp = find(bins(:,1) > metadata.Start(m)*3600,1);
        bins = bins(temp:end,:);
        spectra_raw = spectra_raw(temp:end,:,:);
    end
    
    % If spectra_raw is larger, fix bins (last bin may not be counted if
    % it's too short)
    if size(bins,1) > size(spectra_raw,1)
        bins = bins(1:size(spectra_raw,1),:);
    end
    
    %% Get accelerometer data
    fprintf('Loading accelerometer\n');
    X = nc3chan(0, ceil(session_time), 100, [], 100, fullfile(fpath,[exp,'_AccelX.i16']));
    Y = nc3chan(0, ceil(session_time), 100, [], 100, fullfile(fpath,[exp,'_AccelY.i16']));
    Z = nc3chan(0, ceil(session_time), 100, [], 100, fullfile(fpath,[exp,'_AccelZ.i16']));
    
    dX = diff(X);
    dY = diff(Y);
    dZ = diff(Z);
    
    mfs = 100;
    
    movement = sqrt(dX.^2+dY.^2+dZ.^2);
    
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
    
    move = move(1:size(spectra_raw,1));
    bad = find(isnan(move));
    move(bad) = [];
    spectra_raw(bad,:,:) = [];
    bins(bad,:) = [];
    
    % Normalize
    move = move'; move = zscore(log(move)); move = move./std(move);
    
    %% Normalize spectra
    fprintf('Classifying states\n');
    
    spectra = spectra_raw;
    
    spectra = log(spectra);
    
    spectra = trimmean(spectra,20,3);
    
    spectra = spectra - min(spectra,[],2);
    
    spectra = spectra./sum(spectra,2);
        
    %% PCA
    [coeff,score,~,~,explained,~] = pca(spectra);
    variance = cumsum(explained/sum(explained));
    r = find(variance>0.9,1); % Use number of dimensions that explains more than 90% of variance

    Vr = score(:,1:r);
    Ur = coeff(:,1:r);
    
    % Adjust movement so it affects the sorting but doesn't dominate
    move = move.*(mean(std(Vr(:,1))));
    Vr(:,r+1) = move(1:size(Vr,1));
    
    %% Kmeans with 4 clusters
    % Remove outliers using just the SVD dimensions(5% of furthest points) 
    dist = squareform(pdist(Vr(:,1:r)));
    dist = mean(dist);
    prc = prctile(dist,90);
    tempVr = Vr(dist<prc,:);
    
    % Cluster
    [~, centroids] = kmeans(tempVr,4,'Replicates',50,'MaxIter',1000);
    
    % Find which cluster each point is in
    dist = pdist2(centroids,Vr);
    [~,idx] = min(dist);
    
    %% Determine which cluster is which
    % First find cluster with highest acceleration and denote it as awake
    % and moving. Then find the cluster with highest delta power and denote
    % it to be NREM sleep. Then find the next highest movement cluster and
    % denote it to be awake and at rest. The final remaining cluster is REM
    % sleep.
    % 1 - Awake and moving
    % 2 - Awake and at rest
    % 3 - REM sleep
    % 4 = NREM sleep
    avgspect = [];
    avgmove = [];
    for i = 1:max(idx)
        avgspect(i,:) = mean(spectra(idx==i,:));
        avgmove(i) = mean(move(idx==i));
    end
    
    fixind = zeros(1,4);
    [~,maxmove] = max(avgmove);
    fixind(maxmove) = 5;
    remind = find(fixind==0);
    
    delta = sum(avgspect(remind,1:find(f>5,1))');
    [~,maxdelta] = max(delta);
    fixind(remind(maxdelta)) = 8;
    remind = find(fixind==0);
    
    [~,maxmove] = max(avgmove(remind));
    fixind(remind(maxmove)) = 6;
    remind = find(fixind==0);
    
    fixind(fixind==0) = 7;
    
    idx(idx==1) = fixind(1);
    idx(idx==2) = fixind(2);
    idx(idx==3) = fixind(3);
    idx(idx==4) = fixind(4);
    idx(idx==5) = 1;
    idx(idx==6) = 2;
    idx(idx==7) = 3;
    idx(idx==8) = 4;
    
    % Flipped for some reason for this specific day
    if strcmp(exp,'Jafar_20200217_05')
        temp1 = idx==2;
        temp2 = idx==3;
        idx(temp1) = 3;
        idx(temp2) = 2;
    end
    
    %% Majority filter
    window = 4; % Over ~30 seconds
    smoothidx = majorityFilt(idx, window);
    
    %% Save states
    save(fullfile(fpath,'SortedIdx'),'idx','smoothidx','bins','r');
    
    %% Save low dimensional data for cross validation
    save(fullfile(fpath,'LowDim'), 'Vr', 'centroids', 'r');  
    
    %% Plot and save for debugging
    % Plot the eigenvectors
    figure('visible','off');
    subplot(2,2,1);
    plot(f,Ur(:,1:5),'linewidth',2);
    xlabel('Hz'); ylabel('a.u.'); xlim([f(1),f(end)]);
    legend('EV1','EV2','EV3','EV4');
    title(exp,'Interpreter', 'none');
    
    % Plot the data points in PCA / accelerometer space
    subplot(2,2,2);
    hold off;
    for i = 1:max(smoothidx)
        scatter3(Vr(smoothidx==i,1),Vr(smoothidx==i,2),Vr(smoothidx==i,r+1),'.'); hold on;
    end
    xlabel('PC1'); ylabel('PC2'); zlabel('Accel');
    xticks([]); yticks([]); zticks([]);
    view(-30,15)
    
    % Spectra stacked on top
    subplot(2,2,3);
    colors = get(gca,'colororder'); p = [];
    origspect = trimmean(spectra_raw,20,3); origspect = origspect - min(spectra,[],2);
    origspect = origspect./sum(origspect,2);
    for i = 1:max(smoothidx)
        x = f';
        y = mean(spectra(smoothidx==i,:))';
        dy = std(spectra(smoothidx==i,:))';
        h = fill([x;flipud(x)],[y-dy;flipud(y+dy)],colors(i,:),'edgealpha',0);
        set(h,'facealpha',0.15);
        hold on;
        p(i) = plot(x,mean(spectra(smoothidx==i,:)),'Color',colors(i,:),'linewidth',2);
    end
    xlabel('Frequency (Hz)'); xlim([0,50])
    ylabel('Normalized Spectral Density');
    legend(p,'AM','AR','REM','NREM','Location','northeast');
    
    % States over time
    subplot(2,2,4); ts = mean(bins,2);
    hold on;
    colors = get(gca,'colororder');
    for i = 1:max(smoothidx)
        s(i) = scatter(ts(smoothidx==i)/3600,i*ones(1,sum(smoothidx==i)),20,colors(i,:),'s','filled');
        s(i).MarkerFaceAlpha = 0.05;
        hold on;
    end
    xlabel('Time(h)');
    set(gca, 'YDir','reverse');
    yticks(1:4);
    ylim([0,max(smoothidx)+1]); ylabel('State')
    box off;
    
    % Print to ps file
    print(packet,'-append','-dpsc2','-fillpage');
    close(gcf);
    
end

callps2pdf(packet);
