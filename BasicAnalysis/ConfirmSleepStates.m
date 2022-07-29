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
    
    fprintf('%s - %s\n', animal, exp);
    
    filepath = fullfile('R:\Yun',animal,'Neurochip');
    
    files = dir(fullfile(filepath,exp,'*.mat'));
    filenames = extractfield(files,'name');    
    
    %% Load data
    fprintf('Loading data...'); tic;
    [fpath,fname,Channels,fs,session_time] = getNCData(filepath,exp);
    
    % Load sorted states
    % bins, smoothidx
    load(fullfile(fpath,'SortedIdx'));
    
    % Load low dimensional data
    % Vr, r
    load(fullfile(fpath,'LowDim'));
    
    t = toc;
    fprintf('%f s\n',t);
    
    %%  Cross-validation
    fprintf('Cross Validation...Rep.'); tic
    K = 20; repeat = 50;
    trainerror = zeros(K,repeat);
    testerror = zeros(K,repeat);
    for rep = 1:repeat
        
        fprintf('%d.',rep);
        
        % Find limits of each group
        inds = randperm(size(Vr,1));
        edges = linspace(0,1,K+1); edges = round(edges.*length(inds));
        
        for k = 1:K
            % Get test and train groups
            testinds = inds(edges(k)+1:edges(k+1));
            traininds = 1:length(inds);
            traininds(inds(edges(k)+1:edges(k+1))) = [];
            
            %% Train on train group
            trainVr = Vr(traininds,:);
            
            % Kmeans with 4 clusters
            % Remove outliers using just the SVD dimensions(5% of furthest points)
            dist = squareform(pdist(trainVr(:,1:end-1)));
            dist = mean(dist);
            prc = prctile(dist,90);
            tempVr = trainVr(dist<prc,:);
            
            % Cluster
            [~, centroids] = kmeans(tempVr,4,'Replicates',50,'MaxIter',1000);
            
            % Find which cluster each point is in
            dist = pdist2(centroids,Vr);
            [~,trainidx] = min(dist);
    
            %% Compare with original sorting to figure out which cluster is which
            newidx = zeros(1,length(trainidx));
            for s = 1:4
                state = mode(smoothidx(trainidx==s));
                newidx(trainidx==s) = state;
            end
            trainidx = newidx;
            
            %% Majority filter
            window = 4; % Over ~30 seconds
            smoothtrainidx = majorityFilt(trainidx, window);
            
            %% Calculate errors
            temp = sum(smoothtrainidx(traininds) ~= smoothidx(traininds));
            trainerror(k,rep) = temp ./ length(traininds);
            
            temp = sum(smoothtrainidx(testinds) ~= smoothidx(testinds));
            testerror(k,rep) = temp ./ length(testinds);
            
        end              
    end
    
    save(fullfile(fpath,'k-fold.mat'),'trainerror','testerror','k','repeat');
    
    t = toc;
    fprintf('%f s',t);
    
end






