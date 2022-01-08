function getEOG(fpath,fname,Channels,session_time,dwnfs,stim)

%% Load bins from saveSpectra
% Need to initialize for parfor
bins = [];
load(fullfile(fpath,'Bins.mat'))

%% Get all EOG data 
EOGmean = zeros(size(bins,1),1); 
EOGmed = zeros(size(bins,1),1); 
EOGvar = zeros(size(bins,1),1); 
EOGstd = zeros(size(bins,1),1); 
EOGint = zeros(size(bins,1),1); 
parfor i = 1:size(bins,1)
    t = bins(i,1);
    width = bins(i,2)-bins(i,1);
    if(width < 5)
        continue;
    end
    width = floor(width*(dwnfs/10))/(dwnfs/10); % fixing so it fits sampling rate    
    t = ceil(t*dwnfs)/dwnfs;
    
    try
        [data, ~, ~] = nc3data(Channels, t, width, dwnfs, [1,60], fname);
    catch
        continue;
    end
    
    if~(size(data,2) < 2)
        EOG = data(:,1) - data(:,2);
        EOG = bpfilt(EOG,[1,20],dwnfs,2);
        
        EOGmean(i) = mean(EOG);
        EOGmed(i) = median(EOG);
        EOGvar(i) = var(EOG);
        EOGstd(i) = std(EOG);
        EOGint(i) = sum(abs(EOG));
    end
    
end

save(fullfile(fpath,'EOG.mat'),'EOGmean','EOGmed','EOGvar','EOGstd','EOGint','dwnfs',...
    'bins','-v7.3');

end