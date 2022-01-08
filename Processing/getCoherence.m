function getCoherence(fpath,fname,Channels,session_time,dwnfs,stim)
%% Saves the average of pairwise coherence for each bin
% Bins are 8 seconds wide. If there was stimulation (0.1Hz), ignore
% 1.9 seconds after stimulation and 0.1 seconds before to ensure artifact
% and stimulus response are not included

%% Load bins saved from saveSpectra
% Need to initialize for parfor
bins = [];
load(fullfile(fpath,'Bins.mat'))

% Frequencies of interest (Higher frequencies show coherence)
f = 0:0.2:120;
    
%% Get all coherence 
coherence_raw = [];
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
        Coh = [];
        for a = 1:(size(data,2)-1)
            [temp,~] = mscohere(data(:,a+1:end),data(:,a),[],[],f,dwnfs);
            if(a == (size(data,2)-1))
                temp = temp';
            end
            Coh = [Coh,temp];
        end
    catch
        continue;
    end
    
    coherence_raw(:,i) = mean(Coh,2);
end

save(fullfile(fpath,'Coherence.mat'),'coherence_raw','dwnfs',...
    'bins','f','-v7.3');

end

