function getSpectra(fpath,fname,Channels,session_time,dwnfs,stim)
%% Gets spectra for each bin
% Bins are 8 seconds wide. If there was stimulation (0.1Hz), ignore
% 1.9 seconds after stimulation and 0.1 seconds before to ensure artifact
% and stimulus response are not included

%% Load meta data
events = nc3events(fname);

%% Find all bins to get spectra from
start = 0; 
if stim
    stimtimes = cellfun(@(x) x',events.stim,'uniformoutput',0);
    stimtimes = sort(cell2mat(stimtimes));
    stimtimes(stimtimes>session_time) = [];
    left = [0,stimtimes+1.9];
    right = [stimtimes-0.1,session_time];
else
    width = 8;
    left = start:width:session_time;
    right = [left(2:end),session_time];
end

bins = [left',right'];

% Remove bins that are too small
bins(diff(bins')<5,:) = [];

% Frequencies of interest (higher frequencies are not prominent enough)
f = 0:0.1:50;

%% Get all spectra
spectra_raw = [];
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
        [pxx,~] = pwelch(data,[],[],f,dwnfs);
    catch
        continue;
    end
    
    spectra_raw(i,:,:) = pxx;
end

save(fullfile(fpath,'Bins.mat'),'bins');
save(fullfile(fpath,'Spectra.mat'),'spectra_raw','dwnfs',...
    'bins','f','-v7.3');

end

