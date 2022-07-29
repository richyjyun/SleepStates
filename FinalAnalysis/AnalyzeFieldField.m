%% Get data from excel log
clear; close all;

metafile = 'Experiments.xlsx';
opts = detectImportOptions(metafile);
metadata = readtable(metafile,opts);

bw = [1,0.5,0.5,0.3,0.2,0.1];
MVL = {};
Counts = {};
BinAmp = {};
BinPhase = {};
LFPAmp = {};
MaxMVL = {};
packet = 'CrossFreq.ps';
ind = 1;
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
    
    crossfreqpath = fullfile(fpath,'CrossFreq');
    avgamppath = fullfile(fpath,'AvgAmp');
    
    for c = 1:length(spikechannels)
        
        fprintf('%d\n',spikechannels(c));
        
        if(any(c==skip)), continue; end
        
        spikefile = [num2str(spikechannels(c)),'.mat'];
       
        %,'MVL_PhasetoAmp','Counts_PhasetoAmp','AmpBins','PhaseBins','names','bands','bw'
        load(fullfile(crossfreqpath,spikefile));
        
        MLV{ind} = MVL_PhasetoAmp;
        Counts{ind} = Counts_PhasetoAmp;
        BinAmp{ind} = AmpBins;
        BinPhase{ind} = PhaseBins;
        
        % 'AvgAmp','MVL_Max','bands','names'
        load(fullfile(avgamppath,spikefile),'AvgAmp','MVL_Max');
        
        LFPAmp{ind} = AvgAmp;
        MaxMVL{ind} = MVL_Max;
        
        ind = ind+1;
    end
    
end
temppath = '\_Brain States';
save(fullfile(temppath,'MVL'),'MLV','LFPAmp','MaxMVL','Counts','BinAmp','BinPhase','names','bw','-v7.3');


%% Plot average mean vector per pairing
packetpath = 'Final';
packet = fullfile(packetpath,'AverageMeanVector.ps');
figure; fig = 1;
colors = get(gca,'colororder');
dir = []; mag = [];
for i = 1:5
    for j = (i+1):6
        subplot(5,5,(i-1)*5+j-1)
        for s = 1:4
            vals = cellfun(@(x) x(s,i,j), MLV);
%             norm = cellfun(@(x) x(s,j), LFPAmp); 
            norm = cellfun(@(x) abs(x(s,i,j)), MaxMVL);
%             vals(bad) = [];
%             temp = mean(vals);
%             direction = angle(temp); magnitude = abs(temp);
            direction = circMean(angle(vals)); dir(i,j,s) = direction;
            magnitude = median(abs(vals)./norm); mag(i,j,s) = magnitude;
            polarplot([0,direction],[0,magnitude],'--','color',colors(s,:)); 
            hold on;
            polarscatter(direction,magnitude,50,colors(s,:),'x','linewidth',2);
            set(gca,'Thetaaxisunits','radians');
            set(gca,'Gridlinestyle','--');
            thetaticks(0:pi/4:2*pi);
            set(gca,'FontSize',10);
        end
        title([names{i},' phase to ',names{j},' amp']);
        fig = fig+1;
    end
end
set(gcf,'renderer','painters'); orient(gcf,'landscape')
print(packet,'-append','-dpsc2','-fillpage');
close(gcf);
callps2pdf(packet);

%% Boxplots of MVL
figure; fig = 1;
labels = {[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]};
yl = nan(5,5);
yl(2,2) = .015;
yl(2,3) = .1;
yl(2,4) = .12;
yl(3,3) = .052;
yl(3,4) = .08;
colors = get(gca,'colororder');
pvals = [];
for i = 1:5
    for j = (i+1):6
        subplot(5,5,(i-1)*5+j-1)
        mags = [];
        for s = 1:4
            vals = cellfun(@(x) x(s,i,j), MLV);
            %             norm = cellfun(@(x) x(s,j), LFPAmp);
            norm = cellfun(@(x) abs(x(s,i,j)), MaxMVL);
            mags(:,s) = (abs(vals)./norm);

        end
        plotmags = mags;
        if(~isnan(yl(i,j-1)))
            plotmags(mags>yl(i,j-1)) = nan;
        end
        %         boxplot(log(mags),'notch','on'); hold on;
        violin = violinplot(plotmags);
        for s = 1:4
            violin(s).ScatterPlot.MarkerFaceColor = colors(s,:);
            violin(s).ViolinPlot.FaceColor = colors(s,:);
            violin(s).BoxPlot.Vertices(:,1) = [s-.1;s+.1;s+.1;s-.1];
            violin(s).BoxPlot.FaceColor = [0,0,0];
            violin(s).ScatterPlot.SizeData = 5;
            violin(s).MedianPlot.SizeData = 20;
        end
        
        xticks([]);
        set(gca,'FontSize',10);
        
        
        [~,~,stats] = kruskalwallis(mags,[],'off');
        [c,~,~,~] = multcompare(stats,[],'off');
        
        pvals{i,j} = c(:,[1,2,6]);
        sig = c(:,6) < 0.05;
        sigstar(labels(sig),c(sig,6),0,10,0);
        
    end
end

%% Plot distribution of mean vector lengths
packetpath = 'C:\Users\Richy Yun\Dropbox\Fetz Lab\_Brain States\Figures\Final';
packet = fullfile(packetpath,'MVL Distribution.ps');
figure; fig = 1;
colors = get(gca,'colororder');
yl = [-.5,3; -.5,1.5; -.5,1.5; -.5,1.5; -.5,1.5; -0.02,0.15; ...
    -.05,0.4; -.05,0.4; -.05,0.5; -.05,0.3; -.05,0.3;...
    -.05,0.3; -.05,0.3; -.05,0.3; -.05,0.3];
for i = 1:5
    for j = (i+1):6
        subplot(5,5,(i-1)*5+j-1)
        vals = [];
        for s = 1:4
            vals(s,:) = cellfun(@(x) abs(x(s,i,j)), MLV);
        end
        boxplot(vals','notch','on','symbol','w');
        title([names{i},' phase to ',names{j},' amp']);
        ylim(yl(fig,:));
        fig = fig+1;
    end
end
set(gcf,'renderer','painters'); orient(gcf,'landscape')
print(packet,'-append','-dpsc2','-fillpage');
close(gcf);
callps2pdf(packet);

%% Distribution of MVL 
% Phase and magnitude
figure;
colors = get(gca,'colororder');
widths = [0.25,0.25,0.2,0.2,0.1; 0,0.01,0.025,0.025,0.05; 0,0,0.025,0.025,0.025; 0,0,0,0.025,0.025; 0,0,0,0,0.025];
lims = [5,5,4,2,2; 0,0.3,1,0.5,0.6; 0,0,0.4,0.4,0.6; 0,0,0,0.5,0.5; 0,0,0,0,0.3];
for i = 1:5
    for j = (i+1):6
        subplot(5,5,(i-1)*5+j-1)
        for s = 1:4
            mvl = cellfun(@(x) x(s,i,j), MLV);
%             mvl(bad) = [];
            bins = 0:widths(i,j-1):max(abs(mvl));
            counts = histcounts(abs(mvl),bins);
            counts = [counts,0];
            stairs(bins,counts,'linewidth',1.5);
            hold on; 
        end
        xlim([0,lims(i,j-1)]);
        title([names{i},' phase to ',names{j},' amp']);
        box off; xlabel('MVL'); ylabel('Spikes');
    end
end

figure;
colors = get(gca,'colororder');
for i = 1:5
    for j = (i+1):6
        subplot(5,5,(i-1)*5+j-1)
        for s = 1:4
            mvl = cellfun(@(x) x(s,i,j), MLV);
            polarhistogram(angle(mvl),-pi:pi/8:pi,...
            'displaystyle','stairs','linewidth',1.5)
            hold on; 
        end
        title([names{i},' phase to ',names{j},' amp']);
        set(gca,'Thetaaxisunits','radians');
        set(gca,'Gridlinestyle','--');
        thetaticks(0:pi/4:2*pi);
    end
end


%% Sanity Checks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make a packet for each spike
packetpath = 'C:\Users\Richy Yun\Dropbox\Fetz Lab\_Brain States\Packets\FieldField';
for c = 1:length(Counts)
    disp(c)
    packet = fullfile(packetpath,[num2str(c),'.ps']);
    count = Counts{c};
    for b1 = 1:5
        for b2 = b1+1:6
            figure('visible','off');
            for s = 1:4
                allcounts = count{s,b1,b2};
                allcounts = allcounts./sum(allcounts);
                allcounts(end+1,:) = allcounts(end,:);
                
                subplot(2,2,s)
                polarPcolor(0:bw(b2):(bw(b2)*(size(allcounts,2)-1)),rad2deg(BinPhase{b1}),imgaussfilt(allcounts,2)')
                if(s==1)
                    title([names{b1},' phase to ',names{b2},' amp']);
                end
            end
            %         set(gcf,'renderer','painters');
            orient(gcf,'landscape');
            print(packet,'-append','-dpsc2','-fillpage');
            close(gcf);
        end
    end
    callps2pdf(packet);
end

%% Check to make sure each phase has similar count
% Seems to be generally fine
for c = 1:length(Counts)
    temp = cellfun(@(x) sum(x,2), Counts{c},'uniformoutput',false);
end

%% Sanity check of recalculating MVL
% Seems to be generally fine
temp = Counts{1};
temp = temp{4,1,2};
amp = 0:bw(2):((size(temp,2)-1)*bw(2));
amp = amp+bw(2)/2;
phase = linspace(-pi,pi,360);

mv = 0;
for i = 1:size(temp,1)
    for j = 1:size(temp,2)
        vector = temp(i,j).*amp(j).*exp(1i.*phase(i));
        mv = mv+vector;
    end
end
mv = mv./sum(sum(temp));

