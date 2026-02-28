function [spikecountAve, spikecountraw, SCvar] = ...
    MovieSpikecount(spiketrain, goodcell, imgPara, stimInd, StimTimeIndv,...
    SpikeCount_Twindow, IgnoreTimeWindow)
StimNum = imgPara.stim_time/StimTimeIndv;
SpikeCountbin = SpikeCount_Twindow/0.001;
spikecountAve = []; spikecountraw = []; SCvar = [];
for n = 1:length(goodcell)
    goodn = goodcell(n); 
    tmphist = []; 
    for j1 = 1:length(stimInd)
        j2 = stimInd(j1); % sort respose by temporal freqeuncy 
        for i = 1:imgPara.stimrep
            st = spiketrain(goodn).st{i,j2}*imgPara.dt;
            st(st<0) = []; st(st>imgPara.stim_time) = [];
            stbin = histcounts(st, 0:0.001:imgPara.stim_time); % bin spikecount at 1ms timebin
            stbin = reshape(stbin, StimTimeIndv/0.001, StimNum);
            if IgnoreTimeWindow>0
                ignoreTimebin = IgnoreTimeWindow/0.001;
                stbin(1:ignoreTimebin,:) = [];
            end 
            sctmp = sum(stbin(1:SpikeCountbin,:),1);
%             tmphist1 = histcounts(st, SpikeCount_Twindow);
%             tmphist2 = tmphist1(1:StimTimeIndv/spikecountkeep_Twindow:end);
            tmphist(i,:,j1) = sctmp;
        end
    end
    spikecountraw(n,:,:,:) = tmphist;
    SCvar(n,:,:) = squeeze(nanstd(tmphist,[],1));
    spikecountAve(n,:,:) = squeeze(nanmean(tmphist,1));
end
