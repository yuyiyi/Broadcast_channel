% correlation vs. distance
load('pooldata_SinGrat.mat')
Dir = [0, 135, 270, 45, 180, 315, 90, 225];
[DirSort, DirInd] = sort(Dir);
stimPara = ...
[8	0.16; 
2	0.02;
8   0.02; 
1	0.02; 
1	0.16;
2	0.16;
8	0.04;
2	0.04;
1	0.04];
[stimsort, stimInd] = sortrows(stimPara, [1 2]);
TmFreq = stimsort(:,1);
SpFreq = stimsort(:,2);

%%
SpikeCount_Twindow = 1; % compute spike count in 1 sec time window following IgnoreTimeWindow 
StimTimeIndv = 2; % duration of each stimulus is 2 sec
IgnoreTimeWindow = 0; % ignore first 0 sec response
clear spikecount
for k1 = 1:size(pooldata,1)
    for k2 = 1:4
        if isempty(pooldata(k1,k2).spiketrain)
            continue
        end
        spiketrain = pooldata(k1,k2).spiketrain;
        goodcell = pooldata(k1,k2).goodcell; 
        imgPara = pooldata(k1,k2).imgPara;
        num(k1,k2) = length(goodcell);
        roi = pooldata(k1,k2).Region;
        loc = regionprops(roi,'Centroid');
        [spikecountAve, spikecountraw,SCvar] = SinGratSpikecount(spiketrain,...
            goodcell, imgPara, stimInd, Dir, StimTimeIndv, SpikeCount_Twindow, IgnoreTimeWindow);
        spikecount{k1,k2}.spikecountAve = spikecountAve;
        spikecount{k1,k2}.variance = SCvar;
        spikecount{k1,k2}.spikecountrawSC = spikecountraw;
        spikecount{k1,k2}.neuloc = loc(goodcell);
        spikecount{k1,k2}.ROIloc = ROIloc(k1,(k2-1)*2+1:k2*2)*250;
    end
end

%% noise correlation within ROI
clear nscinROI
d1 = 500/629; d2 = 500/512;
for k1 = 1:size(spikecount,1)
    for k2 = 1:4
        if isempty(spikecount{k1,k2})
            continue
        end
        loc = spikecount{k1,k2}.neuloc;
        scA = spikecount{k1,k2}.spikecountrawSC;
%         scA = zscore(scA,[],2);
        if isempty(scA) || size(scA,1)<3
            continue
        end
        % noise correlation by methods3 (cor(pool spike count from all stimulus))
        sc1 = reshape(scA, size(scA,1), size(scA,2), 72);
        pearsontmp = NeuCor(sc1, sc1, 1);
        pearsonshulf = NeuCor_shulf(sc1, sc1, 1);
        A = ones(size(pearsontmp));
        A = triu(A,1); 
        nscinROI(k1,k2).nscall_SC = pearsontmp(A==1);
        nscinROI(k1,k2).nscall_SCshulf = pearsonshulf(A==1);
        % signal correlation
        FR = spikecount{k1,k2}.spikecountAve;
        FR = reshape(FR, size(FR,1), 72);
        r_sig = corr(FR');
        nscinROI(k1,k2).r_sig = r_sig(A==1); 
        % geometric mean FR
        [meanFR, maxFR] = GeometricMeanFR(FR, FR);
        meanFR1 = [];
        for j = 1:size(meanFR,3)
            tmp = meanFR(:,:,j);
            meanFR1(:,j) = tmp(A==1);
        end
        nscinROI(k1,k2).maxFR = maxFR(A==1);
        nscinROI(k1,k2).meanFR = meanFR1;        
        d = NeuDistance(loc, loc, d1, d2, [0 0], [0 0]);
        nscinROI(k1,k2).d = d(A==1);
    end
end

%% noise correlation between ROI
% clear nscBTWROI 
areapair = [];
k = 0;
for k1 = 1:size(spikecount,1)
    for k2 = 1:3
        for k3 = k2+1:4
            if isempty(spikecount{k1,k2}) || isempty(spikecount{k1,k3})
                continue
            end
            k = k+1;
            scA = spikecount{k1,k2}.spikecountrawSC;
            scB = spikecount{k1,k3}.spikecountrawSC;
            areapair(k,1) = areaID(k1,k2);
            areapair(k,2) = areaID(k1,k3);
            areapair(k,3) = k1;
            areapair(k,4) = k2;
            areapair(k,5) = k3;
            if size(scA,1)<3 || size(scB,1)<3
                continue
            end
            FRA = spikecount{k1,k2}.spikecountAve;
            FRB = spikecount{k1,k3}.spikecountAve;
            FRA = reshape(FRA, size(FRA,1), 72);
            FRB = reshape(FRB, size(FRB,1), 72);
            [meanFR, maxFR] = GeometricMeanFR(FRA, FRB);
            sc1 = reshape(scA, size(scA,1), size(scA,2), 72);
            sc2 = reshape(scB, size(scB,1), size(scA,2), 72);
            pearsontmp = NeuCor(sc1,sc2, 1);
            nscBTWROI(k).nscall_SC = pearsontmp(:);  
            pearsonshulf = NeuCor_shulf(sc1, sc2, 1);
            nscBTWROI(k).nscall_SCshulf = pearsonshulf(:);
            nscBTWROI(k).meanFR = reshape(meanFR, [], 72);
            nscBTWROI(k).maxFR = maxFR(:);
            r_sig = corr(FRA', FRB');
            nscBTWROI(k).r_sig = r_sig(:);
            ROIlocA = spikecount{k1,k2}.ROIloc;
            ROIlocB = spikecount{k1,k3}.ROIloc;
            locA = spikecount{k1,k2}.neuloc;
            locB = spikecount{k1,k3}.neuloc;
            d = NeuDistance(locA, locB, d1, d2, ROIlocA, ROIlocB);
            nscBTWROI(k).d = d(:);
        end
    end
end

%% noise correlation to preferred stimulus
for k1 = 1:size(spikecount,1)
    for k2 = 1:4
        if isempty(spikecount{k1,k2})
            continue
        end
        scA = spikecount{k1,k2}.spikecountrawSC;
%         scA = zscore(scA,[],2);
        if isempty(scA) || size(scA,1)<3
            continue
        end
        meanFR = nscinROI(k1,k2).meanFR;
        [~,stimIDtmp] = max(meanFR,[],2);
        [ind2, ind1] = meshgrid(1:size(scA,1),1:size(scA,1));
        A = ones(size(scA,1), size(scA,1)); A = triu(A,1); 
        ind1 = ind1(A==1); ind2 = ind2(A==1);
        stimID = [ind1,ind2,stimIDtmp];
        % noise correlation to preferred stimulus (max geometric mean FR)
        sc1 = reshape(scA, size(scA,1), size(scA,2), 72);
        pearsontmp = NeuCor(sc1, sc1, 3, stimID);        
        nscinROI(k1,k2).nscpref_SC = pearsontmp';
    end
end

k = 0;
for k1 = 1:size(spikecount,1)
    for k2 = 1:3
        for k3 = k2+1:4
            if isempty(spikecount{k1,k2}) || isempty(spikecount{k1,k3})
                continue
            end
            k = k+1;
            scA = spikecount{k1,k2}.spikecountrawSC;
            scB = spikecount{k1,k3}.spikecountrawSC;
            if size(scA,1)<3 || size(scB,1)<3
                continue
            end
            meanFR = nscBTWROI(k).meanFR;
            [~,stimIDtmp] = max(meanFR,[],2);
            [ind2, ind1] = meshgrid(1:size(scB,1),1:size(scA,1));
            ind1 = ind1(:); ind2 = ind2(:);
            stimID = [ind1,ind2,stimIDtmp];
            sc1 = reshape(scA, size(scA,1), size(scA,2), 72);
            sc2 = reshape(scB, size(scB,1), size(scB,2), 72);
            pearsontmp = NeuCor(sc1,sc2,3,stimID);
            nscBTWROI(k).nscpref_SC = pearsontmp';
        end
    end
end

%% number of NC pairs
load('NSC_dist_SinGrat.mat', 'nscinROI', 'nscBTWROI', 'areaID', 'areapair')
NC_num = [];
for k = 1:50
    n = 0;
    for k2 = 1:4
        n = n + sum(~isnan(nscinROI(k,k2).nscall_SC));
    end
    id = find(areapair(:,3)==k);
    for i = 1:length(id)
        n = n + sum(~isnan(nscBTWROI(id(i)).nscall_SC));
    end
    NC_num(k, 1) = n;
end
areaID([5,35,43],:) = [];
nscinROI([5,35,43],:) = [];
NC_area = [];
for k = 1:6
    id = find(areaID==k-1);
    n = 0;
    for i = 1:length(id)
        n = n + sum(~isnan(nscinROI(id(i)).nscall_SC));
    end
    NC_area(k) = n;
end
%% NSC vs. neuron distance
clear NSCin_pool mean_* se_*
titlelabel = {'V1', 'LM', 'AL', 'PM','LI'};
frthreshold = 0.5; nsc_RsigInROI = []; frackeep_in = []; numinpool = []; 
nsc_distInROI = [];
for k = 1:5
    [id1] = find(areaID==k-1);
    if ~isempty(id1)
        % correlation to preferred stimulus
        pearsontmp = []; r_sigtmp = [];FRtmp = []; 
        dtmp = []; pearsonpref = [];FRmax = []; 
        for k1 = 1:length(id1)
            if ~isempty(nscinROI(id1(k1)).nscall_SC)
                dtmp = cat(1, dtmp, nscinROI(id1(k1)).d);
                r_sigtmp = cat(1,r_sigtmp,nscinROI(id1(k1)).r_sig);
                FRmeantmp = nanmean(nscinROI(id1(k1)).meanFR,2);
                tmp = nscinROI(id1(k1)).nscall_SC;
                tmppref = nscinROI(id1(k1)).nscpref_SC;
                FRmaxtmp = nscinROI(id1(k1)).maxFR;
                FRtmp = cat(1,FRtmp, FRmeantmp);
                FRmax = cat(1,FRmax, FRmaxtmp);
                pearsontmp = cat(1,pearsontmp,tmp);
                pearsonpref = cat(1,pearsonpref,tmppref);
            end
        end
        r_sigtmp(FRmax<frthreshold) = [];
        pearsontmp(FRmax<frthreshold,:) = [];
        dtmp(FRmax<frthreshold) = [];        
        frackeep_in(k) = sum(FRmax<frthreshold)/length(FRmax);
        FRtmp(FRmax<frthreshold) = [];
        FRmax(FRmax<frthreshold) = [];
        mean_pearson(k) = nanmean(nanmean(pearsontmp,2));
        se_pearson(k) = nanstd(nanmean(pearsontmp,2))/sqrt(size(pearsontmp,1));
        numinpool(k,1) = size(pearsontmp,1);
        numinpool(k,2) = length(id1);
        NSCin_pool{k} = [pearsontmp, r_sigtmp, dtmp, FRtmp];
        
        % noise correlation vs. distance 
        nsctmp = pearsontmp; 
        x = 5:80:500; nsc_d = [];
        for j = 1:length(x)-1
            id = intersect(find(dtmp>x(j)), find(dtmp<=x(j+1)));
            if length(id)>50
                nsc_d(j,1) = nanmean(nsctmp(id));
                nsc_d(j,2) = nanstd(nsctmp(id))/sqrt(length(id));
            else
                nsc_d(j,:) = [nan, nan];
            end
        end
        nsc_distInROI(:,k) = nsc_d(:,1);
        figure(5), 
        subplot(3,5,k), shadedErrorBar(x(1:end-1)+40,nsc_d(:,1),nsc_d(:,2))
        xlabel('Neuron distance'); ylabel('Noise correlation')
        title(titlelabel{k}); axis([0 500 0 0.05]);  box off
        figure(6);
        subplot(121), plot(x(1:end-1)+10,nsc_d(:,1)); hold on
        xlabel('Neuron distance'); ylabel('Noise correlation')
        % noise correlation vs. signal correlation
        nsctmp = pearsontmp;
        nsc_d = []; x = -1:0.1:1;
        for j = 1:length(x)-1
            id = intersect(find(r_sigtmp>x(j)), find(r_sigtmp<=x(j+1)));
            if length(id)>20
                nsc_d(j,1) = nanmean(nsctmp(id));
                nsc_d(j,2) = nanstd(nsctmp(id))/sqrt(length(id));
                nsc_d(j,3) = nanmean(FRmax(id));
                nsc_d(j,4) = nanstd(FRmax(id))/sqrt(length(id));
            else
                nsc_d(j,:) = [nan, nan,nan,nan];
            end
        end
        figure(5), 
        subplot(3,5,5+k), shadedErrorBar(x(1:end-1)+0.05,nsc_d(:,1),nsc_d(:,2))
        xlabel('Signal correlation'); ylabel('noise correlation')
        title(titlelabel{k}); axis([-0.5 1 0 0.2]);  box off
        figure(6);
        subplot(122), plot(x(1:end-1)+0.05,nsc_d(:,1)); hold on
        xlabel('Signal correlation'); ylabel('Normalized correlation')
        nsc_RsigInROI(:,k) = nsc_d(:,1);
        mean_pearson(k) = nanmean(nanmean(pearsontmp,2));
        se_pearson(k) = nanstd(nanmean(pearsontmp,2))/sqrt(size(pearsontmp,1));
        numinpool(k) = size(pearsontmp,1);
        mean_r_sig(k) = nanmean(r_sigtmp);
        se_r_sig(k) = nanstd(r_sigtmp)/sqrt(length(pearsontmp));
    end
end
numinpool
figure (10),
subplot(2,2,1), barwitherr(se_r_sig,mean_r_sig)
set(gca,'XTickLabel',titlelabel); title('signal correlation within area')
box off
subplot(2,2,2), barwitherr(se_pearson, mean_pearson)
set(gca,'XTickLabel',titlelabel); title('noise correlation within area')
box off
p_in = [];
for k1 = 1:5
    for k2 = k1+1:5
    [~, p_in(k1,k2)] = ttest2(NSCin_pool{k1}(:,1),NSCin_pool{k2}(:,1));
    end
end

%% NSC vs. max FR
titlelabel = {'V1', 'LM', 'AL', 'PM','LI'};
frackeep_in = []; 
for k = 1:5
    [id1] = find(areaID==k-1);
    if ~isempty(id1)
        % correlation to preferred stimulus
        pearsontmp = []; r_sigtmp = [];FRtmp = []; 
        dtmp = []; pearsonmax = [];FRmax = []; 
        for k1 = 1:length(id1)
            if ~isempty(nscinROI(id1(k1)).nscall_SC)
                dtmp = cat(1, dtmp, nscinROI(id1(k1)).d);
                r_sigtmp = cat(1,r_sigtmp,nscinROI(id1(k1)).r_sig);
                FRmeantmp = nanmean(nscinROI(id1(k1)).meanFR,2);
                tmp = nscinROI(id1(k1)).nscall_SC;
                FRmaxtmp = nscinROI(id1(k1)).maxFR;
                FRmax = cat(1,FRmax, FRmaxtmp);
                FRtmp = cat(1,FRtmp, FRmeantmp);
                pearsontmp = cat(1,pearsontmp,tmp);
            end
        end
        r_sigtmp(FRmax<frthreshold) = [];
        pearsontmp(FRmax<frthreshold,:) = [];
        dtmp(FRmax<frthreshold) = [];        
        frackeep_in(k) = sum(FRmax<frthreshold)/length(FRmax);
        FRtmp(FRmax<frthreshold) = [];
        FRmax(FRmax<frthreshold) = [];
        nsctmp = pearsontmp;
        % noise correlation vs. fr
        nsc_d = []; x = 0:0.05:0.2;
        for j = 1:length(x)-1
            id = intersect(find(FRtmp>x(j)), find(FRtmp<=x(j+1)));
            if length(id)>20
                nsc_d(j,1) = nanmean(nsctmp(id));
                nsc_d(j,2) = nanstd(nsctmp(id))/sqrt(length(id));
            else
                nsc_d(j,:) = [nan, nan];
            end
        end
        figure(7), 
        subplot(1,5,k), shadedErrorBar(x(1:end-1),nsc_d(:,1),nsc_d(:,2))
        xlabel('mean FR'); ylabel('Noise correlation')
        title(titlelabel{k}); % axis([0.4 2 0 0.15]);  box off
    end
end

%% NSC vs. signal correlation binning
mean_r_sig = [];  se_r_sig = []; se_pearsonmax = []; mean_pearsonmax = [];
se_pearson = []; mean_pearson = []; 
clear NSCBTW_pool 
titlelabel = {'V1-V1','V1-LM','V1-AL','V1-PM','V1-LI','LM-LI'};
frthreshold = 0.4; nsc_RsigBtwROI = [];
frackeep_btw = []; numbtwpool = []; nsc_distBtwROI = [];
for k = 1:6
    if k<=5
        id1 = intersect(find(areapair(:,2)==k-1),find(areapair(:,1)==0));
    else
        id2 = intersect(find(areapair(:,2)==4),find(areapair(:,1)==1));
        id3 = intersect(find(areapair(:,2)==1),find(areapair(:,1)==4));        
        id1 = [id2;id3];
    end        
    if ~isempty(id1)
        pearsontmp = []; r_sigtmp = [];FRtmp = []; 
        dtmp = []; pearsonmax = [];FRmax = [];
        for k1 = 1:length(id1)
            if ~isempty(nscBTWROI(id1(k1)).nscall_SC)
                r_sigtmp = cat(1,r_sigtmp,nscBTWROI(id1(k1)).r_sig);
                dtmp = cat(1, dtmp, nscBTWROI(id1(k1)).d);
                tmp = nscBTWROI(id1(k1)).nscall_SC;
                FRmaxtmp = nscBTWROI(id1(k1)).maxFR;
                FRmax = cat(1,FRmax, FRmaxtmp);
                pearsontmp = cat(1,pearsontmp,tmp);
            end
        end
        r_sigtmp(FRmax<frthreshold) = [];
        pearsontmp(FRmax<frthreshold,:) = [];
        frackeep_btw(k) = sum(FRmax<frthreshold)/length(FRmax);
        dtmp(FRmax<frthreshold) = [];        
        FRmax(FRmax<frthreshold) = [];

        NSCBTW_pool{k} = [pearsontmp, r_sigtmp, FRtmp];
        % noise correlation vs. distance 
%         if k==0
%             nsctmp = pearsontmp; 
%         nsctmp = pearsontmp(r_sigtmp>0.5); 
%         nsc_d = [];
%         dtmp = dtmp(r_sigtmp>0.5);
%         x = 200:80:1500;
%         for j = 1:length(x)-1
%             id = intersect(find(dtmp>x(j)), find(dtmp<=x(j+1)));
%             if length(id)>100
%                 nsc_d(j,1) = nanmean(nsctmp(id));
%                 nsc_d(j,2) = nanstd(nsctmp(id))/sqrt(length(id));
%             else
%                 nsc_d(j,:) = [nan, nan];
%             end
%         end
%         nsc_distBtwROI(:,k) = nsc_d(:,1);
%         figure(5), 
%         subplot(3,5,1),hold on, shadedErrorBar(x(1:end-1)+40,nsc_d(:,1),nsc_d(:,2))
%         xlabel('Neuron distance'); ylabel('Noise correlation'); box off
%         figure(6);
%         subplot(121), plot(x(1:end-1)+10,nsc_d(:,1)); hold on
%         xlabel('Neuron distance'); ylabel('Noise correlation to preferred stim')
%         end
        % noise correlation vs. signal correlation
        if k>1
        nsctmp = pearsontmp; nsc_d = []; x = -1:0.1:1;
        for j = 1:length(x)-1
            id = intersect(find(r_sigtmp>x(j)), find(r_sigtmp<=x(j+1)));
            if length(id)>20
                nsc_d(j,1) = nanmean(nsctmp(id));
                nsc_d(j,2) = nanstd(nsctmp(id))/sqrt(length(id));
            else
                nsc_d(j,:) = [nan, nan];
            end
        end
        figure(5),
        subplot(3,5,10+k-1), shadedErrorBar(x(1:end-1)+0.05,nsc_d(:,1),nsc_d(:,2))
        xlabel('Signal correlation'); ylabel('noise correlation'); box off
        title(titlelabel{k}); axis([-0.5 1 0 0.2])
        figure(6);
        subplot(122), plot(x(1:end-1)+0.05,nsc_d(:,1),'--'); hold on
        xlabel('Signal correlation'); ylabel('Normalized correlation')
        end
        mean_pearson(k) = nanmean(nanmean(pearsontmp,2));
        se_pearson(k) = nanstd(nanmean(pearsontmp,2))/sqrt(size(pearsontmp,1));
        numbtwpool(k,1) = size(pearsontmp,1);
        numbtwpool(k,2) = length(id1);
        mean_r_sig(k) = nanmean(r_sigtmp);
        se_r_sig(k) = nanstd(r_sigtmp)/sqrt(numbtwpool(k));
        nsc_RsigBtwROI(:,k) = nsc_d(:,1);
    end
end
numbtwpool
figure(10)
subplot(2,2,3), barwitherr(se_r_sig,mean_r_sig)
set(gca,'XTickLabel',titlelabel); title('signal correlation between area'); box off
subplot(2,2,4), barwitherr(se_pearson, mean_pearson)
set(gca,'XTickLabel',titlelabel); title('noise correlation between area'); box off
for k = 1:5
    for k2 = 1:5
    [~, p_btw(k,k2)] = ttest2(NSCBTW_pool{k}(:,1), NSCBTW_pool{k2}(:,1));
    end
end

%% plot histogram for correlation results
histIn = []; titlelabel = {'V1', 'LM', 'AL', 'PM','LI'};
x1 = -0.16:0.04:0.36;
for k = 1:5
    histIn(:,k) = histcounts(NSCin_pool{k}(:,1),x1)/numinpool(k,1);
end
figure, subplot(2,2,1), plot(x1(1:end-1)+0.02, histIn)
box off; xlabel('Fraction'), ylabel('Noise correlation'); title('Histogram')
subplot(2,2,2), plot(x1(1:end-1)+0.02, cumsum(histIn))
box off; xlabel('Fraction'), ylabel('Noise correlation'); title('Cumulative histogram')
legend(titlelabel)

histbtw = []; titlelabel = {'V1-V1','V1-LM','V1-AL','V1-PM','V1-LI'};
x1 = -0.16:0.04:0.36;
for k = 1:5
    histbtw(:,k) = histcounts(NSCBTW_pool{k}(:,1),x1)/numbtwpool(k,1);
end
subplot(2,2,3), plot(x1(1:end-1)+0.02, histbtw)
box off; xlabel('Fraction'), ylabel('Noise correlation'); title('Histogram')
subplot(2,2,4), plot(x1(1:end-1)+0.02, cumsum(histbtw))
box off; xlabel('Fraction'), ylabel('Noise correlation'); title('Cumulative histogram')
legend(titlelabel)

%% fit exponential function (nsc vs. r_sig)
figure, 
x1 = -1:0.1:1; rmse_in_rsig = [];
for i = 1:5
    x = x1(1:end-1)'+0.05;
    y = nsc_RsigInROI(:,i);
    x(isnan(y)) = [];
    y(isnan(y)) = [];
    f = fittype(@(a,b,x) b.*a.^x,'dependent',{'y'},'independent',{'x'});
    [obj1, gof1] = fit(x, y,f,'StartPoint',[1, 0.1]);
    rmse_in_rsig(i) = gof1.rmse;
    a(i) = obj1.a;
    b(i) = obj1.b;
    hold on, plot(x1,b(i).*a(i).^x1)
end

x1 = -1:0.1:1; rmse_btw_rsig = [];
for i = 2:5
    x = x1(1:end-1)'+0.05;
    y = nsc_RsigBtwROI(:,i);
    x(isnan(y)) = [];
    y(isnan(y)) = [];
    f = fittype(@(a,b,x) b.*a.^x,'dependent',{'y'},'independent',{'x'});
    [obj1, gof1] = fit(x, y,f,'StartPoint',[1, 0.1]);
    rmse_btw_rsig(i) = gof1.rmse;
    a(i) = obj1.a;
    b(i) = obj1.b;
    hold on, plot(x1,b(i).*a(i).^x1,'--')
end
title('Exponential fitting')
xlabel('Signal correlation'); ylabel('Noise correlation')
titlelabel1 = {'V1', 'LM', 'AL', 'PM','LI'};
titlelabel2 = {'V1-LM','V1-AL','V1-PM','V1-LI'};
legend([titlelabel1,titlelabel2])

%%  fit exponential function (nsc vs. neuron distance)
figure, 
x1 = 5:80:500; rmse_in_dist = [];
for i = 1:5
    x = x1(1:end-1)'+20;
    y = nsc_distInROI(:,i);
    x(isnan(y)) = [];
    y(isnan(y)) = [];
    f = fittype(@(a,b,x) b.*a.^x,'dependent',{'y'},'independent',{'x'});
    [obj1, gof1] = fit(x, y,f,'StartPoint',[1, 0.1]);
    rmse_in_dist(i) = gof1.rmse;
    a(i) = obj1.a;
    b(i) = obj1.b;
    hold on, plot(x1,b(i).*a(i).^x1)
end
title('Exponential fitting')
xlabel('Neuron distance (\mum)'); ylabel('Noise correlation')
titlelabel1 = {'V1', 'LM', 'AL', 'PM','LI'};
legend(titlelabel1)
