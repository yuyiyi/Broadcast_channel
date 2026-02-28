load('pooldata_ComboStim.mat')
stimInd = [1,2,3];
SpikeCount_Twindow = 1; % compute spike count in 1 sec time window following IgnoreTimeWindow 
StimTimeIndv = 1; % duration of each stimulus is 2 sec
IgnoreTimeWindow = 0; % ignore first 0 sec response

%% spike counnt in 1 sec bin
clear spikecount num
for k1 = 1:size(pooldata,1)
    for k2 = 1:2
        if isempty(pooldata(k1,k2).spiketrain)
            continue
        end
        spiketrain = pooldata(k1,k2).spiketrain;
        goodcell = pooldata(k1,k2).goodcell; 
        imgPara = pooldata(k1,k2).imgPara;
        num(k1,k2) = length(goodcell);
        roi = pooldata(k1,k2).Region;
        loc = regionprops(roi,'Centroid');
        [spikecountAve, spikecountraw, SCvar] = ...
            MovieSpikecount(spiketrain, goodcell, imgPara, stimInd, StimTimeIndv,...
            SpikeCount_Twindow, IgnoreTimeWindow);
        spikecount{k1,k2}.spikecountAve = spikecountAve;
        spikecount{k1,k2}.variance = SCvar;
        spikecount{k1,k2}.spikecountrawSC = spikecountraw;
        spikecount{k1,k2}.neuloc = loc(goodcell);
        spikecount{k1,k2}.ROIloc = ROIloc(k1,(k2-1)*2+1:k2*2)*250;
    end
end

%% noise correlation vs. group (between area)
d1 = 5000/512; d2 = 5000/512;
clear NC_GratvsNat areapair 
frthreshold = 0.2; 

k = 0; 
for k1 = 1:size(spikecount,1)
    for k2 = 1:2
        for k3 = k2:2
            if isempty(spikecount{k1,k2}) || isempty(spikecount{k1,k3})
                continue
            end
            if isempty(spikecount{k1,k2}.spikecountAve) || isempty(spikecount{k1,k3}.spikecountAve)
                continue
            end
            scA = spikecount{k1,k2}.spikecountrawSC; 
            R = size(scA,2);
            if R<10
                continue
            end
            scB = spikecount{k1,k3}.spikecountrawSC;
            FRA = spikecount{k1,k2}.spikecountAve;
            frm = squeeze(max(FRA, [],2));
            iikeepA = intersect(find(frm(:,1)>frthreshold), find(max(frm(:,2:3),[],2)>frthreshold));            
            FRB = spikecount{k1,k3}.spikecountAve;
            frm = squeeze(max(FRB, [],2));
            iikeepB = intersect(find(frm(:,1)>frthreshold), find(max(frm(:,2:3),[],2)>frthreshold));
            
            locA = spikecount{k1,k2}.neuloc;
            locB = spikecount{k1,k3}.neuloc;
            ROIlocA = spikecount{k1,k2}.ROIloc;
            ROIlocB = spikecount{k1,k3}.ROIloc;

            areaID1 = areaID_correct(k1,k2)+1;
            areaID2 = areaID_correct(k1,k3)+1;
            if areaID1>5 || areaID2>5
                continue
            end         
            if ~isempty(iikeepA) && ~isempty(iikeepB)                
                A = ones(length(iikeepA), length(iikeepB));
                if k2==k3
                    A = triu(A,1);
                end
                FR1 = FRA(iikeepA,:,:);                 
                FR2 = FRB(iikeepB,:,:); 
                meanFR1 = []; meanFR2 = []; meanFR = [];
                for n1 = 1:length(iikeepA)
                    for n2 = 1:length(iikeepB)
                        meanFR1(n1,n2) = sqrt(sum(FR1(n1,:,1),2).*sum(FR2(n2,:,1),2)); 
                        meanFR2(n1,n2) = sqrt(sum(sum(FR1(n1,:,2:3),2),3).*sum(sum(FR2(n2,:,2:3),2),3)); 
                    end
                end
                meanFR = [meanFR1(A==1), meanFR2(A==1)];
                
                % NC to grating
                sc1 = scA(iikeepA,:,:,1);
                sc2 = scB(iikeepB,:,:,1);
                tmp = NeuCor(sc1, sc2, 1);
                pearsongratall= tmp(A==1);
                pearsongrat = [];
                for ri = 1:10
                    ri1 = randperm(R);
                    pearsontmp =  NeuCor(sc1(:,ri1(1:ceil(R/2)),:,:), sc2(:,ri1(1:ceil(R/2)),:,:), 1);
                    pearsongrat(:, ri) = pearsontmp(A==1);
                end
                pearson_shulf = NeuCor_shulf(sc1, sc2, 1);
                pearson_shulf_grat= pearson_shulf(A==1);
                % NC to nature movie
                sc1 = scA(iikeepA,:,:,2:3);
                sc1 = reshape(sc1, length(iikeepA), R, []);
                sc2 = scB(iikeepB,:,:,2:3);
                sc2 = reshape(sc2, length(iikeepB), R, []);
                tmp = NeuCor(sc1, sc2, 1);
                pearsonNatall= tmp(A==1);
                pearsonNat = [];
                for ri = 1:10
                    ri1 = randperm(R);
                    pearsontmp =  NeuCor(sc1(:,ri1(1:ceil(R/2)),:,:), sc2(:,ri1(1:ceil(R/2)),:,:), 1);
                    pearsonNat(:, ri) = pearsontmp(A==1);
                end
                pearson_shulf = NeuCor_shulf(sc1, sc2, 1);
                pearson_shulf_nat= pearson_shulf(A==1);
               
                r_sig1 = corr(FR1(:,:,1)', FR2(:,:,1)');
                r_sig2 = corr(reshape(FR1(:,:,2:3),size(FR1,1),64)', reshape(FR2(:,:,2:3),size(FR2,1),64)');
                
                d = NeuDistance(locA(iikeepA), locB(iikeepB), d1, d2, ROIlocA, ROIlocB);
                if sum(A(:))>5
                    k = k+1;
                    areapair(k,1) = areaID_correct(k1,k2);
                    areapair(k,2) = areaID_correct(k1,k3);
                    areapair(k,3:5) = [k1,k2,k3];
                    NC_GratvsNat(k).nscGrat_SC = pearsongrat;
                    NC_GratvsNat(k).nscNat_SC = pearsonNat;
                    NC_GratvsNat(k).nsc_SC = [pearsongratall, pearsonNatall];
                    NC_GratvsNat(k).nscshulf_SC = [pearson_shulf_grat, pearson_shulf_nat];
                    NC_GratvsNat(k).meanFR = meanFR;
                    NC_GratvsNat(k).r_sig = [r_sig1(A==1), r_sig2(A==1)];
                    NC_GratvsNat(k).neudist = d(A==1);
                end
            end       
        end
    end
end

%% number of NC pairs
NCnum = [];
for k = 1:26
    id = find(areapair(:,3)==k);
    n = 0;
    for i = 1:length(id)
    n = n+size(NC_GratvsNat(id(i)).nsc_SC,1);
    end
    NCnum(k,1) = n;
end

NC_area = [];
for k = 1:6
    id = find(areaID_correct==k-1);
    n = 0;
    for i = 1:length(id)
        n = n + size(NC_GratvsNat(id(i)).nsc_SC,1);
    end
    NC_area(k) = n;
end

%% CI of NC estimation for individual neuron
% bootstraps percentile 
var_explained = []
N = length(NC_GratvsNat);
for k = 1:length(NC_GratvsNat)
    NC_grat = NC_GratvsNat(k).nscGrat_SC;
    a = quantile(NC_grat, [0.025, 0.975], 2);
    % variance explained by random subest of trials 
    mdl = fitlm(NC_grat(:,1), NC_grat(:,2));
    var_explained(k,1) = mdl.Rsquared.Ordinary;

    NC_nat = NC_GratvsNat(k).nscNat_SC;
    b = quantile(NC_nat, [0.025, 0.975], 2);
    % variance explained by random subest of trials 
    mdl = fitlm(NC_nat(:,1), NC_nat(:,2));
    var_explained(k,2) = mdl.Rsquared.Ordinary;    
end


%% correlation(NC_grat, NC_nat) vs. correlation (SC_grat, SC_nat)
c_nc2 = [];
c_nc = []; 
c_rig = [];
N = length(NC_GratvsNat);
for k = 1:length(NC_GratvsNat)
    NC_grat = NC_GratvsNat(k).nscGrat_SC;
    if length(NC_grat)>20
    NC_nat = NC_GratvsNat(k).nscNat_SC;
    NC_shulf = NC_GratvsNat(k).nscshulf_SC;
    nsc_SC = NC_GratvsNat(k).nsc_SC;
    r_sig = NC_GratvsNat(k).r_sig;
    
    mdl = fitlm(zscore(nsc_SC(:,1)), nsc_SC(:,2));
    var_explained(k,1) = mdl.Rsquared.Ordinary;
    mdl = fitlm(zscore(r_sig(:,2)), nsc_SC(:,2));
    var_explained(k,2) = mdl.Rsquared.Ordinary;
    mdl = fitlm([zscore(r_sig(:,2)),zscore(nsc_SC(:,1))], nsc_SC(:,2));
    var_explained(k,3) = mdl.Rsquared.Ordinary;

    c_nc(k) = corr(nanmean(NC_grat,2), nanmean(NC_nat,2));
    c_nc2(k) = corr(nsc_SC(:,1), nsc_SC(:,2));
    c_ncshulf(k) = corr(NC_shulf(:,1), NC_shulf(:,2));
    c_rig(k) = corr(r_sig(:,1), r_sig(:,2));
    else
        c_nc2(k) = nan;
        c_nc(k) = nan; 
        c_ncshulf(k) = nan;
        c_rig(k) = nan; 
        var_explained(k,1:3) = nan;
    end
end
m1 = [nanmean(c_nc), nanstd(c_nc)];
m2 = [nanmean(c_rig), nanstd(c_rig)];
m3 = [nanmean(c_ncshulf), nanstd(c_ncshulf)];

figure(4), clf('reset') 
for k = 1:5
    if k<5
        ii = intersect(find(areapair(:,1)==k-1), find(areapair(:,2)==k-1));
    else
        ii = find(areapair(:,1)-areapair(:,2)~=0);
    end
        hold on, plot(c_rig(ii), c_nc(ii), 'o', 'markersize',4)
end
hold on, plot(c_rig, c_ncshulf, '.', 'markersize', 10, 'color', [0.8 0.8 0.8])
errorbar(m2(1), m1(1), m1(2), m1(2), m2(2), m2(2), '.k', 'markersize', 15)
errorbar(m2(1), m3(1), m3(2), m3(2), m2(2), m2(2), '.', 'markersize', 15, 'color', [0.5 0.5 0.5])
hold on, plot(-0.2:0.01:0.3, -0.2:0.01:0.3, 'k')
box off
xlabel('corr(SCgrat, SCnat)')
ylabel('corr(NCgrat, NCnat)')
legend('V1', 'LM', 'AL', 'PM', 'V1-HVA', 'trialshulf')
axis([-0.2 0.3 -0.2 0.6])


%% show an example
k = 2;
% NC_grat = NC_GratvsNat(k).nscGrat_SC;
% NC_nat = NC_GratvsNat(k).nscNat_SC;
NC_grat = NC_GratvsNat(k).nsc_SC(:,1);
NC_nat = NC_GratvsNat(k).nsc_SC(:,2);
NC_shulf = NC_GratvsNat(k).nscshulf_SC;
r_sig = NC_GratvsNat(k).r_sig;

binwidth = 0.01;
xx = -0.6:binwidth:1; 
[cc, xedges, yedges, xbin, ybin] = histcounts2(nanmean(NC_grat,2), nanmean(NC_nat,2), xx, xx);
cc = cc/size(NC_grat,1);
% figure, imagesc(xx(1:end-1), xx(1:end-1), cc)

hmax = max(cc(:));
colgrad = 10;
c1 = colormap(jet(colgrad));
figure(5), clf('reset')
subplot(121), hold on
for i = 1:size(NC_grat,1)
    a = nanmean(NC_grat(i,:));
    b = nanmean(NC_nat(i,:));
    c2 = ceil(cc(xbin(i), ybin(i))/hmax*colgrad);
    plot(a,b, '.', 'color', c1(c2,:))
end
axis([-0.2 0.4 -0.2 0.4])
xlabel('NC to grating'), ylabel('NC to Nature movie')

subplot(122), hold on,
% plot(r_sig(:,1), r_sig(:,2), '.k')
[cc, xedges, yedges, xbin, ybin] = histcounts2(r_sig(:,1), r_sig(:,2), xx, xx);
cc = cc/size(NC_grat,1);
for i = 1:size(NC_grat,1)
    a = r_sig(i,1);
    b = r_sig(i,2);
    c2 = ceil(cc(xbin(i), ybin(i))/hmax*colgrad);
    plot(a,b, '.', 'color', c1(c2,:))
end
box off
xlabel('SC to grating'), ylabel('SC to Nature movie')

%% NC to natural movie variance explained 
NC_grat = []; 
NC_nat = [];
nsc_SC = [];
r_sig = [];
for k = 1:length(NC_GratvsNat)
    NC_grat = cat(1, NC_grat, NC_GratvsNat(k).nscGrat_SC);
    NC_nat = cat(1, NC_nat, NC_GratvsNat(k).nscNat_SC);
    nsc_SC = cat(1, nsc_SC, NC_GratvsNat(k).nsc_SC);
    r_sig = cat(1, r_sig, NC_GratvsNat(k).r_sig);
end
nsc_SC_norm = [zscore(nsc_SC(:,1)), nsc_SC(:,2)];
r_sig_norm = zscore(r_sig);
x1 = -1:0.1:1;
% plot SC vs. NC
[h1, ~, ind] = histcounts(r_sig_norm(:,2), x1);
m_sc = [];
for i = 1:length(x1)-1
    ii = find(ind == i);
    if ~isempty(ii)
        m_sc(i,1:2) = [nanmean(nsc_SC_norm(ii,2)), nanstd(nsc_SC_norm(ii,2))/sqrt(length(ii))];
        m_sc(i,3:4) = [nanmean(r_sig_norm(ii,2)), nanstd(r_sig_norm(ii,2))/sqrt(length(ii))];
    else
        m_sc(i,1:4) = [nan, nan, nan, nan];
    end
end
figure(6), clf('reset'), 
subplot(121), shadedErrorBar(x1(1:end-1), m_sc(:,1), m_sc(:,2), '-b')
box off, xlabel('normalized SC to natural video'), ylabel('NC to natural video')

[h1, ~, ind] = histcounts(nsc_SC_norm(:,1), x1);
m_nc = [];
for i = 1:length(x1)-1
    ii = find(ind == i);
    if ~isempty(ii)
        m_nc(i,1:2) = [nanmean(nsc_SC_norm(ii,2)), nanstd(nsc_SC_norm(ii,2))/sqrt(length(ii))];
        m_nc(i,3:4) = [nanmean(nsc_SC_norm(ii,1)), nanstd(nsc_SC_norm(ii,1))/sqrt(length(ii))];
    else
        m_nc(i,1:4) = [nan, nan, nan, nan];
    end
end
hold on, shadedErrorBar(x1(1:end-1), m_nc(:,1), m_nc(:,2),'-r')
box off, xlabel('normalized NC to grating'), ylabel('NC to natural video')

N = size(nsc_SC,1);
var_explained = [];
n = 1000;
for perm = 1:500
    ii = randperm(N);
    i1 = ii(1:n);
    mdl = fitlm(nsc_SC_norm(i1,1), nsc_SC_norm(i1,2));
    var_explained(perm,2) = mdl.Rsquared.Ordinary;
    mdl = fitlm(r_sig_norm(i1,2), nsc_SC_norm(i1,2));
    var_explained(perm,1) = mdl.Rsquared.Ordinary;
    mdl = fitlm([r_sig_norm(i1,2), nsc_SC_norm(i1,1)], nsc_SC_norm(i1,2));
    var_explained(perm,3) = mdl.Rsquared.Ordinary;
end
mean(var_explained)
[h,p] = ttest(var_explained(:,1), var_explained(:,2))
figure(6), subplot(122)
boxplot(var_explained*100, 'whisker', 10)
box off, ylabel('variance explained')
set(gca, 'XTickLabel', {'SC', 'NC', 'SC & NC'})
xlabel('Linear predictor')
g = ones(1500,1);
g(501:1000) = 2;
g(1001:1500) = 3;
var_explained = reshape(var_explained, [], 1);
var_explained = cat(2, var_explained*100, g);

%%
ft = fittype('1./(1 + exp(-(a.*x + b)))', 'independent', 'x', 'dependent', 'y');
% ft = fittype('a+c*b.^x', 'independent', 'x', 'dependent', 'y');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [10, 5];

x1 = -1:0.1:1;
% plot SC vs. NC
[h1, ~, ind] = histcounts(r_sig(:,2), x1);
m_sc = [];
for i = 1:length(x1)-1
    ii = find(ind == i);
    if ~isempty(ii)
        m_sc(i,1:2) = [nanmean(nsc_SC(ii,2)), nanstd(nsc_SC(ii,2))/sqrt(length(ii))];
        m_sc(i,3:4) = [nanmean(r_sig(ii,2)), nanstd(r_sig(ii,2))/sqrt(length(ii))];
    else
        m_sc(i,1:4) = [nan, nan, nan, nan];
    end
end
figure(7), clf('reset'), 
subplot(131), 
shadedErrorBar(x1(1:end-1), m_sc(:,1), m_sc(:,2))
m_sc(isnan(m_sc(:,1)),:) = [];
[fitresult1, gof1] = fit(m_sc(:,3), m_sc(:,1), ft, opts);
% [fitresult1, gof1] = fit(r_sig(:,2), nsc_SC(:,2), ft, opts);
hold on, plot(fitresult1)
box off, xlabel('SC to natural video'), ylabel('NC to natural video')
% axis([-0.5 1 -0.05 1])

[h1, ~, ind] = histcounts(nsc_SC(:,1), x1);
m_nc = [];
for i = 1:length(x1)-1
    ii = find(ind == i);
    if ~isempty(ii)
        m_nc(i,1:2) = [nanmean(nsc_SC(ii,2)), nanstd(nsc_SC(ii,2))/sqrt(length(ii))];
        m_nc(i,3:4) = [nanmean(nsc_SC(ii,1)), nanstd(nsc_SC(ii,1))/sqrt(length(ii))];
    else
        m_nc(i,1:4) = [nan, nan, nan, nan];
    end
end
subplot(132), 
shadedErrorBar(x1(1:end-1), m_nc(:,1), m_nc(:,2))
m_nc(isnan(m_nc(:,1)),:) = [];
[fitresult2, gof2] = fit(m_nc(:,3), m_nc(:,1), ft, opts);
% [fitresult2, gof2] = fit( nsc_SC(:,1), nsc_SC(:,2), ft, opts);
hold on, plot(fitresult2)
box off, xlabel('NC to grating'), ylabel('NC to natural video')
% axis([-0.5 1 -0.05 1])
ncpred2 = predint(fitresult2,nsc_SC(:,1)); % 95% confidence interval for prediction

%% fit sigmoidal function for individual experiment
ft = fittype('1./(1 + exp(-(a.*x + b)))', 'independent', 'x', 'dependent', 'y');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [10, 5];

var_explained = [];
for k = 1:length(NC_GratvsNat)
    NC_grat = NC_GratvsNat(k).nscGrat_SC;
    if length(NC_grat)>20
    NC_nat = NC_GratvsNat(k).nscNat_SC;
    NC_shulf = NC_GratvsNat(k).nscshulf_SC;
    nsc_SC = NC_GratvsNat(k).nsc_SC;
    r_sig = NC_GratvsNat(k).r_sig;
    
    [fitresult0, gof1] = fit(r_sig(:,2), nsc_SC(:,2), ft, opts);
    var_explained(k,1) = gof1.rsquare;
    [fitresult0, gof1] = fit(nsc_SC(:,1), nsc_SC(:,2), ft, opts);
    var_explained(k,2) = gof1.rsquare;
    else
        var_explained(k,1:2) = nan;
    end
end
mean(var_explained)
var_explained(isnan(var_explained(:,1)),:) = [];
[h,p] = ttest(var_explained(:,1), var_explained(:,2))

%% sigmoidal fitting for random subset
NC_grat = []; 
NC_nat = [];
nsc_SC = [];
r_sig = [];
for k = 1:length(NC_GratvsNat)
    NC_grat = cat(1, NC_grat, NC_GratvsNat(k).nscGrat_SC);
    NC_nat = cat(1, NC_nat, NC_GratvsNat(k).nscNat_SC);
    nsc_SC = cat(1, nsc_SC, NC_GratvsNat(k).nsc_SC);
    r_sig = cat(1, r_sig, NC_GratvsNat(k).r_sig);
end

ft = fittype('1./(1 + exp(-(a.*x + b)))', 'independent', 'x', 'dependent', 'y');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [10, 5];
N = size(nsc_SC,1);
var_explained = [];
n = 2000;
for perm = 1:100
    ii = randperm(N);
    i1 = ii(1:n);
    [fitresult0, gof1] = fit(r_sig(i1,2), nsc_SC(i1,2), ft, opts);
    var_explained(perm,1) = gof1.rsquare;
    [fitresult0, gof1] = fit(nsc_SC(i1,1), nsc_SC(i1,2), ft, opts);
    var_explained(perm,2) = gof1.rsquare;
end
mean(var_explained)
[h,p] = ttest(var_explained(:,1), var_explained(:,2))
figure(7), subplot(133)
boxplot(var_explained*100, 'whisker', 10)
box off, ylabel('variance explained')
set(gca, 'XTickLabel', {'SC', 'NC', 'SC & NC'})
xlabel('Linear predictor')

%% estimate NC using neural network


