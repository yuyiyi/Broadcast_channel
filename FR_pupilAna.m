load('pooldata_SinGrat.mat', 'pooldata','name1','areaID')
load('pooldata_pupiltrack.mat')
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
stimflag = [repmat(Dir',9,1)/360*2*pi, reshape(repmat(stimPara(:,1),1,8)',[],1),...
    reshape(repmat(stimPara(:,2),1,8)',[],1)];
[stimsort, stimInd] = sortrows(stimPara, [1 2]);
TmFreq = [1,2,8];
SpFreq = [0.02, 0.04, 0.16];
% TmFreq = stimsort(:,1);
% SpFreq = stimsort(:,2);
Twindow = [0 2 4 6 8 10 12 14 16]; % s 
delay = 0;
DirSortVect = [0 1; -sqrt(2)/2 sqrt(2)/2; -1 0;-sqrt(2)/2 -sqrt(2)/2; 0 -1;sqrt(2)/2 -sqrt(2)/2; 1 0;sqrt(2)/2 sqrt(2)/2];

%%
Npc = [];
for k1 = 1:length(poolpupildata)
    if ~isempty(poolpupildata(k1).pupildata)
        pupildata = poolpupildata(k1).pupildata;
        if isfield(pupildata, 'pupilsync')
            Npc(k1,1) = pupildata.video_PCdim;            
            Npc(k1,2) = pupildata.video_PCvar;
        end
    end
end
%%
binWidth = 100; 
smoothSD = 2;
popResp_var = [];
for k1 = 1:length(poolpupildata)
    if ~isempty(poolpupildata(k1).pupildata)
        pupildata = poolpupildata(k1).pupildata;
        if isfield(pupildata, 'pupilsync')
            pupilsync = pupildata.pupilsync;
            sync2p = table2array(poolpupildata(k1).sync2p);

            l0 =  find(pupilsync(:,3)==1,1);            
            l1 =  find(pupilsync(:,3)==1,1,'last');
            pupildyn = [pupildata.eyearea(l0:l1)',pupildata.eyecenter(l0:l1,:),pupildata.videoscore];
            pupildyn(pupildyn(:,1)==0,1:3) = nan;
            pupilPC(k1,1:2) = [pupildata.video_PCdim, pupildata.video_PCvar];
            pupildyn_t = pupilsync(l0:l1,1:2);
            pupildyn_t = bsxfun(@minus, pupildyn_t, pupildyn_t(1,:))+1;
            pupildyn_match2p = [];
            for j = 1:size(pupildyn,2)
                tmp = interp1(pupildyn_t(:,1)/pupildyn_t(end,1), pupildyn(:,j), sync2p(:,1)/sync2p(end,1));
                pupildyn_match2p = cat(2, pupildyn_match2p, tmp);
            end
            
            imgPara = pooldata(k1,1).imgPara;
            % pupil dynamics bined per direction (2 sec window)
            pupildyn_bin = [];
            for i = 1:imgPara.stimrep                
                for j = 1:imgPara.stim_type
                    tmp = pupildyn_match2p(sync2p(:,3)==i & sync2p(:,5)==j,:);
                    tmp1 = imresize(tmp, [8 size(tmp,2)]);                    
                    pupildyn_bin = cat(1, pupildyn_bin, tmp1);
                end
            end
            X_pupildyn = pupildyn_bin;
            X_stim = repmat(stimflag, imgPara.stimrep, 1);
            for k2 = 1:4
                if ~isempty(pooldata(k1,k2).spiketrain) && length(pooldata(k1,k2).goodcell) > 5
                goodcell = pooldata(k1,k2).goodcell; 
                spiketrain = pooldata(k1,k2).spiketrain(goodcell);
                num(k1,k2) = length(goodcell);
                spikecountAve = []; h = []; spikecountraw1 = []; 
                spikecountraw2 = []; snr = [];
                for n = 1:length(goodcell)
                    tmphist = []; 
                    stvec = [];
                    for i = 1:imgPara.stimrep
                        st1 = [];
                        for j = 1:imgPara.stim_type
                            st = spiketrain(n).st{i,j}*imgPara.dt;
                            st(st<0) = [];
                            st(st>imgPara.stim_time) = [];
                            st1 = [st1, st+imgPara.stim_time*(j-1)];
                            tmphist1 = histcounts(st, Twindow);
                            tmphist = cat(1, tmphist, tmphist1');
                        end
                        stvec = vectCat(stvec, st1');
                    end
                    [rate , rateT] = myPSTH(stvec*1000, [0, imgPara.stim_time*imgPara.stim_type]*1000, binWidth, smoothSD);
                    spikecountraw1(n,:) = tmphist;
                    spikecountAve(n,:) = rate;
                end
                
                % total variance of population response
                Y = spikecountraw1';
                Y = bsxfun(@minus, Y, mean(Y,1));                                
                Y_totvar = var(Y(:));
                
                % variance explained by stimulus induced response
                Ym = mean(reshape(Y, 72,imgPara.stimrep,size(Y,2)), 2);
                Ym_var = var(Ym(:));
                Y_stimexp = Ym_var/Y_totvar;       
                
                % residual population response after subtracting stimulus
                % induced portion (estimated by trial averaged response)
                Y = bsxfun(@minus, reshape(Y, 72,imgPara.stimrep,size(Y,2)), Ym);
                Y = reshape(Y,[],size(Y,3));
                Y_resvar = var(Y(:));             
                
                % ordinary regression variance explained by pupil dynamics 
                X = X_pupildyn;
                X(isnan(X)) = 0;
%                 Bols = inv(X'*X)*X'*Y;
%                 Yhat2 = X*Bols;
%                 e = mean((Yhat2-Y).^2);
%                 Y_pupildynexp_OLS = 1-mean(e./var(Y));               
                
                e = []; d = [];
                for k = 1:size(X,2) 
                    [Br, Bols, e(k), d(k)] = RRR_bySVD(X, Y, k);
                end
%                 Y_pupildynexp_RRR = 1 - min(e);
                [~, k] = min(e);
                [Br, Bols, ~, ~] = RRR_bySVD(X, Y, k);
                Yhat2 = X * Br;
                Y_pupildynexp_RRR = var(Yhat2(:))/Y_totvar;
                % residual network response factoring out pupil dynamics
                % correlated portion (graph model)
                Y_res = Y - Yhat2;
                r = []; ypred = [];
                for n = 1:size(Y,2)
                    X = Y_res;
                    X(:,n) = [];
                    [c,s,l] = pca(X);
                    pcscore = s(:,1:find(l<0.1,1));
                    mdl2 = fitlm(pcscore, Y_res(:,n));
%                     mdl2 = fitlm(X, Y_res(:,n));
                    r(n) = mdl2.Rsquared.Ordinary;
                    ypred(:,n) = predict(mdl2,pcscore);
                end
                Y_net = var(ypred(:))/Y_totvar;
%                 Y_net = mean(r);                
                popResp_var = cat(1, popResp_var,[Y_stimexp, Y_pupildynexp_RRR, Y_net,areaID(k1,k2), k]);                
                end
            end 
            
        end
    end
end

figure, boxplot(popResp_var(:,1:3), 'whisker', 10)
box off
set(gca, 'Xticklabels', {'Stimulus','Pupil dynamics','Network'})
ylabel('Variance explained')
title('Partition variance of population neuron response ')

%% 1D pupil dynamics
binWidth = 100; 
smoothSD = 2;
ind = 1;
for k1 = 1:length(poolpupildata)
    if ~isempty(poolpupildata(k1).pupildata)
        pupildata = poolpupildata(k1).pupildata;
        if isfield(pupildata, 'pupilsync')
            pupilsync = pupildata.pupilsync;
            sync2p = table2array(poolpupildata(k1).sync2p);

            l0 =  find(pupilsync(:,3)==1,1);            
            l1 =  find(pupilsync(:,3)==1,1,'last');
            pupildyn = [pupildata.eyearea(l0:l1)',pupildata.eyecenter(l0:l1,:),pupildata.videoscore];
            pupildyn(pupildyn(:,1)==0,1:3) = nan;
            pupilPC(k1,1:2) = [pupildata.video_PCdim, pupildata.video_PCvar];
            pupildyn_t = pupilsync(l0:l1,1:2);
            pupildyn_t = bsxfun(@minus, pupildyn_t, pupildyn_t(1,:))+1;
            pupildyn_match2p = [];
            for j = 1:size(pupildyn,2)
                tmp = interp1(pupildyn_t(:,1)/pupildyn_t(end,1), pupildyn(:,j), sync2p(:,1)/sync2p(end,1));
                pupildyn_match2p = cat(2, pupildyn_match2p, tmp);
            end
            
            imgPara = pooldata(k1,1).imgPara;
            % pupil dynamics bined per direction (2 sec window)
            pupildyn_bin = [];
            for i = 1:imgPara.stimrep                
                for j = 1:imgPara.stim_type
                    tmp = pupildyn_match2p(sync2p(:,3)==i & sync2p(:,5)==j,:);
                    tmp1 = imresize(tmp, [8 size(tmp,2)]);                    
                    pupildyn_bin = cat(1, pupildyn_bin, tmp1);
                end
            end
            X_pupildyn = pupildyn_bin;
            X_stim = repmat(stimflag, imgPara.stimrep, 1);
            for k2 = 1:4
                if ~isempty(pooldata(k1,k2).spiketrain) && length(pooldata(k1,k2).goodcell) > 5
                goodcell = pooldata(k1,k2).goodcell; 
                spiketrain = pooldata(k1,k2).spiketrain(goodcell);
                num(k1,k2) = length(goodcell);
                spikecountAve = []; h = []; spikecountraw1 = []; 
                spikecountraw2 = []; snr = [];
                for n = 1:length(goodcell)
                    tmphist = []; 
                    stvec = [];
                    for i = 1:imgPara.stimrep
                        st1 = [];
                        for j = 1:imgPara.stim_type
                            st = spiketrain(n).st{i,j}*imgPara.dt;
                            st(st<0) = [];
                            st(st>imgPara.stim_time) = [];
                            st1 = [st1, st+imgPara.stim_time*(j-1)];
                            tmphist1 = histcounts(st, Twindow);
                            tmphist = cat(1, tmphist, tmphist1');
                        end
                        stvec = vectCat(stvec, st1');
                    end
                    [rate , rateT] = myPSTH(stvec*1000, [0, imgPara.stim_time*imgPara.stim_type]*1000, binWidth, smoothSD);
                    spikecountraw1(n,:) = tmphist;
                    spikecountAve(n,:) = rate;
                end
                
                % total variance of population response
                Y = spikecountraw1';
                Y = bsxfun(@minus, Y, mean(Y,1));                                
                Y_totvar = var(Y(:));
                
                % variance explained by stimulus induced response
                Ym = mean(reshape(Y, 72,imgPara.stimrep,size(Y,2)), 2);
                Ym_var = var(Ym(:));
                Y_stimexp = Ym_var/Y_totvar;       
                
                % residual population response after subtracting stimulus
                % induced portion (estimated by trial averaged response)
                Y = bsxfun(@minus, reshape(Y, 72,imgPara.stimrep,size(Y,2)), Ym);
                Y = reshape(Y,[],size(Y,3));
                Y_resvar = var(Y(:));             
                
                % ordinary regression variance explained by pupil dynamics 
                X = X_pupildyn;
                X(isnan(X)) = 0;            
                [Br, Bols, ~, ~] = RRR_bySVD(X, Y, 1);
                Yhat2 = X * Br;
                Y_pupildynexp_RRR = var(Yhat2(:))/Y_totvar;
                
                popResp_var(ind, 6) = Y_pupildynexp_RRR;         
                ind = ind+1;
                end
            end 
            
        end
    end
end


figure, boxplot(popResp_var(:,[1,2,6,3]), 'whisker', 10)
box off
set(gca, 'Xticklabels', {'Stimulus','Pupil dynamics','1D Pupil dynamics','Network'})
ylabel('Variance explained')
title('Partition variance of population neuron response ')
