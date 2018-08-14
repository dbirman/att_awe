% load all subject decon.mat files, average everything together, then plot
% them.

%% Load
rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
rmode = {'resp','mresp'};
% subject | roi | model | index | value
time = nan(length(nSIDs),length(rois),2,20,81);
cc = nan(length(nSIDs),length(rois),2,20,81);
for ni = 1:length(nSIDs)
    load(fullfile(datafolder,sprintf('s%04.0f_decon.mat',nSIDs(ni))));
    for ri = 1:length(rois)
        for rj = 1:length(rmode)
            if ~any(ni==[3 4 5])
                time(ni,ri,rj,:,:) = decondata.(rois{ri}).time.(rmode{rj});
            end
            cc(ni,ri,rj,:,:) = decondata.(rois{ri}).cc.(rmode{rj});
        end
    end
end

%% Average over subjects
time_ = nan(2,2,20,81);
time_s = time_;
cc_ = time_;
cc_s = time_;
reps = 250;
disppercent(-1/length(rois));
for ri = 1:length(rois)
    for mi = 1:2
        for ii = 1:20
            tout = bootci(reps,@nanmean,squeeze(time(:,ri,mi,ii,:)));
            toutm = squeeze(mean(tout));
            time_(ri,mi,ii,:) = toutm;
            time_s(ri,mi,ii,:) = tout(2,:)-toutm;
            cout = bootci(reps,@nanmean,squeeze(cc(:,ri,mi,ii,:)));
            coutm = squeeze(mean(cout));
            cc_(ri,mi,ii,:) = coutm;
            cc_s(ri,mi,ii,:) = cout(2,:)-coutm;
        end
    end
    disppercent(ri/length(rois));
end
disppercent(inf);

%%
avg_decon.time = time_;
avg_decon.times = time_s;
avg_decon.cc = cc_;
avg_decon.ccs = cc_s;
avg_decon.conidx = decondata.V1.cc.conidxs;
avg_decon.cohidx = decondata.V1.cc.cohidxs;
avg_decon.timconidx = decondata.V1.time.conidxs;
avg_decon.timcohidx = decondata.V1.time.cohidxs;
avg_decon.timidxs = decondata.V1.time.timidxs;

save(fullfile(datafolder,'avg_decon.mat'),'avg_decon');