%% Load all data and average across shared ROIs
pfxs = {'l','r'};
rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};

allROIs = {};
for ri = 1:length(rois)
    for pi = 1:2
        allROIs{end+1} = sprintf('%s%s',pfxs{pi},rois{ri});
    end
end

% load each subject's relevant data and collapse into a giant monster
% array
allResp = zeros(length(nSIDs),16,32,81);
for ni = 1:length(nSIDs)
    fname = fullfile(datafolder,sprintf('s%04.0f_deconEffects_att.mat',nSIDs(ni)));
    load(fname);
    
    for ri = 1:length(allROIs)
        allResp(ni,ri,:,:) = decondata.(allROIs{ri}).cross.resp;
    end
end
    
% average across hemispheres
allResp_ = zeros(length(nSIDs),8,32,81);
for ri = 1:length(rois)
    roinums = cellfun(@(x) ~isempty(strfind(x,rois{ri})),allROIs,'UniformOutput',false);
    roinums = find([roinums{:}]);
    
    if length(roinums)>2 % V3 being mis-interpreted as V3a/V3b
        roinums = roinums(1:2);
    end
    allResp_(:,ri,:,:) = squeeze(mean(allResp(:,roinums,:,:),2));
end

% average across subjects

mResp = squeeze(bootci(1000,@nanmean,allResp_));
mResp_ = squeeze(mean(mResp));

%% Plot each set of conditions
conidx = decondata.lV1.cross.conidxs;
cohidx = decondata.lV1.cross.cohidxs;
taskidx = decondata.lV1.cross.taskidx;
stim = decondata.lV1.cross.stim;

contrasts = [0.325 0.4 0.55 0.85]-.25;
coherences = [0.15 0.3 0.45 0.6];

%% Prep save
avgdecon.allResp = allResp_;
avgdecon.mResp = mResp;
avgdecon.conidx = conidx;
avgdecon.cohidx = cohidx;
avgdecon.taskidx = taskidx;
avgdecon.stim = stim;
avgdecon.contrasts = contrasts;
avgdecon.coherences = coherences;

%% Remove delta==1
% conidx = conidx(logical(~deltaidx));
% taskidx = taskidx(logical(~deltaidx));
% constim = conStim(logical(~deltaidx));
% cohstim = cohStim(logical(~deltaidx));

% %% Plot V1
% for ri = 8
%     h = figure; hold on
%     for ci = 1:4
%         for mi = 1:4
%             idx = (ci-1)*4+mi;
%             subplot(4,4,idx);
%             plot(squeeze(mcon(ri,idx,:)));
%             axis([0 50 -1 1.5]);
%         end
%     end
% end


%% Save

fname = fullfile(datafolder,'avg_deconatt','avg_deconEffects_att_cross.mat');
save(fname,'avgdecon');