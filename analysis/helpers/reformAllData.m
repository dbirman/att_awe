function allData = reformAllData( decondata, si)
%%
load(fullfile(datafolder,'avg_hrf.mat'));
pull = setdiff(1:11,si);
ahrf = squeeze(mean(hrfs(pull,4,:)));
% ahrf = mean(hrfs,1);
ahrf = ahrf./max(ahrf);
%%
ROIs = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
%%
allData = struct;
% re-organize data into a useful format:
% ROI - IDX - DECON
respOpts = {'resp_prf','resp_all','resp_2','resp_25'};
allData.cc.resp = zeros(8,length(decondata.V1.cc.conidxs),81);
allData.time.resp = zeros(8,length(decondata.V1.time.conidxs),81);
allData.cc.con = decondata.V1.cc.conidxs;
allData.cc.coh = decondata.V1.cc.cohidxs;
allData.cc.time = repmat(5,size(allData.cc.con));
allData.time.con = decondata.V1.time.conidxs;
allData.time.coh = decondata.V1.time.cohidxs;
allData.time.time = decondata.V1.time.timidxs;
allData.basecon = 0.25;
allData.basecoh = 0;

for ro = 1:4
    allData.cc.(respOpts{ro}) = zeros(8,length(decondata.V1.cc.conidxs),81);
    allData.time.(respOpts{ro}) = zeros(8,length(decondata.V1.time.conidxs),81);
    for ri = 1:8
        roi = ROIs{ri};
        decon = decondata.(roi);
        eval(sprintf('allData.cc.%s(ri,:,:) = decon.cc.%s;',respOpts{ro},respOpts{ro}));
        eval(sprintf('allData.time.%s(ri,:,:) = decon.time.%s;',respOpts{ro},respOpts{ro}));
    end
end

allData.ROIs = ROIs;
allData.hrf = ahrf;

