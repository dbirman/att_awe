function allData = reformAllData_cv( decondata, group)
%%
ROIs = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
%%
allData = struct;
% re-organize data into a useful format:
% ROI - IDX - DECON
respOpts = {'resp_prf','resp_all','resp_2','resp_25'};
allData.cc.resp = zeros(8,length(decondata.V1.(group).cc.conidxs),81);
allData.time.resp = zeros(8,length(decondata.V1.(group).time.conidxs),81);
allData.cc.con = decondata.V1.(group).cc.conidxs;
allData.cc.coh = decondata.V1.(group).cc.cohidxs;
allData.cc.time = repmat(5,size(allData.cc.con));
allData.time.con = decondata.V1.(group).time.conidxs;
allData.time.coh = decondata.V1.(group).time.cohidxs;
allData.time.time = decondata.V1.(group).time.timidxs;
allData.basecon = 0.25;
allData.basecoh = 0;

for ro = 1:4
    allData.cc.(respOpts{ro}) = zeros(8,length(decondata.V1.(group).cc.conidxs),81);
    allData.time.(respOpts{ro}) = zeros(8,length(decondata.V1.(group).time.conidxs),81);
    for ri = 1:8
        roi = ROIs{ri};
        decon = decondata.(roi);
        eval(sprintf('allData.(group).cc.%s(ri,:,:) = decon.(group).cc.%s;',respOpts{ro},respOpts{ro}));
        eval(sprintf('allData.(group).time.%s(ri,:,:) = decon.(group).time.%s;',respOpts{ro},respOpts{ro}));
    end
end

allData.ROIs = ROIs;

