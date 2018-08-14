function indiv_analysis(subj,mode)

%% define folders
datafiles = dir(fullfile(datafolder,sprintf('%s*_data.mat',subj)));

%% load data
datas = {};
for fi = 1:length(datafiles)
    load(fullfile(datafolder,datafiles(fi).name));
    % we only want to keep non-fixation task, so check >0
    if all(data.design(:,9)>0)
        datas{end+1}=data;
    end
end

if strfind(mode,'decon') && ~isempty(datas)
    %% Deconvolve the pedestals and increment conditions
    pfxs = {'l','r'};
    rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
    for pi = 1:length(pfxs)
        for ri = 1:length(rois)
            deconvolveEffects_att(datas,sprintf('%s%s',pfxs{pi},rois{ri}),subj);
        end
    end
end