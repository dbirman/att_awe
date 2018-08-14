
% load(fullfile(datafolder,'sigmas.mat'));
r2s = zeros(length(nSIDs),10,2);
allcon = nan(length(nSIDs),101);
allcoh = nan(length(nSIDs),101);
allconfmri = zeros(length(nSIDs),8,101);
allcohfmri = zeros(length(nSIDs),8,101);
for ni = 1:length(nSIDs)
    sid = nSIDs(ni);
    load(fullfile(datafolder,sprintf('s%03.0f_data.mat',sid)));
    load(fullfile(datafolder,sprintf('s%04.0f_fitroi.mat',sid)));
    x = 0:.01:1;
    con = conModel(x,data.fits{1}.params);
    coh = cohModel(x,data.fits{1}.params);
%     if ~isnan(sigmas.(sprintf('s%04.0f',sid))) && sigmas.(sprintf('s%04.0f',sid))<1 && sigmas.(sprintf('s%04.0f',sid))>0
%         allcon(ni,:) = con*sigmas.(sprintf('s%04.0f',sid));
%         allcoh(ni,:) = coh*sigmas.(sprintf('s%04.0f',sid));
%     end
    pos = 1:2:16;
    for ri = 1:8
        lroi = fitroi.roiparams{pos(ri)};
        rroi = fitroi.roiparams{pos(ri)+1};
        confmri = mean([conModel(x,lroi);conModel(x,rroi)]);
        cohfmri = mean([cohModel(x,lroi);cohModel(x,rroi)]);
        allconfmri(ni,ri,:) = confmri-confmri(1);
        allcohfmri(ni,ri,:) = cohfmri-cohfmri(1);
        r2s(ni,ri,1) = corr(coh',cohfmri');
        r2s(ni,ri,2) = corr(con',confmri');
    end
end

fmri.con = allconfmri;
fmri.coh = allcohfmri;

save(fullfile(datafolder,'avg_fmriresponse.mat'),'fmri');
