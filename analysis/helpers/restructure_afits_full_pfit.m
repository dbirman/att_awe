%%
nSIDs = [305 329 43 25 300 346 343 344 338 340 348];
bSIDs = [345 330 337 335 349 354 353 356 334 352]; % behavior only participants
aSIDs = [nSIDs bSIDs];
sids = {};
for ni = 1:length(aSIDs)
    sids{end+1} = sprintf('s%03.0f',aSIDs(ni));
end
rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};

%%
clear afits afits_ ffits
load(fullfile(datafolder,'avg_full_prefits_8.mat'));
ffits_ = afits;
ffits = cell(1,length(aSIDs));
for ai = 1:length(aSIDs)
    ffits{ai} = cell(length(unique(aopts(:,3))),length(unique(aopts(:,4))),length(unique(aopts(:,5))));
end
for ai = 1:size(ffits_,1)
    copt = aopts(ai,:); 
    subj = copt(1);
    shape = copt(2);
    noise = copt(3);
    ropt = copt(4);
    sigma = copt(5);
    ffits{subj}{noise,ropt} = ffits_{ai};
end

return
%% drop si model with lower likelihood

% load(fullfile(datafolder,'avg_indiv_fits.mat'));
afits_ = afits;
numSigmas = size(afits_{1},1);
numModels = length(bmodels);

afits = cell(1,length(aSIDs)); bi = []; sigmas = []; sbs = []; r2s = []; cds = [];
for ai = 1:length(aSIDs)
    for ni = 1:numModels
        for coni = 1:length(rconopts)
            for cohi = 1:length(rcohopts)
                if strfind(bmodels{ni},'roi')
                    afits{ai}{ni,coni,cohi} = afits_{ai}{1,ni,coni,cohi};
                else
                    for si = 1:numSigmas
                        sigmas(end+1) = afits_{ai}{si,ni,coni,cohi}.params.sigma;
                        r2s(end+1) = -sum(afits_{ai}{si,ni,coni,cohi}.cv.like);
                        cds(end+1) = -sum(afits_{ai}{si,ni,coni,cohi}.cv.cd);
                    end
                    bestidx = 1;
                    likelihood = -sum(afits_{ai}{1,ni,coni,cohi}.cv.like);
                    for sigmaopt = 2:numSigmas
                        if -sum(afits_{ai}{sigmaopt,ni,coni,cohi}.cv.like)>likelihood
                            bestidx = sigmaopt;
                        end
                    end
                    afits{ai}{ni,coni,cohi} = afits_{ai}{bestidx,ni,coni,cohi};
                    bi(end+1) = bestidx;
                    sbs(end+1) = afits_{ai}{bestidx,ni,coni,cohi}.params.sigma;
                end
            end
        end
    end
end