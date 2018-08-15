mopts = 1;
models = {'exp'};
% 'sigma','sigma,poisson',
bmodels = {'sigma,roi','sigma,roi,poisson'};%,'doublesigma','doublesigma,poisson'};
% bmodels = {'sigma,roi'};

% options list
aopts = zeros(10000,5);
count = 1;

% build options 

ropts = {[1 8],1:8};

sigmaopts = linspace(.01,.5,6);

for ai = 1:length(aSIDs)
    for mi = 1:length(mopts)
        for ropt = 1:length(ropts)
            for ni = 1:length(bmodels)
                aopts(count,:) = [ai mi ni ropt 1];
                count = count+1;
            end
        end
    end
end
aopts = aopts(1:(count-1),:);

% break into 12*10 size chunks
breaks = 1:240:size(aopts,1);
breaks(end) = size(aopts,1)+1;

if length(breaks)==1
    breaks(2) = breaks(1); breaks(1) = 1;
end
%%
nSIDs = [305 329 43 25 300 346 343 344 338 340 348];
bSIDs = [345 330 337 335 349 354 353 356 334 352]; % behavior only participants
aSIDs = [nSIDs bSIDs];
sids = {};
for ni = 1:length(aSIDs)
    sids{end+1} = sprintf('s%03.0f',aSIDs(ni));
end

%%
clear afits afits_
load(fullfile(datafolder,'avg_indiv_fits.mat'));
afits_ = afits;
afits = cell(1,length(aSIDs));
for ai = 1:length(aSIDs)
    afits{ai} = cell(length(unique(aopts(:,3))),length(unique(aopts(:,4))),length(unique(aopts(:,5))));
end
for ai = 1:size(afits_,1)
    copt = aopts(ai,:); 
    subj = copt(1);
    shape = copt(2);
    noise = copt(3);
    ropt = copt(4);
    sigma = copt(5);
    afits{subj}{noise,ropt} = afits_{ai};
end