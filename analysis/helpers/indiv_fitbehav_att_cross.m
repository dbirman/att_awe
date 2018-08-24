%% Fit the behavioral data (for each of 21 subjects)
% Use the average response functions across all subjects
nSIDs = [305 329 43 25 300 346 343 344 338 340 348];
bSIDs = [345 330 337 335 349 354 353 356 334 352]; % behavior only participants
aSIDs = [nSIDs bSIDs];
sids = {};
for ni = 1:length(aSIDs)
    sids{end+1} = sprintf('s%03.0f',aSIDs(ni));
end
rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
%% Population responses
% load all the linear fits, generate fit data, average across subjects,
% then run in fitCCBehavModel_sigma (to estimate sigma per subject), return
% BIC

load(fullfile(datafolder,'avg_att_cross_fits.mat'));

x = 0:.001:1;
respcon = zeros(11,8,2,length(x));
respcoh = zeros(11,8,2,length(x));

for si = 1:10
    fit = attfits{si}{1}; % 4 refers to resp-25, which is our standard model

    for ri = 1:8
        respcon(si,ri,1,:) = fit.roifit{ri}.conresp_coh;
        respcon(si,ri,2,:) = fit.roifit{ri}.conresp_con;
        respcoh(si,ri,1,:) = fit.roifit{ri}.cohresp_coh;
        respcoh(si,ri,2,:) = fit.roifit{ri}.cohresp_con;
    end
end

respcon_ = squeeze(mean(bootci(10000,@mean,respcon)));
respcoh_ = squeeze(mean(bootci(10000,@mean,respcoh)));

%% Linear functions
respcon_ = respcon_ - repmat(respcon_(:,:,1),1,1,1001);
respcoh_ = respcoh_ - repmat(respcoh_(:,:,1),1,1,1001);
% get linear fits
clear respcon_l respcoh_l
for ri = 1:8
    for ci = 1:2
        respcon_l(ri,ci,:) = x*(x'\squeeze(respcon_(ri,ci,:)));
        respcoh_l(ri,ci,:) = x*(x'\squeeze(respcoh_(ri,ci,:)));
    end
end

%% Check that responses look right
figure; hold on
ro = [1 8];
di = 2;
cons = {'attend coherence','attend contrast'};
cmap = brewermap(7,'PuOr');
for rii = 1:2
    subplot(2,1,rii); hold on
    ri = ro(rii);
    title(sprintf('%s: %s',rois{ri},cons{di}));
%     plot(squeeze(mean(respcon_([1 2 3 4],1,:),1)),'--r');
%     plot(squeeze(mean(respcoh_([5 8],1,:),1)),'--');
    plot(x,squeeze(mean(respcon_l(ri,di,:),1)),'-','Color',cmap(2,:));
    plot(x,squeeze(mean(respcoh_l(ri,di,:),1)),'-','Color',cmap(6,:));
    axis([0 1 0 1.9]);
end

%% Across ROIs
h = figure; hold on
ro = [1 2 3 6 7 8];
for rii = 1:length(ro)
    ri = ro(rii);
    subplot(2,1,1); hold on
    rc = squeeze(respcon_l(ri,1,:));
    plot(x,rc-rc(1),'Color',cmap(2,:));
    subplot(2,1,2); hold on
    rm = squeeze(respcoh_l(ri,1,:));
    plot(x,rm-rm(1),'Color',cmap(6,:));
end

h = figure; hold on
ro = [1 2 3 6 7 8];
for rii = 1:length(ro)
    ri = ro(rii);
    subplot(6,1,rii); hold on
%     subplot(2,1,1); hold on
    rc = squeeze(respcon_l(ri,1,:));
    plot(x,rc-rc(1),'Color',cmap(2,:));
%     subplot(2,1,2); hold on
    rm = squeeze(respcoh_l(ri,1,:));
    plot(x,rm-rm(1),'Color',cmap(6,:));
end

%% Lapse rate calculation

mopts = 1;

lapses = zeros(1,length(aSIDs));
ltrials = lapses;
for ai = 1:length(aSIDs)
    adata = loadadata(sprintf('s%03.0f',aSIDs(ai)));
    % fit the lapse rate parameter
    tdata = adata(adata(:,9)==-1,:);
    tdata(:,13) = abs(tdata(:,5)-tdata(:,4)); % contrast diff
    tdata(:,14) = abs(tdata(:,7)-tdata(:,6)); % coherence diff
    conbins = quantile(tdata(:,13),.95);
    cohbins = quantile(tdata(:,14),.95);
    conidxs = logical((tdata(:,1)==2) .* (tdata(:,13)>conbins));
    cohidxs = logical((tdata(:,1)==1) .* (tdata(:,14)>cohbins));
    lapse = 1-mean([tdata(conidxs,12) ; tdata(cohidxs,12)]);
    ltrials(ai) = size([tdata(conidxs,12) ; tdata(cohidxs,12)],1);
    lapses(ai) = lapse;
end
    
%% Fit to behavior

models = {'exp'};
% 'sigma','sigma,poisson',
bmodels = {'sigma,roi,att','sigma,roi,att,onebeta'};
% bmodels = {'sigma,roi'};

% options list
attopts = zeros(10000,5);
count = 1;
mopts = 1;

% build options 

ropts = {[1 8],[1:8]};
% rconopts = {1, [1 2 3 4], 1:8};
% rcohopts = {8, [5 8], 1:8};

sigmaopts = linspace(.01,.5,6);

for ai = 1:length(aSIDs)
    for mi = 1:length(mopts)
        for ropt = 1:length(ropts)
            for ni = 1:length(bmodels)
                if strfind(bmodels{ni},'roi')
                    % roi models have sigma fixed, no need to do
                    % multiple
                    attopts(count,:) = [ai mi ni ropt 1];
                    count = count+1;
                else
                    for si = 1:length(sigmaopts)
                        attopts(count,:) = [ai mi ni ropt si];
                        count = count+1;
                    end
                end
            end
        end
    end
end
attopts = attopts(1:(count-1),:);

% break into 12*10 size chunks
breaks = 1:240:size(attopts,1);
breaks(end) = size(attopts,1)+1;

if length(breaks)==1
    breaks(2) = breaks(1); breaks(1) = 1;
end
    
%% fit all options
% disppercent(-1/size(attopts,1));

% breaks = [breaks(1) breaks(end)];
attfits = cell(size(attopts,1),1);
wfits = cell(size(attopts,1),1);
for ni = 1:(length(breaks)-1)
    bstart = breaks(ni);
    bend = breaks(ni+1)-1;
    
    parfor ii = bstart:bend
        
        copt = attopts(ii,:);
        
        subj = copt(1);
        adata = loadadata(sprintf('s%03.0f',aSIDs(subj)));
        
        shape = copt(2);
        noise = copt(3);
        ropt = ropts{copt(4)};
        sigma = sigmaopts(copt(5));
        
        info = struct;
        info.sigma = sigma;
        info.model = bmodels{noise};
        info.rois = ropt;
        info.lapse = lapses(subj);
        info.respcon = respcon_l;
        info.respcoh = respcoh_l;
        
        attfits{ii} = fitCCBehavControlModel_fmri(adata,info,1);
        disp(sprintf('Done with %i',subj));
   end
    
%     disppercent(bend/size(attopts,1));
end
% disppercent(inf);

save(fullfile(datafolder,'avg_indiv_fits_att_cross_linear.mat'),'attfits');

% breaks = [breaks(1) breaks(end)];
attfits = cell(size(attopts,1),1);
wfits = cell(size(attopts,1),1);
for ni = 1:(length(breaks)-1)
    bstart = breaks(ni);
    bend = breaks(ni+1)-1;
    
    parfor ii = bstart:bend
        
        copt = attopts(ii,:);
        
        subj = copt(1);
        adata = loadadata(sprintf('s%03.0f',aSIDs(subj)));
        
        shape = copt(2);
        noise = copt(3);
        ropt = ropts{copt(4)};
        sigma = sigmaopts(copt(5));
        
        info = struct;
        info.sigma = sigma;
        info.model = bmodels{noise};
        info.rois = ropt;
        info.lapse = lapses(subj);
        info.respcon = respcon_;
        info.respcoh = respcoh_;
        
        attfits{ii} = fitCCBehavControlModel_fmri(adata,info,1);
        disp(sprintf('Done with %i',subj));
   end
    
%     disppercent(bend/size(attopts,1));
end
% disppercent(inf);

save(fullfile(datafolder,'avg_indiv_fits_att_cross.mat'),'attfits');

%% Restructure attfits
load(fullfile(datafolder,'avg_indiv_fits_att_cross_linear.mat'));
attfits_ = attfits; 
count = 1;
for ai = 1:21
    for mi = 1
        for ropt = 1:2
            for ni = 1:2
                attfits_l{ai,ropt,ni} = attfits_{count};
                count = count+1;
            end
        end
    end
end
disp(count-1)
disp(length(attfits_))
load(fullfile(datafolder,'avg_indiv_fits_att_cross.mat'));
attfits_ = attfits; clear attfits
count = 1;
for ai = 1:21
    for mi = 1
        for ropt = 1:2
            for ni = 1:2
                attfits{ai,ropt,ni} = attfits_{count};
                count = count+1;
            end
        end
    end
end
disp(count-1)
disp(length(attfits_))
%% Compare afits and attfits
afits = restructure_afits('avg_indiv_fits.mat');

clear like cd
for i = 1:21
    % first index is noise add/poi, second index is ropts -- we ignore all
    % poisson models
    for ropt = 1:2
        like(i,ropt) = -sum(afits{i}{1,ropt}.cv.like);
        cd(i,ropt) = afits{i}{1,ropt}.cv.cd;
    end
end

% now pull out like and cd for the cross 

for i = 1:21
    for ropt = 1:2
        for beta = 1:2
            flip = [2 1]; % flip one beta so that 1 is one and 2 is two
            like_l(i,ropt,beta) = -sum(attfits_l{i,ropt,beta}.cv.like);
            cd_l(i,ropt,beta) = attfits_l{i,ropt,beta}.cv.cd;
            
            like_a(i,ropt,beta) = -sum(attfits{i,ropt,beta}.cv.like);
            cd_a(i,ropt,beta) = attfits{i,ropt,beta}.cv.cd;
        end
    end
end

%% Text comparisons for models:
diff_like = like_a(:,:,1)-like_a(:,:,2);

for i = 1:2
    length(ropts{i})
    mean(diff_like(:,i))
    bootci(1000,@mean,diff_like(:,i))
end

diff_cd = cd_a(:,:,1)-cd_a(:,:,2);

for i = 1:2
    length(ropts{i})
    mean(diff_cd(:,i))
    bootci(1000,@mean,diff_cd(:,i))
end

%% Plot model comparison
h = figure; hold on

diff_like = like_a(:,:,1)-like_a(:,:,2);

subplot(121);
barh(diff_like(:,2),'FaceColor',[0.75 0.75 0.75]);

axis([0 300 0 22]);
set(gca,'XTick',0:100:300);
drawPublishAxis('figSize=[6,4.5]');

subplot(122);
barh(diff_like(:,1),'FaceColor',[0.75 0.75 0.75]);

axis([0 300 0 22]);
set(gca,'XTick',[ 0 100 200 300]);
drawPublishAxis('figSize=[4.25,4.25]');

savepdf(h,fullfile(datafolder,'avg_models','flex_inflex_bar.pdf'));

%%
h = figure; hold on

diff_a = cd_a(:,:,1)-cd_a(:,:,2);

subplot(121);
barh(diff_a(:,2),'FaceColor',[0.75 0.75 0.75]);

axis([0 .3 0 22]);
set(gca,'XTick',0:.05:.25);
drawPublishAxis('figSize=[6,4.5]');

subplot(122);
barh(diff_a(:,1),'FaceColor',[0.75 0.75 0.75]);

axis([0 .3 0 22]);
set(gca,'XTick',0:.05:.25);
drawPublishAxis('figSize=[4.25,4.25]');

savepdf(h,fullfile(datafolder,'avg_models','flex_inflex_cd_bar.pdf'));

%% Plot models
roiOpts = {{'V1','MT'},{'V1','V2','V3','V4','V3a','V3b','V7','MT'}};

for ri = 1:2
    ropt = ropts{ri};
    rois = roiOpts{ri};
    
    for beta = 1:2
        if beta==2
            plot_rightchoice_model_att_onebeta(attfits(:,ri,beta),respcon_(ropt,:,:),respcoh_(ropt,:,:),aSIDs,bmodels(beta),rois);
        else
            plot_rightchoice_model_att(attfits(:,ri,beta),respcon_(ropt,:,:),respcoh_(ropt,:,:),aSIDs,bmodels(beta),rois);
        end
        
    end 
end
