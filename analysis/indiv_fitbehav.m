%% Fit the behavioral data (for each of 21 subjects)
% IGNORE the catch trials
%
% Use the average response functions across all subjects
nSIDs = [305 329 43 25 300 346 343 344 338 340 348];
bSIDs = [345 330 337 335 349 354 353 356 334 352]; % behavior only participants
aSIDs = [nSIDs bSIDs];
sids = {};
for ni = 1:length(aSIDs)
    sids{end+1} = sprintf('s%03.0f',aSIDs(ni));
end
rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};

%% Linear + Naka Fit
% load all the linear fits, generate fit data, average across subjects,
% then run in fitCCBehavModel_sigma (to estimate sigma per subject), return
% BIC

load(fullfile(datafolder,'avg_hrffits.mat'));

x = 0:.001:1;
respcon = zeros(11,8,length(x));
respcoh = zeros(11,8,length(x));
derivcon = zeros(11,8,length(x));
derivcoh = zeros(11,8,length(x));

for si = 1:11
    fit = sfits{si}{1,4}; % 4 refers to the non-linear coherence model w/ resp-25, which is our standard model

    for ri = 1:8
        respcon(si,ri,:) = conModel(x,fit.roifit{ri}.params);
        dparams = fit.roifit{ri}.params; dparams.deriv = 1;
        derivcon(si,ri,:) = conModel(x,dparams);
        respcoh(si,ri,:) = cohModel(x,fit.roifit{ri}.params);
        derivcoh(si,ri,:) = cohModel(x,dparams);
    end
end

respcon_ = squeeze(mean(bootci(10000,@mean,respcon)));
respcoh_ = squeeze(mean(bootci(10000,@mean,respcoh)));
derivcon_ = squeeze(mean(bootci(10000,@mean,derivcon)));
derivcoh_ = squeeze(mean(bootci(10000,@mean,derivcoh)));

save(fullfile(datafolder,'fitbehav_resp.mat'),'respcon_','respcoh_','derivcon_','derivcoh_');

%% Or load

x = 0:.001:1;
load(fullfile(datafolder,'fitbehav_resp.mat'));

%% Check that responses look right
figure; hold on
for mi = 1
    plot(squeeze(mean(respcon_([1],:),1)),'-r');
    plot(squeeze(mean(respcoh_([8],:),1)));
end

%% Lapse rate calculation

mopts = 1;

afits = cell(1,length(aSIDs));
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
    cq(ai) = conbins;
    mq(ai) = cohbins;
    conidxs = logical((tdata(:,1)==2) .* (tdata(:,13)>conbins));
    cohidxs = logical((tdata(:,1)==1) .* (tdata(:,14)>cohbins));
    lapse = 1-mean([tdata(conidxs,12) ; tdata(cohidxs,12)]);
    ltrials(ai) = size([tdata(conidxs,12) ; tdata(cohidxs,12)],1);
    lapses(ai) = lapse;
end
    
%% Lapse rate range
mean(cq)
mean(mq)
%% Fit to behavior

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
                % roi models have sigma fixed, no need to do
                % multiple
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
    
%% fit all options
disppercent(-1/size(aopts,1));

% breaks = [breaks(1) breaks(end)];
afits = cell(size(aopts,1),1);
wfits = cell(size(aopts,1),1);
for ni = 1:(length(breaks)-1)
    bstart = breaks(ni);
    bend = breaks(ni+1)-1;
    
    parfor ii = bstart:bend
        copt = aopts(ii,:);
        
        subj = copt(1);
        adata = loadadata(sprintf('s%03.0f',aSIDs(subj)));
        
        shape = copt(2);
        noise = copt(3);
        ropt = ropts{copt(4)};
        
        info = struct;
        info.sigma = 1;
        info.model = bmodels{noise};
        info.rois = ropt;
        info.lapse = lapses(subj);
        info.respcon = respcon_;
        info.respcoh = respcoh_;
        
        afits{ii} = fitCCBehavControlModel_fmri(adata,info,1);
        disp(ii);
   end
    
    disppercent(bend/size(aopts,1));
end
disppercent(inf);

save(fullfile(datafolder,'avg_indiv_fits_fmincon.mat'),'afits');

%% Load all of the data, across all subjects -- compute the response mappings (how each contrast effect changed V1 response)
% and compute the difference scores (R-L)

% load columns
%     1       2         3        4      5      6     7      8       9        10      11        12
%   task - basecon - basecoh - conL - conR - cohL - cohR - resp - catch - pedcon - pedcoh - correct

%  final columns
%     1       2         3        4        5     6     7      8       9       10      11        12       13        14      15
%   subj   - task - basecon - basecoh - conL - conR - cohL - cohR - resp -   catch - pedcon - pedcoh - correct - prevresp - prevcorr

%  response columns
%    1      2         
%    v1L     v1R    ... 

% diff columns
%    1      2     ...    8
%    v1D    v2D    ...  mtD

%% Average fits
plot_rightchoice_model;

%%
plot_indiv_rightchoice_model;

%% Catch trial analysis:
% (1) Attempt to fit catch trials with the existing readout model

%% re-organize afits
afits = restructure_afits('avg_indiv_fits');

% save(fullfile(datafolder,'avg_indiv_fits.mat'),'afits');

%% Fit to behavior: individual fit
% We use the individual responses to see whether they better characterize
% the data than the 
clear afits
parfor ai = 1:length(nSIDs)
    adata = loadadata(sprintf('s%03.0f',nSIDs(ai)));
    
    fits = cell(length(mopts),2);
    for ni = 1:length(bmodels)     
        for ri = 1:length(ropts)
            info = struct;
            info.sigma = 1;
            info.model = bmodels{ni};
            info.rois = ropts{ri};
            info.lapse = lapses(ai);
            info.respcon = squeeze(respcon(ai,:,:));
            info.respcoh = squeeze(respcoh(ai,:,:));

            fits{ni,ri} = fitCCBehavControlModel_fmri(adata,info,1);
        end
    end
    close all
    wfits{ai} = fits;
end

disppercent(-1/size(aopts,1));

save(fullfile(datafolder,'avg_within_fits_fmincon.mat'),'wfits');

%% Permutation
permutation_cohcon;

%% Compare within-subject to across-subject fits
afits = restructure_afits('avg_indiv_fits_fmincon.mat');
allfits = afits;
allfits = allfits(1:11);

clear afits
load(fullfile(datafolder,'avg_within_fits_fmincon.mat'));

withinfits = wfits;

clear ar2 wr2 aaic waic
for ai = 1:11
    for pi = 1:2
        for ri = 1:2
            ar2(ai,pi,ri) = -sum(allfits{ai}{pi,ri}.cv.like);
            wr2(ai,pi,ri) = -sum(withinfits{ai}{pi,ri}.cv.like);
%             aaic(ai,pi,ri) = allfits{ai}{pi,ri}.AIC;
%             waic(ai,pi,ri) = withinfits{ai}{pi,ri}.AIC;

            aaic(ai,pi,ri) = 2*allfits{ai}{pi,ri}.numParams - 2 * ar2(ai,pi,ri);
            waic(ai,pi,ri) = 2*withinfits{ai}{pi,ri}.numParams - 2 * wr2(ai,pi,ri);
        end
    end
end

%%
a_inc = aaic-waic;
a_inc = squeeze(a_inc(:,1,:));

mu = mean(a_inc(:,2));
ci = bootci(10000,@mean,a_inc(:,2));
disp('We compared fitting the linking model on average physiological data with a fully within-subject model for the 11 subjects with matched data.');
disp( 'Fitting on average physiological data compared to within-subject resulted in a change in AIC of');
disp(sprintf( '%1.2f, 95%% CI [%1.2f, %1.2f]',mu, ci(1),ci(2)));

%% Other AIC comparisons

afits = restructure_afits('avg_indiv_fits_fmincon.mat');
for ai = 1:21
    for pi = 1:2
        for ri = 1:2
            ar2(ai,pi,ri) = -sum(afits{ai}{pi,ri}.cv.like);
        end
    end
end

%% Compare 8-area to 2-area for additive
diff_area = ar2(:,1,2)./ar2(:,1,1);

disp('Difference in AIC for 8-area model compared to 2-area model');
ci = bootci(10000,@mean,diff_area);
disp(sprintf('%1.2f, 95%% CI [%1.2f, %1.2f]',mean(diff_area),ci(1),ci(2)));

%% Compare poisson to additive
diff_area = exp(ar2(:,1,2)-ar2(:,2,2));

disp('Difference in AIC for additive vs. poisson, 8-area model');
ci = bootci(10000,@mean,diff_area);
disp(sprintf('%1.2f, 95%% CI [%1.2f, %1.2f]',mean(diff_area),ci(1),ci(2)));

%% Indiv
afits = restructure_afits('avg_indiv_fits_fmincon');

clear r2 sigmas cd
for ai = 1:length(aSIDs)
    for ni = 1:2
        for ropt = 1:2
            cm = afits{ai}{ni,ropt};
%             BIC(ai,ni) = cm.BIC;
%             sigmas(ai,ni,ropt) = cm.params.sigma;
            fr2(ai,ni,ropt) = -sum(cm.cv.like);
            fcd(ai,ni,ropt) = cm.cv.cd;
        end
    end
end

%% get the mean fits for paper
for add = 1
    for ropt = 1:2
        vals = squeeze(fcd(:,add,ropt))*100;
        ci = bootci(1000,@nanmean,vals);
        ci = ci;
        disp(sprintf('%i-area fit, average CD = %1.2f, 95%% CI [%1.2f, %1.2f]',length(ropts{ropt}),mean(vals),ci(1),ci(2)));
    end
end
% average effects across subjects?

% compare 8 roi model to 2 roi model
dr2 = fr2(:,1,2)-fr2(:,1,1);
ci = bootci(10000,@nanmean,dr2);

disp(sprintf('8-area to 2-area cvLL: %1.2f, 95%% CI [%1.2f, %1.2f]',mean(dr2),ci(1),ci(2)));

%%
% add/poiss CD comparison
dcd = fcd(:,1,2)-fcd(:,2,2);
dcd = dcd*100;
ci = bootci(10000,@nanmean,dcd);
disp(sprintf('Additive to poisson CD: %1.2f, 95%% CI [%1.2f, %1.2f]',mean(dcd),ci(1),ci(2)));

% add poiss likelihood
dlike = fr2(:,1,2)-fr2(:,2,2);
ci = bootci(10000,@nanmean,dlike);
disp(sprintf('Additive to poisson cvLL: %1.2f, 95%% CI [%1.2f, %1.2f]',mean(dlike),ci(1),ci(2)));

%% Generate horizontal bar plot of likelihood differences
h = figure; hold on

orig = 1:21;
[d,idx] = sort(dlike,'descend');

for i = 1:21
    if idx(i)<=11
        barh(i,d(i),'FaceColor',[0 0 0]);
    else
        barh(i,d(i),'FaceColor',[0.75 0.75 0.75]);
    end
end

% barh(1:11,d(1:11),'FaceColor',[0 0 0]);
% barh(12:21,d(12:21),'FaceColor',[0.75 0.75 0.75]);
axis([-300 300 0 22]);
drawPublishAxis('figSize=[6,4.5]');

savepdf(h,fullfile(datafolder,'avg_models','add_poiss_bar.pdf'));

%% Generate horizontal bar plot of CD differences
h = figure; hold on

cd_sort = dcd(idx);

for i = 1:21
    if idx(i)<=11
        barh(i,cd_sort(i),'FaceColor',[0 0 0]);
    else
        barh(i,cd_sort(i),'FaceColor',[0.75 0.75 0.75]);
    end
end

axis([-10 10 0 22]);
set(gca,'XTick',[-10:5:10]);

drawPublishAxis('figSize=[6,4.5]');

savepdf(h,fullfile(datafolder,'avg_models','add_poiss_cd.pdf'));


%% Stay switch model (just additive)

% 'sigma','sigma,poisson',
bmodels = {'sigma,roi,stayswitch'};%,'doublesigma','doublesigma,poisson'};
% bmodels = {'sigma,roi'};

% options list
aopts = zeros(10000,5);
count = 1;

% build options 

ropts = {1:8};

sigmaopts = linspace(.01,.5,6);

for ai = 1:length(aSIDs)
    for mi = 1:length(mopts)
        for ropt = 1:length(ropts)
            for ni = 1:length(bmodels)
                % roi models have sigma fixed, no need to do
                % multiple
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
    
%% fit all options
disppercent(-1/size(aopts,1));

% breaks = [breaks(1) breaks(end)];
afits = cell(size(aopts,1),1);
wfits = cell(size(aopts,1),1);
for ni = 1:(length(breaks)-1)
    bstart = breaks(ni);
    bend = breaks(ni+1)-1;
    
    parfor ii = bstart:bend
        copt = aopts(ii,:);
        
        subj = copt(1);
        adata = loadadata(sprintf('s%03.0f',aSIDs(subj)));
        
        shape = copt(2);
        noise = copt(3);
        ropt = ropts{copt(4)};
        
        info = struct;
        info.sigma = 1;
        info.model = bmodels{noise};
        info.rois = ropt;
        info.lapse = lapses(subj);
        info.respcon = respcon_;
        info.respcoh = respcoh_;
        
        afits{ii} = fitCCBehavControlModel_fmri(adata,info,1);
        disp(ii);
   end
    
    disppercent(bend/size(aopts,1));
end
disppercent(inf);

save(fullfile(datafolder,'avg_switch_fits.mat'),'afits');

%% Compare side bias to afits

load(fullfile(datafolder,'avg_switch_fits.mat'));

for ai = 1:length(afits)
    switchfits{ai} = afits{ai};
    like(ai,2) = -sum(switchfits{ai}.cv.like);
end

afits = restructure_afits('avg_indiv_fits_fmincon.mat');

for ai = 1:length(afits)
    indivfits{ai} = afits{ai}{1,2};
    like(ai,1) = -sum(indivfits{ai}.cv.like);
end

diff = like(:,1)-like(:,2);

mean(diff)
bootci(1000,@mean,diff)
