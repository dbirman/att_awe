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
    conidxs = logical((tdata(:,1)==2) .* (tdata(:,13)>conbins));
    cohidxs = logical((tdata(:,1)==1) .* (tdata(:,14)>cohbins));
    lapse = 1-mean([tdata(conidxs,12) ; tdata(cohidxs,12)]);
    ltrials(ai) = size([tdata(conidxs,12) ; tdata(cohidxs,12)],1);
    lapses(ai) = lapse;
end
    
%% Fit to behavior

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
                if strfind(bmodels{ni},'roi')
                    % roi models have sigma fixed, no need to do
                    % multiple
                    aopts(count,:) = [ai mi ni ropt 1];
                    count = count+1;
                else
                    for si = 1:length(sigmaopts)
                        aopts(count,:) = [ai mi ni ropt si];
                        count = count+1;
                    end
                end
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
        sigma = sigmaopts(copt(5));
        
        info = struct;
        info.sigma = sigma;
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

save(fullfile(datafolder,'avg_indiv_fits.mat'),'afits');

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

%% Individual fits

%% Average fits
plot_rightchoice_model;

%%
plot_indiv_rightchoice_model;

%% Catch trial analysis:
% (1) Attempt to fit catch trials with the existing readout model

%% re-organize afits
restructure_afits;

% save(fullfile(datafolder,'avg_indiv_fits.mat'),'afits');

%% Fit to behavior: individual fit
% We use the individual responses to see whether they better characterize
% the data than the 
clear afits
parfor ai = 1:length(nSIDs)
    adata = loadadata(sprintf('s%03.0f',nSIDs(ai)));
    
    fits = cell(length(mopts),2,8,8);
    for mi = 1:length(mopts)
        for ni = 1:length(bmodels)        
            info = struct;
            info.sigma = 1;
            info.model = bmodels{ni};
            info.rois = [1 8];
            info.lapse = lapses(ai);
            info.respcon = squeeze(respcon(ai,:,:));
            info.respcoh = squeeze(respcoh(ai,:,:));

            fits{mi,ni} = fitCCBehavControlModel_fmri(adata,info,1);
        end
        close all
    end
    afits{ai} = fits;
end
save(fullfile(datafolder,'avg_within_fits.mat'),'afits');
%% Permutation test
% Figure out the range of expected R^2 on random data (take the individual
% subject and randomly permute the responses, column 8)
reps = 500; % expect to take about 30 minutes per repetition

permlike = zeros(length(nSIDs),reps,2,2);

tic;
parfor rep = 1:reps
    for ai = 1:length(nSIDs)
        adata = loadadata(sprintf('s%03.0f',nSIDs(ai)));

        localr2 = zeros(reps,2,2);
            pdata = adata;
            pdata(:,8) = pdata(randperm(size(pdata,1)),8);
            for mi = 1:length(mopts)
                for ni = 1:length(bmodels)
                    info = struct;
                    info.sigma = 1;
                    info.model = bmodels{ni};
                    info.rois = [1 8];
                    info.lapse = lapses(ai);
                    info.respcon = squeeze(respcon(ai,:,:));
                    info.respcoh = squeeze(respcoh(ai,:,:));
                    %                     if clapse==0
                    %                         clapse = min(lapses(lapses>0));
                    %                     end
                    fit = fitCCBehavControlModel_fmri(pdata,info,1);
                    localr2(rep,mi,ni) = fit.likelihood;
                end
                close all
            end
    end

    permlike(ai,:,:,:) = localr2;
end
t = toc;
disp(sprintf('Took %01.2f seconds',t));

save(fullfile(datafolder,'permutation.mat','permlike'));

%% Computation permutation R^2 range

quantile(permr2(:),[0.05 0.95 0.99])
hist(permr2(:));


significance = [.05 .01 .001];
pval = [0.0201601085600963 0.0278872183918559 0.0326540667595963];

%% Compute all possible difference scores
pr2 = permr2(:);
diffs = zeros(length(pr2),length(pr2)-1);

for i = 1:length(pr2)
    diffs(i,:) = pr2(setdiff(1:length(pr2),i))-pr2(i);
end
quantile(diffs(:),[0.01 0.05 0.95 0.99])
hist(diffs(:));

%   Columns 1 through 3
% 
%        -0.0276701816421294       -0.0203526997937312        0.0203526997937312
% 
%   Column 4
% 
%         0.0276701816421294

%% Compare within-subject to across-subject fits
restructure_afits;
allfits = afits;
allfits = allfits(1:11);

load(fullfile(datafolder,'avg_within_fits.mat'));
for ai = 1:length(afits)
    temp = afits{ai};
    temp = temp(:,:,1,8);
    afits{ai} = temp;
end
withinfits = afits;

clear ar2 wr2
for ai = 1:11
    for mi = 1:2
        for ni = 1:2
            ar2(ai,mi,ni) = -sum(allfits{ai}{mi,ni}.cv.like);
            wr2(ai,mi,ni) = -sum(withinfits{ai}{mi,ni}.cv.like);
        end
    end
end

%% Use permutation test results to estimate whether there is an improvement within-subject?

%% Indiv
restructure_afits;

clear r2 sigmas cd
for ai = 1:length(aSIDs)
    for ni = 1:2
        for ropt = 1
            cm = afits{ai}{ni};
%             BIC(ai,ni) = cm.BIC;
%             sigmas(ai,ni,ropt) = cm.params.sigma;
            fr2(ai,ni,ropt) = -sum(cm.cv.like);
            fcd(ai,ni,ropt) = cm.cv.cd;
        end
    end
end

% average effects across subjects?

% sigmas(sigmas==1) = NaN;

%% Collect sigmas and indiv r2
restructure_afits;

clear r2 sigmas cd
for ai = 1:length(aSIDs)
    for ni = 1:2
        for ropt = 1:2
            cm = afits{ai}{ni,ropt};
            BIC(ai,ni) = cm.BIC;
            sigmas(ai,ni,ropt) = cm.params.sigma;
            r2(ai,ni,ropt) = -sum(cm.cv.like);
            cd(ai,ni,ropt) = cm.cv.cd;
            % compute 
            if ni ==1
                p = cm.cv.aprobs';
                r = cm.cv.resp;
                t = cm.adata(:,1);
                for ti = 1:2
                    sep_cd(ai,ti) = nanmean(p(logical((t==ti).*(r==1)))) - nanmean(1-p(logical((t==ti).*(r==0))));
                end
            end
        end
    end
end

sigmas(sigmas==1) = NaN;

% r2 = r2(:,:,1,1);
% cd = cd(:,:,1,1);

%% Report R^2 (exp)
add = squeeze(r2(:,1,:));
poi = squeeze(r2(:,2,:));
diffe = add-poi;

aci = bootci(1000,@mean,add);
pci = bootci(1000,@mean,poi);
dci = bootci(1000,@mean,diffe);

figure;
hist(diffe)

%% cd vs r2
add_cd = squeeze(cd(:,1,:));
poi_cd = squeeze(cd(:,2,:));

r2_diff = poi_roi-add_roi;
cd_diff = poi_cd-add_cd;
plot(100*cd_diff,r2_diff,'*k');
xlabel('\Delta CD');
ylabel('\Delta R^2');

%% Generate histogram of R^2 values

% histogram plot: DIFFERENCE FOR NAKA-RUSHTON MODEL OF ADDITIVE VS. POISSON
h = figure; hold on

vline(0,'--k');
[b,i] = hist(diffe,30);
d = diff(i); d = d(1);
for ii = 1:length(b)
    if b(ii)>0
        if abs(i(ii))>2
            rectangle('Position',[i(ii)-d/2 0 d b(ii)],'FaceColor',[0.75 0.75 0.75],'EdgeColor',[0 0 0]);
        else
            rectangle('Position',[i(ii)-d/2 0 d b(ii)],'FaceColor',[0.5 .5 .5],'EdgeColor',[0 0 0]);
        end
    end
end
% a2 = axis;
% mult = max(a2(1)/a(1));
vline(0,'--k');
axis([-125 175 0 3]);
set(gca,'XTick',[-100 -50 0 50 100 150],'XTickLabel',{'Evidence for Poisson','-50', '0', '50', '100', 'Evidence for additive'});
% set(gca,'XTick',[-5 5]*mult,'XTickLabel',{'Naka-Rushton','Linear'});
% set(gca,'XTick',round([-5 5]*mult));
set(gca,'YTick',[0 1 2 3]);
ylabel('Observers');
xlabel('\Delta Likelihood (Additive - Poisson)');
title('Model comparison: additive vs. Poisson noise');
% axis([a(1)*mult a(2)*mult a2(3) a2(4)]);
drawPublishAxis('figSize=[4.5,3.5]');
savepdf(h,fullfile(datafolder,'avg_models','add_poiss_hist.pdf'));

%% CD plot
h = figure; hold on

temp_cd = squeeze(cd(:,:,1));

plot(temp_cd(:,1),temp_cd(:,2),'o','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',5);
plot([min(temp_cd(:)) max(temp_cd(:))],[min(temp_cd(:)) max(temp_cd(:))],'--r');
xlabel('Additive model r^2');
ylabel('Poisson model r^2');
axis([.35 .5 .35 .5]);
set(gca,'XTick',[.35 .4 .45 .5],'YTick',[.35 .4 .45 .5]);
set(gca,'XTickLabel',{'35%','40%','45%','50%'},'YTickLabel',{'35%','40%','45%','50%'});
drawPublishAxis('figSize=[4.5,3.5]');
savepdf(h,fullfile(datafolder,'avg_models','add_poiss_r2.pdf'));

%% Contrast/coherence responses used to fit behavior
h = figure; hold on

cmap = brewermap(7,'PuOr');

x = 0:.001:1;
plot(x,squeeze(respcon_(1,:)),'Color',cmap(2,:));
plot(x,squeeze(respcoh_(8,:)),'Color',cmap(6,:));

l = legend({'V1 Contrast response function','MT Coherence response function'},'FontSize',7,'FontName','Helvetica');
set(l,'box','off');

set(gca,'XTick',[0 1],'XTickLabel',[0 100]);

xlabel('Stimulus strength (%)');
ylabel('\Delta signal (%)');

drawPublishAxis('figSize=[18,14]');

% savepdf(h,fullfile('~/proj/att_awe/talks/data_figures/v1mt.pdf'));




%% indiv_fits_2

% check to see whether the parameters are of similar magnitudes for V1 and
% MT
restructure_afits_2;

rois = {'V1','MT'};
cons = {'coh','con'};
for ni = 1:21
    for ri = 1:2
        for ci = 1:2
            betas(ni,ri,ci) = afits{ni}{1}.params.(sprintf('beta_control_%s_%sw',rois{ri},cons{ci}));
        end
    end
    r2_2(ni) = -sum(afits{ni}{1}.cv.like);
end

v1_con = betas(:,1,2);
mt_coh = betas(:,2,1);

bootci(10000,@mean,v1_con)
mean(ans)
bootci(10000,@mean,mt_coh)
mean(ans)

bootci(10000,@mean,1./[v1_con mt_coh])

%% I have no idea what all thi scode does

dat = 1./[v1_con;mt_coh];
bootci(10000,@mean,dat)

%% Fit ALL data

% aadata = [];
% for ai = 1:length(aSIDs)
%     adata = loadadata(sprintf('s%03.0f',aSIDs(ai)));
%     aadata = [aadata;adata];
% end
%
% jointfit = cell(length(mopts),length(bmodels));
% for mi = 1:length(mopts)
%     for rcon = 1
%         for rcoh = 8
%             for bi = 1:length(bmodels)
%                 clapse = mean(lapses);
%                 jointfit{mi,bi} = fitCCBehavControlModel_fmri(aadata,[],bmodels{bi},squeeze(respcon_(mi,rcon,:)),squeeze(respcoh_(mi,rcoh,:)),clapse);
%             end
%             close all
%         end
%     end
% end

%% Load
load(fullfile(datafolder,'avg_indiv_fits.mat'));

%% Sigmas
clear sigmas
for ai = 1:length(aSIDs)
    for bm = 1:2
        sigmas(ai,bm) = afits{ai}{1,bm,1,8}.params.sigma;
    end
end

%% Probs
clear muprobs
for ai = 1:length(aSIDs)
    for bm = 1:2
        muprobs(ai,bm) = (afits{ai}{1,bm,1,8}.muProb-0.5)*2;
    end
end

%% Indivi plot figures
plot_indiv_rightchoice_model;

%% Figure 7
plot_rightchoice_model;

%% ROI Parameter plot
restructure_afits;

sensitivity = zeros(2,length(aSIDs),8,2);
ROIs = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
cons = {'cohw','conw'};
for bi = 1:2
    for ai = 1:length(aSIDs)
        for ri = 1:8
            for ci = 1:2
                sensitivity(bi,ai,ri,ci) = afits{ai}{bi}.params.(sprintf('beta_control_%s_%s',ROIs{ri},cons{ci}));
            end
        end
    end
end
%%
models = {'exp'};
bmodels_text = {'additive','poisson'};

for bi = 1%:2
    h = figure; clf; hold on

    csensitivity = squeeze(sensitivity(bi,:,:,:));
    
    ci = bootci(1000,@mean,csensitivity);
    s_mean = squeeze(mean(ci));
    s_std = squeeze(ci(2,:,:,:))-s_mean;
        
    tx = [0.75*ones(1,8)];
    ty = [1.5*ones(1,8)];
    conrange = abs([min(s_mean(:,2)) max(s_mean(:,2))]);
    cohrange = abs([min(s_mean(:,1)) max(s_mean(:,1))]);
    for ri = 1:8
        orangeness = [241 163 64] * (s_mean(ri,2)+conrange(1))/sum(conrange);
        purpleness = [153 142 195] * (s_mean(ri,1)+cohrange(1))/sum(cohrange);
        color = orangeness + purpleness / 2;
        if any(color>255)
            color = color / max(color) * 255;
        end
        % plot the horizontal error bar (constant y): CONTRAST
        plot([-1 1]*s_std(ri,2) + s_mean(ri,2),repmat(s_mean(ri,1),1,2),'-','Color',color/255);
        % plot vertical
        plot(repmat(s_mean(ri,2),1,2),s_std(ri,1)*[-1 1] + s_mean(ri,1),'-','Color',color/255);
        plot(s_mean(ri,2),s_mean(ri,1),'o','MarkerFaceColor',color/255,'MarkerEdgeColor','white','MarkerSize',5);
        text(s_mean(ri,2)+tx(ri),s_mean(ri,1)+ty(ri),ROIs{ri},'Color',color/255);
    end
    axis equal
    axis([-13 25 -7 15]);
    axis equal
    set(gca,'XTick',[ -5 0 5 10 15]','YTick',[-5 0 5 10 15]);
    v = hline(0,'--'); set(v,'Color',[0.8 0.8 0.8]);
    v = vline(0,'--'); set(v,'Color',[0.8 0.8 0.8]);
    xlabel('Contrast weight (a.u.)');
    ylabel('Coherence weight (a.u.)');
%     title(sprintf('Weights under %s noise',bmodels_text{bi}));
    drawPublishAxis('figSize=[4.5,4.5]');
%     savepdf(h,fullfile('~/proj/att_awe/talks/data_figures',sprintf('avg_sensitivity_%s.pdf',models{mi})));
    savepdf(h,fullfile(datafolder,'avg_models',sprintf('avg_weights_%s.pdf',bmodels_text{bi})));
end

%% Fit gain parameter 

load(fullfile(datafolder,'avg_indiv_fits.mat'));

clear gfits
for ai = 1:length(aSIDs)
    adata = loadadata(sprintf('s%03.0f',aSIDs(ai)));
    
    rcon = 1; % V1
    rcoh = 8; % MT
    clapse = lapses(ai);
    %                     if clapse==0
    %                         clapse = min(lapses(lapses>0));
    %                     end
    
    % no cross-validation
    fit_g = fitCCBehavControlModel_fmri(adata,afits{ai}{1,1,1,8}.params,'sigma,gain',squeeze(respcon_(1,rcon,:)),squeeze(respcoh_(1,rcoh,:)),clapse,0);
    close all
    gfits{ai} = fit_g;
end

save(fullfile(datafolder,'avg_gain_fits.mat'),'gfits');

%% Pull gain

load(fullfile(datafolder,'avg_gain_fits.mat'));

for ai = 1:length(aSIDs)
    r2_gain(ai) = gfits{ai}.r2;
    gain(ai) = gfits{ai}.params.coh_gain;
end

load(fullfile(datafolder,'avg_indiv_fits.mat'));

for ai = 1:length(aSIDs)
    r2_orig(ai) = afits{ai}{1,1,1,8}.r2;
end

% remove subjects with no r2 improvement

dr2 = r2_gain-r2_orig;
idxs = dr2>.0155;

dci = bootci(10000,@mean,dr2(idxs));
dci = dci*100;

disp(sprintf('R^2 improvement with gain %01.3f 95%% CI [%01.3f %01.3f]',mean(dci),dci(1),dci(2)));

dg = bootci(10000,@mean,gain(idxs));

disp(sprintf('Multiplicative gain %01.3f 95%% CI [%01.3f %01.3f]',mean(dg),dg(1),dg(2)));
%% V1 vs MT vs other areas

parfor ai = 1:length(aSIDs)
    adata = loadadata(sprintf('s%03.0f',aSIDs(ai)));
    
    fits = cell(8,8);
    for mi = 1
        for rcon = 1:8
            for rcoh = 1:8
                for ni = 1
                    clapse = lapses(ai);
                    %                     if clapse==0
                    %                         clapse = min(lapses(lapses>0));
                    %                     end
                    fits{rcon,rcoh} = fitCCBehavControlModel_fmri(adata,[],'freeze',squeeze(respcon_(mi,rcon,:)),squeeze(respcoh_(mi,rcoh,:)),clapse,1);
                end
                close all
            end
        end
    end
    afits{ai} = fits;
end

save(fullfile(datafolder,'avg_indiv_fits_rois.mat'),'afits');

%% Pull CV 

clear r2
for ai = 1:length(aSIDs)
    for rcon = 1:8
        for rcoh = 1:8
            r2(ai,rcon,rcoh) = afits{ai}{rcon,rcoh}.r2;
        end
    end
end

%% Take contrast V1 row and compute error bars
for ri = 1
    r2_v1 = squeeze(r2(:,ri,:));
    r2_v1 = bootci(1000,@mean,r2_v1);

    h= figure; hold on
    errbar(1:8,mean(r2_v1),squeeze(r2_v1(2,:))-mean(r2_v1),'-k');
    plot(1:8,mean(r2_v1),'ok','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',5);
    a = axis;
    axis([1 8 a(3) a(4)]);
    drawPublishAxis('figSize=[8.9 8.9]');
end
%% Take coherence MT row and compute error bars
r2_v1 = squeeze(r2(:,:,5));
r2_v1 = bootci(1000,@mean,r2_v1);

h= figure; hold on
errbar(1:8,mean(r2_v1),squeeze(r2_v1(2,:))-mean(r2_v1),'-k');
plot(1:8,mean(r2_v1),'ok','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',5);
a = axis;
axis([1 8 a(3) a(4)]);
drawPublishAxis('figSize=[8.9 8.9]');
%% plot
h = figure; hold on
r2_ = squeeze(mean(r2));
imagesc(r2_)
colorbar
colormap('gray')
xlabel('Coherence ROI');
set(gca,'XTick',1:8,'XTickLabel',ROIs);
set(gca,'YTick',1:8,'YTickLabel',ROIs);
ylabel('Contrast ROI');
set(gca,'YDir','normal');