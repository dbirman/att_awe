%% This fits the full model INCLUDING the catch trials

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

%%
restructure_afits; % load the prefitted models
global origfits;
origfits = afits;

%% Linear + Naka Fit
% load all the linear fits, generate fit data, average across subjects,
% then run in fitCCBehavModel_sigma (to estimate sigma per subject), return
% BIC

load(fullfile(datafolder,'avg_hrffits.mat'));

x = 0:.001:2;
respcon = zeros(11,8,length(x));
respcoh = zeros(11,8,length(x));

for si = 1:11
    fit = sfits{si}{4}; % 4 refers to resp-25, which is our standard model

    for ri = 1:8
        respcon(si,ri,:) = conModel(x,fit.roifit{ri}.params);
        respcoh(si,ri,:) = cohModel(x,fit.roifit{ri}.params);
    end
end

respcon_ = squeeze(mean(bootci(10000,@mean,respcon)));
respcoh_ = squeeze(mean(bootci(10000,@mean,respcoh)));

%% Check that responses look right
figure; hold on
for mi = 1
    plot(squeeze(mean(respcon_([1 2 3 4],:),1)),'-r');
    plot(squeeze(mean(respcoh_([5 8],:),1)));
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
bmodels = {'sigma,roi,pfit'};%,'doublesigma','doublesigma,poisson'};
% bmodels = {'sigma,roi'};

% options list
aopts = zeros(10000,5);
count = 1;

% build options 

ropts = {[1:8]};
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
        info.subj = subj;
        info.sigma = sigma;
        info.model = bmodels{noise};
        info.rois = ropt;
        info.lapse = lapses(subj);
        info.respcon = respcon_;
        info.respcoh = respcoh_;
        
%         temps = cell{1,4};
%         parfor iii=1:4
%             
%             tinfo = info;
%             tinfo.model = bmodels{iii};
%             temps{iii} = fitCCBehavControlModel_fmri(adata,tinfo,1);
%         end
        afits{ii} = fitCCBehavFullModel_fmri(adata,info,1);
        disp(sprintf('Finished %i',ii));
   end
    
    disppercent(bend/size(aopts,1));
end
disppercent(inf);

save(fullfile(datafolder,'avg_full_prefits_8.mat'),'afits');
% save(fullfile(datafolder,'avg_within_fits.mat'),'wfits');
%     save(fullfile(datafolder,sprintf('avg_indiv_fits_%02.0f.mat',100*sigmaopts(si))),'afits');
%     disp('************************************');
%     disppercent(si/length(sigmaopts));
%     disp('************************************');
% end
% disppercent(inf);
%%
load(fullfile(datafolder,'avg_full_prefits_8.mat'));

% Pull the CD values for the catch trials and regular trials separately
for ni = 1:21
    fit = afits{ni};
    i_cont = fit.adata(:,1)>0;
    i_catc = fit.adata(:,1)<0;
    r = fit.cv.resp;
    p = fit.cv.aprobs;
    cd_co(ni) = nanmean(p(logical(i_cont.*(r==1)))) - nanmean(1-p(logical(i_cont.*(r==0))));
    cd_ca(ni) = nanmean(p(logical(i_catc.*(r==1)))) - nanmean(1-p(logical(i_catc.*(r==0))));
end
ci_co = bootci(10000,@mean,cd_co)*100;
disp(sprintf('%02.2f%%, (95%% CI [%02.2f %02.2f])',mean(ci_co),ci_co(1),ci_co(2)));
ci_ca = bootci(10000,@mean,cd_ca)*100;
disp(sprintf('%02.2f%%, (95%% CI [%02.2f %02.2f])',mean(ci_ca),ci_ca(1),ci_ca(2)));

%% Get relative weights to original model
load(fullfile(datafolder,'avg_full_fits_8.mat'));

attcond = {'control','unatt'};
cons = {'cohw','conw'};
for ai = 1:21
    fit = afits{ai};
    for ri = 1:8
        for atti = 1:2
            for ci = 1:2
                w(ai,ri,atti,ci) = fit.params.(sprintf('beta_%s_%s_%s',attcond{atti},rois{ri},cons{ci}));
            end
        end
    end
    
end

dw = w(:,:,2,:) - w(:,:,1,:);
w_ = squeeze(bootci(1000,@mean,dw));
w = squeeze(mean(w_));

h = figure; hold on

hline(0,'--k');
vline(0,'--k');

for ri = 1:8
    plot([w(ri,2) w(ri,2)],w_(:,ri,1),'-k');
    plot(w_(:,ri,2),[w(ri,1) w(ri,1)],'-k');
end
plot(w(:,2),w(:,1),'o','MarkerFaceColor','k','MarkerEdgeColor','w');
for ri = 1:8
    text(w(ri,2),w(ri,1),rois{ri});
end
hline(0,'--k');
vline(0,'--k');
set(gca,'XTick',-20:5:15,'YTick',-10:5:5);
title('Change in model weight');
xlabel('Discriminating contrast');
ylabel('Discriminating coherence');

drawPublishAxis('figSize=[4.5,4.5]');

savepdf(h,fullfile(datafolder,'avg_models','catch_deltasense.pdf'));

%% Individual fits
plot_indiv_rightchoice_model;

%% Average fits
plot_rightchoice_model_full;
plot_rightchoice_model_full_pfit;

%% re-organize afits
restructure_afits_full;


%% Use permutation test results to estimate whether there is an improvement within-subject?

%% Full
restructure_afits_full;

clear r2 sigmas cd
for ai = 1:length(aSIDs)
    for ni = 1
        for ropt = 1
            cm = ffits{ai}{ni};
%             BIC(ai,ni) = cm.BIC;
%             sigmas(ai,ni,ropt) = cm.params.sigma;
            fr2(ai,ni,ropt) = -sum(cm.cv.like);
            fcd(ai,ni,ropt) = cm.cv.cd;
        end
    end
end

clear b_c b_u
rois = {'V1','MT'};
% pull the parameters
% use rois
cstr = {'cohw','conw'};

for ai = 1:length(aSIDs)
    for ri = 1:2
        for ci = 1:2
            b_c(ai,ri,ci) = ffits{ai}{1}.params.(sprintf('beta_control_%s_%s',rois{ri},cstr{ci}));
            b_u(ai,ri,ci) = ffits{ai}{1}.params.(sprintf('beta_unatt_%s_%s',rois{ri},cstr{ci}));
        end
    end
end

usig = ones(size(b_u,1),1) ./ b_u(:,1,1);
csig = ones(size(b_c,1),1) ./ b_c(:,1,2);
% b_u = b_u ./ repmat(b_u(:,1,1),1,2,2);
% b_c = b_c ./ repmat(b_c(:,1,2),1,2,2);

%% Compute the contrast response from different conditions
% feature | catch | response
x = 0:.001:2;
response = zeros(length(aSIDs),2,2,2001);
rat_con = zeros(length(aSIDs),1);
rat_coh = zeros(length(aSIDs),1);

% catchstr = {'control','unatt'};
roiopts = [1 8];
flip = [2 1];
for ai = 1:length(aSIDs)
    betas_contrast_control = squeeze(b_c(ai,:,2));
    betas_contrast_catch = squeeze(b_u(ai,:,2));
    % contrast
    response(ai,2,1,:) = betas_contrast_control * respcon_(roiopts,:);
    response(ai,2,2,:) = betas_contrast_catch * respcon_(roiopts,:);
    
    rat_con(ai) = squeeze(response(ai,2,2,:))\squeeze(response(ai,2,1,:));
    
    betas_motion_control = squeeze(b_c(ai,:,1));
    betas_motion_catch = squeeze(b_u(ai,:,1));
    % contrast
    response(ai,1,1,:) = betas_motion_control * respcoh_(roiopts,:);
    response(ai,1,2,:) = betas_motion_catch * respcoh_(roiopts,:);
    
    rat_coh(ai) = squeeze(response(ai,1,2,:))\squeeze(response(ai,1,1,:));
end

figure; hist(rat_con);
bootci(1000,@mean,rat_con)
figure; hist(rat_coh);
bootci(1000,@mean,rat_coh)


%% Permute group ratios
for perm = 1:10000
    response_ = squeeze(mean(response(randi(21,1,21),:,:,:)));

    avg_con(perm) = squeeze(response_(2,2,:))\squeeze(response_(2,1,:));
    avg_coh(perm) = squeeze(response_(1,2,:))\squeeze(response_(1,1,:));
end

ci_con = quantile(avg_con,[.025 .5 .975]);
ci_coh = quantile(avg_coh,[.025 .5 .975]);

%%
h = figure; hold on

cmap = brewermap(11,'PuOr');

response_ = squeeze(mean(response));

subplot(211); hold on
plot(x,squeeze(response_(2,1,:)),'-','Color',cmap(3,:),'LineWidth',2);
plot(x,squeeze(response_(2,2,:)),'-','Color',cmap(5,:),'LineWidth',1);
% 
legend({'Contrast - control','Contrast - catch'},'FontSize',7,'FontName','Helvetica');
% 
l = legend(gca,'boxoff');
set(l,'Color','none');
% 
xlabel('Contrast (%)');
ylabel('Response (s.d.)');
axis([0 1 0 60]);
set(gca,'XTick',[0 1],'XTickLabel',{'0%','100%'});
set(gca,'YTick',[0:10:40]);
% 
% % create arrows at the 75% value
con_con = response_(2,1,750);
con_cat = response_(2,2,750);
% 
arrow([.75 con_cat],[.75 con_con],'Length',8);
arrow([.75 con_con],[.75 con_cat],'Length',8);
text(.8,mean([con_cat con_con]),sprintf('Ratio: %01.2f',ci_con(2)));
% 
drawPublishAxis('figSize=[8.9,4]');

subplot(212); hold on
plot(x,squeeze(response_(1,1,:)),'-','Color',cmap(9,:),'LineWidth',2);
plot(x,squeeze(response_(1,2,:)),'-','Color',cmap(8,:),'LineWidth',1);

legend({'Motion - control','Motion - catch'},'FontSize',7,'FontName','Helvetica');
l = legend(gca,'boxoff');
set(l,'Color','none');

xlabel('Motion coherence (%)');
ylabel('Response (s.d.)');
axis([0 1 0 10]);
set(gca,'XTick',[0 1],'XTickLabel',{'0%','100%'});
set(gca,'YTick',[0 10]);


mot_con = response_(1,1,750);
mot_cat = response_(1,2,750);

arrow([.75 mot_cat],[.75 mot_con],'Length',4);
arrow([.75 mot_con],[.75 mot_cat],'Length',4);
text(.8,mean([mot_cat mot_con]),sprintf('Ratio: %01.2f',ci_coh(2)));

drawPublishAxis('figSize=[8.9,4]');

savepdf(h,fullfile(datafolder,'avg_models','catch_resp.pdf'));

%% Compute the various weights
v1_catch_contrast = b_u(:,1,1);
v1_catch_motion = b_u(:,1,2);
mt_catch_contrast = b_u(:,2,1);
mt_catch_motion = b_u(:,2,2);
v1_control_contrast = b_c(:,1,2);
v1_control_motion = b_c(:,1,1);
mt_control_contrast = b_c(:,2,2);
mt_control_motion = b_c(:,2,1);

%%
figure(1);

clf

bins = -1:.1:1;

% plot the effect of normalized MT on contrast readout in control and catch
clf;
hold on
plot([-40 40],[-40 40],'--r');
plot(mt_control_contrast,mt_catch_contrast,'*b');
plot(v1_control_contrast,v1_catch_contrast,'*r');
axis equal

% average effects across subjects?

% sigmas(sigmas==1) = NaN;
use = setdiff(1:21,[3 13]);

contrast_catch_ratio = v1_catch_contrast./mt_catch_contrast;
contrast_control_ratio = v1_control_contrast./mt_control_contrast;
motion_catch_ratio = mt_catch_motion./v1_catch_motion;
motion_control_ratio = mt_control_motion./v1_control_motion;

con_catch = mean(contrast_catch_ratio(use));
con_control = mean(contrast_control_ratio(use));
mot_catch = mean(motion_catch_ratio(use));
mot_control = mean(motion_control_ratio(use));

%% Collect sigmas and indiv r2
restructure_afits;

clear r2 sigmas cd
for ai = 1:length(aSIDs)
    for ni = 1:2
        for ropt = 1
            cm = afits{ai}{ni};
            BIC(ai,ni) = cm.BIC;
            sigmas(ai,ni,ropt) = cm.params.sigma;
            r2(ai,ni,ropt) = -sum(cm.cv.like);
            cd(ai,ni,ropt) = cm.cv.cd;
        end
    end
end

sigmas(sigmas==1) = NaN;

% r2 = r2(:,:,1,1);
% cd = cd(:,:,1,1);

%% Report R^2 (exp)
add = squeeze(r2(:,1));
poi = squeeze(r2(:,2));
diffe = add-poi;

aci = bootci(1000,@mean,add);
pci = bootci(1000,@mean,poi);
dci = bootci(1000,@mean,diffe);

figure;
hist(diffe)

%% cd vs r2
add_cd = squeeze(cd(:,3));
poi_cd = squeeze(cd(:,4));

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
drawPublishAxis('figSize=[8.9,2.5]');
savepdf(h,fullfile(datafolder,'avg_models','add_poiss_hist.pdf'));

%% CD plot
h = figure; hold on

plot(cd(:,1),cd(:,2),'o','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',5);
plot([min(cd(:)) max(cd(:))],[min(cd(:)) max(cd(:))],'--r');
xlabel('Additive model r^2');
ylabel('Poisson model r^2');
axis([.35 .5 .35 .5]);
set(gca,'XTick',[.35 .4 .45 .5],'YTick',[.35 .4 .45 .5]);
drawPublishAxis('figSize=[4.5,4.5]');
savepdf(h,fullfile(datafolder,'avg_models','add_poiss_r2.pdf'));

%% Contrast/coherence responses used to fit behavior
h = figure; hold on

cmap = brewermap(7,'PuOr');

x = 0:.001:1;
plot(x,squeeze(respcon_(1,:)),'Color',cmap(2,:));
plot(x,squeeze(respcoh_(8,:)),'Color',cmap(6,:));

l = legend({'V1 Contrast response function','MT Coherence response function'});
set(l,'box','off');

set(gca,'XTick',[0 1],'XTickLabel',[0 100]);

xlabel('Stimulus strength (%)');
ylabel('\Delta signal (%)');

drawPublishAxis('figSize=[18,14]');

savepdf(h,fullfile('~/proj/att_awe/talks/data_figures/v1mt.pdf'));

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

for bi = 1:length(bmodels)
    h = figure; clf; hold on

    csensitivity = squeeze(sensitivity(bi,:,:,:));
    
    ci = bootci(1000,@mean,csensitivity);
    s_mean = squeeze(mean(ci));
    s_std = squeeze(ci(2,:,:,:))-s_mean;
        
%     tx = [.07 .07 .07 .07 .07 -.1 -.2 -.17];
%     ty = [.07 .07 .07 -.07 .07 -.15 -.13 .07];
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
    axis([-10 20 -5 15]);
    set(gca,'XTick',[-5 0 5 10 15]','YTick',[-5 0 5 10 15]);
    xlabel('Contrast weight (a.u.)');
    ylabel('Coherence weight (a.u.)');
    title(sprintf('Weights under %s noise',bmodels_text{bi}));
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