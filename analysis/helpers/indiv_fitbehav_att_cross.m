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

%% Check that responses look right
figure; hold on
ro = [1 8];
di = 1;
cons = {'attend coherence','attend contrast'};
cmap = brewermap(7,'PuOr');
for rii = 1:2
    subplot(2,1,rii); hold on
    ri = ro(rii);
    title(sprintf('%s: %s',rois{ri},cons{di}));
%     plot(squeeze(mean(respcon_([1 2 3 4],1,:),1)),'--r');
%     plot(squeeze(mean(respcoh_([5 8],1,:),1)),'--');
    plot(x,squeeze(mean(respcon_(ri,di,:),1)),'-','Color',cmap(2,:));
    plot(x,squeeze(mean(respcoh_(ri,di,:),1)),'-','Color',cmap(6,:));
    axis([0 1 0 1.9]);
end

%% Across ROIs
h = figure; hold on
ro = [1 2 3 6 7 8];
for rii = 1:length(ro)
    ri = ro(rii);
    subplot(2,1,1); hold on
    rc = squeeze(respcon_(ri,1,:));
    plot(x,rc-rc(1),'Color',cmap(2,:));
    subplot(2,1,2); hold on
    rm = squeeze(respcoh_(ri,1,:));
    plot(x,rm-rm(1),'Color',cmap(6,:));
end

h = figure; hold on
ro = [1 2 3 6 7 8];
for rii = 1:length(ro)
    ri = ro(rii);
    subplot(6,1,rii); hold on
%     subplot(2,1,1); hold on
    rc = squeeze(respcon_(ri,1,:));
    plot(x,rc-rc(1),'Color',cmap(2,:));
%     subplot(2,1,2); hold on
    rm = squeeze(respcoh_(ri,1,:));
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
bmodels = {'sigma,roi,att','sigma,roi,att,onebeta'};%,'doublesigma','doublesigma,poisson'};
% bmodels = {'sigma,roi'};

% options list
attopts = zeros(10000,5);
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
    
    for ii = bstart:bend
        
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
        
%         temps = cell{1,4};
%         parfor iii=1:4
%             
%             tinfo = info;
%             tinfo.model = bmodels{iii};
%             temps{iii} = fitCCBehavControlModel_fmri(adata,tinfo,1);
%         end
        attfits{ii} = fitCCBehavControlModel_fmri(adata,info,1);
        disp(sprintf('Done with %i',subj));
   end
    
%     disppercent(bend/size(attopts,1));
end
% disppercent(inf);

% save(fullfile(datafolder,'avg_indiv_fits_att_cross.mat'),'attfits');

%% Restructure attfits
load(fullfile(datafolder,'avg_indiv_fits_att_cross_2.mat'));
attfits_ = attfits; clear attfits_2
for ai = 1:21
    for mi = 1:2
        attfits_2{ai,mi} = attfits_{(ai-1)*2+mi};
    end
end
load(fullfile(datafolder,'avg_indiv_fits_att_cross.mat'));
attfits_ = attfits; clear attfits
for ai = 1:21
    for mi = 1:2
        attfits{ai,mi} = attfits_{(ai-1)*2+mi};
    end
end
%% Compare afits and attfits
restructure_afits;

for i = 1:21
    like(i) = -sum(afits{i}{1}.cv.like);
    cd(i) = afits{i}{1}.cv.cd;
end

for i = 1:21
    for mi = 1:2
        cd_like_8(i,mi) = -sum(attfits{i,mi}.cv.like);
        cd_att_8(i,mi) = attfits{i,mi}.cv.cd;
        cd_like(i,mi) = -sum(attfits_2{i,mi}.cv.like);
        cd_att(i,mi) = attfits_2{i,mi}.cv.cd;
    end
end

for i = 1:21
    rois = {'V1','MT'};
    for ri = 1:2
        w(i,ri,1) = attfits_2{i,1}.params.(sprintf('beta_control_%s_cohw',rois{ri}));
        w(i,ri,2) = attfits_2{i,1}.params.(sprintf('beta_control_%s_conw',rois{ri}));
        w1(i,ri) = attfits_2{i,2}.params.(sprintf('beta_control_%s_w',rois{ri}));
    end
    rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
    for ri = 1:8
        w_8(i,ri,1) = attfits{i,1}.params.(sprintf('beta_control_%s_cohw',rois{ri}));
        w_8(i,ri,2) = attfits{i,1}.params.(sprintf('beta_control_%s_conw',rois{ri}));
        w1_8(i,ri) = attfits{i,2}.params.(sprintf('beta_control_%s_w',rois{ri}));
    end
end

%% Get mean weights
w_ci = bootci(10000,@mean,w);
w_ = squeeze(mean(w));

%% Plot likelihood and weights
h = figure;
subplot(211);
[b,x] = hist(cd_like(:,1)-cd_like(:,2));
bar(x,b,'FaceColor',[0.8 0.8 0.8]);
xlabel('Likelihood (Multiple - Single readout)');
set(gca,'XTick',[0 25 50]);
drawPublishAxis('figSize=[3.5,4.5]');

cmap = brewermap(7,'PuOr');
subplot(212); hold on
plot(w_(1,2),w_(1,1),'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','w');
text(w_(1,2),w_(1,1),'V1');
plot(w_(2,2),w_(2,1),'o','MarkerFaceColor',cmap(6,:),'MarkerEdgeColor','w');
text(w_(2,2),w_(2,1),'MT');
v = hline(0,'--'); set(v,'Color',[0.8 0.8 0.8]);
v = vline(0,'--'); set(v,'Color',[0.8 0.8 0.8]);
v = hline(mean(w1(:,1)),'-'); set(v,'Color',cmap(2,:));
v = hline(mean(w1(:,2)),'-'); set(v,'Color',cmap(6,:));
v = vline(mean(w1(:,1)),'-'); set(v,'Color',cmap(2,:));
v = vline(mean(w1(:,2)),'-'); set(v,'Color',cmap(6,:));
axis([-15 25 -10 40]);
set(gca,'XTick',[-10 0 10 20]);
set(gca,'YTick',[-10 0 10 20]);
xlabel('Contrast discrimination');
ylabel('Coherence discrimination');
drawPublishAxis('figSize=[3.5,4.5]');

% plot(cd_att(:,1),cd_att(:,2),'o','MarkerFaceColor','k','MarkerEdgeColor','w');
% x = [min(cd_att(:,1)) max(cd_att(:,1))];
% plot(x,x,'--r');
% xlabel('Multiple readouts');
% ylabel('Single readout');
% title('Variance explained (R^2)');
% axis([0.15 0.4 0.15 0.4]);
% set(gca,'XTick',[0.15 0.25 0.35],'XTickLabel',{'15%','25%','35%'});
% set(gca,'YTick',[0.15 0.25 0.35],'YTickLabel',{'15%','25%','35%'});


savepdf(h,fullfile(datafolder,'avg_models','onebeta_comparison.pdf'));

%% Plot

rois = {'V1','MT'};

% bmodels = {'roi2'};
plot_rightchoice_model_att(attfits_2(:,1),respcon_([1 8],:,:),respcoh_([1 8],:,:),aSIDs,bmodels(1),rois);
plot_rightchoice_model_att_onebeta(attfits_2(:,2),respcon_([1 8],:,:),respcoh_([1 8],:,:),aSIDs,bmodels(2),rois);

    rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
% bmodels = {'roi8'};
plot_rightchoice_model_att(attfits(:,1),respcon_,respcoh_,aSIDs,bmodels(1),rois);
plot_rightchoice_model_att_onebeta(attfits(:,2),respcon_,respcoh_,aSIDs,bmodels(2),rois);


%% Example plots for justin: Readout space
cmap = brewermap(7,'PuOr');
% Compute the readout space under attention for one beta (8 ROIs)
% FEATURE | ATTENTION
clear readout
features = {'respcoh_','respcon_'};
for ci = 1:2
    for di = 1:2
        feat = eval(features{ci});
        readout(ci,di,:) = squeeze(feat(:,di,:))' * squeeze(mean(w1_8))';
    end
end

h = figure; hold on
color = [6 2];
type = {'-','--','--','-'};
for ci = 1:2
    for di = 1:2
        plot(x,squeeze(readout(ci,di,:)),type{(ci-1)*2+di},'Color',cmap(color(ci),:));
    end
end

% Compute the readout space under attention for multiple (8
% ROIs)
clear readout
features = {'respcoh_','respcon_'};
for ci = 1:2
    for di = 1:2
        feat = eval(features{ci});
        readout(ci,di,:) = squeeze(feat(:,di,:))' * squeeze(mean(w_8(:,:,di)))';
    end
end

h = figure; hold on
color = [6 2];
type = {'-','--','--','-'};
for ci = 1:2
    for di = 1:2
        plot(x,squeeze(readout(ci,di,:)),type{(ci-1)*2+di},'Color',cmap(color(ci),:));
    end
end

% Compute the readout space response for the passive responses (8 ROIs)
