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
% load(fullfile(datafolder,'avg_indiv_fits_att_cross_2.mat'));
% attfits_ = attfits; clear attfits_2
% for ai = 1:21
%     for mi = 1:2
%         attfits_2{ai,mi} = attfits_{(ai-1)*2+mi};
%     end
% end
% load(fullfile(datafolder,'avg_indiv_fits_att_cross_68.mat'));
% attfits_ = attfits; 
% count = 1;
% for ai = 1:21
%     for mi = 1:2
% %         for ropt = 1:2
%             attfits_68{ai,mi} = attfits_{count};
%             count = count+1;
% %         end
%     end
% end
load(fullfile(datafolder,'avg_indiv_fits_att_cross_linear.mat'));
attfits_ = attfits; 
count = 1;
for ai = 1:21
    for mi = 1:2
%         for ropt = 1:2
            attfits_l{ai,mi} = attfits_{count};
            count = count+1;
%         end
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
        cd_like_68(i,mi) = -sum(attfits_68{i,mi}.cv.like);
        cd_att_68(i,mi) = attfits_68{i,mi}.cv.cd;
    end
end

for i = 1:21
    rois = {'V1','MT'};
    for ri = 1:2
        w(i,ri,1) = attfits_2{i,1}.params.(sprintf('beta_control_%s_cohw',rois{ri}));
        w(i,ri,2) = attfits_2{i,1}.params.(sprintf('beta_control_%s_conw',rois{ri}));
        w1(i,ri) = attfits_2{i,2}.params.(sprintf('beta_control_%s_w',rois{ri}));
    end
    rois = {'V1','V2','MT'};
    for ri = 1:3
        w_68(i,ri,1) = attfits_68{i,1}.params.(sprintf('beta_control_%s_cohw',rois{ri}));
        w_68(i,ri,2) = attfits_68{i,1}.params.(sprintf('beta_control_%s_conw',rois{ri}));
        w1_68(i,ri) = attfits_68{i,2}.params.(sprintf('beta_control_%s_w',rois{ri}));
    end
    rois = {'V1','V2','MT'};
    for ri = 1:3
        w_l(i,ri,1) = attfits_l{i,1}.params.(sprintf('beta_control_%s_cohw',rois{ri}));
        w_l(i,ri,2) = attfits_l{i,1}.params.(sprintf('beta_control_%s_conw',rois{ri}));
        w1_l(i,ri) = attfits_l{i,2}.params.(sprintf('beta_control_%s_w',rois{ri}));
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

rois = {'V1','V2','MT'};
plot_rightchoice_model_att(attfits_68(:,1),respcon_([1 2 8],:,:),respcoh_([1 2 8],:,:),aSIDs,bmodels(1),rois);
plot_rightchoice_model_att_onebeta(attfits_68(:,2),respcon_([1 2 8],:,:),respcon_([1 2 8],:,:),aSIDs,bmodels(2),rois);

rois = {'V1','V2','MT'};
plot_rightchoice_model_att(attfits_l(:,1),respcon_l([1 2 8],:,:),respcoh_l([1 2 8],:,:),aSIDs,bmodels(1),rois);
plot_rightchoice_model_att_onebeta(attfits_l(:,2),respcon_l([1 2 8],:,:),respcon_l([1 2 8],:,:),aSIDs,bmodels(2),rois);
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

%% Show ideal vs. actual
goalb2 = squeeze(mean(w));
rc = squeeze(respcon_([1 8],2,:));
rm = squeeze(respcoh_([1 8],2,:));
rc = rc - repmat(rc(:,1),1,1001);
rm = rm - repmat(rm(:,1),1,1001);
c2 = rc'*goalb2(:,2);
m2 = rm'*goalb2(:,1);
actualb3 = squeeze(mean(w_68));
rc = squeeze(respcon_([1 2 8],2,:));
rm = squeeze(respcoh_([1 2 8],2,:));
rc = rc - repmat(rc(:,1),1,1001);
rm = rm - repmat(rm(:,1),1,1001);
c3 = rc'*actualb3(:,2);
m3 = rm'*actualb3(:,1);

h = figure;
subplot(211); title('Coherence (attended)'); hold on
plot(x,m2,'-k');
plot(x,m3,'Color',cmap(6,:));
am = axis;
subplot(212); title('Contrast (attended)'); hold on
plot(x,c2,'-k');
plot(x,c3,'Color',cmap(2,:));
ac = axis;

c3 = rc'*actualb3(:,1);
m3 = rm'*actualb3(:,2);
h = figure;
subplot(211); title('Coherence (unatt)'); hold on
hline(0,'-k');
plot(x,m3,'Color',cmap(6,:));
axis([am(1) am(2) min(-1,am(3)) am(4)]);
subplot(212); title('Contrast (unatt)'); hold on
hline(0,'-k');
plot(x,c3,'Color',cmap(2,:));
axis([ac(1) ac(2) min(-1,ac(3)) ac(4)]);
%% Test different beta values (V1/MT)
% just take attend contrast
goalb = squeeze(mean(w));
rc = squeeze(respcon_([1 8],2,:));
rm = squeeze(respcoh_([1 8],2,:));
rc = rc - repmat(rc(:,1),1,1001);
rm = rm - repmat(rm(:,1),1,1001);
gc = rc; gm = rm;

rc = squeeze(respcon_([1 2 8],2,:));
rm = squeeze(respcoh_([1 2 8],2,:));
rc = rc - repmat(rc(:,1),1,1001);
rm = rm - repmat(rm(:,1),1,1001);
cmap = brewermap(7,'PuOr');
bs = -32:8:32;

x = 0:.001:1;
for b3i = 1:length(bs)
    h = figure;
    for b1i = 1:length(bs) % v1 weight
        for b2i = 1:length(bs) % mt weight 
            b1 = bs(b1i);
            b2 = bs(b2i);
            b3 = bs(b3i);
            beta = [b1 b2 b3]';
            subplot(length(bs),length(bs),(b1i-1)*length(bs)+b2i);
            rc_ = rc'*beta;
            rm_ = rm'*beta;
            title(sprintf('V1 %i MT %i R %02.2f',b1,b2,rm_\rc_));
            hold on
            % plot goals
            plot(x,gc'*goalb(:,2),'-k');
            plot(x,gm'*goalb(:,1),'-k');
            % plot actuals
            plot(x,rc_,'Color',cmap(2,:));
            plot(x,rm_,'Color',cmap(6,:));
            set(gca,'XTick',[],'YTick',[]);
        end
    end
end

%% Testing specific values of functions
goalb2 = squeeze(mean(w));
rc = squeeze(respcon_([1 8],2,:));
rm = squeeze(respcoh_([1 8],2,:));
rc = rc - repmat(rc(:,1),1,1001);
rm = rm - repmat(rm(:,1),1,1001);
c2 = rc'*goalb2(:,2);
m2 = rm'*goalb2(:,1);

b_con = [17 0 -2];
b_coh = [-7 0 30];


b_con = [19 0 -7];
b_coh = [-5.5 12.5 20.9];

rc = squeeze(respcon_([1 2 8],2,:));
rm = squeeze(respcoh_([1 2 8],1,:));
rc = rc - repmat(rc(:,1),1,1001);
rm = rm - repmat(rm(:,1),1,1001);
c3 = rc'*b_con';
m3 = rm'*b_coh';

rc = squeeze(respcon_([1 2 8],1,:));
rm = squeeze(respcoh_([1 2 8],2,:));
rc = rc - repmat(rc(:,1),1,1001);
rm = rm - repmat(rm(:,1),1,1001);

h = figure;
subplot(2,2,1); title('Coherence (attended)'); hold on
plot(x,m2,'-k');
plot(x,m3,'Color',cmap(6,:));
am = axis;
axis([am(1) am(2) min(-4,am(3)) am(4)]);
subplot(2,2,2); title('Coherence (unatteded)'); hold on
hline(0,'-k');
plot(x,rm'*b_con','Color',cmap(6,:));
axis([am(1) am(2) min(-4,am(3)) am(4)]);
subplot(2,2,3); title('Contrast (attended)'); hold on
plot(x,c2,'-k');
plot(x,c3,'Color',cmap(2,:));
axis([ac(1) ac(2) min(-4,ac(3)) ac(4)]);
subplot(2,2,4); title('Contrast (unatt)'); hold on
hline(0,'-k');
plot(x,rc'*b_coh','Color',cmap(2,:));
axis([ac(1) ac(2) min(-4,ac(3)) ac(4)]);

%% Linear
% solve the linear system in one move
rc = respcon_;
rm = respcoh_;
areas = 1:8;

% find [a b c] such that [respcon_l ; respcoh_l] * [abc abc] = [ 1 0];
data1 = [squeeze(rc(areas,2,:))' ; squeeze(rm(areas,2,:))'];
out1 = [25*x' ; -.1*x'];

cb = data1\out1;

figure;
subplot(211)
plot(data1*cb);

data2 = [squeeze(rc(areas,1,:))' ; squeeze(rm(areas,1,:))'];
out2 = [6*x' ; .2*x'];

mb = data2\out2;
subplot(212)
plot(data2*mb);

% stack both
data = [data1;data2];
out = [out1;out2];

ab = data\out;
figure;
hold on;
plot(out,'-k');
plot(data*ab,'-r');
axis([0 4000 -5 25]);

%%
afl = attfits_2(:,1);
rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
rois = rois(areas);
cons = {'cohw','conw'};
for ai = 1:21
    for ri = areas
        afl{ai}.params.(sprintf('beta_control_%s_conw',rois{ri})) = cb(ri);
        afl{ai}.params.(sprintf('beta_control_%s_cohw',rois{ri})) = mb(ri);
    end
    afl{ai}.params.bias = 0;
    
end
afl_multi = afl;

plot_rightchoice_model_att(afl,rc(areas,:,:),rm(areas,:,:),aSIDs,bmodels(1),rois);

afl = attfits_2(:,2);
for ai = 1:21
    for ri = areas
        afl{ai}.params.(sprintf('beta_control_%s_w',rois{ri})) = ab(ri);
    end
    afl{ai}.params.bias = 0;
    
end

plot_rightchoice_model_att_onebeta(afl,rc(areas,:,:),rm(areas,:,:),aSIDs,bmodels(1),rois);

%% Get the likelihood for each model

info = struct;
% info.sigma = sigma;
info.model = bmodels{1};
info.rois = 1:8;
info.respcon = respcon_;
info.respcoh = respcoh_;
parfor ai = 1:21
    adata = loadadata(sprintf('s%03.0f',aSIDs(ai)));
    cinfo = info;
    cinfo.lapse = lapses(ai);
    cinfo.fitmodel = afl_multi{ai};
    cinfo.fitmodel.numParams = 16;
    lfit{ai} = fitCCBehavControlModel_fmri(adata,cinfo,1);
end

%%
parfor ai = 1:21
    adata = loadadata(sprintf('s%03.0f',aSIDs(ai)));
    cinfo = info;
    cinfo.lapse = lapses(ai);
    cinfo.model = strcat(cinfo.model,',onebeta');
    cinfo.fitmodel = afl{ai};
    lfit_a{ai} = fitCCBehavControlModel_fmri(adata,cinfo,1);
end

%% compare
for ai = 1:21
    adata = loadadata(sprintf('s%03.0f',aSIDs(ai)));
    n(ai) = size(adata,1);
    l_lin(ai,1) = -lfit{ai}.likelihood;
    bic(ai,1) =  -2*l_lin(ai,1) + 16 * log(size(adata,1));
    l_lin(ai,2) = -lfit_a{ai}.likelihood;
    bic(ai,2) =  -2*l_lin(ai,2) + 8 * log(size(adata,1));
    
    cd(ai,1) = lfit{ai}.cd;
    cd(ai,2) = lfit_a{ai}.cd;
end
d = l_lin(:,1)-l_lin(:,2);
dcd = bic(:,1)-bic(:,2);

dci = bootci(10000,@mean,d);

disp(sprintf('Mean difference %02.2f, 95%% CI [%02.2f %02.2f]',mean(dci),dci(1),dci(2)));
dci = bootci(10000,@mean,dcd);

disp(sprintf('Mean difference %02.2f, 95%% CI [%02.2f %02.2f]',mean(dci),dci(1),dci(2)));