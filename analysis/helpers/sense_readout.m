%% Load all relevant datasets: Sensitivity
nSIDs = [305 329 43 25 300 346 343 344 338 340 348];
bSIDs = [345 330 337 335 349 354 353 356 334 352]; % behavior only participants
aSIDs = [nSIDs bSIDs];
ncorrespond = [1 2 3 4 5 6 7 9 10 11]; % which nSIDs is referenced by each attSID

% Load the sensitivity data for V1->MT during passive viewing
% more specifically: get the contrast/coherence response function averages
load(fullfile(datafolder,'avg_hrffits.mat'));

x = 0:.001:1;
respcon = zeros(11,8,length(x));
respcoh = zeros(11,8,length(x));

for si = 1:11
    fit = sfits{si}{1,4}; % 4 refers to resp-25, which is our standard model

    for ri = 1:8
        respcon(si,ri,:) = conModel(x,fit.roifit{ri}.params)+fit.roifit{ri}.params.offset;
        respcoh(si,ri,:) = cohModel(x,fit.roifit{ri}.params)+fit.roifit{ri}.params.offset;
    end
end

% rc_passive = squeeze(mean(bootci(10000,@mean,respcon)));
% rm_passive = squeeze(mean(bootci(10000,@mean,respcoh)));

% Load the sensitivity data for V1->MT during contrast and coherence
% discrimination
load(fullfile(datafolder,'avg_att_cross_fits_sb.mat'));

x = 0:.001:1;
respcon_ac = zeros(10,8,length(x));
respcon_am = zeros(10,8,length(x));
respcoh_ac = zeros(10,8,length(x));
respcoh_am = zeros(10,8,length(x));

for si = 1:10
    fit = attfits{si}{1}; % 4 refers to resp-25, which is our standard model

    for ri = 1:8
        respcon_am(si,ri,:) = fit.roifit{ri}.conresp_coh;
        respcon_ac(si,ri,:) = fit.roifit{ri}.conresp_con;
        respcoh_am(si,ri,:) = fit.roifit{ri}.cohresp_coh;
        respcoh_ac(si,ri,:) = fit.roifit{ri}.cohresp_con;
    end
end

% rc_attcon = squeeze(mean(bootci(10000,@mean,respcon(:,:,2,:))));
% rm_attcon = squeeze(mean(bootci(10000,@mean,respcoh(:,:,2,:))));
% rc_attcoh = squeeze(mean(bootci(10000,@mean,respcon(:,:,1,:))));
% rm_attcoh = squeeze(mean(bootci(10000,@mean,respcoh(:,:,1,:))));

%% Obtain SnR multipliers for each ROI

snr = zeros(10,2,8);
for si = 1:10
    for ri = 1:8
        pcon = respcon(ncorrespond(si),ri,1);
        pcoh = respcoh(ncorrespond(si),ri,1);
        attcon = respcon_ac(si,ri,1);
        attcoh = respcoh_ac(si,ri,1);
        evals = {'pcon','pcoh','attcon','attcoh'};
        for ei = 1:4
            eval(sprintf('if %s<0, %s=0; end',evals{ei},evals{ei}));
        end
        snr(si,2,ri) = attcon/pcon;
        snr(si,1,ri) = attcoh/pcoh;
    end
end

snr(snr==Inf) = nan;

roi_snr = squeeze(nanmean(snr(:,1,:)));

%% Remove intercepts (we only care about slopes)

datasets = {'respcon','respcoh','respcon_am','respcon_ac','respcoh_am','respcoh_ac'};

for ri = 1:8
    for di = 1:length(datasets)
        
        eval(sprintf('%s(:,ri,:) = %s(:,ri,:) - repmat(%s(:,ri,1),1,1,1001);',datasets{di},datasets{di},datasets{di}));
    end
end

%% Load: weights

% Load the readout weights from the linking model
restructure_afits;

weights_passive = zeros(2,length(aSIDs),8,2);
ROIs = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
cons = {'cohw','conw'};
for bi = 1:2
    for ai = 1:length(aSIDs)
        for ri = 1:8
            for ci = 1:2
                weights_passive(bi,ai,ri,ci) = afits{ai}{bi,2}.params.(sprintf('beta_control_%s_%s',ROIs{ri},cons{ci}));
            end
        end
    end
end

weights_passive = squeeze(weights_passive(1,:,:,:));
w_addpass = squeeze(mean(bootci(10000,@mean,weights_passive)));

%% Load: weights 2

% Load the readout weights from the linking model
restructure_afits;

weights_passive = zeros(2,length(aSIDs),2,2);
ROIs = {'V1','MT'};
cons = {'cohw','conw'};
for bi = 1:2
    for ai = 1:length(aSIDs)
        for ri = 1:2
            for ci = 1:2
                weights_passive(bi,ai,ri,ci) = afits{ai}{bi,1}.params.(sprintf('beta_control_%s_%s',ROIs{ri},cons{ci}));
            end
        end
    end
end

weights_passive = squeeze(weights_passive(1,:,:,:));
w_addpass = squeeze(mean(bootci(10000,@mean,weights_passive)));

%% Load: weights one beta 2

load(fullfile(datafolder,'avg_indiv_fits_onebeta_att_2.mat'));
clear w
for ai = 1:10
    w(ai,1) = attfits{ai}.params.beta_control_V1_w;
    w(ai,2) = attfits{ai}.params.beta_control_MT_w;
end

%% Sensory space sensitivity plots

% normalize functions
norm_con = respcon(:,1,:);
norm_coh = respcoh(:,8,:);

norm_conatt = norm_con;%respcon_ac(:,1,:);
norm_cohatt = norm_coh; %respcoh_am(:,8,:);

clear rc_p rc_c rc_m rm_p rm_c rm_m
ros = [1:8];
x = 0:.001:1;
for rii = 1:length(ros)
    ri = ros(rii);
    
    for ni = 1:10
        rc_p(ni,rii) = x'\squeeze(respcon(ncorrespond(ni),ri,:));
        rc_c(ni,rii) = x'\squeeze(respcon_ac(ni,ri,:));
        rc_m(ni,rii) = x'\squeeze(respcon_am(ni,ri,:));
        
        rm_p(ni,rii) = x'\squeeze(respcoh(ncorrespond(ni),ri,:));
        rm_c(ni,rii) = x'\squeeze(respcoh_ac(ni,ri,:));
        rm_m(ni,rii) = x'\squeeze(respcoh_am(ni,ri,:));
        
%         rc_p(ni,rii) = norm_con(ncorrespond(ni),:)'\squeeze(respcon(ncorrespond(ni),ri,:));
%         rc_c(ni,rii) = norm_conatt(ni,:)'\squeeze(respcon_ac(ni,ri,:));
%         rc_m(ni,rii) = norm_conatt(ni,:)'\squeeze(respcon_am(ni,ri,:));
% 
%         rm_p(ni,rii) = norm_coh(ncorrespond(ni),:)'\squeeze(respcoh(ncorrespond(ni),ri,:));
%         rm_c(ni,rii) = norm_cohatt(ni,:)'\squeeze(respcoh_ac(ni,ri,:));
%         rm_m(ni,rii) = norm_cohatt(ni,:)'\squeeze(respcoh_am(ni,ri,:));
    end
end

%% Correlations

for ni = 1:10
    con(ni) = corr(rc_c(ni,:)',rc_m(ni,:)');
    coh(ni) = corr(rm_c(ni,:)',rm_m(ni,:)');
end

cic = bootci(1000,@mean,con);
cic
mean(cic)
cim = bootci(1000,@mean,coh);
cim
mean(cim)

%% Average across subjects
rc_p_ci = bootci(1000,@median,rc_p);
rc_c_ci = bootci(1000,@median,rc_c);
rc_m_ci = bootci(1000,@median,rc_m);
rm_p_ci = bootci(1000,@median,rm_p);
rm_c_ci = bootci(1000,@median,rm_c);
rm_m_ci = bootci(1000,@median,rm_m);

group = {'rc','rm'};
disc = {'p','c','m'};
for gi = 1:2
    for di = 1:3
        eval(sprintf('%s_%s_i = %s_%s;',group{gi},disc{di},group{gi},disc{di})); % save individual
        eval(sprintf('%s_%s = mean(%s_%s_ci);',group{gi},disc{di},group{gi},disc{di}));
    end
end

%% New plot: dashed line surrounding 
h = figure;  hold on
cmap = brewermap(7,'PuOr');
rois = {'V1','V2','V3','V3a','V3b','V4','V7','MT'};

for ri = 1:8
    p(1) = bar(ri-0.1,rc_c(ri),0.2,'FaceColor',cmap(2,:),'EdgeColor','w');
    p(2) = bar(ri+0.1,rc_m(ri),0.2,'FaceColor',cmap(6,:),'EdgeColor','w');
end
% add error bars
errbar((1:8)-0.1,rc_c,rc_c_ci(2,:)-rc_c,'-k');
errbar((1:8)+0.1,rc_m,rc_m_ci(2,:)-rc_m,'-k');

% add the dashed passive conditions
for ri = 1:8
    plot(ri+[-.2 -.2],rc_p(ri)*[0 1],'--k');
    plot(ri+[.2 .2],rc_p(ri)*[0 1],'--k');
    plot(ri+[-.2 .2],rc_p(ri)*[1 1],'--k');
end
set(gca,'XTick',[1 8],'XTickLabel',rois([1 8]),'YTick',[0 0.5 1]);
ylabel('Contrast sensitivity');
legend(p,{'Discriminating contrast','Discriminating coherence'});
drawPublishAxis('figSize=[8.5, 4.5]');

savepdf(h,fullfile(datafolder,'avg_fitatt','contrast_sensitivity_dashed.pdf'));

h = figure;  hold on
cmap = brewermap(7,'PuOr');
rois = {'V1','V2','V3','V3a','V3b','V4','V7','MT'};

for ri = 1:8
    p(1) = bar(ri-0.1,rm_c(ri),0.2,'FaceColor',cmap(2,:),'EdgeColor','w');
    p(2) = bar(ri+0.1,rm_m(ri),0.2,'FaceColor',cmap(6,:),'EdgeColor','w');
end
% add error bars
errbar((1:8)-0.1,rm_c,rm_c_ci(2,:)-rm_c,'-k');
errbar((1:8)+0.1,rm_m,rm_m_ci(2,:)-rm_m,'-k');

% add the dashed passive conditions

for ri = 1:8
    plot(ri+[-.2 -.2],rm_p(ri)*[0 1],'--k');
    plot(ri+[.2 .2],rm_p(ri)*[0 1],'--k');
    plot(ri+[-.2 .2],rm_p(ri)*[1 1],'--k');
end
a = axis;
axis([a(1) a(2) a(3) max(a(4),0.5)]);
set(gca,'XTick',[1 8],'XTickLabel',rois([1 8]),'YTick',[0 0.5]);
ylabel('Coherence sensitivity');
% legend(p,{'Discriminating contrast','Discriminating coherence'});
drawPublishAxis('figSize=[8.5, 4.5]');

savepdf(h,fullfile(datafolder,'avg_fitatt','coherence_sensitivity_dashed.pdf'));

%% OLD PLOTS

%% Bar plot comparisons
h = figure;
cmap = brewermap(7,'PuOr');
rois = {'V1','V2','V3','V3a','V3b','V4','V7','MT'};

%%%%%%%%%%%%%%%%%%%%%
subplot(211); hold on
bar(1:8,rc_p,'FaceColor',cmap(2,:),'EdgeColor','w');
errbar(1:8,rc_p,abs(rc_p_ci(2,:)-rc_p),'-','Color',[0 0 0]);
a = axis; a(2) = 9;
axis([a(1) a(2) 0 2]);
set(gca,'XTick',[1 8],'XTickLabel',rois([1 8]),'YTick',[0 1 2]);
title('Sensitivity during passive viewing');
ylabel('Contrast');
drawPublishAxis;

%%%%%%%%%%%%%%%%%%%%%
subplot(212); hold on
bar(1:8,rm_p,'FaceColor',cmap(6,:),'EdgeColor','w');
errbar(1:8,rm_p,abs(rm_p_ci(2,:)-rm_p),'-','Color',[0 0 0]);
axis([a(1) a(2) 0 1]);
set(gca,'XTick',[1 8],'XTickLabel',rois([1 8]),'YTick',[0 1]);
ylabel('Coherence');
drawPublishAxis('figSize=[4.5,4.5]');

savepdf(h,fullfile(datafolder,'avg_fitatt','passive_sensitivity_bar.pdf'));
%%

h = figure;
subplot(3,3,[1 2]); hold on
bar(1:8,rc_c,'FaceColor',cmap(2,:),'EdgeColor','w');
errbar(1:8,rc_c,abs(rc_c_ci(2,:)-rc_c),'-','Color',[0 0 0]);
a = axis; a(2) = 9;
axis([a(1) a(2) 0 2]);
set(gca,'XTick',[1 8],'XTickLabel',rois([1 8]),'YTick',[0 1 2]);
drawPublishAxis;

subplot(3,3,[6 9]); hold on

barh(1:8,fliplr(rc_m),'FaceColor',cmap(6,:),'EdgeColor','w');
errbar(fliplr(rc_m),1:8,fliplr(abs(rc_m_ci(2,:)-rc_m)),'-','Color',[0 0 0],'horiz');
% a = axis; a(2) = 9;
% axis([a(1) a(2) 0 2]);
set(gca,'YTick',[1 8],'YTickLabel',fliplr(rois([1 8])),'XTick',[0 1 2]);
drawPublishAxis;

subplot(3,3,[4 5 7 8]); hold on

title('Contrast sensitivity');

plot([0 1],[0 1],'--k');
% error bars
for ri = 1:8
    plot([rc_c(ri) rc_c(ri)],rc_m_ci(:,ri),'-','Color',cmap(2,:));
    plot(rc_c_ci(:,ri),[rc_m(ri) rc_m(ri)],'-','Color',cmap(2,:));
end
% markers
plot(rc_c,rc_m,'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','w','MarkerSize',3);

% lsline
b = [ones(size(rc_c))' rc_c']\rc_m';
x = [min(rc_c) max(rc_c)];
y = b(1) + b(2)*x;
plot(x,y,'-r');

xlabel('Discriminating contrast');
ylabel('Discriminating coherence');
axis equal
axis([0 2 0 1.5]);
set(gca,'Xtick',[0 1],'YTick',[0 1]);

% compute and plot correlation


drawPublishAxis('figSize=[4.5,4.5]');

savepdf(h,fullfile(datafolder,'avg_fitatt','contrast_sensitivity_marginals.pdf'));


%%

h = figure;
subplot(3,3,[1 2]); hold on
bar(1:8,rm_c,'FaceColor',cmap(2,:),'EdgeColor','w');
errbar(1:8,rm_c,abs(rm_c_ci(2,:)-rm_c),'-','Color',[0 0 0]);
a = axis; a(2) = 9;
axis([a(1) a(2) 0 2]);
set(gca,'XTick',[1 8],'XTickLabel',rois([1 8]),'YTick',[0 1]);
drawPublishAxis;

subplot(3,3,[6 9]); hold on

barh(1:8,fliplr(rm_m),'FaceColor',cmap(6,:),'EdgeColor','w');
errbar(fliplr(rm_m),1:8,fliplr(abs(rm_m_ci(2,:)-rm_m)),'-','Color',[0 0 0],'horiz');
% a = axis; a(2) = 9;
% axis([a(1) a(2) 0 2]);
set(gca,'YTick',[1 8],'YTickLabel',fliplr(rois([1 8])),'XTick',[0 1]);
drawPublishAxis;

subplot(3,3,[4 5 7 8]); hold on

title('Contrast sensitivity');

plot([0 1],[0 1],'--k');
% error bars
for ri = 1:8
    plot([rm_c(ri) rm_c(ri)],rm_m_ci(:,ri),'-','Color',cmap(6,:));
    plot(rm_c_ci(:,ri),[rm_m(ri) rm_m(ri)],'-','Color',cmap(6,:));
end
% markers
plot(rm_c,rm_m,'o','MarkerFaceColor',cmap(6,:),'MarkerEdgeColor','w','MarkerSize',3);

% lsline
b = [ones(size(rm_c))' rm_c']\rm_m';
x = [min(rm_c) max(rm_c)];
y = b(1) + b(2)*x;
plot(x,y,'-r');

xlabel('Discriminating contrast');
ylabel('Discriminating coherence');
axis equal
axis([0 2 0 1.5]);
set(gca,'Xtick',[0 1],'YTick',[0 1]);

% compute and plot correlation


drawPublishAxis('figSize=[4.5,4.5]');

savepdf(h,fullfile(datafolder,'avg_fitatt','coherence_sensitivity_marginals.pdf'));



%%%%%%%%%%%%%%
%% OLD CODE %%
%%%%%%%%%%%%%%

%% Plot average sensitivity across ROIs
cmap = brewermap(7,'PuOr');
h = figure; hold on

disc = {'p','c','m'};
col = {[0 0 0],cmap(2,:),cmap(6,:)};

for di = 1:3
    eval(sprintf('vc = rc_%s;',disc{di}));
    eval(sprintf('vm = rm_%s;',disc{di}));
    
    % INDIVIDUAL ROIs
    for ri = 1:8
        plot(vc(ri),vm(ri),'o','MarkerFaceColor',col{di}','MarkerEdgeColor','w','MarkerSize',3);
    end
end

for di = 1:3
    eval(sprintf('vc = rc_%s;',disc{di}));
    eval(sprintf('vm = rm_%s;',disc{di}));
    % ROI AVERAGE
    % plot the mean value
    vc_ci = bootci(1000,@median,vc);
    vm_ci = bootci(1000,@median,vm);
    % compute mu and sd
    vc_mu = mean(vc_ci); vc_sd = vc_ci(2)-vc_mu;
    vm_mu = mean(vm_ci); vm_sd = vm_ci(2)-vm_mu;
    % error bars
    plot([vc_mu vc_mu],[-vm_sd vm_sd]+vm_mu,'-','Color',col{di});
    plot([-vc_sd vc_sd]+vc_mu,[vm_mu vm_mu],'-','Color',col{di});
    % plot
    p(di) = plot(vc_mu,vm_mu,'o','MarkerFaceColor',col{di},'MarkerEdgeColor','w','MarkerSize',6);
end

set(gca,'XTick',0:.5:1.5,'YTick',0:.5:1.5);
% axis([0 1 0 3]);

xlabel('Contrast sensitivity');
ylabel('Coherence sensitivity');

l = legend(p,{'Passive viewing','Discriminate contrast','Discriminate coherence'},'FontSize',7,'FontName','Helvetica');
set(l,'box','off');

drawPublishAxis('figSize=[4.5 4.5]');

savepdf(h,fullfile(datafolder,'avg_fitatt','roi_sensitivity.pdf'));

%% Means and such for paper

%% Compute readout response functions

% RESPONSE TO FEAT | DISCRIMINATING
% First use the weights to run through
clear readout_p readout_c readout_m readout_a readout_l
for ni = 1:10
    for di = 1:2
        readout_p(ni,2,di,:) = squeeze(respcon(ncorrespond(ni),:,:))' * w_addpass(:,di);
        
        readout_p(ni,1,di,:) = squeeze(respcoh(ncorrespond(ni),:,:))' * w_addpass(:,di);
    end
        
%     readout_a(ni,2,2,:) = squeeze(respcon_ac(ni,:,:))' * ones(size(w_addpass(:,di)));%w_addpass(:,di);
%     readout_a(ni,2,1,:) = squeeze(respcon_am(ni,:,:))' * ones(size(w_addpass(:,di)));%w_addpass(:,di);
%     readout_a(ni,1,2,:) = squeeze(respcoh_ac(ni,:,:))' * ones(size(w_addpass(:,di)));%w_addpass(:,di);
%     readout_a(ni,1,1,:) = squeeze(respcoh_am(ni,:,:))' * ones(size(w_addpass(:,di)));%w_addpass(:,di);
%     
%     readout_l(ni,2,2,:) = squeeze(respcon_ac(ni,:,:))' * w_addpass(:,di);
%     readout_l(ni,2,1,:) = squeeze(respcon_am(ni,:,:))' * w_addpass(:,di);
%     readout_l(ni,1,2,:) = squeeze(respcoh_ac(ni,:,:))' * w_addpass(:,di);
%     readout_l(ni,1,1,:) = squeeze(respcoh_am(ni,:,:))' * w_addpass(:,di);
end

% avg
readout_p = squeeze(median(bootci(10000,@mean,readout_p)));
% readout_a = squeeze(median(bootci(1000,@mean,readout_a)));
% readout_l = squeeze(median(bootci(1000,@mean,readout_l)));

readout_p = readout_p - repmat(readout_p(:,:,1),1,1,1001);
% readout_a = readout_a - repmat(readout_a(:,:,1),1,1,1001);
% readout_l = readout_l - repmat(readout_l(:,:,1),1,1,1001);

%% Setup raw data values and multiply through the readout:

fname = fullfile(datafolder,'avg_deconatt','avg_deconEffects_att_cross.mat');
load(fname);

% avgdecon is now available for use

%% Convert allCon to betas

allResp = avgdecon.allResp;
beta = zeros(size(allResp,1),size(allResp,2),size(allResp,3));

load(fullfile(datafolder,'avg_hrf.mat'));

for ni = 1:10
    for ri = 1:8
        for ci = 1:32
            % get the right hrf, i.e. the average of the other observers
            pull = setdiff(1:11,ncorrespond(ni));
            ahrf = squeeze(mean(hrfs(pull,4,:)));
            % ahrf = mean(hrfs,1);
            ahrf = ahrf./max(ahrf);

            beta(ni,ri,ci) = ahrf\squeeze(allResp(ni,ri,ci,:));
        end
    end
end

%% separate data by task, then average over the two dimensions
attcon = beta(ni,ri,taskidx==2);
attcoh = beta(ni,ri,taskidx==1);

% init
betacon_ac = zeros(10,8,4);
betacon_am = betacon_ac;
betacoh_ac = betacon_ac;
betacoh_am = betacon_ac;

% contrast responses
for ci = 1:length(contrasts)
    ccon = contrasts(ci);

    % contrast responses = average across coherences
    idxs_con = logical((conidx==ccon).*(taskidx==2));
    idxs_coh = logical((conidx==ccon).*(taskidx==1));
    betacon_ac(:,:,ci) = mean(beta(:,:,idxs_con),3);
    betacon_am(:,:,ci) = mean(beta(:,:,idxs_coh),3);
end

% coherence responses
for mi = 1:length(coherences)
    ccoh = coherences(mi);
    
    % coherence responses = average across contrasts
    idxs_con = logical((cohidx==ccoh).*(taskidx==2));
    idxs_coh = logical((cohidx==ccoh).*(taskidx==1));
    betacoh_ac(:,:,mi) = mean(beta(:,:,idxs_con),3);
    betacoh_am(:,:,mi) = mean(beta(:,:,idxs_coh),3);
end

%% Multiply through the readout functions
% also remove intercepts
% % bc_ac = squeeze(mean(betacon_ac(:,:,:)));% bc_ac = bc_ac-bc_ac(1);
% % bc_am = squeeze(mean(betacon_am(:,:,:)));% bc_am = bc_am-bc_am(1);
% % bm_ac = squeeze(mean(betacoh_ac(:,:,:)));% bm_ac = bm_ac-bm_ac(1);
% % bm_am = squeeze(mean(betacoh_am(:,:,:)));% bm_am = bm_am-bm_am(1);

% plot model and readout

clear breadc_ac breadm_ac breadc_am breadm_am
for ni = 1:10
    % attend con
    tc = squeeze(betacon_ac(ni,:,:));
    for ri = 1:8
        tc(ri,:) = tc(ri,:) - interp1(contrasts-0.25,tc(ri,:),0,'linear','extrap');
    end
    breadc_ac(ni,:) = tc'*w_addpass(:,2);
    tm = squeeze(betacoh_ac(ni,:,:));
    for ri = 1:8
        tm(ri,:) = tm(ri,:) - interp1(coherences,tm(ri,:),0,'linear','extrap');
    end
    breadm_ac(ni,:) = tm'*w_addpass(:,2);
    % attend coh
    tc = squeeze(betacon_am(ni,:,:));
    for ri = 1:8
        tc(ri,:) = tc(ri,:) - interp1(contrasts-0.25,tc(ri,:),0,'linear','extrap');
    end
    breadc_am(ni,:) = tc'*w_addpass(:,1);
    tm = squeeze(betacoh_am(ni,:,:));
    for ri = 1:8
        tm(ri,:) = tm(ri,:) - interp1(coherences,tm(ri,:),0,'linear','extrap');
    end
    breadm_am(ni,:) = tm'*w_addpass(:,1);
end


%% Averages and SDs

bc_ac_ci = bootci(10000,@mean,breadc_ac);
bc_ac = mean(bc_ac_ci);
bc_am_ci = bootci(10000,@mean,breadc_am);
bc_am = mean(bc_am_ci);
bm_ac_ci = bootci(10000,@mean,breadm_ac);
bm_ac = mean(bm_ac_ci);
bm_am_ci = bootci(10000,@mean,breadm_am);
bm_am = mean(bc_am_ci);

%% Test plot new beta values with SnR readout functions
h = figure; hold on
% CONTRAST
% subplot(121); hold on
x = 0:.001:1;

% estimate scale factor
plot(x,squeeze(readout_p(2,2,:)),'-','Color',cmap(2,:));
plot(x,squeeze(readout_p(2,1,:)),'-','Color',[0.8 0.8 0.8]);
plot(contrasts-0.25,bc_ac,'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','w');
errbar(contrasts-0.25,bc_ac,bc_ac_ci(2,:)-bc_ac,'-','Color',cmap(2,:));
plot(contrasts-0.25+.01,bc_am,'o','MarkerFaceColor',[0.8 0.8 0.8],'MarkerEdgeColor','w');
errbar(contrasts-0.25+.01,bc_am,bc_am_ci(2,:)-bc_am,'-','Color',[0.8 0.8 0.8]);
a = axis;
axis([0 0.75 a(3) a(4)]);
set(gca,'XTick',[0 0.5]','YTick',[-10:10:10]);
xlabel('Contrast (%)');
ylabel('Readout (s.d.)');
l = legend({'Discriminating contrast','Discriminating coherence'},'FontSize',7,'FontName','Helvetica');
set(l,'box','off');

drawPublishAxis('figSize=[4.5 4.5]');

savepdf(h,fullfile(datafolder,'avg_fitatt','readout_betas_con.pdf'));

% COHERENCE
h = figure; hold on

plot(x,squeeze(readout_p(1,1,:)),'-','Color',cmap(6,:));
plot(x,squeeze(readout_p(1,2,:)),'Color',[0.8 0.8 0.8]);
plot(coherences,bm_ac,'o','MarkerFaceColor',[0.8 0.8 0.8],'MarkerEdgeColor','w');
errbar(coherences,bm_ac,bm_ac_ci(2,:)-bm_ac,'-','Color',[0.8 0.8 0.8]);
plot(coherences+.01,bm_am,'o','MarkerFaceColor',cmap(6,:),'MarkerEdgeColor','w');
errbar(coherences+.01,bm_am,bm_am_ci(2,:)-bm_am,'-','Color',cmap(6,:));
legend({'Discriminating coherence','Discriminating contrast'},'FontSize',7,'FontName','Helvetica');
set(l,'box','off');
a = axis;
axis([0 0.75 a(3) a(4)]);
set(gca,'XTick',[0 0.5]','YTick',[-10:10:10]);
xlabel('Coherence (%)');
ylabel('Readout (s.d.)');
drawPublishAxis('figSize=[4.5 4.5]');

savepdf(h,fullfile(datafolder,'avg_fitatt','readout_betas_coh.pdf'));

%% Plot


conds = {'motion','contrast'};
group = {'p','a','a'};
rois = {'V1','MT'};
% split discrimination condition by subplot
h = figure;

subplot(2,2,1); hold on

title('Readout of contrast predicted by linking model');
plot(x,squeeze(readout_p(2,2,:)),'-','Color',cmap(2,:));
plot(x,squeeze(readout_p(2,1,:)),'-','Color',[0.8 0.8 0.8]);
l=legend({'Contrast attended','Contrast unattended'},'FontSize',7,'FontName','Helvetica');
set(l,'box','off');

axis([0 .75 -5 30]);
set(gca,'XTick',[0 .75],'XTickLabel',[0 1],'YTick',0:10:30);
xlabel('Contrast (%)');
ylabel('Readout response (s.d.)');
drawPublishAxis('figSize=[8.9, 8.9]');

% % % % % % % % % % % % % % % % 
subplot(2,2,2); hold on
title('Readout of coherence predicted by linking model');
plot(x,squeeze(readout_p(1,1,:)),'-','Color',cmap(6,:));
plot(x,squeeze(readout_p(1,2,:)),'Color',[0.8 0.8 0.8]);
l=legend({'Coherence attended','Coherence unattended'},'FontSize',7,'FontName','Helvetica');
set(l,'box','off');

axis([0 .75 -5 30]);
set(gca,'XTick',[0 .75],'XTickLabel',[0 1],'YTick',0:10:30);
xlabel('Coherence (%)');
ylabel('Readout response (s.d.)');
drawPublishAxis('figSize=[8.9, 8.9]');

% % % % % % % % % % % % % % % % 
subplot(2,2,3); hold on
title('Readout of contrast (from active viewing)');
% plot(x,squeeze(readout_a(2,2,:)),'-k');
% plot(x,squeeze(readout_l(2,1,:)),'Color',[0.8 0.8 0.8]);
plot(x,squeeze(readout_l(2,2,:)),'-','Color',cmap(2,:));
plot(x,squeeze(readout_l(2,1,:)),'Color',[0.8 0.8 0.8]);

axis([0 .75 -5 30]);
set(gca,'XTick',[0 .75],'XTickLabel',[0 1],'YTick',0:10:30);
xlabel('Contrast (%)');
ylabel('Readout response (s.d.)');
drawPublishAxis('figSize=[8.9, 8.9]');

% % % % % % % % % % % % % % % % 
subplot(2,2,4); hold on
title('Readout of coherence (from active viewing)');
% plot(x,squeeze(readout_a(1,1,:)),'-k');
plot(x,squeeze(readout_l(1,2,:)),'Color',[0.8 0.8 0.8]);
plot(x,squeeze(readout_l(1,1,:)),'-','Color',cmap(6,:));

axis([0 .75 -5 30]);
set(gca,'XTick',[0 .75],'XTickLabel',[0 1],'YTick',0:10:30);
xlabel('Coherence (%)');
ylabel('Readout response (s.d.)');
drawPublishAxis('figSize=[8.9, 8.9]');

savepdf(h,fullfile(datafolder,'avg_fitatt','readout.pdf'));


%% READOUT 2: JUST V1/MT FROM ATTENTION DATA

% RESPONSE TO FEAT | DISCRIMINATING
% First use the weights to run through
clear readout_p readout_c readout_m readout_a readout_l

for ni = 1:10
    for di = 1:2
        readout_p(ni,2,di,:) = squeeze(respcon(ncorrespond(ni),:,:))' * w_addpass(:,di);
        
        readout_p(ni,1,di,:) = squeeze(respcoh(ncorrespond(ni),:,:))' * w_addpass(:,di);
    end
    readout_l(ni,2,2,:) = squeeze(respcon_ac(ni,[1 8],:))' * w(ni,:)';
    readout_l(ni,2,1,:) = squeeze(respcon_am(ni,[1 8],:))' * w(ni,:)';
    readout_l(ni,1,2,:) = squeeze(respcoh_ac(ni,[1 8],:))' * w(ni,:)';
    readout_l(ni,1,1,:) = squeeze(respcoh_am(ni,[1 8],:))' * w(ni,:)';
end

% avg
readout_p = squeeze(median(bootci(1000,@mean,readout_p)));
readout_l = squeeze(median(bootci(1000,@mean,readout_l)));

readout_p = readout_p - repmat(readout_p(:,:,1),1,1,1001);
readout_l = readout_l - repmat(readout_l(:,:,1),1,1,1001);

%% Plot


conds = {'motion','contrast'};
group = {'p','a','a'};
rois = {'V1','MT'};
% split discrimination condition by subplot
h = figure;

subplot(2,2,1); hold on

title('Readout of contrast predicted by linking model');
plot(x,squeeze(readout_p(2,2,:)),'-','Color',cmap(2,:));
plot(x,squeeze(readout_p(2,1,:)),'-','Color',[0.8 0.8 0.8]);
l=legend({'Contrast attended','Contrast unattended'},'FontSize',7,'FontName','Helvetica');
set(l,'box','off');

axis([0 .75 -5 30]);
set(gca,'XTick',[0 .75],'XTickLabel',[0 1],'YTick',0:10:30);
xlabel('Contrast (%)');
ylabel('Readout response (s.d.)');
drawPublishAxis('figSize=[8.9, 8.9]');

% % % % % % % % % % % % % % % % 
subplot(2,2,2); hold on
title('Readout of coherence predicted by linking model');
plot(x,squeeze(readout_p(1,1,:)),'-','Color',cmap(6,:));
plot(x,squeeze(readout_p(1,2,:)),'Color',[0.8 0.8 0.8]);
l=legend({'Coherence attended','Coherence unattended'},'FontSize',7,'FontName','Helvetica');
set(l,'box','off');

axis([0 .75 -5 30]);
set(gca,'XTick',[0 .75],'XTickLabel',[0 1],'YTick',0:10:30);
xlabel('Coherence (%)');
ylabel('Readout response (s.d.)');
drawPublishAxis('figSize=[8.9, 8.9]');

% % % % % % % % % % % % % % % % 
subplot(2,2,3); hold on
title('Readout of contrast (from active viewing)');
% plot(x,squeeze(readout_a(2,2,:)),'-k');
% plot(x,squeeze(readout_l(2,1,:)),'Color',[0.8 0.8 0.8]);
plot(x,squeeze(readout_l(2,2,:)),'-','Color',cmap(2,:));
plot(x,squeeze(readout_l(2,1,:)),'Color',[0.8 0.8 0.8]);

axis([0 .75 -5 30]);
set(gca,'XTick',[0 .75],'XTickLabel',[0 1],'YTick',0:10:30);
xlabel('Contrast (%)');
ylabel('Readout response (s.d.)');
drawPublishAxis('figSize=[8.9, 8.9]');

% % % % % % % % % % % % % % % % 
subplot(2,2,4); hold on
title('Readout of coherence (from active viewing)');
% plot(x,squeeze(readout_a(1,1,:)),'-k');
plot(x,squeeze(readout_l(1,2,:)),'Color',[0.8 0.8 0.8]);
plot(x,squeeze(readout_l(1,1,:)),'-','Color',cmap(6,:));

axis([0 .75 -5 30]);
set(gca,'XTick',[0 .75],'XTickLabel',[0 1],'YTick',0:10:30);
xlabel('Coherence (%)');
ylabel('Readout response (s.d.)');
drawPublishAxis('figSize=[8.9, 8.9]');

savepdf(h,fullfile(datafolder,'avg_fitatt','readout_one2.pdf'));
