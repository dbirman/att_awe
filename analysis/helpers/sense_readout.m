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
    fit = sfits{si}{4}; % 4 refers to resp-25, which is our standard model

    for ri = 1:8
        respcon(si,ri,:) = conModel(x,fit.roifit{ri}.params)+fit.roifit{ri}.params.offset;
        respcoh(si,ri,:) = cohModel(x,fit.roifit{ri}.params)+fit.roifit{ri}.params.offset;
    end
end

% rc_passive = squeeze(mean(bootci(10000,@mean,respcon)));
% rm_passive = squeeze(mean(bootci(10000,@mean,respcoh)));

% Load the sensitivity data for V1->MT during contrast and coherence
% discrimination
load(fullfile(datafolder,'avg_att_fits_sb.mat'));

x = 0:.001:1;
respcon_ac = zeros(10,8,length(x));
respcon_am = zeros(10,8,length(x));
respcoh_ac = zeros(10,8,length(x));
respcoh_am = zeros(10,8,length(x));

for si = 1:10
    fit = afits{si}{1}; % 4 refers to resp-25, which is our standard model

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
                weights_passive(bi,ai,ri,ci) = afits{ai}{bi}.params.(sprintf('beta_control_%s_%s',ROIs{ri},cons{ci}));
            end
        end
    end
end

weights_passive = squeeze(weights_passive(1,:,:,:));
w_addpass = squeeze(mean(bootci(10000,@mean,weights_passive)));

%% Load: weights 2

% Load the readout weights from the linking model
restructure_afits_2;

weights_passive = zeros(2,length(aSIDs),2,2);
ROIs = {'V1','MT'};
cons = {'cohw','conw'};
for bi = 1:2
    for ai = 1:length(aSIDs)
        for ri = 1:2
            for ci = 1:2
                weights_passive(bi,ai,ri,ci) = afits{ai}{bi}.params.(sprintf('beta_control_%s_%s',ROIs{ri},cons{ci}));
            end
        end
    end
end

weights_passive = squeeze(weights_passive(1,:,:,:));
w_addpass = squeeze(mean(bootci(10000,@mean,weights_passive)));

%% Sensory space sensitivity plots

% normalize functions
norm_con = respcon(:,1,:);
norm_coh = respcoh(:,8,:);

norm_conatt = norm_con;%respcon_ac(:,1,:);
norm_cohatt = norm_coh; %respcoh_am(:,8,:);

clear rc_p rc_c rc_m rm_p rm_c rm_m
ros = [1:8];
for rii = 1:length(ros)
    ri = ros(rii);
    
    for ni = 1:10
        rc_p(ni,rii) = norm_con(ncorrespond(ni),:)'\squeeze(respcon(ncorrespond(ni),ri,:));
        rc_c(ni,rii) = norm_conatt(ni,:)'\squeeze(respcon_ac(ni,ri,:));
        rc_m(ni,rii) = norm_conatt(ni,:)'\squeeze(respcon_am(ni,ri,:));

        rm_p(ni,rii) = norm_coh(ncorrespond(ni),:)'\squeeze(respcoh(ncorrespond(ni),ri,:));
        rm_c(ni,rii) = norm_cohatt(ni,:)'\squeeze(respcoh_ac(ni,ri,:));
        rm_m(ni,rii) = norm_cohatt(ni,:)'\squeeze(respcoh_am(ni,ri,:));
    end
end

%% Multiply by SNR

% rc_c = rc_c .* sqrt(squeeze(snr(:,2,:)));
% rc_m = rc_m .* sqrt(squeeze(snr(:,2,:)));
% 
% rm_c = rm_c .* sqrt(squeeze(snr(:,1,:)));
% rm_m = rm_m .* sqrt(squeeze(snr(:,1,:)));

%% Average across subjects
rc_p_ci = bootci(1000,@mean,rc_p);
rc_c_ci = bootci(1000,@mean,rc_c);
rc_m_ci = bootci(1000,@mean,rc_m);
rm_p_ci = bootci(1000,@mean,rm_p);
rm_c_ci = bootci(1000,@mean,rm_c);
rm_m_ci = bootci(1000,@mean,rm_m);

group = {'rc','rm'};
disc = {'p','c','m'};
for gi = 1:2
    for di = 1:3
        eval(sprintf('%s_%s_i = %s_%s;',group{gi},disc{di},group{gi},disc{di})); % save individual
        eval(sprintf('%s_%s = mean(%s_%s_ci);',group{gi},disc{di},group{gi},disc{di}));
    end
end

%% Average over ROIs
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

set(gca,'XTick',[0 1],'YTick',[0 1]);
axis([0 1 0 3]);

xlabel('Normalized contrast sensitivity');
ylabel('Normalized coherence sensitivity');

l = legend(p,{'Passive viewing','Discriminate contrast','Discriminate coherence'},'FontSize',7,'FontName','Helvetica');
set(l,'box','off');

drawPublishAxis('figSize=[4.5 4.5]');

savepdf(h,fullfile(datafolder,'avg_fitatt','roi_sensitivity.pdf'));

%% Compute readout response functions

% RESPONSE TO FEAT | DISCRIMINATING
% First use the weights to run through
clear readout_p readout_c readout_m readout_a readout_l
for ni = 1:10
    for di = 1:2
        readout_p(ni,2,di,:) = squeeze(respcon(ncorrespond(ni),:,:))' * w_addpass(:,di);
        
        readout_p(ni,1,di,:) = squeeze(respcoh(ncorrespond(ni),:,:))' * w_addpass(:,di);
    end
        
    readout_a(ni,2,2,:) = squeeze(respcon_ac(ni,:,:))' * ones(size(w_addpass(:,di)));%w_addpass(:,di);
    readout_a(ni,2,1,:) = squeeze(respcon_am(ni,:,:))' * ones(size(w_addpass(:,di)));%w_addpass(:,di);
    readout_a(ni,1,2,:) = squeeze(respcoh_ac(ni,:,:))' * ones(size(w_addpass(:,di)));%w_addpass(:,di);
    readout_a(ni,1,1,:) = squeeze(respcoh_am(ni,:,:))' * ones(size(w_addpass(:,di)));%w_addpass(:,di);
    
    readout_l(ni,2,2,:) = squeeze(respcon_ac(ni,:,:))' * w_addpass(:,di);
    readout_l(ni,2,1,:) = squeeze(respcon_am(ni,:,:))' * w_addpass(:,di);
    readout_l(ni,1,2,:) = squeeze(respcoh_ac(ni,:,:))' * w_addpass(:,di);
    readout_l(ni,1,1,:) = squeeze(respcoh_am(ni,:,:))' * w_addpass(:,di);
end

% avg
readout_p = squeeze(median(bootci(1000,@mean,readout_p)));
readout_a = squeeze(median(bootci(1000,@mean,readout_a)));
readout_l = squeeze(median(bootci(1000,@mean,readout_l)));

readout_p = readout_p - repmat(readout_p(:,:,1),1,1,1001);
readout_a = readout_a - repmat(readout_a(:,:,1),1,1,1001);
readout_l = readout_l - repmat(readout_l(:,:,1),1,1,1001);
% readout_m = readout_m - repmat(readout_m(:,:,1),1,1,1001);

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
