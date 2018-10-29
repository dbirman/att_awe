
%%
%sids = {'s043','s300','s305','s329','s025'};
nSIDs = [305 329 43 25 300 346 343 344 338 340 348];
bSIDs = [345 330 337 335 349 354 353 356 334 352]; % behavior only participants
aSIDs = [nSIDs bSIDs];
sids = {};
for bi = 1:length(aSIDs)
    sids{end+1} = sprintf('s%03.0f',aSIDs(bi));
end
rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};

%% Ages + gender

% nsids followed by bsids
age = [30 26 27 33 27 26 36 24 25 20 19]-1;
bage = [28 23 29 30 23 37 24 26 56 34]-1;
% f/m
gender = [1 1 1 2 2 1 1 2 1 1 1];
bgender = [2 2 1 2 2 1 2 1 1 1];

%% Statistics

trials = zeros(1,length(sids));
for si = 1:length(sids);
    subj = sids{si};
    
    % Get files
    files = dir(fullfile(datafolder,subj));

    % Load data
    adata = loadadata(subj);
    
    trials(si) = size(adata,1);
end

%% Fit cumulative gaussian functions to data
avg_rightperformancecontrol;


%% Fit behavioral models
parfor si = 1:length(sids)
    subj_behav_analysis(sids{si},'refit');
    close all
end

%% Re-compute thresholds
parfor si = 1:length(sids)
    subj_behav_analysis(sids{si},'thresholds');
    close all
end

%% Plot threshold plots
control = zeros(length(sids),2,4);
attend = zeros(length(sids),2);
unattend = zeros(length(sids),2);
for si = 1:length(sids)
    load(fullfile(datafolder,sprintf('%s_data.mat',sids{si})));
    control(si,:,:) = data.control;
    attend(si,:) = data.attend;
    unattend(si,:) = data.unattend;
end

control(control<=0) = nan;
attend(attend<=0) = nan;
unattend(unattend<=0) = nan;
    
control_ = squeeze(bootci(1000,@nanmedian,control));
control_m = squeeze(mean(control_));
control_s = squeeze(control_(2,:,:))-control_m;

attend_ = squeeze(bootci(1000,@nanmedian,attend));
attend_m = squeeze(mean(attend_));
attend_s = squeeze(attend_(2,:,:))-attend_m;

unattend_ = squeeze(bootci(1000,@nanmedian,unattend));
unattend_m = squeeze(mean(unattend_));
unattend_s = squeeze(unattend_(2,:,:))-unattend_m;

%% Average threshold increase
ratio = unattend ./ attend;
r_mu = nanmean(ratio);
r_ci = bootci(10000,@nanmean,ratio);
disp('Ratio increase');

disp(sprintf('%1.2f, 95%% CI [%1.2f, %1.2f]',r_mu(1),r_ci(1,1),r_ci(2,1)));
disp(sprintf('%1.2f, 95%% CI [%1.2f, %1.2f]',r_mu(2),r_ci(1,2),r_ci(2,2)));


mu = nanmean(ratio(:));
ci = bootci(10000,@nanmean,ratio(:));

disp(sprintf('%1.2f, 95%% CI [%1.2f, %1.2f]',mu,ci(1),ci(2)));

%% Check if the scanning subjects are better than the non-scanning subjects
scan = control(1:11,:,:);
non = control(12:end,:,:);

mu_s = nanmean(scan);
ci_s = bootci(1000,@nanmean,scan);

mu_n = nanmean(non);
ci_n = bootci(1000,@nanmean,non);

clear p h stats
for task=1:2
    for base=1:4
        [p(task,base),h(task,base),stats{task,base}] = ranksum(scan(:,task,base),non(:,task,base));
    end
end

%% Disp info for paper

contrast = [0.325 0.4 0.55 0.85];
coherence = [0.15 0.3 0.45 0.6];

group = {'coherence','contrast'};
for ci = 1:2
    for b = 1:4
        disp(sprintf('%s base %2.0f\\%%, JND = %2.1f\\%% 95\\%% CI [%2.1f, %2.1f]',group{ci},eval(sprintf('%s(b)',group{ci}))*100,control_m(ci,b)*100,control_(1,ci,b)*100,control_(2,ci,b)*100));
    end
end

%% Get weber law slopes

for bi = 1:21
    cont = squeeze(control(bi,:,:));
    
    conslope(bi) = log(contrasts*100)'\log(cont(2,:)*100)';
    cohslope(bi) = log(coherences*100)'\log(cont(1,:)*100)';
end

nanmean(conslope)
bootci(1000,@nanmean,conslope)

nanmean(cohslope)
bootci(1000,@nanmean,cohslope)

%% Compare catch trial performance to control and control to attend
control_attend = control(:,:,2);

diff_attend = control_attend-attend;
diff_unattend = control_attend-unattend;

mean(diff_attend)
bootci(1000,@nanmean,diff_attend)


nanmean(diff_unattend)
bootci(1000,@nanmean,diff_unattend)

%% actual plot (only behavior -- no model fits, note that we don't use this for the paper)

% plot the contrast data
h = figure; hold on
cmap = brewermap(7,'PuOr');

contrast = [0.325 0.4 0.55 0.85];
coherence = [0.15 0.3 0.45 0.6];

errbar(contrast,control_m(2,:),control_s(2,:),'-','Color',cmap(2,:));
plot(contrast,control_m(2,:),'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','w','MarkerSize',7);
errbar(coherence,control_m(1,:),control_s(1,:),'-','Color',cmap(6,:));
plot(coherence,control_m(1,:),'o','MarkerFaceColor',cmap(6,:),'MarkerEdgeColor','w','MarkerSize',7);




% savepdf(h,fullfile(datafolder,'avg_behav','thresholds_coherence.pdf'));

plot(contrast(2),unattend_m(2),'o','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor','w','MarkerSize',7);
errbar(contrast(2),unattend_m(2),unattend_s(2),'-','Color',cmap(1,:));
plot(coherence(2),unattend_m(1),'o','MarkerFaceColor',cmap(7,:),'MarkerEdgeColor','w','MarkerSize',7);
errbar(coherence(2),unattend_m(1),unattend_s(1),'-','Color',cmap(7,:));

% legend('Contrast','Coherence');
xlabel('Base stimulus (%)');
axis([0 1 0 0.6]);
set(gca,'XTick',[0 0.25 0.5 0.75 1]','XTickLabel',[0 25 50 75 100],'YTick',0:0.2:.6,'YTickLabel',0:20:60);
ylabel('Threshold (%)');



drawPublishAxis('figSize=[4.5,4.5]');

savepdf(h,fullfile(datafolder,'avg_behav','thresholds.pdf'));

%% Probit (just a plot, for figures)

h = figure;

x = -2:.01:2;
plot(x,normcdf(x));
xlabel('Trial effect (s.d.)');
ylabel('P(R)');
title('Cumulative gaussian');

savepdf(h,fullfile(datafolder,'avg_behav','probit.pdf'),'figsize=[5,4]');

%% Check for high difficulty confound

for si = 1:length(sids);
    subj = sids{si};
    
    % Load data
    adata = loadadata(subj);
end