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

% Load the sensitivity data for V1->MT during contrast and coherence
% discrimination
load(fullfile(datafolder,'avg_att_cross_fits.mat'));
% single beta version: 
% load(fullfile(datafolder,'avg_att_cross_fits_sb.mat'));

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

%% Separate intercepts (we mostly care about the slope changes)

datasets = {'respcon','respcoh','respcon_am','respcon_ac','respcoh_am','respcoh_ac'};

for ri = 1:8
    for di = 1:length(datasets)
        eval(sprintf('%s_offset(:,ri) = %s(:,ri,1);',datasets{di},datasets{di}));
        eval(sprintf('%s(:,ri,:) = %s(:,ri,:) - repmat(%s(:,ri,1),1,1,1001);',datasets{di},datasets{di},datasets{di}));
    end
end

offset = respcon_offset;
offset_am = respcon_am_offset;
offset_ac = respcon_ac_offset;

%% Sensory space sensitivity plots

% normalize functions
norm_con = respcon(:,1,:);
norm_coh = respcoh(:,8,:);

norm_conatt = norm_con;%respcon_ac(:,1,:);
norm_cohatt = norm_coh; %respcoh_am(:,8,:);

clear rc_p rc_c rc_m rm_p rm_c rm_m
ros = [1:8];
x = 0:.001:1;

X = [ones(size(x')) x'];

for rii = 1:length(ros)
    ri = ros(rii);
    
    for ni = 1:10
        b = X\squeeze(respcon(ncorrespond(ni),ri,:));
        rc_p(ni,rii) = b(2);
        b = X\squeeze(respcon_ac(ni,ri,:));
        rc_c(ni,rii) = b(2);
        b = X\squeeze(respcon_am(ni,ri,:));
        rc_m(ni,rii) = b(2);
        
        b = X\squeeze(respcoh(ncorrespond(ni),ri,:));
        rm_p(ni,rii) = b(2);
        b = X\squeeze(respcoh_ac(ni,ri,:));
        rm_c(ni,rii) = b(2);
        b = X\squeeze(respcoh_am(ni,ri,:));
        rm_m(ni,rii) = b(2);
        
%         rc_p(ni,rii) = norm_con(ncorrespond(ni),:)'\squeeze(respcon(ncorrespond(ni),ri,:));
%         rc_c(ni,rii) = norm_conatt(ni,:)'\squeeze(respcon_ac(ni,ri,:));
%         rc_m(ni,rii) = norm_conatt(ni,:)'\squeeze(respcon_am(ni,ri,:));
% 
%         rm_p(ni,rii) = norm_coh(ncorrespond(ni),:)'\squeeze(respcoh(ncorrespond(ni),ri,:));
%         rm_c(ni,rii) = norm_cohatt(ni,:)'\squeeze(respcoh_ac(ni,ri,:));
%         rm_m(ni,rii) = norm_cohatt(ni,:)'\squeeze(respcoh_am(ni,ri,:));
    end
end

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

%% Text for paper, take differences
adrc = rc_c(:)-rc_m(:);
mu = mean(adrc);
ci = bootci(10000,@mean,adrc);
disp(sprintf('Average difference in contrast sensitivity between tasks: %1.2f, 95%% CI [%1.2f, %1.2f]',mu,ci(1),ci(2)));

adrm = rm_m(:)-rm_c(:);
mu = mean(adrm);
ci = bootci(10000,@mean,adrm);
disp(sprintf('Average difference in coherence sensitivity between tasks: %1.2f, 95%% CI [%1.2f, %1.2f]',mu,ci(1),ci(2)));


%%
drc = rc_c-rc_m;
drm = rm_m-rm_c;
drc_ = mean(drc);
drc_ci = bootci(1000,@mean,drc);
drm_ = mean(drm);
drm_ci = bootci(1000,@mean,drm);

for ri = 1:8
    disp(sprintf('%s: %1.2f (95%% CI [%1.2f, %1.2f]); ',rois{ri},drc_(ri),drc_ci(1,ri),drc_ci(2,ri)));
end

for ri = 1:8
    disp(sprintf('%s: %1.2f (95%% CI [%1.2f, %1.2f]); ',rois{ri},drm_(ri),drm_ci(1,ri),drm_ci(2,ri)));
end


%% offsets
% offset_ac(offset_ac<0) = -10;
% offset_am(offset_am<0) = -10;
o_p_ci = bootci(1000,@nanmedian,offset);
o_p = squeeze(median(offset));
o_c_ci = bootci(1000,@nanmedian,offset_ac);
o_c = squeeze(median(offset_ac));
o_m_ci = bootci(1000,@nanmedian,offset_am);
o_m = squeeze(median(offset_am));

%% Offset text for paper

o_all = cat(3,offset_ac,offset_am);
o_all = mean(o_all,3);
o_a = squeeze(mean(o_all));
o_a_ci = bootci(10000,@nanmean,o_all);

% average offset
for ri = 1:8
    fprintf('%s %1.2f, ',rois{ri},o_a(ri));
%     fprintf('%s %1.2f, 95%% CI [%1.2f, %1.2f], ',rois{ri},o_a(ri),o_a_ci(1,ri),o_a_ci(2,ri));
end
% average difference in offset

o_diff = offset_ac-offset_am;


%% Plot the offset info

bsize = 0.15;

h = figure;  hold on
cmap = brewermap(7,'PuOr');
rois = {'V1','V2','V3','V3A','V3B','V4','V7','MT'};

% add error bars
errbar((1:8)-bsize,o_c,o_c_ci(2,:)-o_c,'-','Color',cmap(2,:));
errbar((1:8)+bsize,o_m,o_m_ci(2,:)-o_m,'-','Color',cmap(6,:));
for ri = 1:8
    p(1) = plot(ri-bsize,o_c(ri),'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','w','MarkerSize',4);
    p(2) = plot(ri+bsize,o_m(ri),'o','MarkerFaceColor',cmap(6,:),'MarkerEdgeColor','w','MarkerSize',4);
%     p(1) = bar(ri-bsize,rc_c(ri),bsize*2,'FaceColor',cmap(2,:),'EdgeColor','w');
%     p(2) = bar(ri+bsize,rc_m(ri),bsize*2,'FaceColor',cmap(6,:),'EdgeColor','w');
end

axis([-0.1 8.1 -1.25 2.25]); 

set(gca,'XTick',[1:8],'XTickLabel',rois([1:8]),'YTick',[0 0.5 1.0]);
ylabel('Offset (% signal change)');
drawPublishAxis('figSize=[5.5, 4.5]');

savepdf(h,fullfile(datafolder,'avg_fitatt','contrast_offset_dashed.pdf'));


%% New plot: dashed line surrounding 
bsize = 0.15;

h = figure;  hold on
cmap = brewermap(7,'PuOr');
rois = {'V1','V2','V3','V3A','V3B','V4','V7','MT'};

% add error bars
errbar((1:8)-bsize,rc_c,rc_c_ci(2,:)-rc_c,'-','Color',cmap(2,:));
errbar((1:8)+bsize,rc_m,rc_m_ci(2,:)-rc_m,'-','Color',cmap(6,:));
for ri = 1:8
    p(1) = plot(ri-bsize,rc_c(ri),'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','w','MarkerSize',4);
    p(2) = plot(ri+bsize,rc_m(ri),'o','MarkerFaceColor',cmap(6,:),'MarkerEdgeColor','w','MarkerSize',4);
%     p(1) = bar(ri-bsize,rc_c(ri),bsize*2,'FaceColor',cmap(2,:),'EdgeColor','w');
%     p(2) = bar(ri+bsize,rc_m(ri),bsize*2,'FaceColor',cmap(6,:),'EdgeColor','w');
end

axis([-0.1 8.1 -1.25 3.25]); 
set(gca,'XTick',[1:8],'XTickLabel',rois([1:8]),'YTick',0:1:2);

ylabel('Relative contrast sensitivity');
% legend(p,{'Discriminating contrast','Discriminating coherence'});
drawPublishAxis('figSize=[5.5, 4.5]');

% savepdf(h,fullfile(datafolder,'avg_fitatt','contrast_sensitivity_dashed.pdf'));

h = figure;  hold on
cmap = brewermap(7,'PuOr');

% add error bars
errbar((1:8)-bsize,rm_c,rm_c_ci(2,:)-rm_c,'-','Color',cmap(2,:));
errbar((1:8)+bsize,rm_m,rm_m_ci(2,:)-rm_m,'-','Color',cmap(6,:));
for ri = 1:8
    p(1) = plot(ri-bsize,rm_c(ri),'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','w','MarkerSize',4);
    p(2) = plot(ri+bsize,rm_m(ri),'o','MarkerFaceColor',cmap(6,:),'MarkerEdgeColor','w','MarkerSize',4);
%     p(1) = bar(ri-bsize,rm_c(ri),bsize*2,'FaceColor',cmap(2,:),'EdgeColor','w');
%     p(2) = bar(ri+bsize,rm_m(ri),bsize*2,'FaceColor',cmap(6,:),'EdgeColor','w');
end

% add the dashed passive conditions

axis([-0.1 8.1 -1.25 3.25]); 
set(gca,'XTick',[1:8],'XTickLabel',rois([1:8]),'YTick',0:1:2);

ylabel('Relative coherence sensitivity');
% legend(p,{'Discriminating contrast','Discriminating coherence'});
drawPublishAxis('figSize=[5.5, 4.5]');

% savepdf(h,fullfile(datafolder,'avg_fitatt','coherence_sensitivity_dashed.pdf'));
