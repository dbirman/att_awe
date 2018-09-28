
afits = restructure_afits('avg_indiv_fits');

%% Actual model fits
% x = 0:.001:1;

coherences = [0.15 0.3 0.45 0.6];
contrasts = [0.325 0.4 0.55 0.85];
% coherences = [0.25 0.5 0.75 1];
% contrasts = [0.25 0.5 0.75 1];

basecoh = 0; basecon = 0.25;
conbins = 0:.01:1;
conbins_ = conbins(1:(length(conbins)-1))+(diff(conbins)/2);
cohbins = 0:.01:1;
cohbins_ = cohbins(1:(length(cohbins)-1))+(diff(cohbins)/2);

% subj | fmri | behav | disc | feat | ped | effect
%   1      2     3       4      5       6     7
model = zeros(length(aSIDs),length(bmodels),length(ropts),2,2,4,length(conbins_));
athresholds = zeros(length(aSIDs),length(bmodels),length(ropts),101,2);

disppercent(-1/length(aSIDs));
for ai = 1:length(aSIDs) % which subject
    %     for mi = 1:length(mopts) % which fMRI model
    
    allROIs = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
    
    cur = model(ai,:,:,:,:,:,:);
    thresholds = athresholds(ai,:,:,:,:);
    
    for bi = 1:length(bmodels)
        for ri = 1:length(ropts)
            rois = allROIs(ropts{ri});

            % get betas
            cons = {'cohw','conw'};
            betas = zeros(length(rois),2);

            for cri = 1:length(rois)
                for ci = 1:2
                    betas(cri,ci) = afits{ai}{bi,ri}.params.(sprintf('beta_control_%s_%s',rois{cri},cons{ci}));
                end
            end

            bias = afits{ai}{bi}.params.bias;
            lapse = afits{ai}{bi}.params.lapse;
            
            threshold = zeros(101,2);
            cthreshold = zeros(101,2);
            
            for ci = 1:2
                
                % compute the threshold at every contrast/coherence
                % increment
                flip = [2 1];
                % first compute the beta-weighted response functions
                tresp_con = abs(betas(:,ci)'*respcon_(ropts{ri},:));
                tresp_coh = abs(betas(:,ci)'*respcoh_(ropts{ri},:));
                
                peds = 0:.01:1;
                for pi = 1:length(peds)
                    ped = peds(pi);
                    % now compute the discrimination threshold for the current
                    % condition at the current pedestal
                    if ci==1
                        baseResp = interp1(x,tresp_coh,ped);
                        thresholdResp = baseResp + 1;
                        % now find the contrast value that corresponds to this
                        threshold(pi,ci) = interp1(tresp_coh,x,thresholdResp) - ped;
                        
                        % also compute the catch threshold for the other
                        % task
                        cthreshold(pi,flip(ci)) = 1 / (x'\tresp_con');
                    else
                        baseResp = interp1(x,tresp_con,ped);
                        thresholdResp = baseResp + 1;
                        % now find the contrast value that corresponds to this
                        threshold(pi,ci) = interp1(tresp_con,x,thresholdResp) - ped;
                        
                        cthreshold(pi,flip(ci)) = 1 / (x'\tresp_coh');
                    end
                end
                
                for ped = 1:4 % which pedestal
                    conped = contrasts(ped)-basecon;
                    cohped = coherences(ped)-basecoh;
                    

                    for bii = 1:length(conbins_) % which bin
                        conbin = conbins_(bii);
                        cohbin = cohbins_(bii);
                        % best for loop chain I've ever written
                        % compute probability of a right response for this
                        % set of conditions

                        % using respcon_ and respcoh_

                        conp = zeros(length(rois),1);
                        coni = conp;
                        cohp = conp;
                        cohi = conp;
                        
                        % get con resp
                        for rii = 1:length(rois)
                            con = respcon_(ropts{ri}(rii),:);
                            
                            conp(rii,1) = interp1(x,con,conped);
                            coni(rii,1) = interp1(x,con,conped+conbin);
                        end
                        
                        conresp = betas(:,ci)' * (coni - conp);
                        coneff = betas(:,ci)'*mean([conp coni],2);
                        % get coh resp
                        for rii = 1:length(rois)
                            coh = respcoh_(ropts{ri}(rii),:);
                            cohp(rii,1) = interp1(x,coh,cohped);
                            cohi(rii,1) = interp1(x,coh,cohped+cohbin);
                        end
                        
                        cohresp = betas(:,ci)' * (cohi - cohp);
                        coheff = betas(:,ci)'*mean([cohp cohi],2);

                        if bi==2
                            % POISSON MODELS: NOISE SCALES
                            % model prediction for CONTRAST
                            cur(1,bi,ri,2,ci,ped,bii) = lapse + (1-2*lapse)*normcdf(0,conresp,sqrt(coneff),'upper');
                            % COHERENCE
                            cur(1,bi,ri,1,ci,ped,bii) = lapse + (1-2*lapse)*normcdf(0,cohresp,sqrt(coheff),'upper');
                        else
                            % model prediction for CONTRAST
                            cur(1,bi,ri,2,ci,ped,bii) = lapse + (1-2*lapse)*normcdf(0,conresp,1,'upper');
                            % COHERENCE
                            cur(1,bi,ri,1,ci,ped,bii) = lapse + (1-2*lapse)*normcdf(0,cohresp,1,'upper');
                        end
                    end
                end
            end
            
            % save threshold somehow
            thresholds(1,bi,ri,:,:) = threshold;
            cthresholds(1,bi,ri,:,:) = cthreshold;
        end
    end
    
    athresholds(ai,:,:,:,:) = thresholds;
    acthresholds(ai,:,:,:,:) = cthresholds;
    model(ai,:,:,:,:,:,:) = cur;

    disppercent(ai/length(aSIDs));
end
disppercent(inf);

%% Check acthresholds to see what it looks like

temp = acthresholds(:,1,2,:,:);
temp = squeeze(temp);

figure; hold on

for i = 1:21
    plot(peds,squeeze(temp(i,:,:)),'*');
end
legend('coh','con');
% 
% figure; hold on
% plot(peds,squeeze(mean(temp(:,1))));
% plot(peds,squeeze(mean(temp(:,2))));

%% Average across subjects

model_ = squeeze(bootci(1000,@nanmean,model));
model_m = squeeze(nanmean(model_));
% model_s = squeeze(model_(2,:,:,:,:,:,:))-model_m;
% model_m = squeeze(mean(model));

thresh__ = squeeze(nanmean(athresholds));
thresh_ = squeeze(bootci(10000,@nanmean,athresholds));
thresh_s = squeeze(thresh_(2,:,:,:,:))-thresh__;

cthresh_ = squeeze(bootci(10000,@nanmean,acthresholds));

cthresh__ = squeeze(nanmean(acthresholds));

%% Average plot across subjects
coherences = [0.15 0.3 0.45 0.6];
contrasts = [0.325 0.4 0.55 0.85];
basecoh = 0; basecon = 0.25;
conbins = 0:.01:1;
conbins_ = conbins(1:(length(conbins)-1))+(diff(conbins)/2);
cohbins = 0:.01:1;
cohbins_ = cohbins(1:(length(cohbins)-1))+(diff(cohbins)/2);
conditions = {'Coherence','Contrast'};
for bm = 1:length(bmodels)
    for ri = 1:length(ropts)
        h = figure;

        amap = brewermap(11,'PuOr');
        maps = [[9 5];[8 3]];
        %  fmri | behav | disc | feat | ped | effect
        %    1       2     3       4      5       6
        % Go through plot now and add model fits

        % split by discrimination conditions

        subopts = [4 1];
        % model_m:
        for ti = 1:2 % which TASK
            subplot(2,1,ti); hold on

            for cii = 1:2 % which FEATURE
                cmap = amap(maps(ti,:),:);
                cmodel = squeeze(mean(model_m(bm,ri,cii,ti,:,:),5));
                plot([-fliplr(conbins_) conbins_],[1-flipud(cmodel) ;cmodel]','-','Color',cmap(cii,:));
            end
        end

        %% Load behavioral data
        load(fullfile(datafolder,'avg_binnedbehav.mat'));
        allData_m = allData.m;
        allData_s = allData.s;
        conbins = allData.conbins;
        cohbins = allData.cohbins;
        %%
        conright = zeros(length(sids),2,4,length(diff(conbins)));

        cohright = zeros(length(sids),2,4,length(diff(cohbins)));

        conditions = {'coherence','contrast'};

        amap = brewermap(11,'PuOr');
        maps = [[9 5];[8 3]];
        bins = [cohbins(1:end-1)+(diff(cohbins)/2);conbins(1:end-1)+(diff(conbins)/2)];

        % average of conditions
        for ti=1:2 % ATTENTION
            subplot(2,1,ti); hold on
            for cii = 1:2 % FEATURE
                cmap = amap(maps(ti,:),:);
                c_m = squeeze(mean(squeeze(allData_m(cii,ti,:,:))));
                c_s = squeeze(mean(squeeze(allData_s(cii,ti,:,:))));
                errbar(bins(ti,:),c_m,c_s,'-','Color',cmap(cii,:));
                pi(ti) = plot(bins(ti,:),c_m,'o','MarkerFaceColor',cmap(cii,:),'MarkerEdgeColor','w','MarkerSize',4);
            end
            
            axis([-0.5 0.5 0 1]);
            set(gca,'XTick',[-0.25 0 0.25],'XTickLabel',[-0.25 0 .25]*100);
            set(gca,'YTick',[0 0.5 1],'YTickLabel',[0 0.5 1]*100);

            xlabel('Stimulus strength (R-L, %)');
            ylabel('Choices to right (%)');
            drawPublishAxis('figSize=[5,5.5]');
        end

%         savepdf(h,fullfile(datafolder,'avg_models',sprintf('right_f2b_%s_%i.pdf',bmodels{bm},length(ropts{ri}))));
    end
end

%% Now plot the threshold plots

% using model_m and conbins_ you can find the threshold necessary to shift
% by X % signal strength

%% Load threshold level data
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

%% Smooth the threshold functions
for b = 1:2
    for r = 1:2
        for t = 6:61
            idxs = (t-5):(t+5);
            thresh_m(b,r,t,:) = squeeze(mean(thresh__(b,r,idxs,:),3));
        end
    end
end
        
thresh_m(:,:,[1:7 62:end],:) = nan;

%% actual plot -- threshold behavior with model fits

% plot the contrast data
h = figure; hold on
cmap = brewermap(7,'PuOr');

contrast = [0.325 0.4 0.55 0.85];

% plot(peds+0.25,squeeze(thresh_m(1,2,:,2)),'-','Color',cmap(2,:));
boundedline(peds+0.25,squeeze(thresh_m(1,2,:,2))',squeeze(thresh_s(1,2,:,2))','-','cmap',cmap(2,:),'nan','remove');
errbar(contrast,control_m(2,:),control_s(2,:),'-','Color',cmap(2,:));
plot(contrast,control_m(2,:),'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','w','MarkerSize',4);
% errbar(coherence,control_m(1,:),control_s(1,:),'-','Color',cmap(6,:));
% plot(coherence,control_m(1,:),'o','MarkerFaceColor',cmap(6,:),'MarkerEdgeColor','w','MarkerSize',7);
% legend('Contrast','Coherence');
% xlabel('Base stimulus (%)');
set(gca,'XTick',[0 0.25 0.5 0.75 1]','XTickLabel',[0 25 50 75 100],'YTick',[0 0.025 0.05 0.075 0.1],'YTickLabel',[0 2.5 5 7.5 10]);
axis([0 1 0 0.1]);
% vline(0.25,'--k');
% ylabel('Threshold (%)');

drawPublishAxis('figSize=[4.25,4.5]');

savepdf(h,fullfile(datafolder,'avg_models','thresholds_contrast.pdf'));
% se[arately plot the coherence data
h = figure; hold on
cmap = brewermap(7,'PuOr');

coherence = [0.15 0.3 0.45 0.6];

% errbar(contrast,control_m(2,:),control_s(2,:),'-','Color',cmap(2,:));
% plot(contrast,control_m(2,:),'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','w','MarkerSize',7);
boundedline(peds,squeeze(thresh_m(1,2,:,1))',squeeze(thresh_s(1,2,:,1))','-','cmap',cmap(6,:),'nan','remove');
errbar(coherence,control_m(1,:),control_s(1,:),'-','Color',cmap(6,:));
% plot(peds,squeeze(thresh_m(1,2,:,1)),'-','Color',cmap(6,:));
plot(coherence,control_m(1,:),'o','MarkerFaceColor',cmap(6,:),'MarkerEdgeColor','w','MarkerSize',4);

% legend('Coherence');
% xlabel('Base stimulus (%)');
set(gca,'XTick',[0 0.25 0.5 0.75]','XTickLabel',[0 25 50 75],'YTick',[0 .05 .1 .15 .2 .25],'YTickLabel',[0:5:25]);
axis([0 .75 0 0.25]);
% vline(0,'--k');
% ylabel('Threshold (%)');

drawPublishAxis('figSize=[4.25,4.5]');

savepdf(h,fullfile(datafolder,'avg_models','thresholds_coherence.pdf'));

%% CATCH trials -- plot the threshold behavior and the model's estimated thresholds


%% actual plot (only behavior -- no model fits, note that we don't use this for the paper)

% plot the contrast data
h = figure; hold on
cmap = brewermap(7,'PuOr');

contrast = [0.325 0.4 0.55 0.85];

errbar(contrast(2),control_m(2,2),control_s(2,2),'-','Color',cmap(2,:));
plot(contrast(2),control_m(2,2),'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','w','MarkerSize',5);
% xlabel('Base stimulus (%)');
axis([0 1 0 1.5]);
set(gca,'XTick',[0 0.25 0.5 0.75 1]','XTickLabel',[0 25 50 75 100],'YTick',0:0.25:1,'YTickLabel',0:25:100);
% ylabel('Threshold (%)');

% contrast catch
errbar(contrast(2),unattend_m(2),unattend_s(2),'-','Color',cmap(1,:));
plot(contrast(2),unattend_m(2),'o','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor','w','MarkerSize',5);

% xlabel('Base stimulus (%)');
axis([0.35 0.45 0 1.5]);
set(gca,'XTick',0.4,'XTickLabel',40,'YTick',0:0.25:1.25,'YTickLabel',{'0','25','50','75','100','Inf'});
% ylabel('Threshold (%)');

noise = 2.1;

% plot the catch model estimate
concatch = cthresh__(1,2,1,2);

% plot the catch model (noise added) estimate
idx = find(peds>=(contrast(2)-0.25),1);
estimate = noise * thresh_m(1,2,idx,2);

plot([0.375 0.425],[concatch concatch],'--','Color',cmap(1,:));
plot([0.375 0.425],[estimate estimate],'-','Color',cmap(1,:));

drawPublishAxis('figSize=[4.5,4.5]');
savepdf(h,fullfile(datafolder,'avg_behav','catch_thresholds_con.pdf'));
%%%%%%%%%%%%%%%%%%%%%%%%

% plot the coherence data
h = figure; hold on

coherence = [0.15 0.3 0.45 0.6];

errbar(coherence(2),control_m(1,2),control_s(1,2),'-','Color',cmap(6,:));
plot(coherence(2),control_m(1,2),'o','MarkerFaceColor',cmap(6,:),'MarkerEdgeColor','w','MarkerSize',5);
% xlabel('Base stimulus (%)');
axis([0.25 0.35 0 1.5]);
set(gca,'XTick',0.3,'XTickLabel',30,'YTick',0:0.25:1.25,'YTickLabel',{'0','25','50','75','100','Inf'});
% ylabel('Threshold (%)');

% coherence catch
errbar(coherence(2),unattend_m(1),unattend_s(1),'-','Color',cmap(7,:));
plot(coherence(2),unattend_m(1),'o','MarkerFaceColor',cmap(7,:),'MarkerEdgeColor','w','MarkerSize',5);


cohcatch = cthresh__(1,2,1,1);
idx = find(peds==coherence(2));
estimate = noise * thresh_m(1,2,idx,1);

plot([0.275 0.325],[cohcatch cohcatch],'--','Color',cmap(7,:));
plot([0.275 0.325],[estimate estimate],'-','Color',cmap(7,:));
drawPublishAxis('figSize=[4.5,4.5]');

savepdf(h,fullfile(datafolder,'avg_behav','catch_thresholds_coh.pdf'));
%%



