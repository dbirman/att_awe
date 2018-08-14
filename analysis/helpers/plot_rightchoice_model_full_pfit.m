
restructure_afits_full_pfit;

%% Subject plot

coherences = [0.15 0.3 0.45 0.6];
contrasts = [0.325 0.4 0.55 0.85];
basecoh = 0; basecon = 0.25;
conbins = 0:.01:1;
conbins_ = conbins(1:(length(conbins)-1))+(diff(conbins)/2);
cohbins = 0:.01:1;
cohbins_ = cohbins(1:(length(cohbins)-1))+(diff(cohbins)/2);

load(fullfile(datafolder,'avg_binnedbehav.mat'));
% allData_m = allData.m;
% allData_s = allData.s;
conbins = allData.conbins;
conbinx = [-inf (conbins + [diff(conbins/2) 0])]; conbinx(end) = inf;
cohbins = allData.cohbins;
cohbinx = [-inf (conbins + [diff(conbins/2) 0])]; cohbinx(end) = inf;

conbins_ca = conbins(1:2:11);
conbinx_ca = [-inf (conbins_ca + [diff(conbins_ca/2) 0])]; conbinx_ca(end) = inf;
cohbins_ca = cohbins(1:2:11);
cohbinx_ca = [-inf (conbins_ca + [diff(conbins_ca/2) 0])]; cohbinx_ca(end) = inf;

adata = zeros(length(aSIDs),length(bmodels),2,2,2,length(conbins));
adata_ca = zeros(length(aSIDs),length(bmodels),2,2,2,length(cohbins_ca));
for ai = 1:length(aSIDs)
    for bi = 1
        dat = ffits{ai}{bi};
        
        % collect data
        ldata = dat.adata;
        probs = dat.probs;
        resp = dat.resp;
        probs(resp==0) = 1-probs(resp==0);
        
        % setup data tracker
        %      task | stimulus | model (1==true, 2==model) | bin
%         idata = zeros(2,2,2,length(conbins));
        n = zeros(2,2,length(conbins));
        
        control_idxs = ldata(:,9)<=0;
        catch_idxs = ldata(:,9)>0;
        
        flip = [2 1];
        % split by attention condition
        for attendCond = 1:2
            catchidxs = catch_idxs .* (ldata(:,1)==-flip(attendCond));
            controlidxs = control_idxs .* (ldata(:,1)==attendCond);

            idxs = logical(controlidxs);
            cdata = ldata(idxs,:);
            cprobs = probs(idxs);
            cresp = resp(idxs);
    
            % CONTROL DATA
            catchCond = 1;
            
            usebins = conbinx;
            dcon = cdata(:,5)-cdata(:,4);
            for ci = 2:length(usebins)
                lidxs = logical((dcon>usebins(ci-1)) .* (dcon<=usebins(ci)));
                if sum(lidxs)>2
                    adata(ai,bi,attendCond,2,1,ci-1) = nanmean(cresp(lidxs));
                    adata(ai,bi,attendCond,2,2,ci-1) = nanmean(cprobs(lidxs));
                else
                    adata(ai,bi,attendCond,2,1,ci-1) = nan;
                    adata(ai,bi,attendCond,2,2,ci-1) = nan;
                end
                n(attendCond,2,ci-1) = sum(lidxs);
            end

            usebins = cohbinx;
            % collect each coherence level
            dcoh = cdata(:,7)-cdata(:,6);
            for mi = 2:length(usebins)
                lidxs = logical((dcoh>usebins(mi-1)) .* (dcoh<=usebins(mi)));
                if sum(lidxs)>2
                    adata(ai,bi,attendCond,1,1,mi-1) = nanmean(cresp(lidxs));
                    adata(ai,bi,attendCond,1,2,mi-1) = nanmean(cprobs(lidxs));
                else
                    adata(ai,bi,attendCond,1,1,mi-1) = nan;
                    adata(ai,bi,attendCond,1,2,mi-1) = nan;
                end
                n(attendCond,1,mi-1) = sum(lidxs);
            end
            
            % CATCH DATA
            
            idxs = logical(catchidxs);
            cdata = ldata(idxs,:);
            cprobs = probs(idxs);
            cresp = resp(idxs);
            
            catchCond = 2;
            
            usebins = conbinx_ca;
            dcon = cdata(:,5)-cdata(:,4);
            for ci = 2:length(usebins)
                lidxs = logical((dcon>usebins(ci-1)) .* (dcon<=usebins(ci)));
                if sum(lidxs)>2
                    adata_ca(ai,bi,attendCond,2,1,ci-1) = nanmean(cresp(lidxs));
                    adata_ca(ai,bi,attendCond,2,2,ci-1) = nanmean(cprobs(lidxs));
                else
                    adata_ca(ai,bi,attendCond,2,1,ci-1) = nan;
                    adata_ca(ai,bi,attendCond,2,2,ci-1) = nan;
                end
                n(attendCond,2,ci-1) = sum(lidxs);
            end

            usebins = cohbinx_ca;
            % collect each coherence level
            dcoh = cdata(:,7)-cdata(:,6);
            for mi = 2:length(usebins)
                lidxs = logical((dcoh>usebins(mi-1)) .* (dcoh<=usebins(mi)));
                if sum(lidxs)>2
                    adata_ca(ai,bi,attendCond,1,1,mi-1) = nanmean(cresp(lidxs));
                    adata_ca(ai,bi,attendCond,1,2,mi-1) = nanmean(cprobs(lidxs));
                else
                    adata_ca(ai,bi,attendCond,1,1,mi-1) = nan;
                    adata_ca(ai,bi,attendCond,1,2,mi-1) = nan;
                end
                n(attendCond,1,mi-1) = sum(lidxs);
            end
        end
    end
end


%%

amap = brewermap(11,'PuOr');

smaps = [[9 8];[3 5]];

bmaps = [[2 3 4 5];[10 9 8 7]];

idata_control = squeeze(nanmean(adata(:,bi,:,:,:,:)));
idata_catch = squeeze(nanmean(adata_ca(:,bi,:,:,:,:)));
% one figure, average across subjects
h = figure;
% MINI PLOTS
subopts = [2;1];
conds = {'Motion','Contrast'};
for attendCond = 1:2
    feat = attendCond;
    subplot(2,1,subopts(attendCond,:)); hold on
    title(sprintf('Discriminating %s',conds{attendCond}));
    cmap = amap(smaps(attendCond,:),:);

    condat = squeeze(idata_control(attendCond,feat,:,:)); % control
    catdat = squeeze(idata_catch(attendCond,feat,:,:)); % catch

    % plot model fit
    plot(conbins,squeeze(condat(2,:)),'-','Color',cmap(1,:),'LineWidth',2);
    plot(conbins_ca,squeeze(catdat(2,:)),'-','Color',cmap(2,:));

    legend({'Control','Catch'},'FontSize',7,'FontName','Helvetica');

    % plot behavioral data
    plot(conbins,squeeze(condat(1,:)),'o','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor',cmap(2,:));
    plot(conbins_ca,squeeze(catdat(1,:)),'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor',cmap(2,:));

    % set axis
    axis([-0.6 0.6 0 1]);
    set(gca,'XTick',[-0.25 0 0.25],'XTickLabel',[-0.25 0 .25]*100);
    set(gca,'YTick',[0 0.5 1],'YTickLabel',[0 0.5 1]*100);

    drawPublishAxis;
end

    savepdf(h,fullfile(datafolder,'avg_models',sprintf('%s_modelfit_full_pfit.pdf',bmodels{bi})));
%     close all

%%

x = 0:.001:2;
if ~(size(respcon_,2)==length(x))
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
end

%% Actual model fits
restructure_afits_full_pfit;

coherences = [0.15 0.3 0.45 0.6];
contrasts = [0.325 0.4 0.55 0.85];
basecoh = 0; basecon = 0.25;
conbins = 0:.01:1;
conbins_ = conbins(1:(length(conbins)-1))+(diff(conbins)/2);
cohbins = 0:.01:1;
cohbins_ = cohbins(1:(length(cohbins)-1))+(diff(cohbins)/2);

% subj | fmri | behav | disc | feat | ped | effect
%   1      2     3       4      5       6     7
model = zeros(length(aSIDs),length(bmodels),2,2,4,length(conbins_));
model_ca = zeros(length(aSIDs),length(bmodels),2,2,4,length(conbins_));

disppercent(-1/length(aSIDs));
for ai = 1:length(aSIDs) % which subject
    %     for mi = 1:length(mopts) % which fMRI model
    
    for bi = 1
        % get betas
        rois = {'V1','MT'};
        cons = {'cohw','conw'};
        betas = zeros(8,2);
        betas_ca = zeros(8,2);
        
        for ri = 1:2
            for ci = 1:2
                betas(ri,ci) = ffits{ai}{bi}.params.(sprintf('beta_control_%s_%s',rois{ri},cons{ci}));
                betas_ca(ri,ci) = ffits{ai}{bi}.params.(sprintf('beta_unatt_%s_%s',rois{ri},cons{ci}));
            end
        end
        
        bias = ffits{ai}{bi}.params.bias;
        lapse = ffits{ai}{bi}.params.lapse;
        lapse = 0;
        
        for ci = 1:2
            for ped = 1:4 % which pedestal
                conped = contrasts(ped);
                cohped = coherences(ped);
                
                for bii = 1:length(conbins_) % which bin
                    conbin = conbins_(bii);
                    cohbin = cohbins_(bii);
                    % best for loop chain I've ever written
                    % compute probability of a right response for this
                    % set of conditions
                    
                    % using respcon_ and respcoh_
                    
                    % get con resp
                    conp = respcon_(:,(find(x>=conped,1)));
                    coni = respcon_(:,find(x>=(conped+conbin),1));
                    conresp = betas(:,ci)' * (coni - conp);
                    conresp_ca = betas_ca(:,ci)' * (coni-conp);
                    coneff = betas(:,ci)'*mean([conp coni],2);
                    % get coh resp
                    cohp = respcoh_(:,(find(x>=cohped,1)));
                    cohi = respcoh_(:,find(x>=(cohped+cohbin),1));
                    cohresp = betas(:,ci)' * (cohi - cohp);
                    cohresp_ca = betas_ca(:,ci)' * (cohi - cohp);
                    coheff = betas(:,ci)'*mean([cohp cohi],2);
                    
                    if bi==2
                        % POISSON MODELS: NOISE SCALES
                        % model prediction for CONTRAST
                        model(ai,bi,ci,2,ped,bii) = lapse + (1-2*lapse)*normcdf(0,conresp,sqrt(coneff),'upper');
                        % COHERENCE
                        model(ai,bi,ci,1,ped,bii) = lapse + (1-2*lapse)*normcdf(0,cohresp,sqrt(coheff),'upper');
                    else
                        % model prediction for CONTRAST
                        model(ai,bi,ci,2,ped,bii) = lapse + (1-2*lapse)*normcdf(0,conresp,1,'upper');
                        % COHERENCE
                        model(ai,bi,ci,1,ped,bii) = lapse + (1-2*lapse)*normcdf(0,cohresp,1,'upper');
                        
                        % catch predictions
                        % model prediction for CONTRAST
                        model_ca(ai,bi,ci,2,ped,bii) = lapse + (1-2*lapse)*normcdf(0,conresp_ca,1,'upper');
                        % COHERENCE
                        model_ca(ai,bi,ci,1,ped,bii) = lapse + (1-2*lapse)*normcdf(0,cohresp_ca,1,'upper');
                    end
                end
            end
        end
    end
    disppercent(ai/length(aSIDs));
end
disppercent(inf);

%% Average across subjects

% model_ = squeeze(bootci(10000,@mean,model));
% model_m = squeeze(mean(model_));
% model_s = squeeze(model_(2,:,:,:,:,:,:))-model_m;
model_m = squeeze(mean(model));
model_ca_m = squeeze(mean(model_ca));
%% Average plot across subjects
coherences = [0.15 0.3 0.45 0.6];
contrasts = [0.325 0.4 0.55 0.85];
basecoh = 0; basecon = 0.25;
conbins2 = 0:.01:1;
conbins_2 = conbins2(1:(length(conbins2)-1))+(diff(conbins2)/2);
cohbins2 = 0:.01:1;
cohbins_2 = cohbins2(1:(length(cohbins2)-1))+(diff(cohbins2)/2);
conditions = {'Coherence','Contrast'};


conbins = allData.conbins;
% conbinx = [-inf (conbins + [diff(conbins/2) 0])]; conbinx(end) = inf;
cohbins = allData.cohbins;
% cohbinx = [-inf (conbins + [diff(conbins/2) 0])]; cohbinx(end) = inf;

conbins_ca = conbins(1:2:11);
% conbinx_ca = [-inf (conbins_ca + [diff(conbins_ca/2) 0])]; conbinx_ca(end) = inf;
cohbins_ca = cohbins(1:2:11);
% cohbinx_ca = [-inf (conbins_ca + [diff(conbins_ca/2) 0])]; cohbinx_ca(end) = inf;

h = figure;
    
amap = brewermap(11,'PuOr');
maps = [[9 8];[3 5]];
%  fmri | behav | disc | feat | ped | effect
%    1       2     3       4      5       6
% Go through plot now and add model fits

% split by discrimination conditions

subopts = [2 1];
flip = [2 1];
% model_m:
feature = {'Motion','Contrast'};
for di = 1:2 % which attention condition
    subplot(2,1,subopts(di)); hold on

    cmap = amap(maps(di,:),:);
    % control
    cmodel = squeeze(mean(model_m(di,di,:,:),3));
    plot([-fliplr(conbins_2) conbins_2],[1-flipud(cmodel) ;cmodel]','-','Color',cmap(1,:));
%     
%     % catch
    cmodel_ca = squeeze(mean(model_ca_m(di,di,:,:),3));
    plot([-fliplr(conbins_2) conbins_2],[1-flipud(cmodel_ca) ;cmodel_ca]','-','Color',cmap(2,:));
    legend({'Control','Catch'},'FontSize',7,'FontName','Helvetica');
% 
    l = legend(gca,'boxoff');
    set(l,'Color','none');
    title(sprintf('%s discrimination',conditions{di}));
    
    condat = squeeze(idata_control(di,di,:,:)); % control
    if di==2
        catdat = squeeze(idata_catch(di,1,:,:)); % catch
    else
        catdat = squeeze(idata_catch(di,di,:,:)); % catch
    end

    % plot behavioral data
    if di==1
        usebins = cohbins;
        usebins_ca = cohbins_ca;
    else
        usebins = conbins;
        usebins_ca = conbins_ca;
    end
    plot(usebins,squeeze(condat(di,:)),'o','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor',cmap(2,:),'MarkerSize',4);
    plot(usebins_ca,squeeze(catdat(di,:)),'^','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor',cmap(2,:),'MarkerSize',4);
    
    
%     axis([-.6 .6 0 1]);
    set(gca,'XTick',[-1 -.5 0 0.5 1],'XTickLabel',{'-100%','-50%','0%','50%','100%'});
    set(gca,'YTick',[0 1],'YTickLabel',{'0%','100%'});
    xlabel(sprintf('\\Delta %s (%)',feature{di}));
    ylabel(sprintf('Right choice probability (%)'));
    
    drawPublishAxis('figSize=[4.5 4.5]');
end

%     savepdf(h,fullfile('~/proj/att_awe/talks/data_figures/',sprintf('right_f2b_%s.pdf',bmodels{bm})));
savepdf(h,fullfile(datafolder,'avg_models','catch_f2b_pfit.pdf'));