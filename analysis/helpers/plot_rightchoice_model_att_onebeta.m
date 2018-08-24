function plot_rightchoice_model_att_onebeta(fitdata,respcon_,respcoh_,aSIDs,bmodels,rois)
% restructure_afits_att;
mopts = 1;

x = 0:.001:1;

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

adata = zeros(length(aSIDs),1,2,2,2,length(conbins));
for ai = 1:length(aSIDs)
    for bi = 1%:length(bmodels)
        dat = fitdata{ai};
        
        % collect data
        ldata = dat.adata;
        probs = dat.probs;
        resp = dat.resp;
        probs(resp==0) = 1-probs(resp==0);
        
        % setup data tracker
        %      task | stimulus | model (1==true, 2==model) | bin
%         idata = zeros(2,2,2,length(conbins));
        n = zeros(2,2,length(conbins));
        
        % split by attention condition
        for attendCond = 1:2
            aidxs = ldata(:,1)==attendCond;
            
            cdata = ldata(aidxs,:);
            cprobs = probs(aidxs);
            cresp = resp(aidxs);
            
            % collect each contrast level
            for ci = 2:length(conbinx)
                dcon = cdata(:,5)-cdata(:,4);
                lidxs = logical((dcon>conbinx(ci-1)) .* (dcon<=conbinx(ci)));
                if sum(lidxs)>2
                    adata(ai,bi,attendCond,2,1,ci-1) = nanmean(cresp(lidxs));
                    adata(ai,bi,attendCond,2,2,ci-1) = nanmean(cprobs(lidxs));
                else
                    adata(ai,bi,attendCond,2,1,ci-1) = nan;
                    adata(ai,bi,attendCond,2,2,ci-1) = nan;
                end
                n(attendCond,2,ci-1) = sum(lidxs);
            end
            
            % collect each coherence level
            for mi = 2:length(cohbinx)
                dcoh = cdata(:,7)-cdata(:,6);
                lidxs = logical((dcoh>cohbinx(mi-1)) .* (dcoh<=cohbinx(mi)));
                if sum(lidxs)>2
                    adata(ai,bi,attendCond,1,1,mi-1) = nanmean(cresp(lidxs));
                    adata(ai,bi,attendCond,1,2,mi-1) = nanmean(cprobs(lidxs));
                else
                    adata(ai,bi,attendCond,1,1,mi-1) = nan;
                    adata(ai,bi,attendCond,1,2,mi-1) = nan;
                end
                n(attendCond,1,mi-1) = sum(lidxs);
            end
        end
    end
end


%% Actual model fits
% x = 0:.001:2;
% restructure_afits;

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

disppercent(-1/length(aSIDs));
for ai = 1:length(aSIDs) % which subject
    %     for mi = 1:length(mopts) % which fMRI model
    
    for bi = 1:length(bmodels)
        % get betas
        cons = {'cohw','conw'};
        
%         rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
        betas = zeros(length(rois),2);
        for ri = 1:length(rois)
            for ci = 1:2

                betas(ri,:) = repmat(fitdata{ai}.params.(sprintf('beta_control_%s_w',rois{ri})),1,2);            
            end
        end
        
        bias = fitdata{ai}.params.bias;
        lapse = fitdata{ai}.params.lapse;
        
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
                    
                        % attending contrast
                    crespcon = squeeze(respcon_(:,ci,:));
                    crespcoh = squeeze(respcoh_(:,ci,:));
                    
                    % get con resp
                    conp = crespcon(:,(find(x>=conped,1)));
                    coni = crespcon(:,find(x>=(conped+conbin),1));
                    if isempty(coni)
                        coni = crespcon(:,end);
                    end
                    conresp = betas(:,ci)' * (coni - conp);
                    coneff = betas(:,ci)'*mean([conp coni],2);
                    % get coh resp
                    cohp = crespcoh(:,(find(x>=cohped,1)));
                    cohi = crespcoh(:,find(x>=(cohped+cohbin),1));
                    if isempty(cohi)
                        cohi = crespcoh(:,end);
                    end
                    cohresp = betas(:,ci)' * (cohi - cohp);
                    coheff = betas(:,ci)'*mean([cohp cohi],2);
                    
                    if bi==2
                        % POISSON MODELS: NOISE SCALES
                        % model prediction for CONTRAST
                        model(ai,bi,2,ci,ped,bii) = lapse + (1-2*lapse)*normcdf(0,conresp,sqrt(coneff),'upper');
                        % COHERENCE
                        model(ai,bi,1,ci,ped,bii) = lapse + (1-2*lapse)*normcdf(0,cohresp,sqrt(coheff),'upper');
                    else
                        % model prediction for CONTRAST
                        model(ai,bi,2,ci,ped,bii) = lapse + (1-2*lapse)*normcdf(0,conresp,1,'upper');
                        % COHERENCE
                        model(ai,bi,1,ci,ped,bii) = lapse + (1-2*lapse)*normcdf(0,cohresp,1,'upper');
                    end
                end
            end
        end
    end
    disppercent(ai/length(aSIDs));
end
disppercent(inf);

%% Average across subjects

model_ = squeeze(bootci(10000,@mean,model));
model_m = squeeze(mean(model_));
model_s = squeeze(model_(2,:,:,:,:,:,:))-model_m;
% model_m = squeeze(mean(model));

%% Average plot across subjects
coherences = [0.15 0.3 0.45 0.6];
contrasts = [0.325 0.4 0.55 0.85];
basecoh = 0; basecon = 0.25;
conbins = 0:.01:1;
conbins_ = conbins(1:(length(conbins)-1))+(diff(conbins)/2);
cohbins = 0:.01:1;
cohbins_ = cohbins(1:(length(cohbins)-1))+(diff(cohbins)/2);
conditions = {'Coherence','Contrast'};
for bm = 1
    h = figure;
    
    amap = brewermap(11,'PuOr');
    maps = [[9 4];[8 3]];
    %  fmri | behav | disc | feat | ped | effect
    %    1       2     3       4      5       6
    % Go through plot now and add model fits
    
    % split by discrimination conditions
    
    subopts = [2 1];
    % model_m:
    for di = 1:2 % which attention condition
        subplot(2,1,di); hold on
        
        for cii = 1:2 % which condition
            cmap = amap(maps(di,:),:);
            cmodel = squeeze(mean(model_m(cii,di,:,:),3));
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
    conright = zeros(length(aSIDs),2,4,length(diff(conbins)));
    
    cohright = zeros(length(aSIDs),2,4,length(diff(cohbins)));
    
    % load(fullfile(datafolder,'avg_binnedbehav.mat'),'
    conditions = {'coherence','contrast'};
    
    amap = brewermap(11,'PuOr');
    maps = [[9 4];[8 3]];
    bins = [cohbins(1:end-1)+(diff(cohbins)/2);conbins(1:end-1)+(diff(conbins)/2)];
    
    % average of conditions
    for di=1:2
        subplot(2,1,di); hold on
        for cii = 1:2
            cmap = amap(maps(di,:),:);
            c_m = squeeze(mean(squeeze(allData_m(cii,di,:,:))));
            c_s = squeeze(mean(squeeze(allData_s(cii,di,:,:))));
            errbar(bins(di,:),c_m,c_s,'-','Color',cmap(cii,:));
            pi(di) = plot(bins(di,:),c_m,'o','MarkerFaceColor',cmap(cii,:),'MarkerEdgeColor','w','MarkerSize',3);
        end
        axis([-0.5 0.5 0 1]);
        set(gca,'XTick',[-0.25 0 0.25],'XTickLabel',[-0.25 0 .25]*100);
        set(gca,'YTick',[0 0.5 1],'YTickLabel',[0 0.5 1]*100);

        xlabel('Stimulus strength (R-L, %)');
        ylabel('Choices to right (%)');
        drawPublishAxis('figSize=[5,5.5]');
    end
    
    savepdf(h,fullfile(datafolder,'avg_models',sprintf('right_f2b_%s_att_onebeta_%i.pdf',bmodels{bm},length(rois))));
end