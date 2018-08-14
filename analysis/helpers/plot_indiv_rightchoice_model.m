
restructure_afits;

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

idatas = cell(length(aSIDs),length(bmodels));
for ai = 1:length(aSIDs)
    for bi = 1:length(bmodels)
        dat = afits{ai}{bi};
        
        % collect data
        adata = dat.adata;
        probs = dat.probs;
        resp = dat.resp;
        probs(resp==0) = 1-probs(resp==0);
        
        % setup data tracker
        %      task | stimulus | model (1==true, 2==model) | bin
        idata = zeros(2,2,2,length(conbins));
        n = zeros(2,2,length(conbins));
        
        % split by attention condition
        for attendCond = 1:2
            aidxs = adata(:,1)==attendCond;
            
            cdata = adata(aidxs,:);
            cprobs = probs(aidxs);
            cresp = resp(aidxs);
            
            % collect each contrast level
            for ci = 2:length(conbinx)
                dcon = cdata(:,5)-cdata(:,4);
                lidxs = logical((dcon>conbinx(ci-1)) .* (dcon<=conbinx(ci)));
                if sum(lidxs)>2
                    idata(attendCond,2,1,ci-1) = nanmean(cresp(lidxs));
                    idata(attendCond,2,2,ci-1) = nanmean(cprobs(lidxs));
                else
                    idata(attendCond,2,1,ci-1) = nan;
                    idata(attendCond,2,2,ci-1) = nan;
                end
                n(attendCond,2,ci-1) = sum(lidxs);
            end
            
            % collect each coherence level
            for mi = 2:length(cohbinx)
                dcoh = cdata(:,7)-cdata(:,6);
                lidxs = logical((dcoh>cohbinx(mi-1)) .* (dcoh<=cohbinx(mi)));
                if sum(lidxs)>2
                    idata(attendCond,1,1,mi-1) = nanmean(cresp(lidxs));
                    idata(attendCond,1,2,mi-1) = nanmean(cprobs(lidxs));
                else
                    idata(attendCond,1,1,mi-1) = nan;
                    idata(attendCond,1,2,mi-1) = nan;
                end
                n(attendCond,1,mi-1) = sum(lidxs);
            end
        end
        idatas{ai,bi} = idata;
    end
end

%% Plotting routines

amap = brewermap(11,'PuOr');

smaps = [[9 5];[8 3]];

bmaps = [[2 3 4 5];[10 9 8 7]];

for ai = 1:length(aSIDs)
    subj = aSIDs(ai);
    for bi = 1:length(bmodels)
        idata = idatas{ai,bi};
        % one figure per subject
        h = figure;
        % MINI PLOTS
        subopts = [2;1];
        conds = {'Motion','Contrast'};
        for attendCond = 1:2
            subplot(2,1,subopts(attendCond,:)); hold on
            title(sprintf('Attending %s',conds{attendCond}));
            cmap = amap(smaps(attendCond,:),:);

            cdat = squeeze(idata(attendCond,2,:,:)); % contrast
            
            % plot model fit
            

            for cii = 1:2 % which FEATURE
                cmap = amap(maps(attendCond,:),:);
                cmodel = squeeze(mean(model(ai,bi,2,cii,attendCond,:,:),6));
                plot([-fliplr(conbins_) conbins_],[1-flipud(cmodel) ;cmodel]','-','Color',cmap(cii,:));
            end
            title(sprintf('%s discrimination: %s',conditions{attendCond},bmodels{bi}));
        
            mdat = squeeze(idata(attendCond,1,:,:)); % contrast

%             plot model fit
            plot(conbins,squeeze(cdat(2,:)),'-','Color',cmap(2,:));
            plot(conbins,squeeze(mdat(2,:)),'-','Color',cmap(1,:));

            % plot behavioral data
            plot(conbins,squeeze(cdat(1,:)),'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor',cmap(2,:));
            plot(conbins,squeeze(mdat(1,:)),'o','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor',cmap(1,:));
            
            % set axis
            axis([-0.6 0.6 0 1]);
            set(gca,'XTick',[-0.25 0 0.25],'XTickLabel',[-0.25 0 .25]*100);
            set(gca,'YTick',[0 0.5 1],'YTickLabel',[0 0.5 1]*100);
            
            drawPublishAxis;
        end

        folder = fullfile(datafolder,'all_models',sprintf('%03.0f',subj));
        if ~isdir(folder)
            mkdir(folder);
        end
        
%         savepdf(h,fullfile(folder,sprintf('%s_modelfit.pdf',bmodels{bi})));
%         close all
    end
end