%%
coherences = [0.15 0.3 0.45 0.6];
contrasts = [0.325 0.4 0.55 0.85];

%% Determine spacing
% load all data for all subjects, find the 10% quantiles for
% contrast/coherence
aadata = zeros(100000,2); count = 1;
for si = 1:length(sids)
    [adata,s] = loadadata(sids{si});
    adata = adata(adata(:,9)==-1,:);
    adata(:,13) = adata(:,5)-adata(:,4); % contrast diff
    adata(:,14) = adata(:,7)-adata(:,6); % coherence diff
    aadata(count:count+size(adata,1)-1,:) = [adata(:,13) adata(:,14)];
    count = count+size(adata,1);
end
aadata = aadata(1:count-1,:);
conbins = [quantile(aadata(:,1),0:.1:1)];
cohbins = [quantile(aadata(:,2),0:.1:1)];
%% rightperformancecontrol
% Load each subjects raw data. Compute for each pedestal the contrast
% response curve for attend vs. unattend. Plot these (2x4 subplots).
% Average across subjects.
conright = zeros(length(sids),2,4,length(diff(conbins)));

cohright = zeros(length(sids),2,4,length(diff(cohbins)));

for si = 1:length(sids)
    [adata,s] = loadadata(sids{si});
    
    % reduce data to control task
    adata = adata(adata(:,9)==-1,:);
    
    % add pedestals
    for ai = 1:size(adata,1)
        pl_con = find(contrasts==adata(ai,4));
        pr_con = find(contrasts==adata(ai,5));
        if ~isempty(pl_con), adata(ai,13)=pl_con; else adata(ai,13)=pr_con; end
        pl_coh = find(coherences==adata(ai,6));
        pr_coh = find(coherences==adata(ai,7));
        if ~isempty(pl_coh), adata(ai,14)=pl_coh; else adata(ai,14)=pr_coh; end
    end
    
    % add R-L differences
    adata(:,15) = adata(:,5)-adata(:,4); % contrast diff
    adata(:,16) = adata(:,7)-adata(:,6); % coherence diff
    
    % bin by pedestal and calculate percent right (not percent correct) as
    % a function of R-L difference for that feature (we assume independence
    % here otherwise this is actually not possible to do). First we do this
    % just for contrast alone
    pright_con = zeros(2,4,length(diff(conbins))); 
    % EFFECT OF CONTRAST ACROSS CONDITIONS
    for di = 1:2
        % 1 = discriminating coherence, 2 = contrast
        adata_con = adata(adata(:,1)==di,:);
        for ped = 1:4
            % get data at current pedestal
            local_con = adata_con(adata_con(:,13)==ped,:);
            for ci = 1:(length(conbins)-1)
                % calculate % right within this bin
                idxs = logical((local_con(:,15)>=conbins(ci)).*(local_con(:,15)<conbins(ci+1)));
                pright_con(di,ped,ci) = nanmean(local_con(idxs,8));
            end
        end
    end
    
    pright_coh = zeros(2,4,length(diff(conbins)));
    % EFFECT OF COHERENCE ACROSS CONDITIONS
    for di = 1:2
        % 1 = discriminating coherence, 2 = cohtrast
        adata_coh = adata(adata(:,1)==di,:);
        for ped = 1:4
            % get data at current pedestal
            local_coh = adata_coh(adata_coh(:,14)==ped,:);
            for ci = 1:(length(cohbins)-1)
                % calculate % right within this bin
                idxs = logical((local_coh(:,16)>=cohbins(ci)).*(local_coh(:,16)<cohbins(ci+1)));
                pright_coh(di,ped,ci) = nanmean(local_coh(idxs,8));
            end
        end
    end

%     figure; hold on
%     cmap = brewermap(7,'PuOr'); cmap = cmap([6 5],:);
%     for di = 1:2
%         for ped = 1:4
%             plot(cohbins(1:end-1)+(diff(cohbins)/2),squeeze(pright_coh(di,ped,:)),'o','MarkerFaceColor',cmap(di,:),'MarkerEdgeColor','w','MarkerSize',10);
%         end
%     end
%     
    conright(si,:,:,:) = pright_con;
    cohright(si,:,:,:) = pright_coh;
end

%% Average
% FEATURE, ATTENTION CONDITION, PEDESTAL, EFFECT SIZE
conright_ = squeeze(bootci(1000,@nanmean,conright));
conright_m = squeeze(mean(conright_));
conright_s = squeeze(conright_(2,:,:,:)) - conright_m;
% coherence
cohright_ = squeeze(bootci(1000,@nanmean,cohright));
cohright_m = squeeze(mean(cohright_));
cohright_s = squeeze(cohright_(2,:,:,:)) - cohright_m;

clear allData_m allData_s
allData_m(1,:,:,:) = cohright_m;
allData_m(2,:,:,:) = conright_m;
allData_s(1,:,:,:) = cohright_s;
allData_s(2,:,:,:) = conright_s;

%% Save
allData.m = allData_m;
allData.s = allData_s;
allData.conbins = conbins;
allData.cohbins = cohbins;
save(fullfile(datafolder,'avg_binnedbehav.mat'),'allData');

%% Fit a cumulative gaussian to allData_m for contrast/coherence and report R^2
% concon
bins = [cohbins(1:end-1)+(diff(cohbins)/2);conbins(1:end-1)+(diff(conbins)/2)];
con_m = squeeze(mean(conright(:,2,:,:),3));
con_x = bins(2,:);
coh_m = squeeze(mean(cohright(:,1,:,:),3));
coh_x = bins(1,:);
figure; hold on
cmap = brewermap(7,'PuOr');
plot(con_x,con_m,'o','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','w');
plot(coh_x,coh_m,'o','MarkerFaceColor',cmap(6,:),'MarkerEdgeColor','w');

%% Fits
stacky = [];
stacky_ = [];
for si = 1:size(con_m,1)
    % drop this person
    test = con_m(si,:);
    train = con_m(setdiff(1:size(con_m,1),si),:);
    
    fitcn = fit_normcdf(con_x,nanmean(train));
    
    ty = (1-2*fitcn.params.lapse)*normcdf(con_x,fitcn.params.mu,fitcn.params.sigma)+fitcn.params.lapse;
    stacky = [stacky test];
    stacky_ = [stacky_ ty];
    
    test = coh_m(si,:);
    train = coh_m(setdiff(1:size(coh_m,1),si),:);
    
    fitcn = fit_normcdf(coh_x,nanmean(train));
    
    ty = (1-2*fitcn.params.lapse)*normcdf(coh_x,fitcn.params.mu,fitcn.params.sigma)+fitcn.params.lapse;
    stacky = [stacky test];
    stacky_ = [stacky_ ty];
end
r2 = 1 - nansum((stacky_-stacky).^2)/nansum((stacky-0.5).^2);

%% Compute statistics (zero conditions and non-zero conditions)
for si = 1:size(con_m,1)
    % compute the sigma (slope) parameters
    vals = con_m(si,:);
    tx = con_x(~isnan(vals));
    vals = vals(~isnan(vals));
    fit = fit_normcdf(tx,vals);
    conslope(si) = (1-2*fit.params.lapse)*normpdf(0,fit.params.mu,fit.params.sigma);
    
    vals = coh_m(si,:);
    tx = coh_x(~isnan(vals));
    vals = vals(~isnan(vals));
    fit = fit_normcdf(tx,vals);
    cohslope(si) = (1-2*fit.params.lapse)*normpdf(0,fit.params.mu,fit.params.sigma);
    
    vals = squeeze(mean(conright(si,1,:,:),3))';
    tx = con_x(~isnan(vals));
    vals = vals(~isnan(vals));
    fit = fit_normcdf(tx,vals);
    conslope_null(si) = (1-2*fit.params.lapse)*normpdf(0,fit.params.mu,fit.params.sigma);
    
    vals = squeeze(mean(cohright(si,2,:,:),3))';
    tx = con_x(~isnan(vals));
    vals = vals(~isnan(vals));
    fit = fit_normcdf(tx,vals);
    cohslope_null(si) = (1-2*fit.params.lapse)*normpdf(0,fit.params.mu,fit.params.sigma);
end

%%
h = figure; hold on

group = {'conslope','cohslope','conslope_null','cohslope_null'};

for gi = 1:4
    dat = eval(group{gi});
    plot(gi*ones(1,length(dat)),dat,'*k');
    plot([-.05 .05]+gi,repmat(mean(bootci(1000,@mean,dat)),1,2),'-k');
    
    ci = bootci(1000,@nanmean,dat);
    disp(sprintf('%s, mean: %01.2f 95%% CI [%01.2f %01.2f]',group{gi},mean(ci),ci(1),ci(2)));
end

ddat = conslope-cohslope;
ci = bootci(1000,@nanmean,ddat);
disp(sprintf('Differences, mean: %01.2f 95%% CI [%01.2f %01.2f]',mean(ci),ci(1),ci(2)))

%% Mini-panels: contrast / coherence discrimination

conditions = {'coherence','contrast'};

amap = brewermap(7,'PuOr');
maps = [[6 5];[3 2]];
conditionOpts = [2 1];
bins = [cohbins(1:end-1)+(diff(cohbins)/2);conbins(1:end-1)+(diff(conbins)/2)];

% contrast effect
con_m = squeeze(allData_m(2,:,:,:));
con_s = squeeze(allData_s(2,:,:,:));
con_con_m = squeeze(mean(con_m(2,:,:),2));
con_coh_m = squeeze(mean(con_m(1,:,:),2));
con_con_s = squeeze(mean(con_s(2,:,:),2));
con_coh_s = squeeze(mean(con_s(1,:,:),2));
% coherence effect
coh_m = squeeze(allData_m(2,:,:,:));
coh_s = squeeze(allData_s(2,:,:,:));
coh_coh_m = squeeze(mean(coh_m(2,:,:),2));
coh_con_m = squeeze(mean(coh_m(1,:,:),2));
coh_coh_s = squeeze(mean(coh_s(2,:,:),2));
coh_con_s = squeeze(mean(coh_s(1,:,:),2));

h = figure; hold on

% coherence effect onc ontrast discrimination
cmc = amap(maps(1,2),:);
errbar(bins(1,:)',coh_con_m,coh_con_s,'-','Color',cmc);
p(1) = plot(bins(1,:),coh_con_m,'o','MarkerFaceColor',cmc,'MarkerEdgeColor','w','MarkerSize',5);
% plot contrast effect on contrast discrimination
ccc = amap(maps(2,2),:);
errbar(bins(2,:)',con_con_m,con_con_s,'-','Color',ccc);
p(2) = plot(bins(2,:),con_con_m,'o','MarkerFaceColor',ccc,'MarkerEdgeColor','w','MarkerSize',5);

title('Contrast discrimination');
xlabel('\Delta contrast / coherence (%)');
ylabel('Proportion right choice (%)');
drawPublishAxis('figSize=[5,3]');

savepdf(h,fullfile(datafolder,'avg_behav','con_disc.pdf'));

h = figure; hold on

% plot contrast effect on coherence discrimination
ccc = amap(maps(2,1),:);
errbar(bins(2,:)',con_coh_m,con_coh_s,'-','Color',ccc);
p(2) = plot(bins(2,:),con_coh_m,'o','MarkerFaceColor',ccc,'MarkerEdgeColor','w','MarkerSize',5);
% coherence effect on coherence discrimination
cmc = amap(maps(1,1),:);
errbar(bins(1,:)',coh_coh_m,coh_coh_s,'-','Color',cmc);
p(1) = plot(bins(1,:),coh_coh_m,'o','MarkerFaceColor',cmc,'MarkerEdgeColor','w','MarkerSize',5);

title('Coherence discrimination');
xlabel('\Delta contrast / coherence (%)');
ylabel('Proportion right choice (%)');
drawPublishAxis('figSize=[5,3]');
savepdf(h,fullfile(datafolder,'avg_behav','coh_disc.pdf'));

%% Big panels concon and cohcoh broken down by pedestal
con_m = squeeze(allData_m(2,2,:,:));
coh_m = squeeze(allData_m(1,1,:,:));
con_s = squeeze(allData_s(2,2,:,:));
coh_s = squeeze(allData_s(1,1,:,:));

% contrast contrast
h = figure; hold on
condition = 2;
amap = brewermap(11,'PuOr');
cmap = amap([2 3 4 5],:);
for ped = 1:4
    x = bins(2,:);
    c_m = con_m(ped,:);

    xshort = x(2:end-1);
    c_mshort = c_m(2:end-1);
    
    b = [ones(size(xshort')) xshort']\c_mshort';
    bs = glmfit(bins(condition,:),squeeze(data_m(condition,ped,:)),'normal','link','probit');
    
    plot(xshort,0.5+xshort*b(2),'-','Color',cmap(ped,:));
    mult = .001;
    errbar(x+mult*ped,squeeze(con_m(ped,:)),squeeze(con_s(ped,:)),'-','Color',cmap(ped,:));
    p(ped) = plot(x+mult*ped,c_m,'o','MarkerFaceColor',cmap(ped,:),'MarkerEdgeColor','w','MarkerSize',5);
%         title(sprintf('%02.1f%%, slope: %01.2f',100*eval(sprintf('%ss(ped)',conditions{condition})),bs(flip(cii),2)));
    legs{ped} = sprintf('Pedestal: %02.1f%%',100*eval(sprintf('%ss(ped)',conditions{condition})));

end

axis([-0.1 0.1 0 1]);
set(gca,'XTick',[-0.1 0 0.1],'XTickLabel',[-0.1 0 .1]*100);
set(gca,'YTick',[0 0.5 1],'YTickLabel',[0 0.5 1]*100);
l = legend(p,legs,'FontSize',7,'FontName','Helvetica');
set(l,'box','off');
xlabel(sprintf('\\Delta %s (R-L, %%)',conditions{condition}));

drawPublishAxis('figSize=[9,3]');
savepdf(h,fullfile(datafolder,'avg_behav','con_ped.pdf'));


% coherence coherence
h = figure; hold on
condition = 1;
amap = brewermap(11,'PuOr');
cmap = amap([10 9 8 7],:);
for ped = 1:4
    x = bins(1,:);
    c_m = coh_m(ped,:);

    xshort = x(2:end-1);
    c_mshort = c_m(2:end-1);
    
    b = [ones(size(xshort')) xshort']\c_mshort';
    bs = glmfit(bins(condition,:),squeeze(data_m(condition,ped,:)),'normal','link','probit');
    
    plot(xshort,0.5+xshort*b(2),'-','Color',cmap(ped,:));
    mult = .001;
    errbar(x+mult*ped,squeeze(coh_m(ped,:)),squeeze(coh_s(ped,:)),'-','Color',cmap(ped,:));
    p(ped) = plot(x+mult*ped,c_m,'o','MarkerFaceColor',cmap(ped,:),'MarkerEdgeColor','w','MarkerSize',5);
%         title(sprintf('%02.1f%%, slope: %01.2f',100*eval(sprintf('%ss(ped)',conditions{condition})),bs(flip(cii),2)));
    legs{ped} = sprintf('Pedestal: %02.1f%%',100*eval(sprintf('%ss(ped)',conditions{condition})));

end

axis([-0.3 0.3 0 1]);
set(gca,'XTick',[-0.3 0 0.3],'XTickLabel',[-0.3 0 .3]*100);
set(gca,'YTick',[0 0.5 1],'YTickLabel',[0 0.5 1]*100);
l = legend(p,legs,'FontSize',7,'FontName','Helvetica');
set(l,'box','off');
xlabel(sprintf('\\Delta %s (R-L, %%)',conditions{condition}));

drawPublishAxis('figSize=[9,3]');
savepdf(h,fullfile(datafolder,'avg_behav','coh_ped.pdf'));