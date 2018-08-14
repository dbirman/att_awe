%%
coherences = [0.15 0.3 0.45 0.6];
thresh_vals = .05:.05:.65;
contrasts = [0.325 0.4 0.55 0.85];
%%
% combine data
control = zeros(length(sids),2,4);
x = 0:.001:1;
respfits = zeros(length(sids),2,8,length(x));
thresholds = zeros(length(sids),2,8,length(thresh_vals));
disppercent(-1/(length(sids)*4));
for si = 1:length(sids)
    load(fullfile(datafolder,sprintf('%s_data.mat',sids{si})));
    control(si,:,:) = data.control;
    
    % for the optimal behavior we need the response functions from each set
    % of fits data.fits{1:4}. Get them and extrapolate
    
    for i = 1:8
        respfits(si,1,i,:) = cohModel(x,data.fits{i}.params);
        respfits(si,2,i,:) = conModel(x,data.fits{i}.params);
        
        maxV = 5;
        con = []; coh = [];
        while isempty(con) || isempty(coh)
            try
                x2 = 0:.001:5;
                con = conModel(x2,data.fits{i}.params);
                coh = cohModel(x2,data.fits{i}.params);
            catch
                maxV = max(maxV-1,1);
            end
        end
        % while we're here, compute the thresholds 
        for ci = 1:length(thresh_vals)
            cvalue = thresh_vals(ci);

            % use the analytical derivative: because these functions are
            % bounded they sometimes can't converge (usually happens for
            % coherence not contrast)
%             ccoh = thresh_vals(ci);
%             if data.fits{i}.params.conmodel==3 && data.fits{i}.params.poissonNoise==0
%                 thresholds(si,2,i,ci) = 1 / (data.fits{i}.params.conalpha*data.fits{i}.params.conkappa*exp(-data.fits{i}.params.conkappa*ccon));
%                 thresholds(si,1,i,ci) = 1 / (data.fits{i}.params.cohalpha*data.fits{i}.params.cohkappa*exp(-data.fits{i}.params.cohkappa*ccoh));
%             elseif data.fits{i}.params.conmodel==3
%                 % poisson noise
%                 thresholds(si,2,i,ci) = sqrt(conModel(ccon,data.fits{i}.params)) / (data.fits{i}.params.conalpha*data.fits{i}.params.conkappa*exp(-data.fits{i}.params.conkappa*ccon));
%                 thresholds(si,1,i,ci) = 2*sqrt(cohModel(ccoh,data.fits{i}.params)) / (data.fits{i}.params.cohalpha*data.fits{i}.params.cohkappa*exp(-data.fits{i}.params.cohkappa*ccoh));
%             elseif data.fits{i}.params.poissonNoise==1
%                 % linear + poisson noise
%                 thresholds(si,2,i,ci) = sqrt(conModel(ccon,data.fits{i}.params)) / data.fits{i}.params.conslope;
%                 thresholds(si,1,i,ci) = 2*sqrt(cohModel(ccoh,data.fits{i}.params)) / data.fits{i}.params.cohslope;
%             else
%                 thresholds(si,2,i,ci) = 1 / data.fits{i}.params.conslope;
%                 thresholds(si,1,i,ci) = 1 / data.fits{i}.params.cohslope;
%             end
            
            if mod(i,2)==1
                poisson = false;
            else
                poisson = true;
            end
            thresholds(si,2,i,ci) = findDprime(cvalue,x2,con,poisson);
            thresholds(si,1,i,ci) = findDprime(cvalue,x2,coh,poisson);
        end
        disppercent(((si-1)*4+i)/(length(sids)*4));
    end
end
disppercent(inf);

control(control<=0) = NaN;
control(control>=0.6) = NaN;

%% Save respfits
save(fullfile(datafolder,'avg_respfits.mat'),'respfits');
%%

controls = squeeze(bootci(10000,@(x) nanmean(x,1),control));
mcontrol = squeeze(mean(controls));
controls = squeeze(controls(2,:,:))-mcontrol;
% mcontrol = squeeze(nanmean(control,1));
% controls = squeeze(nanstd(control,1));

%% add data
    map = brewermap(7,'PuOr');
    h = figure; hold on
avg_performancecontrol_helper

axis([0 1 0 0.3])
set(gca,'XTick',[0 0.125 0.25 0.5 1],'XTickLabel',[0 12.5 25 50 75]);
set(gca,'YTick',[0 0.1 0.2 0.3],'YTickLabel',[0 10 20 30]);
drawPublishAxis('figSize=[8.5,4.25]');

savepdf(h,fullfile(datafolder,'avg_behav','avg_performance.pdf'));

%% BW figure (for talk
    map = brewermap(7,'PuOr');
    map = repmat(0,7,3);
    h = figure; hold on
avg_performancecontrol_helper

axis([0 1 0 0.3])
set(gca,'XTick',[0 0.125 0.25 0.5 1],'XTickLabel',[0 12.5 25 50 75]);
set(gca,'YTick',[0 0.1 0.2 0.3],'YTickLabel',[0 10 20 30]);
xlabel('Base stimulus strength (%)');
ylabel('Just noticeable difference (%)');
drawPublishAxis('figSize=[12,8]');

savepdf(h,fullfile('~/proj/att_awe/talks/data_figures/bw_jnd.pdf'));
% savepdf(h,fullfile(datafolder,'avg_behav','avg_performance.pdf'));
%% Average thresholds
thresholds_ = squeeze(bootci(1000,@nanmean,thresholds));
thresholdz = squeeze(mean(thresholds_));
%% Resp fits
respfits_ = squeeze(bootci(1000,@nanmean,respfits));
respfitz = squeeze(mean(respfits_));
%% Now re-plot the same figure - + optimal behavior functions
    map = brewermap(7,'PuOr');
% but with optimal behavioral functions drawn on
models = {'con-exp,coh-exp','con-exp,coh-exp,poisson','con-linear,coh-linear','con-linear,coh-linear,poisson','con-naka,coh-naka','con-naka,coh-naka,poisson','con-explin,coh-explin','con-explin,coh-explin,poisson'}; %

% plot the average models

h = figure; hold on
% order = [3 4 2 1];
order = 1:2:7;
for ii = 1:4
    i = order(ii);
    subplot(length(order),2,(ii-1)*2+1); hold on
    
    
%     plot(x,squeeze(respfits(2,i,:)),'-','Color',map(2,:));
    if mod(i,2)==0
        boundedline(x,squeeze(respfitz(2,i,:))',sqrt(squeeze(respfitz(2,i,:)))','cmap',map(2,:));
        boundedline(x,squeeze(respfitz(1,i,:))',sqrt(squeeze(respfitz(1,i,:)))','cmap',map(6,:));
        axis([0 0.75 -25 125]);
        set(gca,'YTick',[0 100]);
    else
        boundedline(x,squeeze(respfitz(2,i,:))',ones(size(x)),'cmap',map(2,:));
        boundedline(x,squeeze(respfitz(1,i,:))',ones(size(x)),'cmap',map(6,:));
        axis([0 0.75 -2.5 17.5]);
        set(gca,'YTick',[0 10]);
    end
    
    title(models{i});
    
    set(gca,'XTick',[0 0.75],'XTickLabel',[0 75]);
    if ii==4
        xlabel('Base strength (%)');
        ylabel('Response (s.d.)');
    end
    drawPublishAxis('figSize=[8.5,14]');
end

for ii = 1:4
    i = order(ii);
    subplot(length(order),2,(ii-1)*2+2); hold on
    plot(thresh_vals,squeeze(thresholdz(2,i,:)),'-','Color',map(2,:));
    plot(thresh_vals,squeeze(thresholdz(1,i,:)),'-','Color',map(6,:));
    avg_performancecontrol_helper;

    % add optimal behavior

    axis([0 0.75 0 0.35])
    set(gca,'XTick',[0 0.75],'XTickLabel',[0 75]);
    set(gca,'YTick',[0 0.3],'YTickLabel',[0 30]);
    drawPublishAxis('figSize=[8.5,14]');
    
end
savepdf(h,fullfile(datafolder,'avg_models','optimal_models.pdf'));
