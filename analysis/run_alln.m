

%%
%%%%%%%%%%%%%%%%%%%%% NEURAL %%%%%%%%%%%%%%%%%%%%%%%%%
nSIDs = [305 329 43 25 300 346 343 344 338 340 348];
rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
ROIs = rois;
%% Delete old files
% warning('Are you sure?');
% keyboard
% for ni = 1:length(nSIDs)
%     files = dir(fullfile(datafolder,sprintf('s%04.0f_*.mat',nSIDs(ni))));
%     for fi = 1:length(files)
%         delete(fullfile(datafolder,files(fi).name));
%     end
% end
% delete(fullfile(datafolder,'sigmas.mat'));

%% Deconvolve
clear rr2 asd

parfor ni = 1:length(nSIDs)
    files = dir(fullfile(datafolder,sprintf('s%04.0f*_data.mat',nSIDs(ni))));
    datas = {};
    for fi = 1:length(files)
        loaded = load(fullfile(datafolder,files(fi).name));
        datas{end+1} = loaded.data;
    end
    
%     if 0
    local_rr2 = zeros(length(rois),4,2);
    local_sd = zeros(length(rois),4,2,20);
    for ri = 1:length(rois)
        try
        [local_rr2(ri,:,:), local_sd(ri,:,:,:)] = deconCohCon(datas,rois{ri},sprintf('s%04.0f',nSIDs(ni)));
        catch
            disp(sprintf('Failed for %i',ni));
        end
    end
    rr2(ni,:,:,:) = local_rr2;
    asd(ni,:,:,:,:) = local_sd;
%     end
end

save(fullfile(datafolder,'decon_r2tsnr.mat'),'rr2','asd');
% save(fullfile(datafolder,'decon_r2.mat'),'rr2');
clear ri

%% 
load(fullfile(datafolder,'decon_r2tsnr.mat'));

%% Compute tSNR statistics
asd_ = squeeze(asd([1 2 6:end],:,4,:,:));
asd_ = reshape(asd_,8,8,40); % collapse over cc and time
asd_ = squeeze(mean(asd_,3)); % collapse over conditions

% compute confidence intervals over tSNR
ci = bootci(1000,@mean,asd_);
mu = mean(asd_(:));
disp(sprintf('Mean %1.2f',mu));
for ri = 1:8
    disp(sprintf('%s %1.2f, 95%% CI [%1.2f %1.2f]',rois{ri},mean(ci(:,ri)),ci(1,ri),ci(2,ri)));
end

%% compute linear model analysis of tSNR
load(fullfile(datafolder,'avg_decon_deconv.mat'));
asd_ = squeeze(asd([1 2 6:end],:,4,:,:));

header = {'subj','roi','con','coh','sd'};

data = zeros(10000,5);
count = 1;
for ni = 1:8
    for ri = 1:8
        for cc = 1:20
            data(count,:) = [ni ri avg_decon.conidx(cc) avg_decon.cohidx(cc) asd_(ni,ri,1,cc)];
            count = count+1;
        end
        
        for time = 1:20
            data(count,:) = [ni ri avg_decon.timconidx(cc) avg_decon.timcohidx(cc) asd_(ni,ri,2,time)];
            count = count+1;
        end
    end
end
data = data(1:(count-1),:);

% add to table
ds = mat2dataset(data,'VarNames',header);
ds.subj = nominal(ds.subj);
ds.roi = nominal(ds.roi);

%% Modeling
intercept = fitlme(ds,'sd~1','DummyVarCoding','effects');
con = fitlme(ds,'sd~1+con','DummyVarCoding','effects'); % no improvement
coh = fitlme(ds,'sd~1+coh','DummyVarCoding','effects'); % no improvement
roi = fitlme(ds,'sd~1+roi','DummyVarCoding','effects');
subj = fitlme(ds,'sd~1+(1 | subj)','DummyVarCoding','effects');
full = fitlme(ds,'sd~1+con+coh+roi+(1 | subj)','DummyVarCoding','effects');
concoh = fitlme(ds,'sd~1+con+coh','DummyVarCoding','effects');
%% Check deconvolution r2
load(fullfile(datafolder,'decon_r2.mat'));
csvwriteh(fullfile(datafolder,'decon_r2.csv'),round(squeeze(mean(rr2,3))*1000)/1000,rois);

%% Average deconvolution plots
avg_deconPlot_cc(nSIDs,'nomodel,nozero');
avg_deconPlot_time(nSIDs,'nomodel,nozero');

%% Fit Models
indiv_hrffit_control;

%% Plots
avg_fmri_plots;

%% Indiv plots
indiv_fmri_plots;

%% Fit behavior using fMRI responses
indiv_fitbehav;

indiv_fitbehav_full;

%% Average over contrast/coherence (testing for Justin)

% plot the average across contrast and coherence (instead of individual
% functions)
avg_deconPlot_cc_avg(nSIDs,'nomodel,nozero');
close all

%% Time
avg_deconPlot_time(nSIDs,'nomodel,nozero');


%% Test Linearity (per subject)
% This is actually the exponent testing code, not sure why it's named like
% this
subj_linear_test;

%% Grid and Time Plots
grid_time_plots;
% grid_time_plots_extra;


%% Average HRF plot
t = 0.25:0.5:40.5;
hrfs = zeros(length(nSIDs),length(t));
for ni = 1:length(nSIDs)
    load(fullfile(datafolder,sprintf('s%04.0f_fithrf.mat',nSIDs(ni))));
    hrfs(ni,:) = cc_gamma(t,fithrf.params);
end
save(fullfile(datafolder,'avg_hrf.mat'),'hrfs');

%% test average hrf
% get parameters of avg hrf:
avg_params = fitavghrf(mean(hrfs));
figure; hold on
time = 0.25:0.5:40.5;
plot(mean(hrfs),'b');
plot(cc_gamma(time,avg_params),'r');

%% Plot individual hrfs and average hrf
h = figure; hold on
x = 0.25:.5:40.5;
for ni = 1:length(nSIDs)
    plot(x,hrfs(ni,:),'-k','LineWidth',1);
end
plot(x,cc_gamma(time,avg_params),'--r','LineWidth',2);
xlabel('Time (s)');
ylabel('Normalized signal change (%)');
for hi = 1:length(header)
    text(15,1-(hi*0.15),sprintf('%s: %1.2f',header{hi},avg_params.(header{hi})));
end
axis([0 41 -0.6 1]);
drawPublishAxis('figSize=[8.9,4.5]');
savepdf(h,fullfile(datafolder,'hrf','hrfFigure.pdf'));

%% Build table for paper

load(fullfile(datafolder,'avg_hrf_roidelta.mat'));
header1 = {'tau1','amp2','tau2'};
data1 = zeros(11,8,length(header1));
for ni = 1:length(nSIDs)
%     load(fullfile(datafolder,sprintf('s%04.0f_fithrf.mat',nSIDs(ni))));
    for ri = 1:8
        for hi = 1:length(header1)
            data1(ni,ri,hi) = afits{ni}{ri}.params.(header1{hi});
        end
    end
end

load(fullfile(datafolder,'avg_hrffits.mat'));

header2 = {'conRmax','conc50','cohalpha','cohkappa','offset'};
data2 = zeros(11,8,length(header2));
for ni = 1:length(nSIDs)
    for ri = 1:8
        for hi = 1:length(header2)
            data2(ni,ri,hi) = sfits{ni}{1,4}.roifit{ri}.params.(header2{hi});
        end
    end
end

%% store data in a single dataset to do statistics over parameters

header_ = [header1 header2];
data_ = cat(3,data1,data2);

for hi = 1:length(header_)
    % compute SD over subjects
    for ri = 1:8
        mu = mean(data_(:,ri,hi));
        sd = std(data_(:,ri,hi));
        disp(sprintf('%s %s: %1.2f+-%1.2f',rois{ri},header_{hi},mu,sd));
    end

    % compute SD over areas
%     disp('*********************');
%     for ni = 1:11
%         mu = mean(squeeze(data_(ni,:,hi)));
%         sd = std(squeeze(data_(ni,:,hi)));
%         disp(sprintf('%s %s: %1.2f+-%1.2f',rois{ri},header_{hi},mu,sd));
%     end
end

%% convert the data to longform (1 row = 1 observation)
data = zeros(10000,length(header1)+length(header2)+2);
count = 1;
for ni = 1:length(nSIDs)
    for ri = 1:8
        data(count,:) = [ni ri squeeze(data1(ni,ri,:))' squeeze(data2(ni,ri,:))'];
        count = count+1;
    end
end
data = data(1:(count-1),:);

header = [{'subj' 'roi'} header1 header2];
csvwriteh(fullfile(datafolder,'all_parameters.csv'),data,header);

% add to table
ds = mat2dataset(data,'VarNames',header);
ds.subj = nominal(ds.subj);
ds.roi = nominal(ds.roi);

export(ds,'File',fullfile(datafolder,'dataset.csv'),'Delimiter',',');

%% Mixed effects modeling
ps = zeros(length(header)-2,4);
keep = [2 2 4 4 4 2 2 4];
models = {};
for hi = 3:length(header)
    cprm = header{hi};
    lmes = cell(1,4);
    lmes{1} = fitlme(ds,sprintf('%s~1',cprm),'DummyVarCoding','effects');
    lmes{2} = fitlme(ds,sprintf('%s~subj',cprm),'DummyVarCoding','effects');
    lmes{3} = fitlme(ds,sprintf('%s~roi',cprm),'DummyVarCoding','effects');
    lmes{4} = fitlme(ds,sprintf('%s~subj+roi',cprm),'DummyVarCoding','effects');
%     compares = [[1 2];[1 3];[2 4];[3 4]];
%     for li = 1:4
%         c = compare(lmes{compares(li,1)},lmes{compares(li,2)});
%         p = double(c.pValue);
%         ps(hi-2,li) = p(2);
%     end
%     c = compare(lmes{2},lmes{4});
%     p = double(c.pValue);
%     ps(hi-2,4) = p(2);
    
    models{end+1} = lmes{keep(hi-2)};
end

%% Plot sensitivity
plot_sensitivity;

%% Sensitivity bootstrap for paper
load(fullfile(datafolder,'avg_sensitivity.mat'));

temp_sensitivity = sensitivity;
temp_sensitivity(:,:,2,:) = sense_param(:,:,2,:);
csensitivity = squeeze(temp_sensitivity(:,1,:,:));

% bootstrap 
rois = ROIs;

[ci,stat] = bootci(10000,@mean,csensitivity(:,1,:));

str = '';
for ri = 1:length(rois)
    % significance test
    sigs = [.05 .01 .001 .0001];
    for si = 1:length(sigs)
        if quantile(stat(:,ri),sigs(si))<0
            p(ri) = quantile(stat(:,ri),sigs(si));
            break;
        else
            p(ri) = sigs(si);
        end
    end
    
    str = strcat(str,sprintf(', %s = %1.2f',rois{ri},mean(ci(:,ri))));
end

% coherence
% compare V3A and MT to V1
v3a = csensitivity(:,1,5);
mt = csensitivity(:,1,8);
v1 = csensitivity(:,1,1);

dci = bootci(10000,@mean,v3a-v1)
mean(dci)
dci = bootci(10000,@mean,mt-v1)
mean(dci)

















%% Averaged response plots (per roi)
all_responseplots; % wtf does this do?
%% Average response function plots
avg_responseplots;
avg_responseplotsplot; % actual plot!

%% Parameter plots: for each subject * ROI draw their contrast/coherence function fit and label the parameter values

load(fullfile(datafolder,'avg_hrffits.mat'));
%%
cmap = brewermap(7,'PuOr');

h = figure;
for ni = 1:length(nSIDs)
    for ri = 1:8
        subplot(length(nSIDs),8,(ni-1)*8+ri); hold on
        x = 0:.01:1;
        y = conModel(x,sfits{ni}{4}.roifit{ri}.params);
        rmax = sfits{ni}{4}.roifit{ri}.params.conRmax;
        c50 = sfits{ni}{4}.roifit{ri}.params.conc50;
        plot(x,y,'-','Color',cmap(2,:));
        text(0.2,0.5,num2str(round(c50*100)/100));
        text(0.6,2,num2str(round(rmax*100)/100));
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        axis([0 1 0 3]);
        drawPublishAxis('figSize=[18.3,24]');
    end
end
savepdf(h,fullfile(datafolder,'avg_fmri','params.pdf'));


%% Copy HRF data to data.mat for on OSF repo 
data_repo = '~/proj/att_awe/analysis/osf_predict_func';

load(fullfile(datafolder,'avg_hrffits.mat'));

clear data
for ni = 1:length(nSIDs)
    fit = sfits{ni}{1,4};
    data{ni} = struct;
    data{ni}.id = nSIDs(ni);
    for ri = 1:8
        roi{ri} = struct;
        roi{ri}.cc = struct;
        roi{ri}.cc.resp = fit.roifit{ri}.cc.resp_25;
        roi{ri}.cc.con = fit.roifit{ri}.cc.con;
        roi{ri}.cc.coh = fit.roifit{ri}.cc.coh;
        roi{ri}.cc.time = fit.roifit{ri}.cc.time/2;
        roi{ri}.time = struct;
        roi{ri}.time.resp = fit.roifit{ri}.time.resp_25;
        roi{ri}.time.con = fit.roifit{ri}.time.con;
        roi{ri}.time.coh = fit.roifit{ri}.time.coh;
        roi{ri}.time.time = fit.roifit{ri}.time.time/2;
    end
    data{ni}.areas = ROIs;
    data{ni}.adata = roi;
end

save(fullfile(data_repo,'data.mat'),'data');

%% Copy parameters for prediction function to params.mat on OSF repo
data_repo = '~/proj/att_awe/analysis/osf_predict_func';

load(fullfile(datafolder,'avg_hrffits.mat'));

for ni = 1:length(nSIDs)
    fit = sfits{ni}{1,4};
    
    param = struct;
    
    for ri = 1:8
        param.(ROIs{ri}).conalpha = fit.roifit{ri}.params.conRmax;
        param.(ROIs{ri}).consigma = fit.roifit{ri}.params.conc50;
        param.(ROIs{ri}).conp = fit.roifit{ri}.params.conp;
        param.(ROIs{ri}).conq = fit.roifit{ri}.params.conq;
        param.(ROIs{ri}).cohalpha = fit.roifit{ri}.params.cohalpha;
        param.(ROIs{ri}).cohkappa = fit.roifit{ri}.params.cohkappa;
    end
    
    params{ni} = param;
    
end

save(fullfile(data_repo,'params.mat'),'params');