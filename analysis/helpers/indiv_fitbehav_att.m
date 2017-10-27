%% Fit the behavioral data (for each of 21 subjects)
% Use the average response functions across all subjects
nSIDs = [305 329 43 25 300 346 343 344 338 340 348];
bSIDs = [345 330 337 335 349 354 353 356 334 352]; % behavior only participants
aSIDs = [nSIDs bSIDs];
sids = {};
for ni = 1:length(aSIDs)
    sids{end+1} = sprintf('s%03.0f',aSIDs(ni));
end
rois = {'V1','V2','V3','V4','V3a','V3b','V7','MT'};
%% Population responses
% load all the linear fits, generate fit data, average across subjects,
% then run in fitCCBehavModel_sigma (to estimate sigma per subject), return
% BIC

load(fullfile(datafolder,'avg_att_fits.mat'));

x = 0:.001:1;
respcon = zeros(11,8,2,length(x));
respcoh = zeros(11,8,2,length(x));

for si = 1:10
    fit = afits{si}{1}; % 4 refers to resp-25, which is our standard model

    for ri = 1:8
        respcon(si,ri,1,:) = fit.roifit{ri}.conresp_coh;
        respcon(si,ri,2,:) = fit.roifit{ri}.conresp_con;
        respcoh(si,ri,1,:) = fit.roifit{ri}.cohresp_coh;
        respcoh(si,ri,2,:) = fit.roifit{ri}.cohresp_con;
    end
end

respcon_ = squeeze(mean(bootci(10000,@mean,respcon)));
respcoh_ = squeeze(mean(bootci(10000,@mean,respcoh)));

%% Check that responses look right
figure; hold on
for mi = 1
%     plot(squeeze(mean(respcon_([1 2 3 4],1,:),1)),'--r');
%     plot(squeeze(mean(respcoh_([5 8],1,:),1)),'--');
    plot(squeeze(mean(respcon_([1 2 3 4],1,:),1)),'-r');
    plot(squeeze(mean(respcoh_([5 8],1,:),1)),'-');
end

%% Lapse rate calculation

mopts = 1;

lapses = zeros(1,length(aSIDs));
ltrials = lapses;
for ai = 1:length(aSIDs)
    adata = loadadata(sprintf('s%03.0f',aSIDs(ai)));
    % fit the lapse rate parameter
    tdata = adata(adata(:,9)==-1,:);
    tdata(:,13) = abs(tdata(:,5)-tdata(:,4)); % contrast diff
    tdata(:,14) = abs(tdata(:,7)-tdata(:,6)); % coherence diff
    conbins = quantile(tdata(:,13),.95);
    cohbins = quantile(tdata(:,14),.95);
    conidxs = logical((tdata(:,1)==2) .* (tdata(:,13)>conbins));
    cohidxs = logical((tdata(:,1)==1) .* (tdata(:,14)>cohbins));
    lapse = 1-mean([tdata(conidxs,12) ; tdata(cohidxs,12)]);
    ltrials(ai) = size([tdata(conidxs,12) ; tdata(cohidxs,12)],1);
    lapses(ai) = lapse;
end
    
%% Fit to behavior

models = {'exp'};
% 'sigma','sigma,poisson',
bmodels = {'sigma,roi,att'};%,'doublesigma','doublesigma,poisson'};
% bmodels = {'sigma,roi'};

% options list
attopts = zeros(10000,5);
count = 1;

% build options 

ropts = {[1:8]};
% rconopts = {1, [1 2 3 4], 1:8};
% rcohopts = {8, [5 8], 1:8};

sigmaopts = linspace(.01,.5,6);

for ai = 1:length(aSIDs)
    for mi = 1:length(mopts)
        for ropt = 1:length(ropts)
            for ni = 1:length(bmodels)
                if strfind(bmodels{ni},'roi')
                    % roi models have sigma fixed, no need to do
                    % multiple
                    attopts(count,:) = [ai mi ni ropt 1];
                    count = count+1;
                else
                    for si = 1:length(sigmaopts)
                        attopts(count,:) = [ai mi ni ropt si];
                        count = count+1;
                    end
                end
            end
        end
    end
end
attopts = attopts(1:(count-1),:);

% break into 12*10 size chunks
breaks = 1:240:size(attopts,1);
breaks(end) = size(attopts,1)+1;

if length(breaks)==1
    breaks(2) = breaks(1); breaks(1) = 1;
end
    
%% fit all options
disppercent(-1/size(attopts,1));

% breaks = [breaks(1) breaks(end)];
attfits = cell(size(attopts,1),1);
wfits = cell(size(attopts,1),1);
for ni = 1:(length(breaks)-1)
    bstart = breaks(ni);
    bend = breaks(ni+1)-1;
    
    parfor ii = bstart:bend
        copt = attopts(ii,:);
        
        subj = copt(1);
        adata = loadadata(sprintf('s%03.0f',aSIDs(subj)));
        
        shape = copt(2);
        noise = copt(3);
        ropt = ropts{copt(4)};
        sigma = sigmaopts(copt(5));
        
        info = struct;
        info.sigma = sigma;
        info.model = bmodels{noise};
        info.rois = ropt;
        info.lapse = lapses(subj);
        info.respcon = respcon_;
        info.respcoh = respcoh_;
        
%         temps = cell{1,4};
%         parfor iii=1:4
%             
%             tinfo = info;
%             tinfo.model = bmodels{iii};
%             temps{iii} = fitCCBehavControlModel_fmri(adata,tinfo,1);
%         end
        attfits{ii} = fitCCBehavControlModel_fmri(adata,info,1);
   end
    
    disppercent(bend/size(attopts,1));
end
disppercent(inf);

save(fullfile(datafolder,'avg_indiv_fits_att.mat'),'attfits');
% save(fullfile(datafolder,'avg_within_fits.mat'),'wfits');
%     save(fullfile(datafolder,sprintf('avg_indiv_fits_%02.0f.mat',100*sigmaopts(si))),'afits');
%     disp('************************************');
%     disppercent(si/length(sigmaopts));
%     disp('************************************');
% end
% disppercent(inf);

%% Plot

plot_rightchoice_model_full;


%% Compare afits and attfits
restructure_afits;

for i = 1:21
    cd(i) = afits{i}{1}.cv.cd;
    cd_att(i) = attfits{i}.cv.cd;
end
