%% run_all_time
nSIDs = [305 329 43 25 300 346 343 344 338 340 348];

tSIDs = [300,343,346,043,338];

tcorrespond = [5 7 6 3 9];
sids = {};
for ti = 1:length(tSIDs)
    sids{end+1} = sprintf('s%03.0f',tSIDs(ti));
end

%% Load response functions
load(fullfile(datafolder,'avg_hrffits.mat'));

x = 0:.001:1;
respcon = zeros(11,length(mopts),8,length(x));
respcoh = zeros(11,length(mopts),8,length(x));

disppercent(-1/11);
for si = 1:11
    for mi = 1:length(mopts)
        fit = sfits{si}{mopts(mi),4}; % 4 refers to resp-25, which is our standard model
        
        for ri = 1:8
            respcon(si,mi,ri,:) = conModel(x,fit.roifit{ri}.params);
            respcoh(si,mi,ri,:) = cohModel(x,fit.roifit{ri}.params);
        end
    end
    disppercent(si/11);
end
disppercent(inf);

respcon_ = squeeze(mean(bootci(1000,@mean,respcon)));
respcoh_ = squeeze(mean(bootci(1000,@mean,respcoh)));

%% Run fitting procedure
for si = 1:5
    adata = loadadata_time(sids{si});

    fits{si} = fitCCBehavControlModel_time(adata,0,'sigma',squeeze(respcon_(1,1,:)),squeeze(respcoh_(1,8,:)),lapses(tcorrespond(si)),1);
%     fits{si} = fitCCBehavControlModel_time(adata,0,'sigma',squeeze(respcon(tcorrespond(si),2,1,:)),squeeze(respcoh(tcorrespond(si),1,8,:)),0,1);
end

save(fullfile(datafolder,'avg_time_fits.mat'),'fits');

%%
for si = 1:5
    h = figure;
    hold on
    plot(squeeze(respcon(tcorrespond(si),2,1,:)));
    plot(squeeze(respcoh(tcorrespond(si),1,8,:)));
end
%% Get R^2
% load(fullfile(datafolder,'avg_indiv_fits.mat'));
% for ai = 1:21
% end

load(fullfile(datafolder,'avg_time_fits.mat'));
clear r2_
for si = 1:5
    r2_(si,:,:) = fits{si}.r2_;
end

out = squeeze(mean(r2_));

out = out*100;

times = [250,500,1000];
for i = 1:3
    disp(sprintf('%i ms: %01.2f [%01.2f %01.2f]',times(i),mean(out(i,:)),out(i,1),out(i,2)));
end

