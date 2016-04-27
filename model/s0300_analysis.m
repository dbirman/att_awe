%% Joint Analysis for s0025

%% Load Files
fdir = '~/data/cohcon_localizer'; 
files = dir(fullfile(fdir,'s0300*.mat'));

ROIs = {'V1','V2','V3','V4','V3a','V3b','V7','LO1','LO2','MT'};

analyses = {'CohxCon','Timing'};
%% Load each file
datas = {};

for fi = 1:length(files)
    load(fullfile(fdir,files(fi).name));
    datas{fi} = data;
end

time = datas{fi}.deconvo{1}.time;
time2 = datas{fi}.deconvo{1}.time2;

%% Move everything into these convenient hodlers
deconvo = {};
basecon = [];
basecoh = [];
con = [];
coh = [];
timing = [];
idx_groups = {};
count = 1;

for fi = 1:length(datas)
    for ri = 1:length(ROIs)
        if fi == 1
            deconvo{ri}.ehdr = [];
            deconvo{ri}.ehdrste = [];
            deconvo{ri}.ehdr2 = [];
        end
        deconvo{ri}.ehdr = [deconvo{ri}.ehdr ;datas{fi}.deconvo{ri}.ehdr];
        deconvo{ri}.ehdrste = [deconvo{ri}.ehdrste ;datas{fi}.deconvo{ri}.ehdrste];
        deconvo{ri}.ehdr2 = [deconvo{ri}.ehdr2 ;datas{fi}.deconvo{ri}.ehdr2];
    end
    idx_groups{fi} = count:(count-1+size(datas{fi}.deconvo{1}.ehdr,1));
    count = count + size(datas{fi}.deconvo{1}.ehdr,1);
    if length(datas{fi}.basecon)==1
        basecon = [basecon ;repmat(datas{fi}.basecon,size(datas{fi}.deconvo{1}.ehdr,1),1)];
    else
        basecon = [basecon;datas{fi}.basecon];
    end
    if length(datas{fi}.basecoh)==1
        basecoh = [basecoh ;repmat(datas{fi}.basecoh,size(datas{fi}.deconvo{1}.ehdr,1),1)];
    else
        basecoh = [basecoh;datas{fi}.basecoh];
    end
    if length(datas{fi}.con)==1
        con = [con ;repmat(datas{fi}.con,size(datas{fi}.deconvo{1}.ehdr,1),1)];
    else
        con = [con;datas{fi}.con];
    end
    if length(datas{fi}.coh)==1
        coh = [coh ;repmat(datas{fi}.coh,size(datas{fi}.deconvo{1}.ehdr,1),1)];
    else
        coh = [coh;datas{fi}.coh];
    end
    if length(datas{fi}.timing)==1
        timing = [timing ;repmat(datas{fi}.timing,size(datas{fi}.deconvo{1}.ehdr,1),1)];
    else
        timing = [timing;datas{fi}.timing];
    end
end

%% Remove conditions with con=0 coh=0
remove_idxs = logical(logical((basecon-con)==0) .* logical((basecoh-coh)==0));
if any(remove_idxs)
    warning('You included conditions where con=0 and coh=0 (i.e. no change occurred), the model can''t account for these so they will be removed');  
    hrf = hrf(~remove_idxs,:);
    con = con(~remove_idxs);
    coh = coh(~remove_idxs);
    basecon = basecon(~remove_idxs);
    basecoh = basecoh(~remove_idxs);
    timing = timing(~remove_idxs);
end

%% New analysis: fit HRF

% we will model an HRF (lag/peak/etc) using a single gamma function and
% then fit it across all the deconvolutions, using the model:
% HRF = conv(baseHRF, filter);
% Where filter is an exponential weighted (e^-lambda*x) design matrix.
%
% Need to fit:
% lambda: exponent constant
% [amp tau timelag offset exponent]

% Do this per ROI
bestfits = {};
for ri = 1:length(deconvo)
    %%
    decon = deconvo{ri};
    
    ehdr = decon.ehdr2;
    
    % build conv_timing
%     conv_timing = [1 2 4 8 16 1 2 4 8 16]';
%     basecon = 0.25;
%     basecoh = 0;
%     con = 0.25;
%     coh = [0.25 0.25 0.25 0.25 0.25 1 1 1 1 1]';

    incl = 1:2;
    idxs = [];
    for i = incl
        idxs = [idxs idx_groups{i}];
    end
        
    bestfits{ri} = ccVoxelModel(ehdr(idxs,:),basecon(idxs),basecoh(idxs),con(idxs),coh(idxs),timing(idxs),time2);
    
end

%% Plot
for ri = 1:length(bestfits)
    bestfit = bestfits{ri};
    %%
    figure, hold on
    plot(bestfit.full.fcon,bestfit.full.fconr,'r');
    plot(bestfit.full.fcoh,bestfit.full.fcohr,'b');
    
    %%
    clist = brewermap(20,'PuOr');
    res = [0.25 0.5 1 2 4];
%     if ~isdir(fullfile(pwd,'Figures/linear')), mkdir(fullfile(pwd,'Figures/linear')); end

    figure, hold on
    
    count = 1;
    for i = incl
        subplot(length(datas),1,i), hold on
        for j = 1:size(datas{i}.deconvo{ri}.ehdr,1)
            % plot original
            plot(time,datas{i}.deconvo{ri}.ehdr(j,:),'o','MarkerFaceColor',clist(j,:),'MarkerEdgeColor',[1 1 1]);
            h = errbar(time,datas{i}.deconvo{ri}.ehdr(j,:),datas{i}.deconvo{ri}.ehdrste(j,:),'Color',clist(j,:));
        end
        
        for k = 1:size(datas{i}.deconvo{ri}.ehdr,1)
            % plot the fit
            plot(time2,bestfit.out(count,:),'-','Color',clist(k,:));
            count = count + 1;
        end
        title(analyses{i})
    end
    axis([0 30 -0.5 1.5]);
    xlabel('Time (s)');
    ylabel('BOLD Amplitude (%s Signal Change)');
    title(sprintf('25%% Coherence, R^2: %0.2f',bestfit.r2));
    %     drawPublishAxis;

    %%
    if ~isdir(fullfile(pwd,'Figures')), mkdir(fullfile(pwd,'Figures')); end
    if ~isdir(fullfile(pwd,'Figures/linear')), mkdir(fullfile(pwd,'Figures/linear')); end
    fname = fullfile(pwd,'Figures/linear',sprintf('%s_lintiming.pdf',ROIs{ri}));
    print(fname,'-dpdf');
end