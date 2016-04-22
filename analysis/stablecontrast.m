%% Analyze the data from s030020160309: stablecon=1 run
cd('~/data/cohcon/s030020160309');
% ROIs = {'lV1','lV2v','lV2d','lV3v','lV3d','lMT','rV1','rV2v','rV2d','rV3v','rV3d','rMT','lV3a','lV3b','rV3a','rV3b','lLO1','lLO2','lV4','rLO1','rLO2','rV4'};
% shortROIs = {'V3a','V3b','V1','V2','V3','MT','LO1','LO2','V4'};

ROIs = {'lV1','lMT','rV1','rMT'};
shortROIs = {'V1','MT'};


%% Setup a view, don't load ROIs
view = newView();
view = viewSet(view,'curGroup','Concatenation');
view = viewSet(view,'curScan',1);

%% Generate SCM
[stimvolL, stimNamesL, ~] = getStimvol(view,'lCoh');
[stimvolR, stimNamesR, ~] = getStimvol(view,'rCoh');

%% Load data
view = loadAnalysis(view,sprintf('erAnal/%s','both_ER'));
analysis = viewGet(view,'analysis');
d = analysis.d{1};
d.scanNum = 1;
d.groupNum = view.curGroup;
thresh = 0.12;
d = loadroi(d,ROIs,'keepNAN=true');
concatInfo = viewGet(view,'concatInfo');
r2 = viewGet(view,'overlayData',d.scanNum);
d.roi = getSortIndex(view,d.roi,r2);

tSeries = cell(1,length(d.roi)); % meaned tSeries

scanDims = viewGet(view,'scanDims');

for ri = 1:length(d.roi)
    r = d.roi{ri};
    r.linearScanCoords = sub2ind(scanDims,r.scanCoords(1,:),r.scanCoords(2,:),r.scanCoords(3,:));
    
    r.r2 = r2(r.linearScanCoords);    

    tSeries{ri} = nanmean(r.tSeries(r.r2>thresh,:));
end

%%
deconvo = cell(1,length(d.roi));
for ri = 1:length(d.roi)
    if d.roi{ri}.name(1)=='l'
        cd = constructD(tSeries{ri},stimvolR,0.5,20,concatInfo,'none','deconv',0);
    else
        cd = constructD(tSeries{ri},stimvolL,0.5,20,concatInfo,'none','deconv',0);
    end
    decon = getr2timecourse(cd.timecourse,cd.nhdr,cd.hdrlenTR,cd.scm,cd.framePeriod,cd.verbose);
    decon = rmfield(decon,'scm');
    decon = rmfield(decon,'covar');
    
    deconvo{ri} = decon;
end

%% Non-linear fit

fitter = 'fitType=nonlin';
amper = 'amplitudeType=fit2';
fit = cell(1,length(d.roi));
for ri = 1:length(d.roi)
    if d.roi{ri}.name(1)=='l'
        fit{ri} = fitTimecourse(tSeries{ri},stimvolR,.5,'concatInfo',concatInfo,fitter,amper);
    else
        fit{ri} = fitTimecourse(tSeries{ri},stimvolL,.5,'concatInfo',concatInfo,fitter,amper);
    end
end

%% Convert to amps
amps = cell(1,length(d.roi));
% either:
for ri = 1:length(d.roi)
    amps{ri} = mean(deconvo{ri}.ehdr(:,10:16),2);
end
% or
for ri = 1:length(d.roi)
    amps{ri} = fit{ri}.amps;
end

%%
[con, coh, task] = parseNamesSimple(stimNamesR);

%% Plots?
close all
for ri = 1:length(d.roi)
%     if strfind(d.roi{ri}.name,'V1')
        figure
        decon= deconvo{ri};
        amps = mean(decon.ehdr(:,10:15),2);
%         amps = amps{ri};
%         amps = amps';
        cohs = coh';
        cohs(:,2) = 1;
%         ab = cohs\amps;
        hold on
        x = [0 1];
%         plot(x,x*ab(1)+ab(2));
        plot(coh,amps,'*');
        title(d.roi{ri}.name);
%     end
    drawPublishAxis
end

%% Get left/right ROIs
lindxs = []; rindxs = [];
for ri = 1:length(d.roi)
    if strcmp(d.roi{ri}.name(1),'l')
        lindxs = [lindxs ri];
    else
        rindxs = [rindxs ri];
    end
end

lROIs = d.roi(lindxs);
rROIs = d.roi(rindxs);
%% Try Decoding?
% We will take all voxels > r^2 cutoff, and then we will do a regularized
% regression fit to see if we can CV out the true coherence.

lInst = getInstances(view,lROIs,stimvolR,'startLag=8','blockLen=12','n=100');
rInst = getInstances(view,rROIs,stimvolL,'startLag=8','blockLen=12','n=100');

%% Have instances

% header: ROI voxel amplitude repeat coherence
data = inst2long(lInst, coh);

%% Cross-val regression

rois = unique(data(:,1));
for ri = 6
    maxV = max(data(:,2));

    data_r = sel(data,1,ri);
    % generate data of the form: coh ~ amp_v1 amp_v2 amp_v3 ...
    % we will subsample to train and cross-validate to see performance
    % get one repeat (should have one coherence)
    Y = [];
    X = zeros(10000,maxV); count = 1;
    % get one coherence
    cohs = unique(data_r(:,5));
    for cohi = 1:length(cohs)
        ccoh = cohs(cohi);
        data_c = sel(data_r,5,ccoh);

        % get one repetition
        reps = unique(data_c(:,4));
        for repi = 1:length(reps)
            rep = reps(repi);
            % looking at data from one trial
            data_t = sel(data_c,4,rep);
            % convert to long form
            Y(count) = ccoh;
            X(count,:) = data_t(:,3)';
            count = count + 1;
        end
    end
    X = X(1:count-1,:);
    CV = 5;
    for cv = 1:CV
        hold_idxs = rand(1,length(Y))<1/CV;
        tr_idxs = ~hold_idxs;
        
        [b,stats] = lasso(X(tr_idxs,:),Y(tr_idxs)','CV',3);
        b_ = b(:,stats.IndexMinMSE);
        % use lasso 
        % make predictions
        pred_Y = X(hold_idxs,:)*b_+stats.Intercept(stats.IndexMinMSE);
        % compute R^2
        r2 = 1-sum((Y(hold_idxs)-pred_Y').^2)/sum((Y(hold_idxs)'-mean(Y(hold_idxs))).^2);
        disp(r2)
    end
end