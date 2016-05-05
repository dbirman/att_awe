function savedata_contrastxcoherence( cfolder, thresh, analysis)



%% Pick ROIs we will generate figures for
mrQuit
cd(fullfile('~/data/cohcon_localizer/',cfolder));
if ~isdir(fullfile(pwd,'Figures')), mkdir(fullfile(pwd,'Figures')); end
% ROIs = {'lV1','lV2v','lV2d','lV3v','lV3d','lMT','rV1','rV2v','rV2d','rV3v','rV3d','rMT','lV3a','lV3b','rV3a','rV3b','lLO1','lLO2','lV4','rLO1','rLO2','rV4'};
% shortROIs = {'V3a','V3b','V1','V2','V3','MT','LO1','LO2','V4'};

ROIs = {'V1','V2','V3','V4','V3a','V3b','V7','LO1','LO2','MT'};
% ROIs = {'V1','MT'};

%% Setup a view
view = newView();
view = viewSet(view,'curGroup','Concatenation');
view = viewSet(view,'curScan',1); % make sure scan # is correct
%% Generate SCM
[stimvol, stimNames, ~] = getStimvol(view,'contrast_x_coherence','taskNum=1','phaseNum=2');

%% Load data and threshold
view = loadAnalysis(view,sprintf('erAnal/%s',analysis)); % check analysis name!
analysis = viewGet(view,'analysis');
d = analysis.d{1};
d.scanNum = 1;
d.groupNum = view.curGroup;
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

%% Get deconvolution
deconvo = cell(1,length(d.roi));
for ri = 1:length(d.roi)
    curd = constructD(tSeries{ri},stimvol,0.5,15,concatInfo,'none','deconv',0);
    decon = getr2timecourse(curd.timecourse,curd.nhdr,curd.hdrlenTR,curd.scm,curd.framePeriod,curd.verbose);
    decon = rmfield(decon,'scm');
    decon = rmfield(decon,'covar');
    
    deconvo{ri} = decon;
end

%% Compute rdata
[con,coh] = parseNames(stimNames,'contrast=','coherence=',' and '); % get individual contrast/coherence values
cons = unique(con); cohs = unique(coh);

rdata = {};
rdatas = {};
for ri = 1:length(ROIs)
    % get amplitude data
    rdata_ = zeros(length(cons),length(cohs));
    rdata_ste = zeros(length(cons),length(cohs));
    for coni = 1:length(cons)
        for cohi = 1:length(cohs)
            ccon = cons(coni);
            ccoh = cohs(cohi);
            
            idx = logical(logical(con==ccon).*logical(coh==ccoh));
            rdata_(coni,cohi) = mean(deconvo{ri}.ehdr(idx,8:12));
            rdata_ste(coni,cohi) = mean(deconvo{ri}.ehdrste(idx,8:12));
        end
    end
    rdata{ri} = rdata_;
    rdatas{ri} = rdata_ste;
end

%% Interpolate deconvolutions into 0.25s space (easier for convolutions)

for di = 1:length(deconvo)
    decon = deconvo{di};
    
    decon.time2 = decon.time(1):.25:decon.time(end);
    decon.ehdr2 = zeros(size(decon.ehdr,1),size(decon.time2,2));
    for wi = 1:size(decon.ehdr,1)
        decon.ehdr2(wi,:) = interp1(decon.time,decon.ehdr(wi,:),decon.time2,'linear');
    end
    deconvo{di} = decon;
end

%% Save data
savefolder = '~/data/cohcon_localizer';
fname = sprintf('%s.mat',cfolder);
fname = fullfile(savefolder,fname);

[con,coh] = parseNames(stimNames,'contrast=','coherence=',' and '); % get individual contrast/coherence values

% you need: deconvo, basecon, basecoh, con, coh, timing
data.pre.tSeries = tSeries;
data.pre.stimvol = stimvol;
data.pre.stimNames = stimNames;
data.pre.concatInfo = concatInfo;
data.pre.TR = 0.5;
data.ROIs = ROIs;
data.deconvo = deconvo;
data.basecon = 0.25;
data.basecoh = 0;
data.con = con';
data.coh = coh';
data.timing = repmat(10,size(data.con));
save(fname,'data');
