function neural = cc_decoding(neural,sid)

stop = 1;


%% Step 1: Load Data
folders = neural.folders;

if ~isfield(neural,'tSeries')
    neural.tSeries = struct;
end

preF = fullfile('~/data/cohcon/');
cdir = pwd;
name = 'full';
ds = {};
for fi = 1:length(folders)
    %%
    folder = folders{fi};
    folderz = sprintf('f%s',folder);
    fullFolder = fullfile(preF,sprintf('%s%s',sid,folder));
    
    mrQuit;
    cd(fullFolder);
    
    %% Setup a view, don't load ROIs
    view = newView();
    view = viewSet(view,'curGroup','Concatenation');
    view = viewSet(view,'curScan',neural.concatScan);
    view = loadAnalysis(view,sprintf('erAnal/%s','both_ER'));
    analysis = viewGet(view,'analysis');
    d = analysis.d{neural.concatScan};
    d.scanNum = neural.concatScan;
    d.groupNum = view.curGroup;
    
    
    d = loadroi(d,neural.ROIs,'keepNAN=true');
    
    r2 = viewGet(view,'overlayData',d.scanNum);
    
    d.roi = getSortIndex(view,d.roi,r2);
    
    ds{fi} = d;
    % 

    
end
cd(cdir);

%% Step 2: Get Instances
% now we want to get the average response to every single condition that
% exists in our experiment. We use 'get instances' to do this. Then we want
% to combine this into a giant long-form dataset.
%
% Header
% roi - voxel - con - coh - task - amplitude - instance
header = {'roi','voxel','con','coh','task','amplitude','instance'};
data = []; count = 1;
prefixes = {{'rCon','rCoh'},{'lCon','lCoh'}};
for fi = 1:length(folders)
    folderz = sprintf('f%s',folders{fi});
    d = ds{fi};
    
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
%     
%     lInst = getInstances(view,lROIs,neural.SCM.(folderz).lStim.stimVol,'startLag=6','blockLen=10','minResponseLen=4');
%     rInst = getInstances(view,rROIs,neural.SCM.(folderz).rStim.stimVol,'startLag=6','blockLen=10','minResponseLen=4');
    lInst = getInstances(view,lROIs,neural.SCM_f.(folderz).lStim.stimVol,'type','glm','r2cutoff=0.07','hdrlen=15');
    rInst = getInstances(view,rROIs,neural.SCM_f.(folderz).rStim.stimVol,'type','glm','r2cutoff=0.07','hdrlen=15');

    %% calculate left/right data
    disp('Starting left...');
    disppercent(-1/length(lInst));
    for li = 1:length(lInst)
        inst = lInst{li};
        
        sN = neural.SCM_f.(folderz).lStim.stimNames;
        sV = neural.SCM_f.(folderz).lStim.stimVol;
        tsN = neural.SCM_f.(folderz).lStim.taskNames;
        tsV = neural.SCM_f.(folderz).lStim.taskSV;
        iV = inst.classify.instanceVol;
        amps = inst.classify.instances;
        [con,coh,task, sv] = sv2long(sV,sN,tsV,tsN);
        
        roi_pos = find(cellfun(@(x) strcmp(x,inst.name),neural.ROIs));
        croi = neural.ROIs{roi_pos};
        [s, r] = parseROI(croi,neural.shortROIs);
        
        
        for ci=1:length(amps)
            if ~isempty(amps)
                camps = amps{ci};
                for i = 1:size(camps,1)
                    camp = camps(i,:);
                    vol = iV{ci}(i);
                    
                    sv_pos = find(sv==vol);
                    if ~isempty(sv_pos)
                        ccon = con(sv_pos);
                        ccoh = coh(sv_pos);
                        ctask = task(sv_pos);
                        
                        for ai = 1:length(camp) % for each voxel
                            
                            a = camp(ai);
                            
                            % roi - voxel - con - coh - task - amplitude - instance
                            % - folder - side
                            dat = [r ai ccon ccoh ctask a i fi s];
                            data(count,:) = dat;
                            count = count+1;
                        end
                    end
                end
                
            end
        end
        
        disppercent(li/length(lInst));
    end
    disppercent(inf);
    
    disp('Starting right...');
    disppercent(-1/length(rInst));
    for li = 1:length(rInst)
        inst = rInst{li};
        
        sN = neural.SCM_f.(folderz).rStim.stimNames;
        sV = neural.SCM_f.(folderz).rStim.stimVol;
        tsN = neural.SCM_f.(folderz).rStim.taskNames;
        tsV = neural.SCM_f.(folderz).rStim.taskSV;
        iV = inst.classify.instanceVol;
        amps = inst.classify.instances;
        [con,coh,task, sv] = sv2long(sV,sN,tsV,tsN);
        
        roi_pos = find(cellfun(@(x) strcmp(x,inst.name),neural.ROIs));
        croi = neural.ROIs{roi_pos};
        [s, r] = parseROI(croi,neural.shortROIs);
        
        
        for ci=1:length(amps)
            if ~isempty(amps)
                camps = amps{ci};
                for i = 1:size(camps,1)
                    camp = camps(i,:);
                    vol = iV{ci}(i);
                    
                    sv_pos = find(sv==vol);
                    if ~isempty(sv_pos)
                        ccon = con(sv_pos);
                        ccoh = coh(sv_pos);
                        ctask = task(sv_pos);
                        
                        for ai = 1:length(camp) % for each voxel
                            
                            a = camp(ai);
                            
                            % roi - voxel - con - coh - task - amplitude - instance
                            % - folder - side
                            dat = [r ai ccon ccoh ctask a i fi s];
                            data(count,:) = dat;
                            count = count+1;
                        end
                    end
                end
                
            end
        end
    end
    disppercent(li/length(rInst));
end
disppercent(inf);

%% Write CSV (then can be used for R analyses)

                            % roi - voxel - con - coh - task - amplitude - instance
                            % - folder - side
header = {'roi','voxel','con','coh','task','amplitude','instance','folder','side'};
fname = sprintf('~/proj/att_awe/analysis/data/%s_voxels.csv',sid);
csvwriteh(fname,data,header);

%% Step 3: Lasso Decoding

% data is obs x 8
% roi - voxel - con - coh - task - amplitude - instance - folder - side

rois = neural.shortROIs;

figure
% first do attending contrast
for ri = 1:length(rois)
    subplot(1,length(rois),ri), hold on
    bcon = [];
    bcoh = [];
    pvcon = [];
    pvcoh = [];
    disp(sprintf('Current ROI: %s',rois{ri}));
    dat = select(data,1,ri);
    dat = select(dat,5,1); % contrast task
    voxels = unique(dat(:,2));
    for vi = 1:length(voxels)
        voxel = voxels(vi);
        dat2 = select(dat,2,voxel);
        X = dat2(:,3:4);
%         X = [ones(size(X,1),1),X]; % intercept
%         Y = dat2(:,6)*100-100; % reduce range to percentage
%         increase/decrease (only if using mean)
        Y = dat2(:,6);

        if ~any(isnan(Y))
            if any(Y>10) || any(Y<-10)
                X = X(Y<10,:);
                Y = Y(Y<10,:);
                X = X(Y>-10,:);
                Y = Y(Y>-10,:);
            end
            lm = fitlm(X,Y);
    %         B = X\Y;
            int = lm.Coefficients.Estimate(1);
            Bcon = lm.Coefficients.Estimate(2);
            Bcoh = lm.Coefficients.Estimate(3);
            bcon = [bcon Bcon];
            pvcon = [pvcon lm.Coefficients.pValue(2)];
            bcoh = [bcoh Bcoh];
            pvcoh = [pvcoh lm.Coefficients.pValue(3)];
        end
    end
%     X = [ones(size(bcon,2),1),bcon']; Y = bcoh';
%     Broi = X\Y;
    lm = fitlm(bcon',bcoh');
    x = -5:5;
    for i = 1:length(bcon)
        if pvcon(i) < .05 && pvcoh(i) < .05
            plot(bcon(i),bcoh(i),'*');
        else
            plot(bcon(i),bcoh(i),'o');
        end
    end
    text(-0.5,1,sprintf('sd: %1.2f',std(bcoh)));
    text(0.5,-1,sprintf('sd: %1.2f',std(bcon)));
    plot(x,x,'-.k');
    axis([-2 2 -2 2]);
    plot(x,lm.Coefficients.Estimate(1)+lm.Coefficients.Estimate(2)*x,'-.r');
    title(sprintf('%s: B = %1.2f, p = %1.3f',rois{ri},lm.Coefficients.Estimate(2),lm.Coefficients.pValue(2)));
end


function data = select(data,col,value)

scol = data(:,col);
data = data(scol==value,:);