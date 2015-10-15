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
    
    ds{fi} = loadroi(d,neural.ROIs,'keepNAN=true');

%     concatInfo = viewGet(view,'concatInfo');
% 
% %     tSeries = {}; % meaned tSeries
%     fTSeries = {}; % all voxels surviving r^2 cutoff
% 
%     scanDims = viewGet(view,'scanDims');
end
cd(cdir);

%% Step 2: Concatenate all the ROIs
% this combines the voxels across sessions
% neural2 = cc_concatER(neural,name,'SCM_f',false);

%% Step 2: Get Instances
% now we want to get the average response to every single condition that
% exists in our experiment. We use 'get instances' to do this. Then we want
% to combine this into a giant fucking long-form dataset. This shit is
% GIANT! FUCKING HUGE! mostly because there are a billion voxels.
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
    
    lInst = getInstances(view,lROIs,neural.SCM_f.(folderz).lStim.stimVol,'startLag=9','blockLen=10');
    rInst = getInstances(view,rROIs,neural.SCM_f.(folderz).rStim.stimVol,'startLag=9','blockLen=10');

    for li = 1:length(lInst)
        inst = lInst{li};
        
        prefixs = prefixes{1};
        sN = neural.SCM_f.(folderz).lStim.stimNames;
        sV = neural.SCM_f.(folderz).lStim.stimVol;
        tsN = neural.SCM_f.(folderz).lStim.taskNames;
        tsV = neural.SCM_f.(folderz).lStim.taskSV;
        [con,coh,task,~] = parseSCM(sN,sV,tsN,tsV,prefixs,true);
        
        roi_pos = find(cellfun(@(x) strcmp(x,inst.name),neural.ROIs));
        
        sv = inst.classify.stimvol;
        amps = inst.classify.instances;
        
        for ci=1:length(sv)
            csv = sv{ci};
            amp = amps{ci};
            ccon = con(ci);
            ccoh = coh(ci);
            ctask = task(ci);
            for ni = 1:length(csv)
                if ~isempty(ni)
                    stim = csv(ni);
                    camp = amp(ni,:);
                    for ai = 1:length(camp)
                        a = amp(ai);

                        % roi - voxel - con - coh - task - amplitude - instance
                        % - folder
                        dat = [roi_pos ai ccon ccoh ctask a ni fi];
                        data(count,:) = dat;
                        count = count+1;

                    end
                end
            end
        end
    end
end

%% Step 3: Cross-Validate and Decode
% we are going to subset our data by taking only some of the instances for
% each roi, then we are going to build a linear regression model that
% predicts amplitude from con and coh. We will then use that model to
% predict the left-out dataset, our prediction accuracy (




