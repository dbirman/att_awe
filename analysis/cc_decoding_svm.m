function neural = cc_decoding_svm(neural,sid)

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

%% Combine l/r ROIs using stackROIs

roi_idxs = cell(1,length(ds));

for fi = 1:length(ds)
    d = ds{fi};
    cf_idxs = cell(1,40);
    maxidx = 0;
    for ri = 1:length(d.roi)
        roi_pos = find(cellfun(@(x) strcmp(x,d.roi{ri}.name),neural.ROIs));
        croi = neural.ROIs{roi_pos};
        [s, r] = parseROI(croi,neural.shortROIs);
        cf_idxs{r} = [cf_idxs{r} roi_pos];
        if r > maxidx, maxidx=r; end
    end
    cf_idxs = cf_idxs(1:maxidx);
    roi_idxs{fi} = cf_idxs;
end

%% Get ROI Stacks

aroi = {}; % this is where we will store the stacked rois

for ri = 1:length(neural.shortROIs)
    % iterate by ROIs on the outside, this is easier to deal with
    roi_stack = {}; % we will store all the ROIs that need to be stacked here
    for fi = 1:length(ds)
        d = ds{fi};
        % for each folder
        folder_idxs = roi_idxs{fi};
        if ~isempty(folder_idxs{ri})
            crois = folder_idxs{ri};
            for iii = 1:length(crois)
                roi_stack{end+1} = d.roi{crois(iii)};
            end
        end
    end
    % collapse roi_stack
    nroi = roi_stack{1};
    for i = 2:length(roi_stack)
        nroi = stackROI(nroi,roi_stack{i});
    end
    aroi{end+1} = nroi;
    aroi{end}.name = neural.shortROIs{ri};
end

%% Task: Coherence
inst = getInstances(view,aroi,neural.SCM_st.(folderz).stimVol,'type','glm','n=inf','r2cutoff=0.07','hdrlen=15');

%% Task: Contrast
% two sides, so chance is 50%
inst_con = inst;
for ii = 1:length(inst_con)
    inst_con{ii}.classify.instances = inst_con{ii}.classify.instances(1:2);
end
out_con = leaveOneOut(inst_con);

%% Task: Coherence
inst_coh = inst;
for ii = 1:length(inst_coh)
    inst_coh{ii}.classify.instances = inst_coh{ii}.classify.instances(3:4);
end
out_coh = leaveOneOut(inst_coh);

%% Pull out Data and Graph

con_ROI = {''};
coh_ROI = {};
con_dec = []; con_decCI = [];
coh_dec = []; coh_decCI = [];
for ii = 1:length(out_con)
    dec = out_con{ii}.classify.leaveOneOut.correct;
    decs = out_con{ii}.classify.leaveOneOut.correctSTE*1.96;
    con_dec = [con_dec dec];
    con_decCI = [con_decCI decs];
    con_ROI{end+1} = out_con{ii}.name;
end
for ii = 1:length(out_coh)
    dec = out_coh{ii}.classify.leaveOneOut.correct;
    decs = out_coh{ii}.classify.leaveOneOut.correctSTE*1.96;
    coh_dec = [coh_dec dec];
    coh_decCI = [coh_decCI decs];
    coh_ROI{end+1} = out_coh{ii}.name;
end
con_ROI{end+1} = '';

decode_data = struct;

decode_data.con_ROI = con_ROI;
decode_data.coh_ROI = coh_ROI;
decode_data.con_dec = con_dec;
decode_data.con_decCI = con_decCI;
decode_data.coh_dec = coh_dec;
decode_data.coh_decCI = coh_decCI;

dir = fullfile('~/proj/att_awe/analysis/data/decode');
if ~isdir(dir), mkdir(dir); end
fname = fullfile(dir,sprintf('%s_ddata.mat',sid));
save(fname,'decode_data');

%%
Y = [con_dec',coh_dec'];
figure, hold on
bar(1:size(Y,1),Y,1);
% drawPublishAxis
set(gca,'XtickL',con_ROI);
errorbar([(1:9)-0.15, (1:9)+0.15],[con_dec,coh_dec],[con_decCI,coh_decCI],'ok');
hline(0.5);
axis([0 10 0.4 1.0]);
legend('Attend Contrast','Attend Motion');
drawPublishAxis

%% Print

dir = fullfile('~/proj/att_awe/analysis/figures',sid);
fname = fullfile(dir,'decode_ROIs');
print(fname,'-dpdf');

