function neural = cc_decoding_svm(neural,sid,varargin)

zscore = []; type = []; folders = [];
getArgs(varargin,{'zscore=0','type=svm','folders=0','all=1'});


folders = neural.folders;
% folders = folders(1:2);

if ~isfield(neural,'tSeries')
    neural.tSeries = struct;
end

preF = fullfile('~/data/cohcon/');
cdir = pwd;
name = 'full';
doload = 1;

if ~doload
    %% Step 1: Load Data
    ds = {};
    for fi = 1:length(folders)
        %%
        folder = folders{fi};
        folderz = sprintf('f%s',folder);
        fullFolder = fullfile(preF,sprintf('%s%s',sid,folder));

        mrQuit;
        cd(fullFolder);

        %% Setup a view, load ROIs
        view = newView();
        view = viewSet(view,'curGroup','Concatenation');
        view = viewSet(view,'curScan',neural.concatScan);
        view = loadAnalysis(view,sprintf('erAnal/%s','both_ER'));    

        analysis = viewGet(view,'analysis');
        d = analysis.d{neural.concatScan};
        d.scanNum = neural.concatScan;
        d.groupNum = view.curGroup;
        d.concatInfo = viewGet(view,'concatInfo');


        d = loadroi(d,neural.ROIs,'keepNAN=true');

        r2 = viewGet(view,'overlayData',d.scanNum);

        d.roi = getSortIndex(view,d.roi,r2);
    %     
        ds{end+1} = d;
    end
    %% Clear stuff
    clear analysis

    %% Collapse ROIs by folder
    dAll = {}; % new struct to hold collapse folders

    for fi = 1:length(folders)
        % setup
        cf_idxs = cell(1,40);
        d = ds{fi};
        % find collapses
        maxidx = 0;
        for ri = 1:length(d.roi)
            roi_pos = find(cellfun(@(x) strcmp(x,d.roi{ri}.name),neural.ROIs));
            croi = neural.ROIs{roi_pos};
            [s, r] = parseROI(croi,neural.shortROIs);
            cf_idxs{r} = [cf_idxs{r} roi_pos];
            if r > maxidx, maxidx=r; end
        end
        cf_idxs = cf_idxs(1:maxidx);
        % actually do the collapsing
        aroi = {};
        for ri = 1:length(neural.shortROIs)
            % iterate by ROIs on the outside, this is easier to deal with
            roi_stack = {}; % we will store all the ROIs that need to be stacked here
            % for each folder
            if ~isempty(cf_idxs{ri})
                crois = cf_idxs{ri};
                for ci = 1:length(crois)
                    roi_stack{end+1} = d.roi{crois(ci)};
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
        % save data
        dAll{end+1} = struct;
        dAll{end}.roi = aroi;
        dAll{end}.concatInfo = d.concatInfo;
    end

    %% Clear space
    clear ds d

    %% Instances across ROIs
    insts = {};
    for fi = 1:length(folders)
        folder = folders{fi};
        folderz = sprintf('f%s',folder);
        insts{end+1} = getInstances(view,dAll{fi}.roi,neural.SCM_st.(folderz).stimVol,'type','glm','n=inf','r2cutoff=0.04','hdrlen=15','concatInfo',dAll{fi}.concatInfo);
    end
    %% Save

    dir = fullfile('~/proj/att_awe/analysis/data/decode');
    fname = fullfile(dir,sprintf('%s_inst.mat',sid));
    save(fname,'insts');
    
else
    %% Load
    disp('(cc_decode) Loading instances file');
    dir = fullfile('~/proj/att_awe/analysis/data/decode');
    fname = fullfile(dir,sprintf('%s_inst.mat',sid));
    load(fname);
    
end
%% Clear space??

    
%% Pre-process

if zscore
    for fi = 1:length(folders)
        inst = insts{fi};
        for iii = 1:length(inst)
            cinst = inst{iii};
            cinst.classify.instances = preprocessInstances(cinst.classify.instances);
            inst{iii} = cinst;
        end
        insts{fi} = inst;
    end
end

%% Full graphs by folder
if folders 
    for fi = 1:length(folders)
        folder = folders{fi};
        folderz = sprintf('f%s',folder);
        inst = insts{fi};
        %% Task: Contrast
        % two sides, so chance is 50%
        inst_con = inst;
        for ii = 1:length(inst_con)
            inst_con{ii}.classify.instances = inst_con{ii}.classify.instances(1:2);
        end
    %     out_con = leaveOneOut(inst_con,'permutation=1','balancByBootSt=1');
%         out_con = leaveOneOut(inst_con,'balancByRemovI=1');
        out_con = leaveOneOut(inst_con);
        clear inst_con

        %% Task: Coherence
        inst_coh = inst;
        for ii = 1:length(inst_coh)
            inst_coh{ii}.classify.instances = inst_coh{ii}.classify.instances(3:4);
        end
    %     out_coh = leaveOneOut(inst_coh,'permutation=1','balancByBootSt=1');
%         out_coh = leaveOneOut(inst_coh,'balancByRemovI=1');
        out_coh = leaveOneOut(inst_coh);
        clear inst_coh

        %% Clear
        clear inst

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

    %     dir = fullfile('~/proj/att_awe/analysis/data/decode');
    %     if ~isdir(dir), mkdir(dir); end
    %     fname = fullfile(dir,sprintf('%s_ddata.mat',sid));
    %     save(fname,'decode_data');

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
        fname = fullfile(dir,sprintf('decode_ROIs_%s',folderz));
        print(fname,'-dpdf');
    end
end

%% Stack folders and make a single graph for all folders

ainst = {};
for ri = 1:length(neural.shortROIs)
    ninst = struct;
    ninst.classify.instances = insts{1}{ri}.classify.instances;
    for fi = 2:length(folders)
        for ii = 1:length(insts{1}{ri}.classify.instances)
            ninst.classify.instances{ii} = [ninst.classify.instances{ii} ; insts{fi}{ri}.classify.instances{ii}];
        end
    end
    ninst.name = insts{1}{ri}.name;
    ainst{end+1} = ninst;
end

%% Run
if all
    %% Task: Contrast
    % two sides, so chance is 50%
    inst_con = ainst;
    for ii = 1:length(inst_con)
        inst_con{ii}.classify.instances = inst_con{ii}.classify.instances(1:2);
    end
    %     out_con = leaveOneOut(inst_con,'permutation=1','balancByBootSt=1');
%     out_con = leaveOneOut(inst_con,'balancByRemovI=1');
        out_con = leaveOneOut(inst_con,'uselibsvm=1');
    clear inst_con

    %% Task: Coherence
    inst_coh = ainst;
    for ii = 1:length(inst_coh)
        inst_coh{ii}.classify.instances = inst_coh{ii}.classify.instances(3:4);
    end
    %     out_coh = leaveOneOut(inst_coh,'permutation=1','balancByBootSt=1');
%     out_coh = leaveOneOut(inst_coh,'balancByRemovI=1');
        out_coh = leaveOneOut(inst_coh,'uselibsvm=1');
    clear inst_coh

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

    %     dir = fullfile('~/proj/att_awe/analysis/data/decode');
    %     if ~isdir(dir), mkdir(dir); end
    %     fname = fullfile(dir,sprintf('%s_ddata.mat',sid));
    %     save(fname,'decode_data');

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
    fname = fullfile(dir,'decode_ROIs_z');
    print(fname,'-dpdf');
end