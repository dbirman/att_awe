function deconroidata( data,subj )

warning('This is old code, you use deconCohCon now');

%DECONCOHCON Deconvolve the cohxcon experiment
decondata = {};
for roinum = 1:20

    %% Concatenate Sessions if Necessary
    if iscell(data)
        disp('(roimodel) Concatenating what appear to be different sessions...');
        % to run the model we really only need:
        % data.tSeries
        % data.design
        % data.runtrans
        % but we need these to be concatenated across the different sessions.
        % We have to be a little careful here that we make sure we add the
        % lengths of each run correctly or we'll screw ourselves over.
        data_old = data;
        data = struct;
        data.ROIs = data_old{1}.ROIs;
        % (we can do this the slow way, perf doesn't really matter
        data.tSeries = cell(1,length(data_old{1}.ROIs));
        length_sofar = 0;
        for ri = 1:length(data.tSeries), data.tSeries{ri} = []; end
        data.design = [];
        data.runtrans = [];
        for si = 1:length(data_old)
            for ri = 1:length(data_old{si}.tSeries)
                data.tSeries{ri} = [data.tSeries{ri} data_old{si}.rtSeries{ri}];
            end
            % tweak the SV by adding 
            cdes = data_old{si}.design;
            cdes(:,1) = cdes(:,1) + length_sofar;
            data.design = [data.design ; cdes];
            data.runtrans = [data.runtrans ; data_old{si}.runtrans+length_sofar];
            length_sofar = length_sofar + length(data_old{si}.tSeries{1});
        end
    end

    %%
    % Split the design into TIMING, CONTRAST_x_COHERENCE, and TASK
    % This assumes that the design is already concatenated (i.e. can't deal
    % with {datas} structure

    task_idxs = logical(~(data.design(:,9)==0));
    cohxcon_idxs = logical(logical(data.design(:,8)==5).*logical(~task_idxs));
    timing_idxs = logical((data.design(:,8)~=5).*~task_idxs);

    taskdesign = data.design(task_idxs,:);
    cohxcondesign = data.design(cohxcon_idxs,:);
    timingdesign = data.design(timing_idxs,:);

    designs = {cohxcondesign,timingdesign,taskdesign};

    %% reduce space of cohxcon values (for task)

    %% deconvolve for each timecourse

    %% Timing
    %  - Ignore l vs. right since identical
    contrast = timingdesign(:,3); ucon = unique(contrast);
    coherence = timingdesign(:,6); ucoh = unique(coherence);
    timing = timingdesign(:,8); ut = unique(timing);
    sv = timingdesign(:,1);
    timing_sv = {};
    conidxs = [];
    cohidxs = [];
    timidxs = [];
    for coni = 1:length(ucon)
        for cohi = 1:length(ucoh)
            for ti = 1:length(ut)
                timing_sv{end+1} = sv(logical((contrast==ucon(coni)).*(coherence==ucoh(cohi)).*(timing==ut(ti))));
                conidxs(end+1) = ucon(coni);
                cohidxs(end+1) = ucoh(cohi);
                timidxs(end+1) = ut(ti);
            end
        end
    end
    concatInfo.runTransition = data.runtrans;
    curd = constructD(data.tSeries{roinum},timing_sv,0.5,30,concatInfo,'none','deconv',0);
    decon = getr2timecourse(curd.timecourse,curd.nhdr,curd.hdrlenTR,curd.scm,curd.framePeriod,curd.verbose);
    decon = rmfield(decon,'scm');
    decon = rmfield(decon,'covar');

    %% Save timing data
    time = logical((decon.time>4.5).*(decon.time<7.5));
    decondata{roinum}.time.conidxs = conidxs;
    decondata{roinum}.time.cohidxs = cohidxs;
    decondata{roinum}.time.timidxs = timidxs;
    decondata{roinum}.time.resp = decon.ehdr(:,time);

    %% CohxCon
    contrast = cohxcondesign(:,3); ucon = unique(contrast);
    coherence = cohxcondesign(:,6); ucoh = unique(coherence);
    sv = cohxcondesign(:,1);
    cohxcon_sv = {};
    stimNames = {};
    conidx = [];
    cohidx = [];
    for coni = 1:length(ucon)
        for cohi = 1:length(ucoh)
            cohxcon_sv{end+1} = sv(logical((contrast==ucon(coni)).*(coherence==ucoh(cohi))));
            stimNames{end+1} = sprintf('Contrast=%0.2f Coherence=%0.2f',ucon(coni),ucoh(cohi));
            conidx(end+1) = ucon(coni);
            cohidx(end+1) = ucoh(cohi);
        end
    end
    concatInfo.runTransition = data.runtrans;
    curd = constructD(data.tSeries{roinum},cohxcon_sv,0.5,30,concatInfo,'none','deconv',0);
    decon = getr2timecourse(curd.timecourse,curd.nhdr,curd.hdrlenTR,curd.scm,curd.framePeriod,curd.verbose);
    decon = rmfield(decon,'scm');
    decon = rmfield(decon,'covar');

    %% save cohxcon data
    decondata{roinum}.cc.conidxs = conidx;
    decondata{roinum}.cc.cohidxs = cohidx;
    decondata{roinum}.cc.resp = decon.ehdr(:,time);
end

save(fullfile(datafolder,sprintf('%sdecon.mat',subj)),'decondata');