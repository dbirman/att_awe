function deconvolveEffects_att( data, roiname ,subj)
%DECONCOHCON Helper function designed to pull out the different conditions
%in the attention experiment as a way of seeing whether there is any
%reliable effect of attention on the contrast response functions. So we
%pull out the following conditions:
%
% Keep ROIs lateralized (1)
% Pull out pedestal contrast and coherence responses, by condition
% Pull out binned increments (assuming some average effect)


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
        if ~isempty(data_old{si}.tSeries)
            for ri = 1:length(data_old{si}.tSeries)
                data.tSeries{ri} = [data.tSeries{ri} data_old{si}.rtSeries25{ri}];
            end
            % tweak the SV by adding 
            cdes = data_old{si}.design;
            cdes(:,1) = cdes(:,1) + length_sofar;
            data.design = [data.design ; cdes];
            data.runtrans = [data.runtrans ; data_old{si}.runtrans+length_sofar];
            length_sofar = length_sofar + length(data_old{si}.tSeries{1});
        end
    end
end

%%
% if ~isempty(strfind(fit.mode,'droptiming'))
%     disp('%% All Timing Responses Dropped %%');
%     idxs = data.design(:,8)==5;
%     data.design = data.design(idxs,:);
% end

% % %% Concatenate multiple ROIs when requested
roinums = cellfun(@(x) ~isempty(strfind(x,roiname)),data.ROIs,'UniformOutput',false);
roinums = find([roinums{:}]);
% % 
% % if length(roinums)>2
% %     disp('Found more than two. Probably V3 found V3a/V3b. Taking first two.');
% %     roinums = roinums(1:2);
% % end
% % 
% % if length(roinums)==2
% %     disp('Using two ROIs concatenated');
% %     tSeries = [data.tSeries{roinums(1)} data.tSeries{roinums(2)}];
% %     % double the design
% %     design2 = data.design;
% %     design2(:,1) = design2(:,1)+length(data.tSeries{1}); % add to stimvols
% %     design = [data.design ; design2];
% %     % run transitions
% %     runtrans2 = data.runtrans+length(data.tSeries{1});
% %     runtrans = [data.runtrans ; runtrans2];
% % else
% %     roinum = roinums(1);
% %     tSeries = data.tSeries{roinum};
% %     design = data.design;
% %     runtrans = data.runtrans;
% % end

roinum = roinums(1);
tSeries = data.tSeries{roinum};
design = data.design;
runtrans = data.runtrans;

%%
% Split the design into TIMING, CONTRAST_x_COHERENCE, and TASK
% This assumes that the design is already concatenated (i.e. can't deal
% with {datas} structure

% [csv basecon lCon rCon basecoh lCoh rCoh tim task];
% for task, remove everything where no task was given (fixation task)
task_idxs = logical(~(design(:,9)==0));

taskdesign = design(task_idxs,:);

%% Remove base contrast effect
taskdesign(:,3) = taskdesign(:,3)-taskdesign(:,2);
taskdesign(:,4) = taskdesign(:,4)-taskdesign(:,2);

%% Use task to run the deconvolution model
%  1     2       3     4       5       6     7      8        9 
% SV  basecon   lcon  rcon  basecoh  lcoh  rcoh   timing   task

% these are the pedestal values
contrasts = [0.325 0.4 0.55 0.85]-taskdesign(1,2);
coherences = [0.15 0.3 0.45 0.6];
sv = taskdesign(:,1);

if roiname(1)=='l'
    stimopts = [4 7];
else
    stimopts = [3 6];
end

contrast = taskdesign(:,stimopts(1));
coherence = taskdesign(:,stimopts(2));
task = taskdesign(:,9);

% add delta flag (for contrast)
taskdesign(:,10) = zeros(size(taskdesign,1),1);
for ti = 1:size(taskdesign,1)
    taskdesign(ti,10) = ~any(contrasts==contrast(ti));
end
taskdesign(:,11) = zeros(size(taskdesign,1),1);
for ti = 1:size(taskdesign,1)
    taskdesign(ti,11) = ~any(coherences==coherence(ti));
end
% we need to set up a set of stimvols, stimnames which we'll use to drag
% out the right data. We're going to do this for contrast alone, and then
% for coherence alone, because otherwise we have a bajillion things to
% estimate which would be bad. 
% So we'll just get the contrast pedestals in each condition, and then the
% coherence pedestals in each condition, and see where that gets us.

%% Contrast
con_sv = {};
stimNames = {};
conidx = [];
taskidx = [];
deltaidx = [];

delta = taskdesign(:,10);
contrast_ = roundnearest(contrast,contrasts);

for taski = 1:2
    for di = 0:1
        for coni = 1:length(contrasts)
            con_sv{end+1} = sv(logical((delta==di) .* (task==taski) .* (contrasts(coni)==contrast_)));
            stimNames{end+1} = sprintf('Task=%i Delta? %i Contrast=%0.2f',taski,di,contrasts(coni));
            conidx(end+1) = contrasts(coni);
            taskidx(end+1) = taski;
            deltaidx(end+1) = di;
        end
    end
end

concatInfo.runTransition = runtrans;
curd = constructD(tSeries,con_sv,0.5,40,concatInfo,'none','deconv',0);
decon = getr2timecourse(curd.timecourse,curd.nhdr,curd.hdrlenTR,curd.scm,curd.framePeriod,curd.verbose);
decon = rmfield(decon,'scm');
decon = rmfield(decon,'covar');

%% Save
fname = fullfile(datafolder,sprintf('%s_deconEffects_att.mat',subj));
if exist(fname,'file')==2, load(fname); end
decondata.(roiname).cc.conidxs = conidx;
decondata.(roiname).cc.conStim = stimNames;
decondata.(roiname).cc.taskidx = taskidx;
decondata.(roiname).cc.deltaidx = deltaidx;
decondata.(roiname).cc.conresp = decon.ehdr;
save(fname,'decondata');

%% Plot Contrast
% h = figure; hold on
% 
% taskdisp = {'-o','--o'};
% concolor = brewermap(11,'PuOr');
% concolor = flipud(concolor(1:4,:));
% 
% ehdr = decon.ehdr(logical(deltaidx),:);
% conidx = conidx(logical(deltaidx));
% taskidx = taskidx(logical(deltaidx));
% stimNames = stimNames(logical(deltaidx));
% 
% for si = 1:length(stimNames)
%     subplot(1,2,taskidx(si)); hold on
%     plot(decon.time,ehdr(si,:),'o','MarkerSize',10,'MarkerFaceColor',concolor(find(conidx(si)==contrasts,1),:),'MarkerEdgeColor',[1 1 1]);
% end
% 
% tasks = {'Motion','Contrast'};
% for i = 1:2
%     subplot(1,2,i);
%     title(sprintf('Task %s',tasks{i}));
%     a = axis;
%     axis([0 15 a(3) a(4)]);
%     
%     if i==2, legend(stimNames); end
%     drawPublishAxis
% end
   

%% Coherence
coh_sv = {};
stimNames = {};
cohidx = [];
taskidx = [];
deltaidx = [];

delta = taskdesign(:,11);

coherence_ = roundnearest(coherence,coherences);

for taski = 1:2
    for di = 0:1
        for cohi = 1:length(coherences)
            coh_sv{end+1} = sv(logical((delta==di) .* (task==taski) .* (coherences(cohi)==coherence_)));
            stimNames{end+1} = sprintf('Task=%i Delta? %i Contrast=%0.2f',taski,di,coherences(cohi));
            cohidx(end+1) = coherences(cohi);
            taskidx(end+1) = taski;
            deltaidx(end+1) = di;
        end
    end
end

concatInfo.runTransition = runtrans;
curd = constructD(tSeries,coh_sv,0.5,40,concatInfo,'none','deconv',0);
decon = getr2timecourse(curd.timecourse,curd.nhdr,curd.hdrlenTR,curd.scm,curd.framePeriod,curd.verbose);
decon = rmfield(decon,'scm');
decon = rmfield(decon,'covar');

%% Save
fname = fullfile(datafolder,sprintf('%s_deconEffects_att.mat',subj));
if exist(fname,'file')==2, load(fname); end
decondata.(roiname).cc.cohidxs = cohidx;
decondata.(roiname).cc.cohStim = stimNames;
decondata.(roiname).cc.taskidx = taskidx;
decondata.(roiname).cc.deltaidx = deltaidx;
decondata.(roiname).cc.cohresp = decon.ehdr;
save(fname,'decondata');


%% Plot Coherence
% h = figure; hold on
% 
% taskdisp = {'-o','--o'};
% cohcolor = brewermap(11,'PuOr');
% cohcolor = cohcolor(8:11,:);
% 
% ehdr = decon.ehdr(logical(~deltaidx),:);
% cohidx = cohidx(logical(~deltaidx));
% taskidx = taskidx(logical(~deltaidx));
% stimNames = stimNames(logical(~deltaidx));
% 
% for si = 1:length(stimNames)
%     subplot(1,2,taskidx(si)); hold on
%     plot(decon.time,ehdr(si,:),'o','MarkerSize',10,'MarkerFaceColor',cohcolor(find(cohidx(si)==coherences,1),:),'MarkerEdgeColor',[1 1 1]);
% end
% 
% tasks = {'Motion','Contrast'};
% for i = 1:2
%     subplot(1,2,i);
%     title(sprintf('Task %s',tasks{i}));
%     a = axis;
%     axis([0 15 a(3) a(4)]);
%     
%     if i==2, legend(stimNames); end
%     drawPublishAxis
% end